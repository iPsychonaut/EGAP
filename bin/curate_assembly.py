#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
curate_assembly.py

Post-polishing curation: purge duplicates, reference correction (RagTag),
and gap closing (ONT/PacBio via TGS-GapCloser or Illumina via ABySS-Sealer).

Created on Wed Aug 16 2023

Updated on Wed Sept 3 2025

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""

import os
import sys
import pandas as pd
import shutil
from pathlib import Path
from Bio import SeqIO
from utilities import run_subprocess_cmd, get_current_row_data


# --------------------------------------------------------------
# Validate FASTA file
# --------------------------------------------------------------
def validate_fasta(file_path):
    if not file_path:
        print("ERROR:\tFASTA path is empty/None")
        return False
    if not os.path.exists(file_path):
        print(f"ERROR:\tFASTA file not found: {file_path}")
        return False
    if os.path.getsize(file_path) < 100:
        print(f"ERROR:\tFASTA file is suspiciously small: {file_path}")
        return False
    try:
        with open(file_path, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                if not record.seq:
                    print(f"ERROR:\tFASTA file contains empty sequences: {file_path}")
                    return False
                if not all(c.upper() in "ATCGN" for c in record.seq):
                    print(f"ERROR:\tFASTA file contains non-nucleotide sequences: {file_path}")
                    return False
                return True
    except Exception as e:
        print(f"ERROR:\tInvalid FASTA format in {file_path}: {str(e)}")
        return False
    return False


# --------------------------------------------------------------
# Purge duplicates using long reads
# --------------------------------------------------------------
def long_reads_purge_dups(curation_out_dir, polished_assembly, ont_raw_reads,
                          illu_raw_f_reads, illu_raw_r_reads, pacbio_raw_reads,
                          highest_mean_qual_long_reads, sample_id, cpu_threads):
    cwd = os.getcwd()
    os.makedirs(curation_out_dir, exist_ok=True)
    os.chdir(curation_out_dir)

    dup_purged_assembly = polished_assembly

    if pd.notna(ont_raw_reads) and pd.isna(pacbio_raw_reads):
        pd_fofn = os.path.join(curation_out_dir, "ont_reads.fofn")
        with open(pd_fofn, "w") as fofn_out:
            fofn_out.write(highest_mean_qual_long_reads + "\n")
    elif pd.notna(pacbio_raw_reads) and pd.isna(ont_raw_reads):
        pd_fofn = os.path.join(curation_out_dir, "pacbio_reads.fofn")
        with open(pd_fofn, "w") as fofn_out:
            fofn_out.write(highest_mean_qual_long_reads + "\n")
    else:
        print("WARN:\tNo valid long reads (ONT or PacBio) provided for purge_dups. Skipping.")
        os.chdir(cwd)
        return dup_purged_assembly

    pd_json = os.path.join(curation_out_dir, "purge_dups_config.json")
    purged_output = os.path.join(curation_out_dir, f"{sample_id}_final_polish_assembly", "seqs",
                                 f"{sample_id}_final_polish_assembly.purged.fa")
    dup_purged_assembly = os.path.join(curation_out_dir, f"{sample_id}_purged.fasta")

    if not validate_fasta(polished_assembly):
        print(f"ERROR:\tInvalid polished assembly: {polished_assembly}")
        os.chdir(cwd)
        return polished_assembly
    if not os.path.exists(highest_mean_qual_long_reads):
        print(f"ERROR:\tLong reads not found: {highest_mean_qual_long_reads}")
        os.chdir(cwd)
        return polished_assembly

    print(f"DEBUG - Using polished assembly: {polished_assembly}")
    print(f"DEBUG - Using long reads: {highest_mean_qual_long_reads}")
    print(f"DEBUG - Writing FOFN to: {pd_fofn}")

    if os.path.exists(pd_json):
        print(f"SKIP:\tPurge Dupes JSON already exists: {pd_json}. Overwriting to ensure fresh run.")
        os.remove(pd_json)

    pd_config_path = shutil.which("pd_config.py")
    if not pd_config_path:
        print("WARN:\tpd_config.py not found on PATH. Skipping purge_dups.")
        os.chdir(cwd)
        return dup_purged_assembly

    purge_dupes_config_cmd = ["python3", pd_config_path, polished_assembly,
                              pd_fofn, "-l", curation_out_dir, "-n", pd_json]
    print(f"DEBUG - Running purge_dups config: {' '.join(purge_dupes_config_cmd)}")
    result = run_subprocess_cmd(purge_dupes_config_cmd, shell_check=False)
    if result != 0:
        print(f"WARN:\tFailed to generate purge_dups config: pd_config.py returned {result}")
        os.chdir(cwd)
        return dup_purged_assembly

    if os.path.exists(dup_purged_assembly):
        print(f"SKIP:\tDuplicate Purged Assembly already exists: {dup_purged_assembly}. Overwriting to ensure fresh run.")
        os.remove(dup_purged_assembly)

    run_purge_dups_path = shutil.which("run_purge_dups.py")
    if not run_purge_dups_path:
        print("WARN:\trun_purge_dups.py not found on PATH. Skipping purge_dups.")
        os.chdir(cwd)
        return dup_purged_assembly

    purge_dupes_cmd = ["python3", run_purge_dups_path, pd_json, os.path.dirname(run_purge_dups_path),
                       sample_id, "-p", "bash"]
    print(f"DEBUG - Running purge_dups: {' '.join(purge_dupes_cmd)}")
    result = run_subprocess_cmd(purge_dupes_cmd, shell_check=False)
    if result != 0:
        print(f"WARN:\tFailed to run purge_dups: run_purge_dups.py returned {result}")

    if os.path.exists(purged_output):
        print(f"DEBUG - Copying purge_dups output: {purged_output} to {dup_purged_assembly}")
        shutil.copy(purged_output, dup_purged_assembly)
    else:
        print(f"WARN:\tPurge_dups output not found: {purged_output}")
        os.chdir(cwd)
        return dup_purged_assembly

    os.chdir(cwd)
    print(f"PASS:\tCompleted purge_dups curation: {dup_purged_assembly}.")
    return dup_purged_assembly


# --------------------------------------------------------------
# Correct assembly with reference sequence
# --------------------------------------------------------------
def ref_seq_ragtag(dup_purged_assembly, ref_seq, curation_out_dir, sample_id, cpu_threads, ram_gb):
    final_assembly = dup_purged_assembly

    if not validate_fasta(dup_purged_assembly):
        print(f"ERROR:\tInvalid duplicate-purged assembly: {dup_purged_assembly}")
        return final_assembly
    if not (ref_seq and os.path.exists(ref_seq)):
        print(f"ERROR:\tReference sequence not found: {ref_seq}")
        return final_assembly

    ragtag_output_dir = os.path.join(curation_out_dir, "ragtag_work")
    ragtag_ref_assembly = os.path.join(ragtag_output_dir, f"{sample_id}_ragtag_final.fasta")

    ragtag_scaff_output_dir = os.path.join(ragtag_output_dir, f"{sample_id}_scaffolded")
    scaff_fasta = os.path.join(ragtag_scaff_output_dir, "ragtag.scaffold.fasta")
    os.makedirs(ragtag_scaff_output_dir, exist_ok=True)

    ragtag_corr_output_dir = os.path.join(ragtag_output_dir, f"{sample_id}_corrected")
    corr_fasta = os.path.join(ragtag_corr_output_dir, "ragtag.correct.fasta")
    os.makedirs(ragtag_corr_output_dir, exist_ok=True)

    ragtag_patched_output_dir = os.path.join(ragtag_output_dir, f"{sample_id}_patched")
    patch_fasta = os.path.join(ragtag_patched_output_dir, "ragtag.patch.fasta")
    os.makedirs(ragtag_patched_output_dir, exist_ok=True)

    if os.path.exists(ragtag_ref_assembly):
        print(f"SKIP:\tRagTag Reference-Corrected Assembly already exists: {ragtag_ref_assembly}.")
    else:
        if os.path.exists(scaff_fasta):
            print(f"SKIP:\tRagTag Scaffolded Assembly already exists: {scaff_fasta}.")
        else:
            ragscaf_cmd = ["ragtag.py", "scaffold", ref_seq, dup_purged_assembly,
                           "-o", ragtag_scaff_output_dir,
                           "-t", str(cpu_threads), "-C", "-u"]
            print(f"DEBUG - Running RagTag scaffold: {' '.join(ragscaf_cmd)}")
            result = run_subprocess_cmd(ragscaf_cmd, shell_check=False)
            if result != 0 or not os.path.exists(scaff_fasta):
                print(f"WARN:\tRagTag scaffold failed: ragtag.py returned {result}")
                return final_assembly

        if os.path.exists(corr_fasta):
            print(f"SKIP:\tRagTag Corrected Assembly already exists: {corr_fasta}.")
        else:
            ragcorr_cmd = ["ragtag.py", "correct", ref_seq, scaff_fasta,
                           "-o", ragtag_corr_output_dir,
                           "-t", str(cpu_threads), "-u"]
            print(f"DEBUG - Running RagTag correct: {' '.join(ragcorr_cmd)}")
            result = run_subprocess_cmd(ragcorr_cmd, shell_check=False)
            if result != 0 or not os.path.exists(corr_fasta):
                print(f"WARN:\tRagTag correct failed: ragtag.py returned {result}")
                return final_assembly

        if os.path.exists(patch_fasta):
            print(f"SKIP:\tRagTag Patched Assembly already exists: {patch_fasta}.")
        else:
            ragpatch_cmd = ["ragtag.py", "patch", ref_seq, corr_fasta,
                            "-o", ragtag_patched_output_dir,
                            "-t", str(cpu_threads), "-u"]
            print(f"DEBUG - Running RagTag patch: {' '.join(ragpatch_cmd)}")
            result = run_subprocess_cmd(ragpatch_cmd, shell_check=False)
            if result != 0 or not os.path.exists(patch_fasta):
                print(f"WARN:\tRagTag patch failed: ragtag.py returned {result}")
                return final_assembly

        if os.path.exists(patch_fasta):
            shutil.copy(patch_fasta, ragtag_ref_assembly)
            print(f"DEBUG - Copied RagTag output: {patch_fasta} to {ragtag_ref_assembly}")
        else:
            print(f"ERROR:\tRagTag patch output not found: {patch_fasta}")
            return final_assembly

    print(f"PASS:\tCompleted RagTag reference curation: {ragtag_ref_assembly}.")
    return ragtag_ref_assembly


# --------------------------------------------------------------
# Close gaps with long reads
# --------------------------------------------------------------
def long_reads_tgs_gapcloser(curation_out_dir, ragtag_ref_assembly, highest_mean_qual_long_reads, sample_id, cpu_threads):
    final_assembly = ragtag_ref_assembly

    if not validate_fasta(ragtag_ref_assembly):
        print(f"ERROR:\tInvalid RagTag-corrected assembly: {ragtag_ref_assembly}")
        return final_assembly
    if not (highest_mean_qual_long_reads and os.path.exists(highest_mean_qual_long_reads)):
        print(f"ERROR:\tLong reads not found: {highest_mean_qual_long_reads}")
        return final_assembly

    tgs_gapcloser_dir = os.path.join(curation_out_dir, "tgs_gapcloser")
    os.makedirs(tgs_gapcloser_dir, exist_ok=True)
    original_cwd = os.getcwd()
    os.chdir(tgs_gapcloser_dir)

    tgs_gapcloser_prefix = os.path.join(tgs_gapcloser_dir, "purged_gapclosed")
    tgs_gapcloser_output = os.path.join(tgs_gapcloser_dir, "purged_gapclosed.scaff_seqs")
    gap_filled_assembly = os.path.join(curation_out_dir, f"{sample_id}_EGAP_final_curated.fasta")

    if os.path.exists(tgs_gapcloser_output):
        print(f"SKIP:\tTGS-GapCloser Assembly already exists: {tgs_gapcloser_output}. Overwriting to ensure fresh run.")
        os.remove(tgs_gapcloser_output)

    tgs_gapcloser_cmd = ["tgsgapcloser",
                         "--scaff", ragtag_ref_assembly,
                         "--reads", highest_mean_qual_long_reads,
                         "--output", tgs_gapcloser_prefix,
                         "--thread", str(cpu_threads),
                         "--ne"]
    print(f"DEBUG - Running TGS-GapCloser: {' '.join(tgs_gapcloser_cmd)}")
    result = run_subprocess_cmd(tgs_gapcloser_cmd, shell_check=False)
    if result != 0:
        print(f"WARN:\tTGS-GapCloser failed: tgsgapcloser returned {result}")
        os.chdir(original_cwd)
        return final_assembly

    os.chdir(original_cwd)

    if os.path.exists(gap_filled_assembly):
        print(f"SKIP:\tGap-Filled Assembly already exists: {gap_filled_assembly}. Overwriting to ensure fresh run.")
        os.remove(gap_filled_assembly)
    if os.path.exists(tgs_gapcloser_output):
        shutil.copy(tgs_gapcloser_output, gap_filled_assembly)
        print(f"DEBUG - Copied TGS-GapCloser output: {tgs_gapcloser_output} to {gap_filled_assembly}")
    else:
        print(f"WARN:\tTGS-GapCloser output not found: {tgs_gapcloser_output}")
        return final_assembly

    print(f"PASS:\tCompleted TGS-GapCloser curation: {gap_filled_assembly}.")
    return gap_filled_assembly


# --------------------------------------------------------------
# Close gaps with Illumina reads
# --------------------------------------------------------------
def illu_abyss_sealer(curation_dir, ragtag_ref_assembly, illu_f_dedup, illu_r_dedup, cpu_threads):
    abyss_output_dir = os.path.join(curation_dir, "abyss_sealer")
    os.makedirs(abyss_output_dir, exist_ok=True)
    output_prefix = os.path.join(abyss_output_dir, "abyss_sealer_work")
    sealer_output_file = f"{output_prefix}_scaffold.fa"
    renamed_sealer_output_file = os.path.join(curation_dir, "abyss_sealer_final_curated.fasta")

    print(f"DEBUG - ABySS-Sealer output paths: {sealer_output_file}, {renamed_sealer_output_file}")

    if os.path.isfile(renamed_sealer_output_file):
        print(f"SKIP:\tABySS Sealer output file already exists: {renamed_sealer_output_file}.")
        return renamed_sealer_output_file
    if os.path.isfile(sealer_output_file):
        print(f"SKIP:\tABySS Sealer output file already exists: {sealer_output_file}.")
        if not os.path.exists(renamed_sealer_output_file):
            print(f"DEBUG - Copying {sealer_output_file} to {renamed_sealer_output_file}")
            shutil.copy(sealer_output_file, renamed_sealer_output_file)
    else:
        print(f"DEBUG - Checking inputs: ragtag_ref_assembly={ragtag_ref_assembly}, illu_f_dedup={illu_f_dedup}, illu_r_dedup={illu_r_dedup}")
        if not validate_fasta(ragtag_ref_assembly):
            print(f"ERROR:\tInvalid RagTag-corrected assembly: {ragtag_ref_assembly}")
            return ragtag_ref_assembly
        if not (illu_f_dedup and illu_r_dedup and os.path.exists(illu_f_dedup) and os.path.exists(illu_r_dedup)):
            print(f"ERROR:\tIllumina deduplicated reads not found: {illu_f_dedup}, {illu_r_dedup}")
            return ragtag_ref_assembly

        kmer_sizes = [55, 75, 95]
        abyss_sealer_cmd = ["abyss-sealer",
                            "-o", output_prefix,
                            "-S", ragtag_ref_assembly,
                            "-L", "400",
                            "-j", str(cpu_threads),
                            "-b", "500M"]
        for k in kmer_sizes:
            abyss_sealer_cmd.extend(["-k", str(k)])
        abyss_sealer_cmd.extend([illu_f_dedup, illu_r_dedup])

        print(f"DEBUG - Running ABySS-Sealer command: {' '.join(abyss_sealer_cmd)}")
        result = run_subprocess_cmd(abyss_sealer_cmd, shell_check=False)
        if result != 0 or not os.path.exists(sealer_output_file):
            print(f"WARN:\tABySS Sealer failed: abyss-sealer returned {result}")
            return ragtag_ref_assembly
        print(f"DEBUG - ABySS-Sealer produced output: {sealer_output_file}")
        print(f"DEBUG - Copying {sealer_output_file} to {renamed_sealer_output_file}")
        shutil.copy(sealer_output_file, renamed_sealer_output_file)

    print(f"PASS:\tCompleted Abyss-Sealer curation: {renamed_sealer_output_file}.")
    return renamed_sealer_output_file


# --------------------------------------------------------------
# Orchestrate assembly curation
# --------------------------------------------------------------
def curate_assembly(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    input_df = pd.read_csv(input_csv)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]

    illumina_sra = current_series["ILLUMINA_SRA"]
    illumina_f_raw_reads = current_series["ILLUMINA_RAW_F_READS"]
    illumina_r_raw_reads = current_series["ILLUMINA_RAW_R_READS"]
    ont_sra = current_series["ONT_SRA"]
    ont_raw_reads = current_series["ONT_RAW_READS"]
    pacbio_sra = current_series["PACBIO_SRA"]
    pacbio_raw_reads = current_series["PACBIO_RAW_READS"]
    ref_seq_gca = current_series["REF_SEQ_GCA"]
    ref_seq = current_series["REF_SEQ"]
    species_id = current_series["SPECIES_ID"]
    est_size = current_series["EST_SIZE"]

    # --- ABSOLUTIZE EVERYTHING ---
    output_dir_abs = str(Path(output_dir).resolve())
    species_dir = os.path.join(output_dir_abs, species_id)
    sample_dir = os.path.join(species_dir, sample_id)
    curation_out_dir = os.path.join(sample_dir, "curated_assembly")

    if pd.notna(ont_sra) and pd.isna(ont_raw_reads):
        ont_raw_reads = os.path.join(species_dir, "ONT", f"{ont_sra}.fastq")
    if pd.notna(illumina_sra) and pd.isna(illumina_f_raw_reads) and pd.isna(illumina_r_raw_reads):
        illumina_f_raw_reads = os.path.join(species_dir, "Illumina", f"{illumina_sra}_1.fastq")
        illumina_r_raw_reads = os.path.join(species_dir, "Illumina", f"{illumina_sra}_2.fastq")
    if pd.notna(pacbio_sra) and pd.isna(pacbio_raw_reads):
        pacbio_raw_reads = os.path.join(species_dir, "PacBio", f"{pacbio_sra}.fastq")
    if pd.notna(ref_seq_gca) and pd.isna(ref_seq):
        ref_seq = os.path.join(species_dir, "RefSeq", f"{species_id}_{ref_seq_gca}_RefSeq.fasta")

    # Resolve to absolute strings if present
    def _abs(p): return str(Path(p).resolve()) if isinstance(p, str) else p
    species_dir = _abs(species_dir)
    sample_dir = _abs(sample_dir)
    curation_out_dir = _abs(curation_out_dir)
    illumina_f_raw_reads = _abs(illumina_f_raw_reads) if pd.notna(illumina_f_raw_reads) else illumina_f_raw_reads
    illumina_r_raw_reads = _abs(illumina_r_raw_reads) if pd.notna(illumina_r_raw_reads) else illumina_r_raw_reads
    ont_raw_reads = _abs(ont_raw_reads) if pd.notna(ont_raw_reads) else ont_raw_reads
    pacbio_raw_reads = _abs(pacbio_raw_reads) if pd.notna(pacbio_raw_reads) else pacbio_raw_reads
    ref_seq = _abs(ref_seq) if pd.notna(ref_seq) else ref_seq

    print(f"DEBUG - illumina_sra - {illumina_sra}")
    print(f"DEBUG - illumina_f_raw_reads - {illumina_f_raw_reads}")
    print(f"DEBUG - illumina_r_raw_reads - {illumina_r_raw_reads}")
    print(f"DEBUG - ont_sra - {ont_sra}")
    print(f"DEBUG - ont_raw_reads - {ont_raw_reads}")
    print(f"DEBUG - pacbio_sra - {pacbio_sra}")
    print(f"DEBUG - pacbio_raw_reads - {pacbio_raw_reads}")
    print(f"DEBUG - ref_seq_gca - {ref_seq_gca}")
    print(f"DEBUG - ref_seq - {ref_seq}")   
    print(f"DEBUG - species_id - {species_id}")
    print(f"DEBUG - est_size - {est_size}")

    if (pd.isna(illumina_f_raw_reads) and pd.isna(illumina_sra)) and \
       (pd.isna(illumina_r_raw_reads) and pd.isna(illumina_sra)) and \
       (pd.isna(ont_raw_reads) and pd.isna(ont_sra)) and \
       (pd.isna(pacbio_raw_reads) and pd.isna(pacbio_sra)):
        print("SKIP:\tNo valid reads provided, required for assembly comparison.")
        return None

    # Polished assembly (absolute)
    polished_assembly = os.path.join(sample_dir, f"{sample_id}_final_polish_assembly.fasta")
    if not validate_fasta(polished_assembly):
        print(f"ERROR:\tInvalid polished assembly: {polished_assembly}")
        return None
    print(f"DEBUG - polished_assembly - {polished_assembly}")

    # Illumina dedup (absolute if present)
    illu_dedup_f_reads = None
    illu_dedup_r_reads = None
    if pd.notna(illumina_f_raw_reads) and pd.notna(illumina_r_raw_reads):
        illu_dedup_f_reads = os.path.join(species_dir, "Illumina", f"{species_id}_illu_forward_dedup.fastq")
        illu_dedup_r_reads = os.path.join(species_dir, "Illumina", f"{species_id}_illu_reverse_dedup.fastq")
        print(f"DEBUG - Checking Illumina deduplicated reads: {illu_dedup_f_reads}, {illu_dedup_r_reads}")
        if not (os.path.exists(illu_dedup_f_reads) and os.path.exists(illu_dedup_r_reads)):
            print(f"ERROR:\tIllumina deduplicated reads not found: {illu_dedup_f_reads}, {illu_dedup_r_reads}")
            return None

    print(f"DEBUG - illu_dedup_f_reads - {illu_dedup_f_reads}")
    print(f"DEBUG - illu_dedup_r_reads - {illu_dedup_r_reads}")

    # Long-read pick (absolute)
    highest_mean_qual_long_reads = None
    if pd.notna(ont_raw_reads):
        print("DEBUG - ONT RAW READS EXIST!")
        candidate = os.path.join(species_dir, "ONT", f"{species_id}_ONT_highest_mean_qual_long_reads.fastq")
        highest_mean_qual_long_reads = candidate if os.path.exists(candidate) else ont_raw_reads
    elif pd.notna(pacbio_raw_reads):
        print("DEBUG - PACBIO RAW READS EXIST!")
        candidate = os.path.join(species_dir, "PacBio", f"{species_id}_PacBio_highest_mean_qual_long_reads.fastq")
        highest_mean_qual_long_reads = candidate if os.path.exists(candidate) else pacbio_raw_reads
    print(f"DEBUG - highest_mean_qual_long_reads - {highest_mean_qual_long_reads}")

    # Work in absolute curation dir
    os.makedirs(curation_out_dir, exist_ok=True)
    os.chdir(curation_out_dir)
    print(f"DEBUG - CWD (curation_out_dir) - {os.getcwd()}")

    # 1) Purge duplicates (if long reads)
    print(f"DEBUG - Starting purge duplicates for {sample_id}")
    if pd.notna(ont_raw_reads) or pd.notna(pacbio_raw_reads):
        print("Purging Haplotigs using Long Reads (ONT or PacBio)...")
        dup_purged_assembly = long_reads_purge_dups(curation_out_dir,
                                                    polished_assembly,
                                                    ont_raw_reads,
                                                    illumina_f_raw_reads,
                                                    illumina_r_raw_reads,
                                                    pacbio_raw_reads,
                                                    highest_mean_qual_long_reads,
                                                    sample_id, cpu_threads)
        if not os.path.exists(dup_purged_assembly):
            print(f"WARN:\tPurge duplicates failed to produce output: {dup_purged_assembly}. Falling back to polished assembly.")
            dup_purged_assembly = polished_assembly
    else:
        print("SKIP:\tPurge Duplicates; No Long Reads (ONT or PacBio) Provided.")
        dup_purged_assembly = polished_assembly

    print(f"DEBUG - Purge duplicates output: {dup_purged_assembly}")

    # 2) RagTag reference correction (if ref provided)
    print(f"DEBUG - Starting RagTag correction for {sample_id}")
    if pd.notna(ref_seq) and os.path.exists(ref_seq):
        print("Running RagTag with Reference Sequence...")
        ragtag_ref_assembly = ref_seq_ragtag(dup_purged_assembly, ref_seq, curation_out_dir,
                                             sample_id, cpu_threads, ram_gb)
        if not os.path.exists(ragtag_ref_assembly):
            print(f"WARN:\tRagTag correction failed to produce output: {ragtag_ref_assembly}. Falling back to purged assembly.")
            ragtag_ref_assembly = dup_purged_assembly
    else:
        print(f"SKIP:\tRagTag; No valid Reference Sequence Provided: {ref_seq}")
        ragtag_ref_assembly = dup_purged_assembly

    print(f"DEBUG - RagTag output: {ragtag_ref_assembly}")

    # 3) Gap closing
    print(f"DEBUG - Starting gap closing for {sample_id}")
    if pd.notna(ont_raw_reads) or pd.notna(pacbio_raw_reads):
        print("Running TGS-GapCloser with Long Reads (ONT or PacBio)...")
        curated_assembly = long_reads_tgs_gapcloser(curation_out_dir, ragtag_ref_assembly,
                                                    highest_mean_qual_long_reads,
                                                    sample_id, cpu_threads)
    elif pd.notna(illumina_f_raw_reads) and pd.notna(illumina_r_raw_reads):
        print("Running ABySS-Sealer with Illumina Reads...")
        curated_assembly = illu_abyss_sealer(curation_out_dir, ragtag_ref_assembly,
                                             illu_dedup_f_reads, illu_dedup_r_reads,
                                             cpu_threads)
    else:
        print("SKIP:\tNo available read type for gap closing.")
        curated_assembly = ragtag_ref_assembly

    # Validate curated assembly; fall back to polished
    if not validate_fasta(curated_assembly):
        print(f"WARN:\tInvalid curated assembly: {curated_assembly}. Falling back to polished assembly.")
        curated_assembly = polished_assembly

    final_curated_assembly = os.path.join(sample_dir, f"{sample_id}_final_curated.fasta")
    if os.path.exists(final_curated_assembly):
        print(f"SKIP:\tFinal curated assembly already exists: {final_curated_assembly}. Overwriting to ensure fresh run.")
        os.remove(final_curated_assembly)
    if os.path.exists(curated_assembly):
        shutil.copy(curated_assembly, final_curated_assembly)
        print(f"PASS:\tCopied curated assembly to {final_curated_assembly}")
    else:
        print(f"WARN:\tCannot copy {curated_assembly} to {final_curated_assembly}: curated assembly does not exist. Copying polished assembly.")
        if os.path.exists(polished_assembly):
            shutil.copy(polished_assembly, final_curated_assembly)
            print(f"PASS:\tCopied polished assembly to {final_curated_assembly}")
        else:
            print(f"ERROR:\tPolished assembly not found: {polished_assembly}")
            return None

    print(f"NOTE:\tKeeping polished assembly for debugging: {polished_assembly}")
    print(f"PASS:\tCuration complete for {polished_assembly} -> {final_curated_assembly}")
    os.chdir(sample_dir)
    return final_curated_assembly


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python3 curate_assembly.py <sample_id> <input_csv> <output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)

    final_curated_assembly = curate_assembly(sys.argv[1], sys.argv[2], sys.argv[3], str(sys.argv[4]), str(sys.argv[5]))
