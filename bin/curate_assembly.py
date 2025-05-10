#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
curate_assembly.py

This module provides functions for post-polishing steps of a genomic assembly,
including purging duplicate/haplotig sequences, reference-based correction,
and gap filling (for both long-read and short-read scenarios).

Updated on Mon May 12 2025

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""

import os
import sys
import pandas as pd
import shutil
from Bio import SeqIO
from utilities import run_subprocess_cmd, get_current_row_data


# --------------------------------------------------------------
# Validate FASTA file
# --------------------------------------------------------------
def validate_fasta(file_path):
    """Validate that a FASTA file exists, is non-empty, and contains valid nucleotide sequences.

    Args:
        file_path (str): Path to the FASTA file.

    Returns:
        bool: True if valid, False otherwise.
    """
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
def long_reads_purge_dups(curation_out_dir, polished_assembly, ont_raw_reads, illu_raw_f_reads, illu_raw_r_reads, pacbio_raw_reads, highest_mean_qual_long_reads, sample_id, cpu_threads):
    """Purge duplicates and haplotigs from an assembly using purge_dups.

    Uses long reads (ONT or PacBio) to identify and remove haplotigs and duplicates,
    producing a refined assembly.

    Args:
        curation_out_dir (str): Directory for curation output files.
        polished_assembly (str): Path to the polished assembly FASTA.
        ont_raw_reads (str or float): Path to ONT raw reads or NaN if absent.
        illu_raw_f_reads (str or float): Path to Illumina forward raw reads or NaN.
        illu_raw_r_reads (str or float): Path to Illumina reverse raw reads or NaN.
        pacbio_raw_reads (str or float): Path to PacBio raw reads or NaN.
        highest_mean_qual_long_reads (str): Path to high-quality long reads.
        sample_id (str): Sample identifier.
        cpu_threads (str): Number of CPU threads to use.

    Returns:
        str: Path to the purged assembly FASTA or input assembly if purge_dups fails/skipped.
    """
    # Define the working directory for purge_dups and change into it
    cwd = os.getcwd()
    os.makedirs(curation_out_dir, exist_ok=True)
    os.chdir(curation_out_dir)

    # Default to returning the input assembly if purge_dups cannot proceed
    dup_purged_assembly = polished_assembly

    # Create a FOFN file referencing the (ONT or PacBio) long reads
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

    # Define output paths
    pd_json = os.path.join(curation_out_dir, "purge_dups_config.json")
    purged_output = os.path.join(curation_out_dir, f"{sample_id}_final_polish_assembly", "seqs", f"{sample_id}_final_polish_assembly.purged.fa")
    dup_purged_assembly = os.path.join(curation_out_dir, f"{sample_id}_purged.fasta")

    # Verify input files exist and are valid
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

    # Step 1: Generate a purge_dups configuration JSON if not already present
    if os.path.exists(pd_json):
        print(f"SKIP:\tPurge Dupes JSON already exists: {pd_json}. Overwriting to ensure fresh run.")
        os.remove(pd_json)
    # Try to locate pd_config.py dynamically
    pd_config_path = shutil.which("pd_config.py") or "/opt/conda/envs/EGAP_env/bin/pd_config.py"
    if not os.path.exists(pd_config_path):
        print(f"WARN:\tpd_config.py not found at {pd_config_path}. Skipping purge_dups.")
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

    # Step 2: Run the purge_dups pipeline if the final output does not exist
    if os.path.exists(dup_purged_assembly):
        print(f"SKIP:\tDuplicate Purged Assembly already exists: {dup_purged_assembly}. Overwriting to ensure fresh run.")
        os.remove(dup_purged_assembly)
    run_purge_dups_path = shutil.which("run_purge_dups.py") or "/opt/conda/envs/EGAP_env/bin/run_purge_dups.py"
    if not os.path.exists(run_purge_dups_path):
        print(f"WARN:\trun_purge_dups.py not found at {run_purge_dups_path}. Skipping purge_dups.")
        os.chdir(cwd)
        return dup_purged_assembly
    purge_dupes_cmd = ["python3", run_purge_dups_path, pd_json, "/opt/conda/envs/EGAP_env/bin", sample_id, "-p", "bash"]
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

    # Restore the original working directory
    os.chdir(cwd)
    print(f"PASS:\tCompleted purge_dups curation: {dup_purged_assembly}.")
    return dup_purged_assembly


# --------------------------------------------------------------
# Correct assembly with reference sequence
# --------------------------------------------------------------
def ref_seq_ragtag(dup_purged_assembly, ref_seq, curation_out_dir, sample_id, cpu_threads, ram_gb):
    """Perform reference-guided correction using RagTag.

    Scaffolds, corrects, and patches the assembly against a reference sequence
    using RagTag.

    Args:
        dup_purged_assembly (str): Path to the duplicate-purged assembly FASTA.
        ref_seq (str or float): Path to the reference FASTA or NaN if absent.
        curation_out_dir (str): Directory for curation output files.
        sample_id (str): Sample identifier.
        cpu_threads (str): Number of CPU threads to use.
        ram_gb (str): RAM in GB (not always used by RagTag).

    Returns:
        str: Path to the RagTag-corrected assembly FASTA or input assembly if RagTag fails.
    """
    # Default to returning input assembly if RagTag fails
    final_assembly = dup_purged_assembly

    # Validate input files
    if not validate_fasta(dup_purged_assembly):
        print(f"ERROR:\tInvalid duplicate-purged assembly: {dup_purged_assembly}")
        return final_assembly
    if not os.path.exists(ref_seq):
        print(f"ERROR:\tReference sequence not found: {ref_seq}")
        return final_assembly

    # Define directories and filenames for RagTag output
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

    # If final output already exists, skip
    if os.path.exists(ragtag_ref_assembly):
        print(f"SKIP:\tRagTag Reference-Corrected Assembly already exists: {ragtag_ref_assembly}.")
    else:
        # Scaffold step
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

        # Correct step
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

        # Patch step
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

        # Copy the final patched FASTA to the final assembly name
        if os.path.exists(patch_fasta):
            shutil.copy(patch_fasta, ragtag_ref_assembly)
            print(f"DEBUG - Copied RagTag output: {patch_fasta} to {ragtag_ref_assembly}")
        else:
            print(f"ERROR:\tRagTag patch output not found: {patch_fasta}")
            return final_assembly

    print(f"PASS:\tCompleted RagTag reference curation: {ragtag_ref_assembly}.")
    return ragtag_ref_assembly


# --------------------------------------------------------------
# Close gaps using ONT or PacBio reads
# --------------------------------------------------------------
def long_reads_tgs_gapcloser(curation_out_dir, ragtag_ref_assembly, highest_mean_qual_long_reads, sample_id, cpu_threads):
    """Close gaps in an assembly using TGS-GapCloser with ONT or PacBio reads.

    Applies TGS-GapCloser to fill gaps in the assembly and produces a scaffolded output.

    Args:
        curation_out_dir (str): Directory for curation output files.
        ragtag_ref_assembly (str): Path to the RagTag-corrected assembly FASTA.
        highest_mean_qual_long_reads (str): Path to high-quality ONT or PacBio reads.
        sample_id (str): Sample identifier.
        cpu_threads (str): Number of CPU threads to use.

    Returns:
        str: Path to the gap-filled assembly FASTA or input assembly if gap closing fails.
    """
    # Default to returning input assembly if TGS-GapCloser fails
    final_assembly = ragtag_ref_assembly

    # Validate input files
    if not validate_fasta(ragtag_ref_assembly):
        print(f"ERROR:\tInvalid RagTag-corrected assembly: {ragtag_ref_assembly}")
        return final_assembly
    if not os.path.exists(highest_mean_qual_long_reads):
        print(f"ERROR:\tLong reads not found: {highest_mean_qual_long_reads}")
        return final_assembly
        
    # Prepare TGS-GapCloser directory and output paths
    tgs_gapcloser_dir = os.path.join(curation_out_dir, "tgs_gapcloser")
    os.makedirs(tgs_gapcloser_dir, exist_ok=True)
    original_cwd = os.getcwd()
    os.chdir(tgs_gapcloser_dir)

    tgs_gapcloser_prefix = os.path.join(tgs_gapcloser_dir, "purged_gapclosed")
    tgs_gapcloser_output = os.path.join(tgs_gapcloser_dir, "purged_gapclosed.scaff_seqs")
    gap_filled_assembly = os.path.join(curation_out_dir, f"{sample_id}_EGAP_final_curated.fasta")

    # Run TGS-GapCloser if the output doesn't exist
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

    # Check if a final gap-filled assembly already exists; if not, copy over
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
# Close gaps using Illumina reads
# --------------------------------------------------------------
def illu_abyss_sealer(curation_dir, ragtag_ref_assembly, illu_f_dedup, illu_r_dedup, cpu_threads):
    """Close gaps in an assembly using ABySS Sealer with Illumina reads.

    Applies ABySS Sealer with multiple k-mer sizes to fill gaps in the assembly.

    Args:
        curation_dir (str): Directory for curation output files.
        ragtag_ref_assembly (str): Path to the RagTag-corrected assembly FASTA.
        illu_f_dedup (str): Path to deduplicated forward Illumina reads.
        illu_r_dedup (str): Path to deduplicated reverse Illumina reads.
        cpu_threads (str): Number of CPU threads to use.

    Returns:
        str: Path to the gap-filled assembly FASTA or input assembly if ABySS fails.
    """
    # Prepare output paths
    abyss_output_dir = os.path.join(curation_dir, "abyss_sealer")
    os.makedirs(abyss_output_dir, exist_ok=True)
    output_prefix = os.path.join(abyss_output_dir, "abyss_sealer_work")
    sealer_output_file = f"{output_prefix}_scaffold.fa"
    renamed_sealer_output_file = os.path.join(curation_dir, "abyss_sealer_final_curated.fasta")

    print(f"DEBUG - ABySS-Sealer output paths: {sealer_output_file}, {renamed_sealer_output_file}")

    # If the final output already exists, skip
    if os.path.isfile(renamed_sealer_output_file):
        print(f"SKIP:\tABySS Sealer output file already exists: {renamed_sealer_output_file}.")
        return renamed_sealer_output_file
    if os.path.isfile(sealer_output_file):
        print(f"SKIP:\tABySS Sealer output file already exists: {sealer_output_file}.")
        if not os.path.exists(renamed_sealer_output_file):
            print(f"DEBUG - Copying {sealer_output_file} to {renamed_sealer_output_file}")
            shutil.copy(sealer_output_file, renamed_sealer_output_file)
    else:
        # Validate input files
        print(f"DEBUG - Checking inputs: ragtag_ref_assembly={ragtag_ref_assembly}, illu_f_dedup={illu_f_dedup}, illu_r_dedup={illu_r_dedup}")
        if not validate_fasta(ragtag_ref_assembly):
            print(f"ERROR:\tInvalid RagTag-corrected assembly: {ragtag_ref_assembly}")
            return ragtag_ref_assembly
        if not os.path.exists(illu_f_dedup) or not os.path.exists(illu_r_dedup):
            print(f"ERROR:\tIllumina deduplicated reads not found: {illu_f_dedup}, {illu_r_dedup}")
            return ragtag_ref_assembly

        # Example k-mer sizes for gap closing
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

        # Run ABySS Sealer
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
    """Orchestrate the curation of a polished genomic assembly.

    Performs duplicate purging, reference-guided correction, and gap closing using
    available read types (ONT, PacBio, or Illumina).

    Args:
        sample_id (str): Sample identifier.
        input_csv (str): Path to metadata CSV file.
        output_dir (str): Directory for output files.
        cpu_threads (int or str): Number of CPU threads to use.
        ram_gb (int or str): Available RAM in GB.

    Returns:
        str: Path to the final curated assembly FASTA or None if curation fails.
    """
    # Read the CSV file and filter to the row corresponding to the sample of interest
    input_df = pd.read_csv(input_csv)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]  # Convert to Series (single row)

    # Identify read paths, reference, and BUSCO lineage info from CSV
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

    species_dir = os.path.join(output_dir, species_id)
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
    
    # Check if only reference is provided and skip
    if (pd.isna(illumina_f_raw_reads) and pd.isna(illumina_sra)) and \
       (pd.isna(illumina_r_raw_reads) and pd.isna(illumina_sra)) and \
       (pd.isna(ont_raw_reads) and pd.isna(ont_sra)) and \
       (pd.isna(pacbio_raw_reads) and pd.isna(pacbio_sra)):
        print("SKIP:\tNo valid reads provided, required for assembly comparison.")
        return None
    
    # Set and validate polished assembly path
    polished_assembly = os.path.join(sample_dir, f"{sample_id}_final_polish_assembly.fasta")
    if not validate_fasta(polished_assembly):
        print(f"ERROR:\tInvalid polished assembly: {polished_assembly}")
        return None

    print(f"DEBUG - polished_assembly - {polished_assembly}")

    # Set Illumina deduplicated read paths only if Illumina reads are present
    illu_dedup_f_reads = None
    illu_dedup_r_reads = None
    if pd.notna(illumina_f_raw_reads) and pd.notna(illumina_r_raw_reads):
        illu_dedup_f_reads = os.path.join(species_dir, "Illumina", f"{species_id}_illu_forward_dedup.fastq")
        illu_dedup_r_reads = os.path.join(species_dir, "Illumina", f"{species_id}_illu_reverse_dedup.fastq")
        print(f"DEBUG - Checking Illumina deduplicated reads: {illu_dedup_f_reads}, {illu_dedup_r_reads}")
        if not os.path.exists(illu_dedup_f_reads) or not os.path.exists(illu_dedup_r_reads):
            print(f"ERROR:\tIllumina deduplicated reads not found: {illu_dedup_f_reads}, {illu_dedup_r_reads}")
            return None

    print(f"DEBUG - illu_dedup_f_reads - {illu_dedup_f_reads}")
    print(f"DEBUG - illu_dedup_r_reads - {illu_dedup_r_reads}")

    # Set long-read paths (ONT or PacBio), prefer prefiltered, fallback to raw
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

    # Ensure work directory output
    starting_work_dir = os.getcwd()
    if "work" not in starting_work_dir:
        current_work_dir = curation_out_dir
    else:
        current_work_dir = starting_work_dir
    os.chdir(current_work_dir)

    # -------------------------------------------------------------------------
    # Step 1: Purge duplicates (haplotigs) if ONT or PacBio reads are present
    # -------------------------------------------------------------------------
    print(f"DEBUG - Starting purge duplicates for {sample_id}")
    if pd.notna(ont_raw_reads) or pd.notna(pacbio_raw_reads):
        print("Purging Haplotigs using Long Reads (ONT or PacBio)...")
        dup_purged_assembly = long_reads_purge_dups(current_work_dir,
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

    # -------------------------------------------------------------------------
    # Step 2: Perform RagTag correction/scaffolding if a reference sequence is provided
    # -------------------------------------------------------------------------
    print(f"DEBUG - Starting RagTag correction for {sample_id}")
    if pd.notna(ref_seq) and os.path.exists(ref_seq):
        print("Running RagTag with Reference Sequence...")
        ragtag_ref_assembly = ref_seq_ragtag(dup_purged_assembly, ref_seq, current_work_dir,
                                             sample_id, cpu_threads, ram_gb)
        if not os.path.exists(ragtag_ref_assembly):
            print(f"WARN:\tRagTag correction failed to produce output: {ragtag_ref_assembly}. Falling back to purged assembly.")
            ragtag_ref_assembly = dup_purged_assembly
    else:
        print(f"SKIP:\tRagTag; No valid Reference Sequence Provided: {ref_seq}")
        ragtag_ref_assembly = dup_purged_assembly

    print(f"DEBUG - RagTag output: {ragtag_ref_assembly}")

    # -------------------------------------------------------------------------
    # Step 3: Gap closing with ONT, PacBio, or Illumina reads (if available)
    # -------------------------------------------------------------------------
    print(f"DEBUG - Starting gap closing for {sample_id}")
    if pd.notna(ont_raw_reads) or pd.notna(pacbio_raw_reads):
        print("Running TGS-GapCloser with Long Reads (ONT or PacBio)...")
        curated_assembly = long_reads_tgs_gapcloser(current_work_dir, ragtag_ref_assembly,
                                                   highest_mean_qual_long_reads,
                                                   sample_id, cpu_threads)
    elif pd.notna(illumina_f_raw_reads) and pd.notna(illumina_r_raw_reads):
        print("Running ABySS-Sealer with Illumina Reads...")
        curated_assembly = illu_abyss_sealer(current_work_dir, ragtag_ref_assembly,
                                             illu_dedup_f_reads, illu_dedup_r_reads,
                                             cpu_threads)
    else:
        print("SKIP:\tNo available read type for gap closing.")
        curated_assembly = ragtag_ref_assembly

    # Validate curated assembly
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

    # Avoid removing polished assembly to preserve it for debugging
    print(f"NOTE:\tKeeping polished assembly for debugging: {polished_assembly}")

    print(f"PASS:\tCuration complete for {polished_assembly} -> {final_curated_assembly}")
    os.chdir(sample_dir)
    
    return final_curated_assembly


if __name__ == "__main__":
    # Usage information if incorrect arguments are provided
    if len(sys.argv) != 6:
        print(
            "Usage: python3 curate_assembly.py <sample_id> <input_csv> "
            "<output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)

    # Run the curate_assembly function with CLI arguments
    final_curated_assembly = curate_assembly(sys.argv[1],       # sample_id
                                            sys.argv[2],       # input_csv
                                            sys.argv[3],       # output_dir
                                            str(sys.argv[4]),  # cpu_threads
                                            str(sys.argv[5]))  # ram_gb
