#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
curate_assembly.py

This module provides functions for post-polishing steps of a genomic assembly,
including purging duplicate/haplotig sequences, reference-based correction,
and gap filling (for both long-read and short-read scenarios).

Updated on Sat Apr 11 2025

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""

import os
import sys
import pandas as pd
import shutil
from utilities import run_subprocess_cmd, pigz_compress, pigz_decompress, get_current_row_data


# --------------------------------------------------------------
# Purge duplicates using long reads
# --------------------------------------------------------------
def long_reads_purge_dups(curation_out_dir, polished_assembly, ont_raw_reads, illu_raw_f_reads, illu_raw_r_reads, pacbio_raw_reads, highest_mean_qual_long_reads_gz, sample_id, cpu_threads):
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
        highest_mean_qual_long_reads_gz (str): Path to gzipped high-quality long reads.
        sample_id (str): Sample identifier.
        cpu_threads (str): Number of CPU threads to use.

    Returns:
        str: Path to the purged assembly FASTA or input assembly if purge_dups fails/skipped.
    """
    # Define the working directory for purge_dups and change into it.
    cwd = os.getcwd()
    os.makedirs(curation_out_dir, exist_ok=True)
    os.chdir(curation_out_dir)

    # Default to returning the input assembly if purge_dups cannot proceed
    dup_purged_assembly = polished_assembly

    # Create a FOFN file referencing the (ONT or PacBio) long reads.
    if pd.notna(ont_raw_reads) and pd.isna(pacbio_raw_reads):
        pd_fofn = os.path.join(curation_out_dir, "ont_reads.fofn")
        with open(pd_fofn, "w") as fofn_out:
            fofn_out.write(highest_mean_qual_long_reads_gz + "\n")
    elif pd.notna(pacbio_raw_reads) and pd.isna(ont_raw_reads):
        pd_fofn = os.path.join(curation_out_dir, "pacbio_reads.fofn")
        with open(pd_fofn, "w") as fofn_out:
            fofn_out.write(highest_mean_qual_long_reads_gz + "\n")
    else:
        print("WARN:\tNo valid long reads (ONT or PacBio) provided for purge_dups. Skipping.")
        os.chdir(cwd)
        return dup_purged_assembly

    # Define output paths
    pd_json = os.path.join(curation_out_dir, "purge_dups_config.json")
    purged_output = os.path.join(curation_out_dir, f"{sample_id}_final_polish", "seqs", f"{sample_id}_final_polish.purged.fa")
    dup_purged_assembly = os.path.join(curation_out_dir, f"{sample_id}_purged.fasta")

    # Verify input files exist
    if not os.path.exists(polished_assembly):
        raise FileNotFoundError(f"ERROR:\tPolished assembly not found: {polished_assembly}")
    if not os.path.exists(highest_mean_qual_long_reads_gz):
        raise FileNotFoundError(f"ERROR:\tLong reads not found: {highest_mean_qual_long_reads_gz}")

    print(f"DEBUG - Using polished assembly: {polished_assembly}")
    print(f"DEBUG - Using long reads: {highest_mean_qual_long_reads_gz}")
    print(f"DEBUG - Writing FOFN to: {pd_fofn}")

    # Step 1: Generate a purge_dups configuration JSON if not already present.
    if os.path.exists(pd_json):
        print(f"SKIP:\tPurge Dupes JSON already exists: {pd_json}.")
    else:
        purge_dupes_config_cmd = ["python3", "/opt/conda/envs/EGAP_env/bin/pd_config.py", polished_assembly,
                                 pd_fofn, "-l", curation_out_dir, "-n", pd_json]
        result = run_subprocess_cmd(purge_dupes_config_cmd, shell_check=False)
        if result != 0:
            print(f"WARN:\tFailed to generate purge_dups config: pd_config.py returned {result}")
            os.chdir(cwd)
            return dup_purged_assembly

    # Step 2: Run the purge_dups pipeline if the final output does not exist.
    if os.path.exists(dup_purged_assembly):
        print(f"SKIP:\tDuplicate Purged Assembly already exists: {dup_purged_assembly}.")
    else:
        purge_dupes_cmd = ["python3", "/opt/conda/envs/EGAP_env/bin/run_purge_dups.py", pd_json, "/opt/conda/envs/EGAP_env/bin", sample_id, "-p", "bash"]
        result = run_subprocess_cmd(purge_dupes_cmd, shell_check=False)
        if result != 0:
            print(f"WARN:\tFailed to run purge_dups: run_purge_dups.py returned {result}")
            os.chdir(cwd)
            return dup_purged_assembly

        # Check if purge_dups produced an output (assuming it creates purged.fa)
        if os.path.exists(purged_output):
            shutil.move(purged_output, dup_purged_assembly)
        else:
            raise FileNotFoundError (f"WARN:\tPurge_dups output not found: {purged_output}")

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
        str: Path to the gzipped RagTag-corrected assembly FASTA or input assembly if RagTag fails.
    """
    # Default to returning input assembly if RagTag fails
    final_assembly = dup_purged_assembly

    # Decompress inputs if necessary and validate
    if ".gz" in dup_purged_assembly:
        dup_purged_assembly = pigz_decompress(dup_purged_assembly, cpu_threads)
        if not os.path.exists(dup_purged_assembly) or dup_purged_assembly.endswith(".gz"):
            print(f"ERROR:\tFailed to decompress duplicate-purged assembly: {dup_purged_assembly}")
            return final_assembly
    if ".gz" in ref_seq:
        ref_seq_unzipped = ref_seq.replace(".gz", "")
        if not os.path.exists(ref_seq_unzipped):
            print(f"NOTE:\tUnzipping reference sequence: {ref_seq}")
            ref_seq_unzipped = pigz_decompress(ref_seq, cpu_threads)
            if not os.path.exists(ref_seq_unzipped):
                print(f"ERROR:\tFailed to decompress reference sequence: {ref_seq}")
                return final_assembly
        ref_seq = ref_seq_unzipped

    # Validate input files
    if not os.path.exists(dup_purged_assembly):
        raise FileNotFoundError(f"ERROR:\tDuplicate-purged assembly not found: {dup_purged_assembly}")
    if not os.path.exists(ref_seq):
        raise FileNotFoundError(f"ERROR:\tReference sequence not found: {ref_seq}")

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
            result = run_subprocess_cmd(ragpatch_cmd, shell_check=False)
            if result != 0 or not os.path.exists(patch_fasta):
                print(f"WARN:\tRagTag patch failed: ragtag.py returned {result}")
                return final_assembly

        # Copy the final patched FASTA to the final assembly name
        if os.path.exists(patch_fasta):
            shutil.copyfile(patch_fasta, ragtag_ref_assembly)
        else:
            raise FileNotFoundError(f"ERROR:\tRagTag patch output not found: {patch_fasta}")

    # Compress the final ragtag'ed assembly
    if ".gz" not in ragtag_ref_assembly:
        ragtag_ref_assembly_gz = pigz_compress(ragtag_ref_assembly, cpu_threads)
        if not os.path.exists(ragtag_ref_assembly_gz):
            raise FileNotFoundError(f"ERROR:\tFailed to compress RagTag output: {ragtag_ref_assembly_gz}")
    else:
        ragtag_ref_assembly_gz = ragtag_ref_assembly

    # Clean up temporary unzipped reference file
    if not ref_seq.endswith(".gz") and os.path.exists(ref_seq):
        print(f"NOTE:\tRemoving temporary reference sequence: {ref_seq}")
        os.remove(ref_seq)

    print(f"PASS:\tCompleted RagTag reference curation: {ragtag_ref_assembly_gz}.")
    return ragtag_ref_assembly_gz


# --------------------------------------------------------------
# Close gaps using ONT reads
# --------------------------------------------------------------
def ont_tgs_gapcloser(curation_out_dir, ragtag_ref_assembly, highest_mean_qual_long_reads_gz, sample_id, cpu_threads):
    """Close gaps in an assembly using TGS-GapCloser with ONT reads.

    Applies TGS-GapCloser to fill gaps in the assembly and produces a scaffolded output.

    Args:
        curation_out_dir (str): Directory for curation output files.
        ragtag_ref_assembly (str): Path to the RagTag-corrected assembly FASTA.
        highest_mean_qual_long_reads_gz (str): Path to gzipped high-quality ONT reads.
        sample_id (str): Sample identifier.
        cpu_threads (str): Number of CPU threads to use.

    Returns:
        str: Path to the gzipped gap-filled assembly FASTA or input assembly if gap closing fails.
    """
    # Default to returning input assembly if TGS-GapCloser fails
    final_assembly = ragtag_ref_assembly

    # Validate input files
    if not os.path.exists(ragtag_ref_assembly):
        print(f"ERROR:\tRagTag-corrected assembly not found: {ragtag_ref_assembly}")
        return final_assembly
    if not os.path.exists(highest_mean_qual_long_reads_gz):
        raise FileNotFoundError (f"ERROR:\tONT reads not found: {highest_mean_qual_long_reads_gz}")

    # Decompress ragtag_ref_assembly if necessary
    assembly_input = ragtag_ref_assembly
    if ".gz" in ragtag_ref_assembly:
        assembly_input = pigz_decompress(ragtag_ref_assembly, cpu_threads)
        if not os.path.exists(assembly_input) or assembly_input.endswith(".gz"):
            raise FileNotFoundError(f"ERROR:\tFailed to decompress RagTag assembly: {ragtag_ref_assembly}")

    # Prepare TGS-GapCloser directory and output paths.
    tgs_gapcloser_dir = os.path.join(curation_out_dir, "tgs_gapcloser")
    os.makedirs(tgs_gapcloser_dir, exist_ok=True)
    original_cwd = os.getcwd()
    os.chdir(tgs_gapcloser_dir)

    tgs_gapcloser_prefix = os.path.join(tgs_gapcloser_dir, "purged_gapclosed")
    tgs_gapcloser_output = os.path.join(tgs_gapcloser_dir, "purged_gapclosed.scaff_seqs")
    gap_filled_assembly = os.path.join(curation_out_dir, f"{sample_id}_EGAP_final_curated.fasta")

    # Run TGS-GapCloser if the output doesn't exist.
    if os.path.exists(tgs_gapcloser_output):
        print(f"SKIP:\tTGS-GapCloser Assembly already exists: {tgs_gapcloser_output}.")
    else:
        tgs_gapcloser_cmd = ["tgsgapcloser",
                             "--scaff", assembly_input,
                             "--reads", highest_mean_qual_long_reads_gz,
                             "--output", tgs_gapcloser_prefix,
                             "--thread", str(cpu_threads),
                             "--ne"]
        result = run_subprocess_cmd(tgs_gapcloser_cmd, shell_check=False)
        if result != 0:
            print(f"WARN:\tTGS-GapCloser failed: tgsgapcloser returned {result}")
            os.chdir(original_cwd)
            return final_assembly

    os.chdir(original_cwd)

    # Check if a final gap-filled assembly already exists; if not, copy over.
    if os.path.exists(gap_filled_assembly):
        print(f"SKIP:\tGap-Filled Assembly already exists: {gap_filled_assembly}.")
    else:
        if os.path.exists(tgs_gapcloser_output):
            shutil.copyfile(tgs_gapcloser_output, gap_filled_assembly)
        else:
            raise FileNotFoundError(f"ERROR:\tTGS-GapCloser output not found: {tgs_gapcloser_output}")
            return final_assembly

    # Compress the final gap-filled assembly.
    gap_filled_assembly_gz = gap_filled_assembly + ".gz"
    if not os.path.exists(gap_filled_assembly_gz):
        pigz_compress(gap_filled_assembly, cpu_threads)
        if not os.path.exists(gap_filled_assembly_gz):
            raise FileNotFoundError(f"ERROR:\tFailed to compress gap-filled assembly: {gap_filled_assembly_gz}")
            return final_assembly
    print(f"PASS:\tCompleted TGS-Gapcloser curation: {gap_filled_assembly_gz}.")
    return gap_filled_assembly_gz


# --------------------------------------------------------------
# Close gaps using Illumina reads
# --------------------------------------------------------------
def illu_abyss_sealer(curation_dir, ragtag_ref_assembly, illu_f_dedup_gz, illu_r_dedup_gz, cpu_threads):
    """Close gaps in an assembly using ABySS Sealer with Illumina reads.

    Applies ABySS Sealer with multiple k-mer sizes to fill gaps in the assembly.

    Args:
        curation_dir (str): Directory for curation output files.
        ragtag_ref_assembly (str): Path to the RagTag-corrected assembly FASTA.
        illu_f_dedup_gz (str): Path to gzipped deduplicated forward Illumina reads.
        illu_r_dedup_gz (str): Path to gzipped deduplicated reverse Illumina reads.
        cpu_threads (str): Number of CPU threads to use.

    Returns:
        str: Path to the gzipped gap-filled assembly FASTA.
    """
    # Prepare output paths.
    abyss_output_dir = os.path.join(curation_dir, "abyss_sealer")
    os.makedirs(abyss_output_dir, exist_ok=True)
    output_prefix = os.path.join(abyss_output_dir, "abyss_sealer_work")
    sealer_output_file = f"{output_prefix}-scaffold.fa"
    renamed_sealer_output_file = os.path.join(curation_dir, "abyss_sealer_final_curated.fasta")
    renamed_sealer_output_file_gz = renamed_sealer_output_file + ".gz"

    # If the final output already exists, skip.
    if os.path.isfile(renamed_sealer_output_file_gz):
        print(f"SKIP:\tABySS Sealer output file already exists: {renamed_sealer_output_file_gz}.")
        return renamed_sealer_output_file_gz
    if os.path.isfile(sealer_output_file):
        print(f"SKIP:\tABySS Sealer output file already exists: {sealer_output_file}.")
        if not os.path.exists(renamed_sealer_output_file):
            shutil.move(sealer_output_file, renamed_sealer_output_file)
    else:
        # Validate input files
        if not os.path.exists(ragtag_ref_assembly):
            print(f"ERROR:\tRagTag-corrected assembly not found: {ragtag_ref_assembly}")
            return ragtag_ref_assembly
        if not os.path.exists(illu_f_dedup_gz) or not os.path.exists(illu_r_dedup_gz):
            print(f"ERROR:\tIllumina deduplicated reads not found: {illu_f_dedup_gz}, {illu_r_dedup_gz}")
            return ragtag_ref_assembly

        # Example k-mer sizes for gap closing.
        kmer_sizes = [55, 75, 95]
        abyss_sealer_cmd = ["abyss-sealer",
                            "-o", output_prefix,
                            "-S", ragtag_ref_assembly,
                            "-L", "400",
                            "-j", str(cpu_threads),
                            "-b", "500M"]
        # Add each k-mer size.
        for k in kmer_sizes:
            abyss_sealer_cmd.extend(["-k", str(k)])
        # Add Illumina read paths.
        abyss_sealer_cmd.extend([illu_f_dedup_gz, illu_r_dedup_gz])

        # Run ABySS Sealer.
        result = run_subprocess_cmd(abyss_sealer_cmd, shell_check=False)
        if result != 0 or not os.path.exists(sealer_output_file):
            print(f"ERROR:\tABySS Sealer failed: abyss-sealer returned {result}")
            return ragtag_ref_assembly
        shutil.move(sealer_output_file, renamed_sealer_output_file)

    # Compress the final gap-filled assembly.
    if not os.path.exists(renamed_sealer_output_file_gz):
        gap_filled_assembly_gz = pigz_compress(renamed_sealer_output_file, cpu_threads)
        if not os.path.exists(gap_filled_assembly_gz):
            print(f"ERROR:\tFailed to compress gap-filled assembly: {gap_filled_assembly_gz}")
            return ragtag_ref_assembly
    else:
        gap_filled_assembly_gz = renamed_sealer_output_file_gz

    print(f"PASS:\tCompleted Abyss-Sealer curation: {gap_filled_assembly_gz}.")
    return gap_filled_assembly_gz


# --------------------------------------------------------------
# Orchestrate assembly curation
# --------------------------------------------------------------
def curate_assembly(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    # Read the CSV file and filter to the row corresponding to the sample of interest.
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

    if pd.notna(ont_sra) and pd.isna(ont_raw_reads):
        ont_raw_reads = os.path.join(species_dir, "ONT", f"{ont_sra}.fastq.gz")
    if pd.notna(illumina_sra) and pd.isna(illumina_f_raw_reads) and pd.isna(illumina_r_raw_reads):
        illumina_f_raw_reads = os.path.join(species_dir, "Illumina", f"{illumina_sra}_1.fastq.gz")
        illumina_r_raw_reads = os.path.join(species_dir, "Illumina", f"{illumina_sra}_2.fastq.gz")    
    if pd.notna(pacbio_sra) and pd.isna(pacbio_raw_reads):
        pacbio_raw_reads = os.path.join(species_dir, "PacBio", f"{pacbio_sra}.fastq.gz")
    if pd.notna(ref_seq_gca) and pd.isna(ref_seq):
        ref_seq = os.path.join(species_dir, "RefSeq", f"{species_id}_{ref_seq_gca}_RefSeq.fasta.gz")

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
    if (pd.isna(illumina_f_raw_reads) or pd.isna(illumina_sra)) and (pd.isna(illumina_r_raw_reads) or pd.isna(illumina_sra)) and (pd.isna(ont_raw_reads) or pd.isna(ont_sra)) and (pd.isna(pacbio_raw_reads) or pd.isna(pacbio_sra)):
        print("SKIP:\tNo valid reads provided, required for assembly curation.")
        return None
    
    # Create output directory for curation
    species_dir = os.path.join(output_dir, species_id)    
    sample_dir = os.path.join(species_dir, sample_id)    
    curation_out_dir = os.path.join(sample_dir, "curated_assembly")
    os.makedirs(curation_out_dir, exist_ok=True)
    os.chdir(curation_out_dir)
    
    # Set and validate polished assembly path
    polished_assembly_gz = os.path.join(sample_dir, f"{sample_id}_final_polish.fasta.gz")
    polished_assembly = polished_assembly_gz.replace(".gz", "")
    if not os.path.exists(polished_assembly_gz) and not os.path.exists(polished_assembly):
        print(f"ERROR:\tPolished assembly not found at: {polished_assembly_gz} or {polished_assembly}")
        return None
    if os.path.exists(polished_assembly_gz) and not os.path.exists(polished_assembly):
        print(f"NOTE:\tUnzipping {polished_assembly_gz}")
        polished_assembly = pigz_decompress(polished_assembly_gz, cpu_threads)
        if not os.path.exists(polished_assembly):
            print(f"ERROR:\tFailed to decompress polished assembly: {polished_assembly_gz}")
            return None

    # Check for and remove any existing polished assembly copy in curated_assembly
    curated_polished_copy = os.path.join(curation_out_dir, f"{sample_id}_final_polish.fasta")
    if os.path.exists(curated_polished_copy):
        print(f"NOTE:\tRemoving existing polished assembly copy: {curated_polished_copy}")
        os.remove(curated_polished_copy)

    print(f"DEBUG - polished_assembly_gz - {polished_assembly_gz}")
    print(f"DEBUG - polished_assembly - {polished_assembly}")

    # Set Illumina deduplicated read paths only if Illumina reads are present
    illu_dedup_f_reads = None
    illu_dedup_r_reads = None
    if pd.notna(illumina_f_raw_reads) and pd.notna(illumina_r_raw_reads):
        illu_dedup_f_reads = os.path.join(species_dir, "Illumina", f"{species_id}_illu_forward_dedup.fastq.gz")
        illu_dedup_r_reads = os.path.join(species_dir, "Illumina", f"{species_id}_illu_reverse_dedup.fastq.gz")
        if not os.path.exists(illu_dedup_f_reads) or not os.path.exists(illu_dedup_r_reads):
            print(f"ERROR:\tIllumina deduplicated reads not found: {illu_dedup_f_reads}, {illu_dedup_r_reads}")
            return None

    print(f"DEBUG - illu_dedup_f_reads - {illu_dedup_f_reads}")
    print(f"DEBUG - illu_dedup_r_reads - {illu_dedup_r_reads}")

    # Set long-read paths
    highest_mean_qual_long_reads_gz = None
    highest_mean_qual_long_reads = None
    if pd.notna(ont_raw_reads):
        print("DEBUG - ONT RAW READS EXIST!")
        highest_mean_qual_long_reads_gz = os.path.join(species_dir, "ONT", f"{species_id}_ONT_highest_mean_qual_long_reads.fastq.gz")
        highest_mean_qual_long_reads = highest_mean_qual_long_reads_gz.replace(".gz", "")
    elif pd.notna(pacbio_raw_reads):
        print("DEBUG - PACBIO RAW READS EXIST!")
        highest_mean_qual_long_reads_gz = os.path.join(species_dir, "PacBio", f"{species_id}_PacBio_highest_mean_qual_long_reads.fastq.gz")
        highest_mean_qual_long_reads = highest_mean_qual_long_reads_gz.replace(".gz", "")

    print(f"DEBUG - highest_mean_qual_long_reads_gz - {highest_mean_qual_long_reads_gz}")

    # Validate and ensure long-read .gz file exists
    if (pd.notna(pacbio_raw_reads) or pd.notna(ont_raw_reads)) and highest_mean_qual_long_reads_gz:
        if not os.path.exists(highest_mean_qual_long_reads_gz):
            if os.path.exists(highest_mean_qual_long_reads):
                print(f"NOTE:\tCompressing {highest_mean_qual_long_reads} to {highest_mean_qual_long_reads_gz}")
                highest_mean_qual_long_reads_gz = pigz_compress(highest_mean_qual_long_reads, cpu_threads)
            else:
                print(f"ERROR:\tHighest mean quality long reads not found: {highest_mean_qual_long_reads_gz} or {highest_mean_qual_long_reads}")
                return None
        print(f"DEBUG - Using long reads .gz: {highest_mean_qual_long_reads_gz}")

    print(f"DEBUG - highest_mean_qual_long_reads - {highest_mean_qual_long_reads}")

    # -------------------------------------------------------------------------
    # Step 1: Purge duplicates (haplotigs) if ONT or PacBio reads are present.
    # -------------------------------------------------------------------------
    if pd.notna(ont_raw_reads) or pd.notna(pacbio_raw_reads):
        print("Purging Haplotigs using Long Reads (ONT or PacBio)...")
        dup_purged_assembly = long_reads_purge_dups(curation_out_dir,
                                                    polished_assembly,
                                                    ont_raw_reads,
                                                    illumina_f_raw_reads,
                                                    illumina_r_raw_reads,
                                                    pacbio_raw_reads,
                                                    highest_mean_qual_long_reads_gz,
                                                    sample_id, cpu_threads)
    else:
        print("SKIP:\tPurge Duplicates; No Long Reads (ONT or PacBio) Provided.")
        dup_purged_assembly = polished_assembly

    # -------------------------------------------------------------------------
    # Step 2: Perform RagTag correction/scaffolding if a reference sequence is provided.
    # -------------------------------------------------------------------------
    if pd.notna(ref_seq):
        print("Running RagTag with Reference Sequence...")
        ragtag_ref_assembly = ref_seq_ragtag(dup_purged_assembly, ref_seq, curation_out_dir,
                                             sample_id, cpu_threads, ram_gb)
    else:
        print("SKIP:\tRagTag; No Reference Sequence Provided.")
        ragtag_ref_assembly = dup_purged_assembly

    # -------------------------------------------------------------------------
    # Step 3: Gap closing with ONT or Illumina reads (if available).
    # -------------------------------------------------------------------------
    if pd.notna(ont_raw_reads):
        print("Running TGS-GapCloser with ONT Reads...")
        curated_assembly_gz = ont_tgs_gapcloser(curation_out_dir, ragtag_ref_assembly,
                                                highest_mean_qual_long_reads_gz,
                                                sample_id, cpu_threads)
    elif pd.notna(illumina_f_raw_reads) and pd.notna(illumina_r_raw_reads):
        print("Running ABySS-Sealer with Illumina Reads...")
        curated_assembly_gz = illu_abyss_sealer(curation_out_dir, ragtag_ref_assembly,
                                                illu_dedup_f_reads, illu_dedup_r_reads,
                                                cpu_threads)
    elif pd.notna(pacbio_raw_reads):
        print("SKIP:\tGap Closing; PacBio reads were used or no other option provided.")
        curated_assembly_gz = ragtag_ref_assembly
    else:
        # If no reads for gap closing, we just keep the RagTag or purged assembly.
        print("SKIP:\tNo available read type for gap closing.")
        curated_assembly_gz = ragtag_ref_assembly

    # Validate curated assembly
    if not os.path.exists(curated_assembly_gz):
        print(f"ERROR:\tCurated assembly not produced: {curated_assembly_gz}")
        return None

    final_curated_assembly_gz = os.path.join(sample_dir, f"{sample_id}_final_curated.fasta.gz")

    if os.path.exists(final_curated_assembly_gz):
        print(f"SKIP:\tFinal curated assembly already exists: {final_curated_assembly_gz}")
    elif os.path.exists(curated_assembly_gz):
        shutil.move(curated_assembly_gz, final_curated_assembly_gz)
        print(f"PASS:\tMoved curated assembly to {final_curated_assembly_gz}")
    else:
        print(f"ERROR:\tCannot move {curated_assembly_gz} to {final_curated_assembly_gz}: curated assembly does not exist")
        return None

    # Clean up temporary unzipped polished assembly file in sample_dir
    if os.path.exists(polished_assembly):
        print(f"NOTE:\tRemoving temporary polished assembly: {polished_assembly}")
        os.remove(polished_assembly)

    print(f"PASS:\tCuration complete for {polished_assembly_gz} -> {final_curated_assembly_gz}")
    os.chdir(sample_dir)
    
    return final_curated_assembly_gz


if __name__ == "__main__":
    # Usage information if incorrect arguments are provided.
    if len(sys.argv) != 6:
        print(
            "Usage: python3 curate_assembly.py <sample_id> <input_csv> "
            "<output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)

    # Run the polish_assembly function with CLI arguments.
    gap_filled_assembly_gz = curate_assembly(sys.argv[1],       # sample_id
                                             sys.argv[2],       # input_csv
                                             sys.argv[3],       # output_dir
                                             str(sys.argv[4]),  # cpu_threads
                                             str(sys.argv[5]))   # ram_gb