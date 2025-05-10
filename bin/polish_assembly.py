#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polish_assembly.py

This module processes genomic assembly data by applying polishing steps using Racon
for long reads (ONT or PacBio) and Pilon for Illumina reads. It handles input validation,
subprocess execution, and file management for assembly refinement.

Updated on Sun May 11 2025

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os
import sys
import shutil
import pandas as pd
from utilities import run_subprocess_cmd, get_current_row_data


# --------------------------------------------------------------
# Polish assembly with Racon using long reads
# --------------------------------------------------------------
def racon_polish_assembly(input_assembly, long_reads, racon_out_dir, sample_id, cpu_threads, iteration_count):
    """Polish a genomic assembly with Racon using long reads.

    Aligns long reads (ONT or PacBio) to the assembly using minimap2 and polishes
    with Racon for the specified iteration.

    Args:
        input_assembly (str): Path to the input assembly FASTA.
        long_reads (str): Path to the long reads FASTQ.
        racon_out_dir (str): Directory for polishing output.
        sample_id (str): Sample identifier.
        cpu_threads (str): Number of CPU threads to use.
        iteration_count (int): Current polishing iteration (e.g., 1 or 2).

    Returns:
        str: Path to the Racon-polished assembly FASTA or input assembly if Racon fails.
    """
    # Validate inputs
    if not os.path.exists(input_assembly):
        print(f"ERROR:\tInput assembly not found: {input_assembly}")
        return input_assembly
    if not os.path.exists(long_reads):
        print(f"ERROR:\tLong reads not found: {long_reads}")
        return input_assembly

    # Align reads to the assembly using minimap2, generating a PAF file
    racon_paf = os.path.join(racon_out_dir, f"racon_round{iteration_count}.paf")
    if os.path.exists(racon_paf):
        print(f"SKIP:\tRacon PAF {iteration_count} already exists: {racon_paf}.")
    else:
        minimap2_cmd = f"minimap2 -t {cpu_threads} -x map-ont {input_assembly} {long_reads} > {racon_paf}"
        print(f"DEBUG - Running minimap2: {minimap2_cmd}")
        result = run_subprocess_cmd(minimap2_cmd, shell_check=True)
        if result != 0 or not os.path.exists(racon_paf):
            print(f"WARN:\tminimap2 failed for iteration {iteration_count}. Skipping Racon.")
            return input_assembly

    # Perform polishing with Racon using the PAF alignment
    racon_assembly = os.path.join(racon_out_dir, f"{sample_id}_racon_polish_{iteration_count}.fasta")
    if os.path.exists(racon_assembly):
        print(f"SKIP:\tRacon Assembly {iteration_count} already exists: {racon_assembly}.")
    else:
        racon_cmd = f"racon -t {cpu_threads} {long_reads} {racon_paf} {input_assembly} > {racon_assembly}"
        print(f"DEBUG - Running Racon: {racon_cmd}")
        result = run_subprocess_cmd(racon_cmd, shell_check=True)
        if result != 0 or not os.path.exists(racon_assembly):
            print(f"WARN:\tRacon failed for iteration {iteration_count}. Returning input assembly.")
            return input_assembly

    return racon_assembly


# --------------------------------------------------------------
# Prepare alignment for Pilon polishing
# --------------------------------------------------------------
def pilon_prep(input_assembly, illu_f_dedup, illu_r_dedup, assembly_out_dir, cpu_threads):
    """Prepare BAM file for Pilon polishing using Illumina reads.

    Aligns deduplicated Illumina reads to the assembly using bwa-mem2, converts to BAM,
    sorts, and indexes for Pilon.

    Args:
        input_assembly (str): Path to the assembly FASTA.
        illu_f_dedup (str): Path to deduplicated forward Illumina FASTQ.
        illu_r_dedup (str): Path to deduplicated reverse Illumina FASTQ.
        assembly_out_dir (str): Directory for output files.
        cpu_threads (str): Number of CPU threads to use.

    Returns:
        str: Path to the sorted, indexed BAM file or None if preparation fails.
    """
    # Validate inputs
    if not os.path.exists(input_assembly):
        print(f"ERROR:\tInput assembly not found: {input_assembly}")
        return None
    if not os.path.exists(illu_f_dedup) or not os.path.exists(illu_r_dedup):
        print(f"ERROR:\tIllumina reads not found: {illu_f_dedup}, {illu_r_dedup}")
        return None

    # Construct output file names for SAM, unsorted BAM, and final sorted BAM
    pilon_bam = os.path.join(assembly_out_dir, os.path.basename(input_assembly).replace(".fasta", ".bam"))
    output_sam = pilon_bam.replace(".bam", ".sam")
    sorted_bam = pilon_bam.replace(".bam", "_sorted.bam")

    os.makedirs(assembly_out_dir, exist_ok=True)

    # Build a BWA index and map Illumina reads to the assembly (SAM)
    if not os.path.exists(output_sam):
        bwa_index_cmd = ["bwa-mem2", "index", input_assembly]
        print(f"DEBUG - Running BWA index: {' '.join(bwa_index_cmd)}")
        result = run_subprocess_cmd(bwa_index_cmd, shell_check=False)
        if result != 0:
            print(f"WARN:\tBWA index failed. Skipping Pilon prep.")
            return None

        bwa_cmd = f"bwa-mem2 mem -t {cpu_threads} {input_assembly} {illu_f_dedup} {illu_r_dedup} > {output_sam}"
        print(f"DEBUG - Running BWA: {bwa_cmd}")
        result = run_subprocess_cmd(bwa_cmd, shell_check=True)
        if result != 0 or not os.path.exists(output_sam):
            print(f"WARN:\tBWA mapping failed. Skipping Pilon prep.")
            return None

    if os.path.exists(output_sam) and os.path.getsize(output_sam) < 1000:
        print(f"WARN:\tGenerated SAM file is suspiciously small: {output_sam}")
        return None

    # Convert SAM to an unsorted BAM file
    if not os.path.exists(sorted_bam):
        samview_cmd = f"samtools view -@ {cpu_threads} -S -b {output_sam} > {sorted_bam}"
        print(f"DEBUG - Running samtools view: {samview_cmd}")
        result = run_subprocess_cmd(samview_cmd, shell_check=True)
        if result != 0 or not os.path.exists(sorted_bam):
            print(f"WARN:\tsamtools view failed. Skipping Pilon prep.")
            return None

    # Sort and index the BAM file
    if not os.path.exists(pilon_bam):
        bamsort_cmd = ["bamtools", "sort", "-in", sorted_bam, "-out", pilon_bam]
        print(f"DEBUG - Running bamtools sort: {' '.join(bamsort_cmd)}")
        result = run_subprocess_cmd(bamsort_cmd, shell_check=False)
        if result != 0 or not os.path.exists(pilon_bam):
            print(f"WARN:\tbamtools sort failed. Skipping Pilon prep.")
            return None

        samtools_index_cmd = ["samtools", "index", pilon_bam]
        print(f"DEBUG - Running samtools index: {' '.join(samtools_index_cmd)}")
        result = run_subprocess_cmd(samtools_index_cmd, shell_check=False)
        if result != 0 or not os.path.exists(pilon_bam + ".bai"):
            print(f"WARN:\tsamtools index failed. Skipping Pilon prep.")
            return None

    if os.path.exists(pilon_bam) and os.path.getsize(pilon_bam) < 1000:
        print(f"WARN:\tSorted BAM file is suspiciously small: {pilon_bam}")
        return None

    return pilon_bam


# --------------------------------------------------------------
# Polish assembly with Pilon using Illumina reads
# --------------------------------------------------------------
def pilon_polish(best_assembly, second_racon_assembly, pilon_bam, assembly_out_dir, sample_id, cpu_threads, ram_gb):
    """Polish an assembly with Pilon using aligned Illumina reads.

    Runs Pilon to refine the assembly based on Illumina read alignments, producing
    a polished FASTA file.

    Args:
        best_assembly (str): Path to the initial assembly FASTA.
        second_racon_assembly (str): Path to the Racon-polished assembly.
        pilon_bam (str): Path to the sorted, indexed BAM file.
        assembly_out_dir (str): Directory for output files.
        sample_id (str): Sample identifier.
        cpu_threads (str): Number of CPU threads to use.
        ram_gb (str): RAM in GB for Pilon's Java process.

    Returns:
        str: Path to the Pilon-polished assembly FASTA or Racon assembly if Pilon fails.
    """
    # Validate inputs
    if not os.path.exists(second_racon_assembly):
        print(f"ERROR:\tRacon assembly not found: {second_racon_assembly}")
        return second_racon_assembly
    if not os.path.exists(pilon_bam):
        print(f"ERROR:\tPilon BAM not found: {pilon_bam}")
        return second_racon_assembly

    # Prepare output filenames for Pilon run
    pilon_out_prefix = f"{sample_id}_pilon_assembly"
    pilon_out_dir = assembly_out_dir
    os.makedirs(pilon_out_dir, exist_ok=True)

    pilon_assembly = os.path.join(pilon_out_dir, f"{sample_id}_pilon_assembly.fasta")
    if os.path.exists(pilon_assembly):
        print(f"SKIP:\tPilon Polished Assembly already exists: {pilon_assembly}.")
    else:
        # Construct Pilon command and run
        pilon_cmd = [
            "pilon", f"-Xmx{ram_gb}g",
            "--genome", second_racon_assembly,
            "--frags", pilon_bam,
            "--output", pilon_out_prefix,
            "--outdir", pilon_out_dir,
            "--changes", "--vcf", "--tracks",
            "--chunksize", str(5000000),
            "--fix", ",".join(["indels", "local", "snps"])
        ]
        print(f"DEBUG - Running Pilon: {' '.join(pilon_cmd)}")
        result = run_subprocess_cmd(pilon_cmd, shell_check=False)
        if result != 0 or not os.path.exists(pilon_assembly):
            print(f"WARN:\tPilon failed: pilon returned {result}. Returning Racon assembly.")
            return second_racon_assembly

    if os.path.exists(pilon_assembly) and os.path.getsize(pilon_assembly) < 1000:
        print(f"WARN:\tPilon assembly is suspiciously small: {pilon_assembly}. Returning Racon assembly.")
        return second_racon_assembly

    print(f"DEBUG - pilon_assembly - {pilon_assembly}")
    return pilon_assembly


# --------------------------------------------------------------
# Perform full assembly polishing
# --------------------------------------------------------------
def polish_assembly(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """Polish an assembly by running Racon (two rounds with long reads)
    and then Pilon (with Illumina reads).

    Args:
        sample_id (str): Identifier for the sample.
        input_csv (str): Path to the CSV file with sample metadata.
        output_dir (str): Base directory for all outputs.
        cpu_threads (int): Number of CPU threads to use.
        ram_gb (int): Amount of RAM (in GB) to allocate for Pilon.

    Returns:
        str: Path to the final polished assembly or None if polishing fails.
    """
    # Read the CSV file and filter to the row corresponding to the sample of interest
    input_df = pd.read_csv(input_csv)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]

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
        ont_raw_reads = os.path.join(species_dir, "ONT", f"{ont_sra}.fastq")
    if pd.notna(illumina_sra) and pd.isna(illumina_f_raw_reads) and pd.isna(illumina_r_raw_reads):
        illumina_f_raw_reads = os.path.join(species_dir, "Illumina", f"{illumina_sra}_1.fastq")
        illumina_r_raw_reads = os.path.join(species_dir, "Illumina", f"{illumina_sra}_2.fastq")
    if pd.notna(pacbio_sra) and pd.isna(pacbio_raw_reads):
        pacbio_raw_reads = os.path.join(species_dir, "PacBio", f"{pacbio_sra}.fastq")
    if pd.notna(ref_seq_gca) and pd.isna(ref_seq):
        ref_seq = os.path.join(species_dir, "RefSeq", f"{species_id}_{ref_seq_gca}_RefSeq.fasta")

    print(f"DEBUG - illumina_sra        - {illumina_sra}")
    print(f"DEBUG - illumina_f_raw_reads- {illumina_f_raw_reads}")
    print(f"DEBUG - illumina_r_raw_reads- {illumina_r_raw_reads}")
    print(f"DEBUG - ont_sra             - {ont_sra}")
    print(f"DEBUG - ont_raw_reads       - {ont_raw_reads}")
    print(f"DEBUG - pacbio_sra          - {pacbio_sra}")
    print(f"DEBUG - pacbio_raw_reads    - {pacbio_raw_reads}")
    print(f"DEBUG - ref_seq_gca         - {ref_seq_gca}")
    print(f"DEBUG - ref_seq             - {ref_seq}")
    print(f"DEBUG - species_id          - {species_id}")
    print(f"DEBUG - est_size            - {est_size}")

    # Check if only reference is provided and skip
    if (pd.isna(illumina_f_raw_reads) or pd.isna(illumina_sra)) and \
       (pd.isna(illumina_r_raw_reads) or pd.isna(illumina_sra)) and \
       (pd.isna(ont_raw_reads) or pd.isna(ont_sra)) and \
       (pd.isna(pacbio_raw_reads) or pd.isna(pacbio_sra)):
        print("SKIP:\tNo valid reads provided, required for assembly polishing.")
        return None

    # Create output directory for polishing
    species_dir = os.path.join(output_dir, species_id)
    sample_dir = os.path.join(species_dir, sample_id)
    polish_out_dir = os.path.join(sample_dir, "polished_assembly")
    os.makedirs(polish_out_dir, exist_ok=True)

    best_assembly = os.path.join(sample_dir, f"{species_id}_best_assembly.fasta")
    print(f"DEBUG - best_assembly - {best_assembly}")
    if not os.path.exists(best_assembly):
        print(f"ERROR:\tBest assembly not found: {best_assembly}")
        return None

    # Set Illumina deduplicated read paths only if Illumina reads are present
    illu_dedup_f_reads = None
    illu_dedup_r_reads = None
    if pd.notna(illumina_sra) or (pd.notna(illumina_f_raw_reads) and pd.notna(illumina_r_raw_reads)):
        illu_dedup_f_reads = os.path.join(species_dir, "Illumina", f"{species_id}_illu_forward_dedup.fastq")
        illu_dedup_r_reads = os.path.join(species_dir, "Illumina", f"{species_id}_illu_reverse_dedup.fastq")

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
        current_work_dir = polish_out_dir
    else:
        current_work_dir = starting_work_dir
    os.chdir(current_work_dir)

    # -------------------------------------------------------------------------
    # Step 1: Two rounds of Racon polishing if ONT or PacBio reads exist
    # -------------------------------------------------------------------------
    racon_work_dir = os.path.join(current_work_dir, "racon_polish")
    os.makedirs(racon_work_dir, exist_ok=True)
    racon_final = os.path.join(polish_out_dir, f"{sample_id}_racon.fasta")

    if pd.notna(ont_raw_reads) or pd.notna(pacbio_raw_reads):
        print("2x Racon Polishing Long (ONT or PacBio) Reads...")
        initial_racon_assembly = racon_polish_assembly(best_assembly,
                                                      highest_mean_qual_long_reads,
                                                      racon_work_dir,
                                                      sample_id,
                                                      cpu_threads,
                                                      1)
        racon_assembly = racon_polish_assembly(initial_racon_assembly,
                                               highest_mean_qual_long_reads,
                                               racon_work_dir,
                                               sample_id,
                                               cpu_threads,
                                               2)
    else:
        print("SKIP:\t2x-Racon Polish; No Long Reads Provided.")
        racon_assembly = best_assembly

    # Safe-copy Racon
    if not os.path.exists(racon_assembly):
        print(f"ERROR:\tRacon output not found: {racon_assembly}. Falling back to best assembly.")
        racon_assembly = best_assembly
    if os.path.exists(racon_final):
        print(f"SKIP:\tRacon final already exists: {racon_final}")
    else:
        print(f"DEBUG - Copying Racon assembly: {racon_assembly} to {racon_final}")
        shutil.copy(racon_assembly, racon_final)
        if not os.path.exists(racon_final):
            print(f"ERROR:\tFailed to copy Racon assembly to {racon_final}")
            return None

    # -------------------------------------------------------------------------
    # Step 2: Pilon polishing if Illumina paired-end reads exist
    # -------------------------------------------------------------------------
    pilon_work_dir = os.path.join(polish_out_dir, "pilon_polish")
    os.makedirs(pilon_work_dir, exist_ok=True)
    pilon_final = os.path.join(polish_out_dir, f"{sample_id}_pilon.fasta")

    if pd.notna(illumina_f_raw_reads) and pd.notna(illumina_r_raw_reads):
        print("Pilon Polishing with Illumina Reads...")
        pilon_bam = pilon_prep(racon_final,
                               illu_dedup_f_reads, illu_dedup_r_reads,
                               pilon_work_dir, cpu_threads)
        if pilon_bam is None:
            print(f"WARN:\tPilon prep failed. Falling back to Racon assembly.")
            polished_assembly = racon_final
        else:
            polished_assembly = pilon_polish(best_assembly, racon_final,
                                             pilon_bam, pilon_work_dir, sample_id,
                                             cpu_threads, ram_gb)
    else:
        print("SKIP:\tPilon Polish; No Illumina Reads Provided.")
        polished_assembly = racon_final

    print(f"DEBUG - polished_assembly - {polished_assembly}")

    if not os.path.exists(polished_assembly):
        print(f"WARN:\tPolished assembly not found: {polished_assembly}. Falling back to Racon assembly.")
        polished_assembly = racon_final
        if not os.path.exists(polished_assembly):
            print(f"ERROR:\tRacon assembly not found: {polished_assembly}")
            return None

    if os.path.getsize(polished_assembly) < 1000:
        print(f"WARN:\tPolished assembly is suspiciously small: {polished_assembly}. Falling back to Racon assembly.")
        polished_assembly = racon_final

    # Copy to pilon_final
    if os.path.exists(pilon_final):
        print(f"SKIP:\tPilon final already exists: {pilon_final}")
    else:
        print(f"DEBUG - Copying polished assembly: {polished_assembly} to {pilon_final}")
        shutil.copy(polished_assembly, pilon_final)
        if not os.path.exists(pilon_final):
            print(f"ERROR:\tFailed to copy polished assembly to {pilon_final}")
            return None

    print(f"DEBUG - pilon_final - {pilon_final}")

    # Prepare final polished assembly path
    final_polished_assembly = os.path.join(
        sample_dir,
        f"{sample_id}_final_polish_assembly.fasta")

    print(f"DEBUG - final_polished_assembly - {final_polished_assembly}")

    # Copy to final_polished_assembly
    if os.path.exists(final_polished_assembly):
        print(f"SKIP:\tFinal polished assembly already exists: {final_polished_assembly}")
    else:
        if os.path.exists(pilon_final):
            print(f"DEBUG - Copying pilon final: {pilon_final} to {final_polished_assembly}")
            shutil.copy(pilon_final, final_polished_assembly)
            if not os.path.exists(final_polished_assembly):
                print(f"ERROR:\tFailed to copy pilon final to {final_polished_assembly}")
                return None
        else:
            print(f"WARN:\tPilon final not found: {pilon_final}. Falling back to Racon assembly.")
            if os.path.exists(racon_final):
                print(f"DEBUG - Copying Racon final: {racon_final} to {final_polished_assembly}")
                shutil.copy(racon_final, final_polished_assembly)
                if not os.path.exists(final_polished_assembly):
                    print(f"ERROR:\tFailed to copy Racon final to {final_polished_assembly}")
                    return None
            else:
                print(f"ERROR:\tRacon final not found: {racon_final}")
                return None

    print(f"PASS:\tPolishing complete for {best_assembly} -> {final_polished_assembly}")
    return final_polished_assembly


if __name__ == "__main__":
    # Provide a usage message if the correct number of arguments is not supplied
    if len(sys.argv) != 6:
        print("Usage: python3 polish_assembly.py <sample_id> <input_csv> "
              "<output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)

    # Run the polish_assembly function with CLI arguments
    final_polished_assembly = polish_assembly(sys.argv[1],       # sample_id
                                             sys.argv[2],       # input_csv
                                             sys.argv[3],       # output_dir
                                             str(sys.argv[4]),  # cpu_threads
                                             str(sys.argv[5]))  # ram_gb