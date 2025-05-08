#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polish_assembly.py

This module processes genomic assembly data by applying polishing steps using Racon
for long reads (ONT or PacBio) and Pilon for Illumina reads. It handles input validation,
subprocess execution, and file management for assembly refinement.

Updated on Sat Apr 11 2025

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""

import os, sys, shutil
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
        str: Path to the Racon-polished assembly FASTA.
    """
    # Align reads to the assembly using minimap2, generating a PAF file.
    racon_paf = os.path.join(racon_out_dir, f"racon_round{iteration_count}.paf")
    if os.path.exists(racon_paf):
        print(f"SKIP:\tRacon PAF {iteration_count} already exists: {racon_paf}.")
    else:
        minimap2_cmd = (f"minimap2 -t {cpu_threads} -x map-ont {input_assembly} {long_reads} "
                        f"> {racon_paf}")
        _ = run_subprocess_cmd(minimap2_cmd, shell_check=True)

    # Perform polishing with Racon using the PAF alignment.
    racon_assembly = os.path.join(racon_out_dir, f"{sample_id}_racon_polish_{iteration_count}.fasta")
    if os.path.exists(racon_assembly):
        print(f"SKIP:\tRacon Assembly {iteration_count} already exists: {racon_assembly}.")
    else:
        racon_cmd = (f"racon -t {cpu_threads} {long_reads} {racon_paf} {input_assembly} "
                     f"> {racon_assembly}")
        _ = run_subprocess_cmd(racon_cmd, shell_check=True)

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
        illu_f_dedup (str): Path to gzipped deduplicated forward Illumina FASTQ.
        illu_r_dedup (str): Path to gzipped deduplicated reverse Illumina FASTQ.
        assembly_out_dir (str): Directory for output files.
        cpu_threads (str): Number of CPU threads to use.

    Returns:
        str: Path to the sorted, indexed BAM file.
    """
    # Construct output file names for SAM, unsorted BAM, and final sorted BAM.
    pilon_bam = os.path.join(assembly_out_dir, os.path.basename(input_assembly).replace(".fasta", ".bam"))
    output_sam = pilon_bam.replace(".bam", ".sam")
    sorted_bam = pilon_bam.replace(".bam", "_sorted.bam")

    # Build a BWA index if necessary, then map Illumina reads to the assembly (SAM).
    if not os.path.exists(output_sam):
        bwa_index_cmd = ["bwa-mem2", "index", input_assembly]
        _ = run_subprocess_cmd(bwa_index_cmd, shell_check=False)

        bwa_cmd = (
            f"bwa-mem2 mem -t {cpu_threads} {input_assembly} {illu_f_dedup} {illu_r_dedup} "
            f"> {output_sam}"
        )
        _ = run_subprocess_cmd(bwa_cmd, shell_check=True)

    if os.path.exists(output_sam) and os.path.getsize(output_sam) < 1000:
        raise ValueError(f"Generated SAM file is suspiciously small: {output_sam}")

    # Convert SAM to an unsorted BAM file.
    if not os.path.exists(sorted_bam):
        samview_cmd = f"samtools view -@ {cpu_threads} -S -b {output_sam} > {sorted_bam}"
        _ = run_subprocess_cmd(samview_cmd, shell_check=True)

    # Sort and index the BAM file if final BAM doesn't exist.
    if not os.path.exists(pilon_bam):
        bamsort_cmd = ["bamtools", "sort", "-in", sorted_bam, "-out", pilon_bam]
        _ = run_subprocess_cmd(bamsort_cmd, shell_check=False)

        samtools_index_cmd = ["samtools", "index", pilon_bam]
        _ = run_subprocess_cmd(samtools_index_cmd, shell_check=False)

    if os.path.exists(sorted_bam) and os.path.getsize(sorted_bam) < 1000:
        raise ValueError(f"Sorted BAM file is suspiciously small: {sorted_bam}")

    return pilon_bam


# --------------------------------------------------------------
# Polish assembly with Pilon using Illumina reads
# --------------------------------------------------------------
def pilon_polish(best_assembly, second_racon_assembly, pilon_bam, assembly_out_dir, sample_id, cpu_threads, ram_gb):
    """Polish an assembly with Pilon using aligned Illumina reads.

    Runs Pilon to refine the assembly based on Illumina read alignments, producing
    a polished FASTA file.

    Args:
        best_assembly (str): Path to the initial gzipped assembly FASTA.
        second_racon_assembly (str): Path to the Racon-polished assembly.
        pilon_bam (str): Path to the sorted, indexed BAM file.
        assembly_out_dir (str): Directory for output files.
        sample_id (str): Sample identifier.
        cpu_threads (str): Number of CPU threads to use.
        ram_gb (str): RAM in GB for Pilon's Java process.

    Returns:
        str: Path to the gzipped Pilon-polished assembly FASTA.
    """
    # Prepare output filenames for Pilon run.
    pilon_out_prefix = f"{sample_id}_pilon_assembly"
    pilon_out_dir = assembly_out_dir
    os.makedirs(pilon_out_dir, exist_ok=True)

    pilon_assembly = os.path.join(pilon_out_dir, f"{sample_id}_pilon_assembly.fasta")
    pilon_ext_code = 0
    if os.path.exists(pilon_assembly):
        print(f"SKIP:\tPilon Polished Assembly already exists: {pilon_assembly}.")
    else:
        # Check if a partial Pilon assembly already exists.
        if os.path.exists(pilon_assembly):
            print(f"SKIP:\tPilon Polished Assembly already exists: {pilon_assembly}.")
        else:
            # Construct Pilon command and run.
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
            pilon_ext_code = run_subprocess_cmd(pilon_cmd, shell_check=False)

            # If Pilon fails (often memory issues), skip moving the Pilon output.
            if pilon_ext_code != 0:
                print("WARN:\tPilon was not able to finish; attempting to continue...")

    print(f"DEBUG - pilon_assembly - {pilon_assembly}")

    return pilon_assembly


# --------------------------------------------------------------
# Perform full assembly polishing
# --------------------------------------------------------------
def polish_assembly(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """
    Polish an assembly by running Racon (two rounds with long reads)
    and then Pilon (with Illumina reads), preserving your original
    variable names, flow, and debug statements.

    Args:
        sample_id (str): Identifier for the sample.
        input_csv (str): Path to the CSV file with sample metadata.
        output_dir (str): Base directory for all outputs.
        cpu_threads (int): Number of CPU threads to use.
        ram_gb (int): Amount of RAM (in GB) to allocate for Pilon.

    Returns:
        str: Path to the final gzipped polished assembly.
    """
    # Read the CSV file and filter to the row corresponding to the sample of interest.
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

    # Create output directory for polishing.
    species_dir = os.path.join(output_dir, species_id)
    sample_dir = os.path.join(species_dir, sample_id)
    polish_out_dir = os.path.join(sample_dir, "polished_assembly")
    os.makedirs(polish_out_dir, exist_ok=True)
    cwd = os.getcwd()

    best_assembly = os.path.join(sample_dir, f"{species_id}_best_assembly.fasta")
    print(f"DEBUG - best_assembly - {best_assembly}")

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

    print(f"DEBUG - highest_mean_qual_long_reads    - {highest_mean_qual_long_reads}")

    # Ensure work directory output
    starting_work_dir = os.getcwd()
    if "work" not in starting_work_dir:
        current_work_dir = polish_out_dir
    else:
        current_work_dir = starting_work_dir
    os.chdir(current_work_dir)

    # -------------------------------------------------------------------------
    # Step 1: Two rounds of Racon polishing if ONT or PacBio reads exist.
    # -------------------------------------------------------------------------
    racon_work_dir = os.path.join(cwd, "racon_polish")
    os.makedirs(racon_work_dir, exist_ok=True)
    racon_final = os.path.join(current_work_dir, f"{sample_id}_racon.fasta")

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
        raise FileNotFoundError(f"Racon output not found: {racon_assembly}")
    if os.path.exists(racon_final):
        print(f"SKIP:\tRacon final already exists: {racon_final}")
    else:
        shutil.copy(racon_assembly, racon_final)

    # -------------------------------------------------------------------------
    # Step 2: Pilon polishing if Illumina paired-end reads exist.
    # -------------------------------------------------------------------------
    pilon_work_dir = os.path.join(cwd, "pilon_polish")
    os.makedirs(pilon_work_dir, exist_ok=True)
    pilon_final = os.path.join(current_work_dir, f"{sample_id}_pilon.fasta")

    if pd.notna(illumina_f_raw_reads) and pd.notna(illumina_r_raw_reads):
        print("Pilon Polishing with Illumina Reads...")
        pilon_bam = pilon_prep(racon_final,
                               illu_dedup_f_reads, illu_dedup_r_reads,
                               pilon_work_dir, cpu_threads)
        polished_assembly = pilon_polish(best_assembly, racon_final,
                                            pilon_bam, pilon_work_dir, sample_id,
                                            cpu_threads, ram_gb)
    else:
        print("SKIP:\tPilon Polish; No Illumina Reads Provided.")
        polished_assembly = racon_final

    print(f"DEBUG - polished_assembly - {polished_assembly}")

    if not os.path.exists(polished_assembly):
        raise FileNotFoundError(f"ERROR:\t{polished_assembly} does not exist.")

    if os.path.exists(pilon_final):
        print(f"SKIP:\tPilon final already exists: {pilon_final}")
    else:
        shutil.copy(polished_assembly, pilon_final)

    print(f"DEBUG - pilon_final - {pilon_final}")

    final_polished_assembly = os.path.join(
        sample_dir,
        os.path.basename(polished_assembly)
            .replace("_best_assembly", "_final_polish_assembly")
            .replace("_pilon", "_final_polish")
            .replace("_racon", "_final_polish_assembly"))

    print(f"DEBUG - final_polished_assembly - {final_polished_assembly}")

    if os.path.exists(pilon_final) and not os.path.exists(final_polished_assembly):
        shutil.move(pilon_final, final_polished_assembly)
    else:
        raise FileNotFoundError(f"ERROR:\t{polished_assembly} cannot be moved to {final_polished_assembly}")

    print(f"PASS:\tPolishing complete for {best_assembly} -> {final_polished_assembly}")

    return final_polished_assembly


if __name__ == "__main__":
    # Provide a usage message if the correct number of arguments is not supplied.
    if len(sys.argv) != 6:
        print("Usage: python3 polish_assembly.py <sample_id> <input_csv> "
            "<output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)

    # Run the polish_assembly function with CLI arguments.
    final_polished_assembly = polish_assembly(sys.argv[1],       # sample_id
                                              sys.argv[2],       # input_csv
                                              sys.argv[3],       # output_dir
                                              str(sys.argv[4]),  # cpu_threads
                                              str(sys.argv[5]))  # ram_gb
