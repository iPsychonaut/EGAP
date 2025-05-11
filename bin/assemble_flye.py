#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
assemble_flye.py

This script performs genome assembly on long reads (ONT or PacBio) using the Flye assembler.
It supports customizable genome size estimation, CPU thread allocation, and haplotype retention.

Updated on Sat Mar 29 2025

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os, sys, shutil
import pandas as pd
from utilities import run_subprocess_cmd, get_current_row_data
from qc_assessment import qc_assessment


# --------------------------------------------------------------
# Run Flye assembly with long reads
# --------------------------------------------------------------
def assemble_flye(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """Assemble genomic data using Flye with ONT or PacBio reads.

    Executes Flye assembly on long reads and performs quality control on the
    resulting assembly.

    Args:
        sample_id (str): Sample identifier.
        input_csv (str): Path to metadata CSV file.
        output_dir (str): Directory for output files.
        cpu_threads (int): Number of CPU threads to use.
        ram_gb (int): Available RAM in GB.

    Returns:
        str or None: Path to the gzipped Flye assembly FASTA, or None if no long reads are provided.
    """
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

    if pd.isna(ont_raw_reads) and pd.isna(pacbio_raw_reads):
        print("SKIP:\tNo reads available for processing")
        return None

    sample_dir = os.path.join(species_dir, sample_id)
    os.makedirs(sample_dir, exist_ok=True)    
    flye_out_dir = os.path.join(sample_dir, "flye_assembly")
    os.makedirs(flye_out_dir, exist_ok=True)
    egap_flye_assembly_path = os.path.join(flye_out_dir, f"{sample_id}_flye.fasta")

    # Ensure work directory output
    starting_work_dir = os.getcwd()
    if "work" not in starting_work_dir:
        current_work_dir = flye_out_dir
    else:
        current_work_dir = starting_work_dir
    os.chdir(current_work_dir)

    flye_path = os.path.join(current_work_dir, "assembly.fasta")
    if os.path.exists(egap_flye_assembly_path):
        egap_flye_assembly_path, flye_stats_list, _ = qc_assessment("flye", input_csv, sample_id, output_dir, cpu_threads, ram_gb)
        print(f"SKIP:\tFinal Flye Assembly already exists: {egap_flye_assembly_path}.")
    else:        
        if pd.notna(ont_raw_reads):
            flye_cmd = ["flye",
                        "--nano-corr", highest_mean_qual_long_reads,
                        "--out-dir", current_work_dir, 
                        "--genome-size", str(est_size), "--threads", str(cpu_threads),
                        "--iterations", "3", "--keep-haplotypes"]
        elif pd.notna(pacbio_raw_reads):
            flye_cmd = ["flye",
                        "--pacbio-corr", highest_mean_qual_long_reads,
                        "--out-dir", current_work_dir, 
                        "--genome-size", str(est_size), "--threads", str(cpu_threads),
                        "--iterations", "3", "--keep-haplotypes"]
        _ = run_subprocess_cmd(flye_cmd, shell_check = False)

    if os.path.exists(flye_path) and not os.path.exists(egap_flye_assembly_path):    
        shutil.move(flye_path, egap_flye_assembly_path)
            
    egap_flye_assembly_path, flye_stats_list, _ = qc_assessment("flye", input_csv, sample_id, output_dir, cpu_threads, ram_gb)

    return egap_flye_assembly_path


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python3 assemble_flye.py <sample_id> <input_csv> "
            "<output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)
        
    egap_flye_assembly_path = assemble_flye(sys.argv[1],       # sample_id
                                               sys.argv[2],       # input_csv
                                               sys.argv[3],       # output_dir
                                               str(sys.argv[4]),  # cpu_threads
                                               str(sys.argv[5]))  # ram_gb
