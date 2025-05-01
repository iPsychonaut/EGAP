#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
preprocess_ont.py

Updated on Sat Mar 29 2025

This script preprocesses ONT reads with NanoPlot, Filtlong, and Ratatosk.
Handles raw directory concatenation and SRA downloads.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import os, shutil, sys, subprocess, glob
import pandas as pd
from utilities import run_subprocess_cmd, pigz_compress, get_current_row_data, select_long_reads
from qc_assessment import nanoplot_qc_reads


# --------------------------------------------------------------
# Preprocess ONT sequencing reads
# --------------------------------------------------------------
def preprocess_ont(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """Preprocess ONT reads with NanoPlot, Filtlong, and Ratatosk.

    Handles concatenation of raw directory files, SRA downloads, quality control with
    NanoPlot, filtering with Filtlong, and correction with Ratatosk (if Illumina reads
    are available), then selects the highest quality reads.

    Args:
        sample_id (str): Sample identifier.
        input_csv (str): Path to metadata CSV file.
        output_dir (str): Directory for output files.
        cpu_threads (int or str): Number of CPU threads to use.
        ram_gb (int or str): Available RAM in GB.

    Returns:
        str or None: Path to the highest quality compressed ONT reads file, or None if no reads are available.
    """
    print(f"Preprocessing ONT reads for {sample_id}...")
    
    # Parse the CSV and retrieve relevant row data
    input_df = pd.read_csv(input_csv)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]  # Convert to Series (single row)

    # Identify read paths, reference, and BUSCO lineage info from CSV
    ont_raw_reads = current_series["ONT_RAW_READS"]
    ont_raw_dir = current_series["ONT_RAW_DIR"]
    ont_sra = current_series["ONT_SRA"]
    species_id = current_series["SPECIES_ID"]
    est_size = current_series["EST_SIZE"]
   
    print(f"DEBUG - ont_sra - {ont_sra}")
    print(f"DEBUG - ont_raw_reads - {ont_raw_reads}")
    print(f"DEBUG - ont_raw_dir - {ont_raw_dir}")
    if pd.isna(ont_sra) and pd.isna(ont_raw_dir) and pd.isna(ont_raw_reads):
        print("SKIP:\tONT preprocessing; no reads provided")
        return None

    species_dir = os.path.join(output_dir, species_id)
    os.makedirs(species_dir, exist_ok=True)
    sample_dir = os.path.join(species_dir, sample_id)
    os.makedirs(sample_dir, exist_ok=True)
    ont_dir = os.path.join(species_dir, "ONT")
    os.makedirs(ont_dir, exist_ok=True)
    os.chdir(ont_dir)
    ont_raw_reads = os.path.join(ont_dir, f"{ont_sra}.fastq")
    ont_raw_reads_gz = ont_raw_reads + ".gz"

    illu_dedup_f_reads = os.path.join(species_dir, "Illumina", f"{species_id}_illu_forward_dedup.fastq.gz")
    illu_dedup_r_reads = os.path.join(species_dir, "Illumina", f"{species_id}_illu_reverse_dedup.fastq.gz")

    print(f"DEBUG - illu_dedup_f_reads - {illu_dedup_f_reads}")
    print(f"DEBUG - illu_dedup_r_reads - {illu_dedup_r_reads}")

    # Handle Raw Data Directory    
    if isinstance(ont_raw_dir, str):
        print(f"NOTE:\tConcatenating ONT files from {ont_raw_dir}")
        ont_files = sorted(glob.glob(f"{ont_raw_dir}/*.fastq.gz"))
        if not ont_files:
            print(f"ERROR:\tNo ONT files found in {ont_raw_dir}")
            sys.exit(1)
        ont_raw_reads = f"{species_id}_ont_combined.fastq.gz"
        run_subprocess_cmd(["cat"] + ont_files + [">", ont_raw_reads], True)
        ont_raw_reads_gz = pigz_compress(ont_raw_reads, cpu_threads)

    # Handle SRA download
    if isinstance(ont_sra, str):
        print(f"Downloading SRA {ont_sra} from GenBank...")
        if os.path.exists(ont_raw_reads) :
            print(f"SKIP:\tSRA download already exists: {ont_raw_reads}.")
            ont_raw_reads_gz = pigz_compress(ont_raw_reads, cpu_threads)
        elif os.path.exists(ont_raw_reads_gz):
            print(f"SKIP:\tSRA download already exists: {ont_raw_reads_gz}.")
        else:
            run_subprocess_cmd(["fasterq-dump", "--threads", str(cpu_threads), ont_sra], False)
            ont_raw_reads_gz = pigz_compress(ont_raw_reads, cpu_threads)
    
    if pd.isna(ont_raw_reads) and not os.path.exists(ont_raw_reads_gz):
        print(f"SKIP:\tInput files not found: {ont_raw_reads_gz}.")
        return None
    
    if os.path.exists(ont_raw_reads):
        ont_raw_reads_gz = pigz_compress(ont_raw_reads, cpu_threads)

    # Parse estimated genome size
    if isinstance(est_size, float):
        estimated_size = est_size
    else:
        estimated_size = est_size.lower().strip() if est_size != "None" else "25m"
    multipliers = {'m': 10**6, 'g': 10**9}
    if estimated_size in multipliers:
        est_size_bp = int(float(estimated_size) * multipliers[estimated_size])
    else:
        print(f"NOTE:\tUnable to parse input estimated size {est_size}, using default: 25000000")
        est_size_bp = 25000000

    # NanoPlot Raw Reads
    sample_stats_dict = nanoplot_qc_reads(ont_raw_reads_gz, "Raw_ONT_", cpu_threads, sample_stats_dict)

    # Filtlong
    filtered_ont = os.path.join(ont_dir, os.path.basename(ont_raw_reads_gz.replace(f"{ont_sra}", f"{species_id}_ont_filtered")))
    coverage = 75
    target_bases = est_size_bp * coverage

    if not os.path.exists(filtered_ont):
        illumina_opt = f"-1 {illu_dedup_f_reads} -2 {illu_dedup_r_reads}" if illu_dedup_f_reads != "None" and illu_dedup_r_reads != "None" else ""
        filtlong_cmd = f"filtlong {illumina_opt} --trim --min_length 1000 --min_mean_q 8 --keep_percent 90 --target_bases {target_bases} {ont_raw_reads_gz} | pigz -p {cpu_threads} > {filtered_ont}"
        _ = run_subprocess_cmd(filtlong_cmd, True)
    else:
        print(f"SKIP\tFiltlong filtered reads exist: {filtered_ont}.")

    # NanoPlot Filtered Reads
    sample_stats_dict = nanoplot_qc_reads(filtered_ont, "Filt_ONT_", cpu_threads, sample_stats_dict)
    os.chdir(ont_dir)

    # Ratatosk (only if Illumina reads available)
    corrected_out = os.path.join(ont_dir, "ratatosk_corrected.fastq")
    corrected_ont = os.path.join(ont_dir, os.path.basename(ont_raw_reads_gz).replace(f"{ont_sra}", f"{species_id}_ont_corrected"))
    if illu_dedup_f_reads == "None" or illu_dedup_r_reads == "None":
        print("SKIP:\tRatatosk correction; no Illumina reads provided")
        if not os.path.exists(corrected_ont):
            subprocess.run(["ln", "-sf", filtered_ont, corrected_ont], check=True)
    else:
        if os.path.exists(corrected_out + ".fastq.gz"):
            shutil.move(corrected_out + ".fastq.gz", corrected_ont)
        if not os.path.exists(corrected_ont):
            ratatosk_cmd = ["Ratatosk", "correct", "-s", illu_dedup_f_reads, "-s", illu_dedup_r_reads,
                            "-l", filtered_ont, "-o", corrected_out, "-c", str(cpu_threads), "-G", "-v"]
            _ = run_subprocess_cmd(ratatosk_cmd, False)
            shutil.move(corrected_out + ".fastq.gz", corrected_ont)
        else:
            print("SKIP:\tRatatosk correct reads exist: {corrected_ont}.")

    # NanoPlot Corrected Reads
    sample_stats_dict = nanoplot_qc_reads(corrected_ont, "Corr_ONT_", cpu_threads, sample_stats_dict)    
    os.chdir(ont_dir)
    
    highest_mean_qual_long_reads_gz = select_long_reads(output_dir, input_csv, sample_id, cpu_threads)
    
    print(f"PASS:\tPreprocessed Raw ONT Reads for {sample_id}: {highest_mean_qual_long_reads_gz}.")
    
    return highest_mean_qual_long_reads_gz


if __name__ == "__main__":
    # Log raw sys.argv immediately
    print(f"DEBUG: Raw sys.argv = {sys.argv}")
    print(f"DEBUG: Length of sys.argv = {len(sys.argv)}")
    
    # Check argument count
    if len(sys.argv) != 6:
        print(f"ERROR: Expected 5 arguments (plus script name), got {len(sys.argv)-1}: {sys.argv[1:]}", 
              file=sys.stderr)
        print("Usage: python3 preprocess_ont.py <sample_id> <input_csv> <output_dir> <cpu_threads> <ram_gb>", 
              file=sys.stderr)
        sys.exit(1)
    
    # Log each argument
    for i, arg in enumerate(sys.argv):
        print(f"DEBUG: sys.argv[{i}] = '{arg}'")
    
    sample_id = sys.argv[1]
    input_csv = sys.argv[2]
    output_dir = sys.argv[3]
    cpu_threads = sys.argv[4]
    ram_gb = int(sys.argv[5]) if sys.argv[5] != " " else 8
    
    print(f"DEBUG: Parsed sample_id = '{sample_id}'")
    print(f"DEBUG: Parsed input_csv = '{input_csv}'")
    print(f"DEBUG: Parsed output_dir = '{output_dir}'")
    print(f"DEBUG: Parsed cpu_threads = '{sys.argv[4]}' (converted to {cpu_threads})")
    print(f"DEBUG: Parsed ram_gb = '{sys.argv[5]}' (converted to {ram_gb})")
    
    highest_mean_qual_long_reads_gz = preprocess_ont(sample_id, input_csv, output_dir, cpu_threads, ram_gb)