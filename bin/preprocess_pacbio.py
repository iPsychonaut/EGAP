#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
preprocess_pacbio.py

Updated on Sat Mar 29 2025

This script preprocesses PacBio reads with NanoPlot and Filtlong.
Handles raw directory concatenation and SRA downloads.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import sys, os, glob, shutil
import pandas as pd
from utilities import run_subprocess_cmd, get_current_row_data, select_long_reads
from qc_assessment import nanoplot_qc_reads


# --------------------------------------------------------------
# Preprocess PacBio sequencing reads
# --------------------------------------------------------------
def preprocess_pacbio(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """Preprocess PacBio reads with NanoPlot and Filtlong.

    Handles concatenation of raw directory files, SRA downloads, quality control with
    NanoPlot, and filtering with Filtlong, then selects the highest quality reads.

    Args:
        sample_id (str): Sample identifier.
        input_csv (str): Path to metadata CSV file.
        output_dir (str): Directory for output files.
        cpu_threads (int or str): Number of CPU threads to use.
        ram_gb (int or str): Available RAM in GB.

    Returns:
        str or None: Path to the highest quality PacBio reads file, or None if no reads are available.
    """
    print(f"Preprocessing PacBio reads for {sample_id}...")
    # Parse the CSV and retrieve relevant row data
    input_df = pd.read_csv(input_csv)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]  # Convert to Series (single row)

    # Identify read paths, reference, and BUSCO lineage info from CSV
    pacbio_sra = current_series["PACBIO_SRA"]
    pacbio_raw_reads = current_series["PACBIO_RAW_READS"]
    pacbio_raw_dir = current_series["PACBIO_RAW_DIR"]
    species_id = current_series["SPECIES_ID"]
    est_size = current_series["EST_SIZE"]

    print(f"DEBUG - pacbio_sra - {pacbio_sra}")
    print(f"DEBUG - pacbio_raw_reads - {pacbio_raw_reads}")
    print(f"DEBUG - pacbio_raw_dir - {pacbio_raw_dir}")
    if pd.isna(pacbio_sra) and pd.isna(pacbio_raw_reads) and pd.isna(pacbio_raw_dir):
        print("SKIP:\tPacBio preprocessing; no reads provided")
        return None

    species_dir = os.path.join(output_dir, species_id)
    os.makedirs(species_dir, exist_ok=True)
    pacbio_dir = os.path.join(species_dir, "PacBio")
    os.makedirs(pacbio_dir, exist_ok=True)
    os.chdir(pacbio_dir)
    pacbio_raw_reads = os.path.join(pacbio_dir, f"{pacbio_sra}.fastq")

    # Handle Raw Data Directory    
    if pacbio_raw_dir != "None" and not pacbio_raw_reads != "None":
        print(f"NOTE:\tConcatenating PacBio files from {pacbio_raw_dir}")
        pacbio_files = sorted(glob.glob(f"{pacbio_raw_dir}/*.fastq"))
        if not pacbio_files:
            print(f"ERROR:\tNo PacBio files found in {pacbio_raw_dir}")
            sys.exit(1)
        pacbio_raw_reads = f"{species_id}_pacbio_combined.fastq"
        run_subprocess_cmd(["cat"] + pacbio_files + [">", pacbio_raw_reads], True)

    if not pacbio_raw_reads or pacbio_raw_reads == "None":
        print("SKIP:\tPacBio preprocessing; no valid reads")
        return None

    # Handle SRA download
    if isinstance(pacbio_sra, str):
        print(f"Downloading SRA {pacbio_sra} from GenBank...")
        if os.path.exists(pacbio_raw_reads) :
            print(f"SKIP:\tSRA download already exists: {pacbio_raw_reads}.")
        else:
            run_subprocess_cmd(["fasterq-dump", "--threads", str(cpu_threads), pacbio_sra], False)
    
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
    sample_stats_dict = nanoplot_qc_reads(pacbio_raw_reads, "Raw_PacBio_", cpu_threads, sample_stats_dict)    

    # Filtlong
    filtered_pacbio = os.path.join(pacbio_dir, os.path.basename(pacbio_raw_reads.replace(f"{pacbio_sra}", f"{species_id}_filtered")))
    coverage = 75
    target_bases = est_size_bp * coverage
    filt_min_trim_length = 1000
    filt_min_mean_q = 8
    filt_keep_percent = 90

    if not os.path.exists(filtered_pacbio):
        filtlong_cmd = f"filtlong --min_length {filt_min_trim_length} --min_mean_q {filt_min_mean_q} --keep_percent {filt_keep_percent} --target_bases {target_bases} {pacbio_raw_reads} > {filtered_pacbio}"
        _ = run_subprocess_cmd(filtlong_cmd, True)
    else:
        print("SKIP\tFiltlong filtered reads exist: {filtered_pacbio}.")

    # NanoPlot Filtered Reads
    sample_stats_dict = nanoplot_qc_reads(filtered_pacbio, "Filt_PacBio_", cpu_threads, sample_stats_dict)    
    os.chdir(pacbio_dir)

    highest_mean_qual_long_reads = select_long_reads(output_dir, input_csv, sample_id, cpu_threads)

    print(f"PASS:\tPreprocessed Raw PacBio Reads for {sample_id}: {highest_mean_qual_long_reads}.")
    
    # Force-renaming in case select_long_reads() didn't
    canonical_path = os.path.join(pacbio_dir, f"{species_id}_PacBio_highest_mean_qual_long_reads.fastq")
    if not os.path.exists(canonical_path):
        print(f"FALLBACK:\tCopying highest-mean-qual long reads to canonical path: {canonical_path}")
        shutil.copy(highest_mean_qual_long_reads, canonical_path)
    
    return canonical_path


if __name__ == "__main__":
    # Log raw sys.argv immediately
    print(f"DEBUG: Raw sys.argv = {sys.argv}")
    print(f"DEBUG: Length of sys.argv = {len(sys.argv)}")
    
    # Check argument count
    if len(sys.argv) != 6:
        print(f"ERROR: Expected 5 arguments (plus script name), got {len(sys.argv)-1}: {sys.argv[1:]}", 
              file=sys.stderr)
        print("Usage: python3 preprocess_pacbio.py <sample_id> <input_csv> <output_dir> <cpu_threads> <ram_gb>", 
              file=sys.stderr)
        sys.exit(1)
    
    # Log each argument
    for i, arg in enumerate(sys.argv):
        print(f"DEBUG: sys.argv[{i}] = '{arg}'")
    
    sample_id = sys.argv[1]
    input_csv = sys.argv[2]
    output_dir = sys.argv[3]
    cpu_threads = sys.argv[4] if isinstance(sys.argv[4], int) else 1
    ram_gb = sys.argv[5] if isinstance(sys.argv[5], int) else 8
    
    print(f"DEBUG: Parsed sample_id = '{sample_id}'")
    print(f"DEBUG: Parsed input_csv = '{input_csv}'")
    print(f"DEBUG: Parsed output_dir = '{output_dir}'")
    print(f"DEBUG: Parsed cpu_threads = '{sys.argv[4]}' (converted to {cpu_threads})")
    print(f"DEBUG: Parsed ram_gb = '{sys.argv[5]}' (converted to {ram_gb})")
    
    highest_mean_qual_long_reads = preprocess_pacbio(sample_id, input_csv, output_dir, cpu_threads, ram_gb)
