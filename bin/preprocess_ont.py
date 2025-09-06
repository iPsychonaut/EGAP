#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
preprocess_ont.py

This script preprocesses ONT reads with NanoPlot, Filtlong, and Ratatosk.
Handles raw directory concatenation and SRA downloads.

Created on Wed Aug 16 2023

Updated on Wed Sept 3 2025

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import os, shutil, sys, subprocess, glob, re
import pandas as pd
from utilities import run_subprocess_cmd, get_current_row_data, select_long_reads
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
        str or None: Path to the highest quality ONT reads file, or None if no reads are available.
    """
    print(f"Preprocessing ONT reads for {sample_id}...")

    # Always anchor everything to absolute output_dir to avoid nested paths after chdir
    output_dir_abs = os.path.abspath(output_dir)
    input_csv_abs = os.path.abspath(input_csv)
    # Keep your original species_id source of truth:
    input_df = pd.read_csv(input_csv_abs)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]  # Convert to Series (single row)

    # Identify read paths, reference, and BUSCO lineage info from CSV
    ont_raw_reads = current_series["ONT_RAW_READS"]
    ont_raw_dir = current_series["ONT_RAW_DIR"]
    ont_sra = current_series["ONT_SRA"]
    species_id = current_series["SPECIES_ID"]
    est_size = current_series["EST_SIZE"]

    # Recompute species_dir_abs now that species_id is known
    species_dir_abs = os.path.join(output_dir_abs, species_id)
    ont_dir_abs = os.path.join(species_dir_abs, "ONT")
    os.makedirs(ont_dir_abs, exist_ok=True)

    # Normalize ONT inputs BEFORE any chdir
    if pd.notna(ont_sra) and pd.isna(ont_raw_reads):
        # Expected dump location for SRA
        ont_raw_reads = os.path.join(ont_dir_abs, f"{ont_sra}.fastq")
    elif isinstance(ont_raw_reads, str) and not os.path.isabs(ont_raw_reads):
        # If given a relative path, make it absolute relative to output_dir_abs (project root),
        # not the current cwd (prevents nested repetition)
        ont_raw_reads = os.path.abspath(os.path.join(output_dir_abs, ont_raw_reads))

    # Illumina dedup (absolute paths so downstream tools can find them regardless of cwd)
    illu_dedup_f_reads = os.path.join(species_dir_abs, "Illumina", f"{species_id}_illu_forward_dedup.fastq")
    illu_dedup_r_reads = os.path.join(species_dir_abs, "Illumina", f"{species_id}_illu_reverse_dedup.fastq")

    print(f"DEBUG - ont_sra - {ont_sra}")
    print(f"DEBUG - ont_raw_reads - {ont_raw_reads}")
    print(f"DEBUG - ont_raw_dir - {ont_raw_dir}")
    print(f"DEBUG - illu_dedup_f_reads - {illu_dedup_f_reads}")
    print(f"DEBUG - illu_dedup_r_reads - {illu_dedup_r_reads}")

    if pd.isna(ont_sra) and pd.isna(ont_raw_dir) and pd.isna(ont_raw_reads):
        print("SKIP:\tONT preprocessing; no reads provided")
        return None

    # Work within the ONT directory
    prev_cwd = os.getcwd()
    os.chdir(ont_dir_abs)
    
    try:
        # Handle Raw Data Directory (concatenate .fastq files)
        if isinstance(ont_raw_dir, str):
            print(f"NOTE:\tConcatenating ONT files from {ont_raw_dir}")
            ont_raw_dir_abs = ont_raw_dir if os.path.isabs(ont_raw_dir) else os.path.join(output_dir_abs, ont_raw_dir)
            ont_files = sorted(glob.glob(os.path.join(ont_raw_dir_abs, "*.fastq")))
            if not ont_files:
                print(f"ERROR:\tNo ONT files found in {ont_raw_dir_abs}")
                sys.exit(1)
            ont_raw_reads = os.path.join(ont_dir_abs, f"{species_id}_ont_combined.fastq")
            with open(ont_raw_reads, "wb") as w:
                for f in ont_files:
                    with open(f, "rb") as r:
                        shutil.copyfileobj(r, w)
    
        # SRA download to the ONT directory (no relative -O that repeats the path)
        if isinstance(ont_sra, str):
            if not ont_raw_reads or not os.path.exists(ont_raw_reads):
                print(f"Downloading SRA {ont_sra} from GenBank...")
                _ = run_subprocess_cmd(["prefetch", "--force", "yes", ont_sra], False)
                _ = run_subprocess_cmd(["fasterq-dump", "-e", str(cpu_threads), "-O", ont_dir_abs, ont_sra], False)
                dumped = os.path.join(ont_dir_abs, f"{ont_sra}.fastq")
                if not os.path.exists(dumped):
                    print("ERROR:\tExpected FASTQ not found after fasterq-dump.")
                    return None
                if os.path.getsize(dumped) == 0:
                    print(f"ERROR:\tDownloaded FASTQ is empty: {dumped}")
                    return None
                ont_raw_reads = dumped
                print(f"PASS:\tSRA converted to FASTQ: {ont_raw_reads}")
            else:
                print(f"SKIP:\tSRA already exists: {ont_raw_reads}")
    
        # Parse estimated genome size
        m = re.match(r"^(\d+(?:\.\d+)?)(\D+)$", str(est_size))
        if m:
            est_size_numb, est_size_mult = m.group(1), m.group(2)
            multipliers = {'m': 10**6, 'g': 10**9}
            est_size_bp = int(float(est_size_numb) * multipliers.get(est_size_mult.lower(), 25_000_000))
        else:
            print(f"NOTE:\tUnable to parse input estimated size {est_size}, using default: 25000000")
            est_size_bp = 25_000_000
    
        if not (isinstance(ont_raw_reads, str) and os.path.exists(ont_raw_reads)):
            print(f"ERROR:\tONT reads not found: {ont_raw_reads}")
            return None
    
        # NanoPlot Raw Reads (soft-guard)
        try:
            sample_stats_dict = nanoplot_qc_reads(ont_raw_reads, "Raw_ONT_", cpu_threads, sample_stats_dict)
        except Exception as e:
            print(f"WARN:\tNanoPlot failed on raw ONT reads ({e}); continuing without NanoStats.")
    
        # Filtlong (only include Illumina if those files actually exist)
        filtered_ont = os.path.join(ont_dir_abs,
                                    f"{species_id}_ont_filtered.fastq" if not isinstance(ont_sra, str)
                                    else os.path.basename(ont_raw_reads).replace(f"{ont_sra}", f"{species_id}_ont_filtered"))
        coverage = 75
        target_bases = est_size_bp * coverage
        use_illumina = os.path.exists(illu_dedup_f_reads) and os.path.exists(illu_dedup_r_reads)
        illumina_opt = f"-1 {illu_dedup_f_reads} -2 {illu_dedup_r_reads}" if use_illumina else ""
    
        if not os.path.exists(filtered_ont):
            filtlong_cmd = f"filtlong {illumina_opt} --trim --min_length 1000 --min_mean_q 8 --keep_percent 90 --target_bases {target_bases} {ont_raw_reads} > {filtered_ont}"
            _ = run_subprocess_cmd(filtlong_cmd, True)
            if (not os.path.exists(filtered_ont)) or os.path.getsize(filtered_ont) == 0:
                print(f"ERROR:\tFiltlong did not produce reads at {filtered_ont}")
                return None
        else:
            print(f"SKIP\tFiltlong filtered reads exist: {filtered_ont}.")
    
        # NanoPlot Filtered Reads (soft-guard)
        try:
            sample_stats_dict = nanoplot_qc_reads(filtered_ont, "Filt_ONT_", cpu_threads, sample_stats_dict)
        except Exception as e:
            print(f"WARN:\tNanoPlot failed on filtered ONT reads ({e}); continuing.")
    
        # Ratatosk (only if Illumina reads available)
        corrected_out = os.path.join(ont_dir_abs, "ratatosk_corrected")
        if (not use_illumina):
            print("SKIP:\tRatatosk correction; no usable Illumina reads provided")
            final_corrected_ont = filtered_ont
            safe_alias = os.path.join(ont_dir_abs, f"{species_id}_ont_corrected.fastq")
            if not os.path.exists(safe_alias):
                subprocess.run(["ln", "-sf", filtered_ont, safe_alias], check=True)
        else:
            final_corrected_ont = os.path.join(ont_dir_abs, f"{species_id}_ont_corrected.fastq")
            if os.path.exists(corrected_out + ".fastq"):
                shutil.move(corrected_out + ".fastq", final_corrected_ont)
            if not os.path.exists(final_corrected_ont):
                ratatosk_cmd = ["Ratatosk", "correct", "-s", illu_dedup_f_reads, "-s", illu_dedup_r_reads,
                                "-l", filtered_ont, "-o", corrected_out, "-c", str(cpu_threads), "-v"]
                _ = run_subprocess_cmd(ratatosk_cmd, False)
                if os.path.exists(corrected_out + ".fastq"):
                    shutil.move(corrected_out + ".fastq", final_corrected_ont)
                else:
                    print(f"WARN:\tRatatosk did not produce {corrected_out}.fastq; using filtered reads.")
                    final_corrected_ont = filtered_ont
    
        # NanoPlot Corrected Reads (soft-guard)
        try:
            sample_stats_dict = nanoplot_qc_reads(final_corrected_ont, "Corr_ONT_", cpu_threads, sample_stats_dict)
        except Exception as e:
            print(f"WARN:\tNanoPlot failed on corrected ONT reads ({e}); continuing.")
    
        # Select best long reads (your helper uses output_dir/input_csv paths, unchanged)
        highest_mean_qual_long_reads = select_long_reads(output_dir_abs, input_csv_abs, sample_id, cpu_threads)
        highest = select_long_reads(output_dir_abs, input_csv_abs, sample_id, cpu_threads)
        if not highest:
            highest = final_corrected_ont
    
    finally:
        os.chdir(prev_cwd)

    print(f"PASS:\tPreprocessed Raw ONT Reads for {sample_id}: {highest_mean_qual_long_reads}.")
    
    return highest_mean_qual_long_reads


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
    
    highest_mean_qual_long_reads = preprocess_ont(sample_id, input_csv, output_dir, cpu_threads, ram_gb)
