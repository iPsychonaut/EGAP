#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
preprocess_illumina.py

Updated on Sat Mar 29 2025

This script preprocesses Illumina reads with FastQC, Trimmomatic, BBDuk, and Clumpify.
Reads input from a CSV file, processes all Illumina datasets, and updates the CSV
with declumpified FASTQ file paths in ${params.output_dir}/${sample_prefix}/Illumina/.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import os, sys
import pandas as pd
from utilities import run_subprocess_cmd, get_current_row_data, md5_check


# --------------------------------------------------------------
# Combine and verify Illumina reads
# --------------------------------------------------------------
def illumina_extract_and_check(folder_name, SAMPLE_ID):
    """Combine paired-end Illumina reads after MD5 verification.

    Verifies MD5 checksums, concatenates forward and reverse FASTQ files.

    Args:
        folder_name (str): Directory containing Illumina FASTQ files and MD5.txt.
        SAMPLE_ID (str): Sample identifier for naming output files.

    Returns:
        list or None: Paths to combined forward and reverse files, or None if failed.
    """
    print(f"Running MD5 Checksum Analysis on Raw Illumina FASTQ files in {folder_name}...")
    illumina_df = pd.DataFrame(columns=["MD5", "Filename"])
    base_folder = SAMPLE_ID.split("-")[0]  # e.g., "Es_coli"
    illumina_dir = os.path.join(os.getcwd(), base_folder, "Illumina")  # Use current dir, not hardcoded EGAP_Test_Data
    os.makedirs(illumina_dir, exist_ok=True)
    combined_1_file = os.path.join(illumina_dir, f"{SAMPLE_ID}_combined_1.fq")
    combined_2_file = os.path.join(illumina_dir, f"{SAMPLE_ID}_combined_2.fq")
    combined_list = [combined_1_file, combined_2_file]
    
    if not os.path.isfile(combined_list[0]) or not os.path.isfile(combined_list[1]):
        if not os.path.isfile(combined_1_file) or not os.path.isfile(combined_2_file):
            md5_check(folder_name, illumina_df)
            raw_1_list = []
            raw_2_list = []
            for filename in os.listdir(folder_name):
                if "_1.fastq" in filename:
                    raw_1_list.append(os.path.join(folder_name, filename))
                elif "_2.fastq" in filename:
                    raw_2_list.append(os.path.join(folder_name, filename))
            if not raw_1_list or not raw_2_list:
                print(f"ERROR:\tNo paired Illumina files found in {folder_name}")
                return None
            fwd_cat_cmd = f"cat {' '.join(raw_1_list)} > {combined_1_file}"
            _ = run_subprocess_cmd(fwd_cat_cmd, shell_check=True)
            rev_cat_cmd = f"cat {' '.join(raw_2_list)} > {combined_2_file}"
            _ = run_subprocess_cmd(rev_cat_cmd, shell_check=True)
        else:
            print(f"SKIP:\tCombined FASTQ files already exist: {combined_list[0]}; {combined_list[1]}.")
    else:
        print(f"SKIP:\tGzipped Combined FASTQ files already exist: {combined_list[0]}; {combined_list[1]}.")

    return combined_list if os.path.exists(combined_list[0]) and os.path.exists(combined_list[1]) else None


# --------------------------------------------------------------
# Preprocess Illumina sequencing reads
# --------------------------------------------------------------
def preprocess_illumina(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """Preprocess Illumina reads with FastQC, Trimmomatic, BBDuk, and Clumpify.

    Processes Illumina paired-end reads from SRA or raw files, performs quality control,
    trimming, filtering, and deduplication, and saves processed reads.

    Args:
        sample_id (str): Sample identifier.
        input_csv (str): Path to metadata CSV file.
        output_dir (str): Directory for output files.
        cpu_threads (int): Number of CPU threads to use.
        ram_gb (int): Available RAM in GB.

    Returns:
        tuple or (None, None): Paths to deduplicated forward and reverse FASTQ files,
                               or (None, None) if no reads are available or processing is skipped.
    """ 
    print(f"Preprocessing Illumina reads for {sample_id.split('-')[0]}...")
    # Parse the CSV and retrieve relevant row data
    input_df = pd.read_csv(input_csv)
    print(f"DEBUG - input_df - {input_df}")
    
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]  # Convert to Series (single row)

    print(f"DEBUG - current_series - {current_series}")

    # Identify read paths, reference, and BUSCO lineage info from CSV
    illu_sra = current_series["ILLUMINA_SRA"]
    illu_raw_f_reads = current_series["ILLUMINA_RAW_F_READS"]
    illu_raw_r_reads = current_series["ILLUMINA_RAW_R_READS"]
    illu_raw_dir = current_series["ILLUMINA_RAW_DIR"]
    species_id = current_series["SPECIES_ID"]

    print(f"DEBUG - illu_sra - {illu_sra}")   
    print(f"DEBUG - illu_raw_f_reads - {illu_raw_f_reads}")
    print(f"DEBUG - illu_raw_r_reads - {illu_raw_f_reads}")
    print(f"DEBUG - illu_raw_dir - {illu_raw_dir}")
    
    if pd.isna(illu_sra) and pd.isna(illu_raw_f_reads) and pd.isna(illu_raw_r_reads) and pd.isna(illu_raw_dir): 
        print(f"SKIP:\tSample does not include Illumina Reads: {sample_id}.")
        return None, None

    species_dir = os.path.join(output_dir, species_id)
    os.makedirs(species_dir, exist_ok=True)
    illumina_dir = os.path.join(species_dir, "Illumina")
    os.makedirs(illumina_dir, exist_ok=True)
    os.chdir(illumina_dir)

    if pd.notna(illu_sra) and pd.isna(illu_raw_r_reads) and pd.isna(illu_raw_f_reads):
        if pd.isna(illu_raw_dir):
            illu_raw_f_reads = os.path.join(illumina_dir, f"{illu_sra}_1.fastq")
            illu_raw_r_reads = os.path.join(illumina_dir, f"{illu_sra}_2.fastq")
            if not os.path.exists(illu_raw_f_reads) and not os.path.exists(illu_raw_r_reads):
                print(f"Downloading Illumina SRA: {illu_sra}...")
                os.chdir(illumina_dir)
                illu_cmd = f"prefetch --force yes {illu_sra} && fasterq-dump -e {cpu_threads} -O {illumina_dir} {illu_sra}"
                if run_subprocess_cmd(illu_cmd, True) != 0:
                    print(f"ERROR:\tFailed to download Illumina SRA {illu_sra}")
                    illu_raw_f_reads = None
                    illu_raw_r_reads = None
            else:
                print(f"PASS:\tIllumina SRA processed: {illu_raw_f_reads}, {illu_raw_r_reads}")
        else:
            print(f"Process Illumina Raw Directory: {illu_raw_dir}...")
            combined_list = illumina_extract_and_check(illu_raw_dir, sample_id)
            if combined_list:
                print("PASS:\tSucessfully processed Illumina Raw Directory.")
                illu_raw_f_reads = combined_list[0]
                illu_raw_r_reads = combined_list[1]
                illu_sra = f"{combined_list[0]},{combined_list[1]}"  # Update ILLUMINA_SRA
            else:
                print(f"ERROR:\tFailed to process Illumina Raw Directory for {sample_id}")
                return None, None
    
    illu_dedup_f_reads = os.path.join(illumina_dir, f"{species_id}_illu_forward_dedup.fastq")
    illu_dedup_r_reads = os.path.join(illumina_dir, f"{species_id}_illu_reverse_dedup.fastq")

    if os.path.exists(illu_dedup_f_reads) and os.path.exists(illu_dedup_r_reads):
        print(f"SKIP:\tIllumina preprocessing already completed: {illu_dedup_f_reads}, {illu_dedup_r_reads}.")
        return illu_dedup_f_reads, illu_dedup_r_reads

    # FastQC
    run_subprocess_cmd(["fastqc", "-t", str(cpu_threads), "-o", "fastqc_results", illu_raw_f_reads, illu_raw_r_reads], False)

    # Trimmomatic
    trimmo_f_pair = illu_raw_f_reads.replace("_1","_forward_paired") 
    trimmo_r_pair = illu_raw_r_reads.replace("_2","_reverse_paired")
    trimmo_f_unpair = illu_raw_f_reads.replace("_1","_forward_unpaired")
    trimmo_r_unpair = illu_raw_r_reads.replace("_2","_reverse_unpaired")
    
    if os.path.exists(trimmo_f_pair) and os.path.exists(trimmo_r_pair):
        print(f"SKIP:\tTrimmomatic files exist: {trimmo_f_pair} & {trimmo_r_pair}.")
    else:
        run_subprocess_cmd(["trimmomatic", "PE",
                            "-threads", str(cpu_threads),
                            "-phred33",
                            illu_raw_f_reads,           # Input R1
                            illu_raw_r_reads,           # Input R2
                            trimmo_f_pair,              # R1 paired output
                            trimmo_f_unpair,            # R1 unpaired output
                            trimmo_r_pair,              # R2 paired output
                            trimmo_r_unpair,            # R2 unpaired output
                            "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:11",
                            "HEADCROP:10",
                            "CROP:145",
                            "SLIDINGWINDOW:50:25",
                            "MINLEN:125"],
                           False)
    # BBDuk
    bbduk_f_map = illu_raw_f_reads.replace("_1","_forward_mapped") 
    bbduk_r_map = illu_raw_r_reads.replace("_2","_reverse_mapped") 
    if os.path.exists(bbduk_f_map) and os.path.exists(bbduk_r_map):
        print(f"SKIP:\tbbduk mapped files exist: {bbduk_f_map} & {bbduk_r_map}.")
    else:
        run_subprocess_cmd(["bbduk.sh", f"in1={trimmo_f_pair}", f"in2={trimmo_r_pair}",
                            f"out1={bbduk_f_map}", f"out2={bbduk_r_map}",
                            "ktrim=r", "k=23", "mink=11", "hdist=1",
                            "tpe", "tbo", "qtrim=rl", "trimq=20"], False)

    # Clumpify
    if os.path.exists(illu_dedup_f_reads) and os.path.exists(illu_dedup_r_reads):
        print(f"SKIP:\tClumpify deduplicated files exist: {illu_dedup_f_reads} & {illu_dedup_r_reads}.")
    else:
        run_subprocess_cmd(["clumpify.sh", f"in={bbduk_f_map}", f"in2={bbduk_r_map}",
                        f"out={illu_dedup_f_reads}", f"out2={illu_dedup_r_reads}", "dedupe"], False)

    print(f"PASS:\tPreprocessed Raw Illumina Reads for {species_id}: {illu_dedup_f_reads}, {illu_dedup_r_reads}.")
    
    # FastQC
    run_subprocess_cmd(["fastqc", "-t", str(cpu_threads), "-o", "fastqc_results", illu_dedup_f_reads, illu_dedup_r_reads], False)
    
    return illu_dedup_f_reads, illu_dedup_r_reads


if __name__ == "__main__":
    # Log raw sys.argv immediately
    print(f"DEBUG: Raw sys.argv = {sys.argv}")
    print(f"DEBUG: Length of sys.argv = {len(sys.argv)}")
    
    # Check argument count
    if len(sys.argv) != 6:
        print(f"ERROR: Expected 5 arguments (plus script name), got {len(sys.argv)-1}: {sys.argv[1:]}", 
              file=sys.stderr)
        print("Usage: python3 preprocess_illumina.py <sample_id> <input_csv> <output_dir> <cpu_threads> <ram_gb>", 
              file=sys.stderr)
        sys.exit(1)
    
    # Log each argument
    for i, arg in enumerate(sys.argv):
        print(f"DEBUG: sys.argv[{i}] = '{arg}'")
    
    sample_id = sys.argv[1]
    input_csv = sys.argv[2]
    output_dir = sys.argv[3]
    cpu_threads = int(sys.argv[4])
    ram_gb = int(sys.argv[5]) if sys.argv[5] != " " else 8
    
    print(f"DEBUG: Parsed sample_id = '{sample_id}' {type(sample_id)}")
    print(f"DEBUG: Parsed input_csv = '{input_csv}' {type(input_csv)}")
    print(f"DEBUG: Parsed output_dir = '{output_dir}' {type(output_dir)}")
    print(f"DEBUG: Parsed cpu_threads = '{sys.argv[4]}' {sys.argv[4]} (converted to {cpu_threads}) {type(cpu_threads)}")
    print(f"DEBUG: Parsed ram_gb = '{sys.argv[5]}' {sys.argv[5]} (converted to {ram_gb}) {type(ram_gb)}")
    
    illu_dedup_f_reads, illu_dedup_r_reads = preprocess_illumina(sample_id, input_csv, output_dir, cpu_threads, ram_gb)
