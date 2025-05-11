#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
preprocess_refseq.py

Updated on Sat Mar 29 2025


@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import os, sys, glob, shutil
import pandas as pd
from utilities import run_subprocess_cmd, get_current_row_data

# --------------------------------------------------------------
# Preprocess reference sequence data
# --------------------------------------------------------------
def preprocess_refseq(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """Preprocess reference sequence data for the assembly pipeline.

    Downloads and verifies reference sequence data (GCA or FASTA)
    based on metadata from a CSV file.

    Args:
        sample_id (str): Sample identifier.
        input_csv (str): Path to metadata CSV file.
        output_dir (str): Directory for output files.
        cpu_threads (int): Number of CPU threads to use.
        ram_gb (int): Available RAM in GB.

    Returns:
        str or None: Path to reference sequence file, or None if no reference is available.
    """  
    print(f"Preprocessing Reference Sequence assembly for {sample_id.split('-')[0]}...")
    # Parse the CSV and retrieve relevant row data
    input_df = pd.read_csv(input_csv)
    print(f"DEBUG - input_df - {input_df}")
    
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]  # Convert to Series (single row)

    # Identify read paths, reference, and BUSCO lineage info from CSV
    ref_seq_gca = current_series["REF_SEQ_GCA"]
    ref_seq = current_series["REF_SEQ"]
    species_id = current_series["SPECIES_ID"]

    species_dir = os.path.join(output_dir, species_id)
    os.makedirs(species_dir, exist_ok=True)
    refseq_dir = os.path.join(species_dir, "RefSeq")
    os.makedirs(refseq_dir, exist_ok=True)
    os.chdir(refseq_dir)
    
    if pd.notna(ref_seq_gca) and pd.isna(ref_seq):
        ref_seq = os.path.join(species_dir, "RefSeq", f"{ref_seq_gca}.fasta")  

    print(f"DEBUG - ref_seq_gca - {ref_seq_gca}")   
    print(f"DEBUG - ref_seq - {ref_seq}")
    
    if pd.isna(ref_seq_gca) and pd.isna(ref_seq): 
        print(f"SKIP:\tSample does not include Reference Sequence: {sample_id}.")
        return None
    
    if pd.notna(ref_seq_gca) and "." not in ref_seq_gca:
        print(f"ERROR:\tReference Sequence GCA requires version number: {ref_seq_gca} has no '.#' ")
    renamed_gca = os.path.join(refseq_dir, f"{species_id}_{ref_seq_gca}_RefSeq.fasta")
    if os.path.exists(ref_seq_gca):
        shutil.move(ref_seq_gca, renamed_gca)
        print(f"Found existing Reference Sequence: {renamed_gca}.")  
        return renamed_gca
    elif os.path.exists(renamed_gca):
        print(f"Found existing Reference Sequence: {renamed_gca}.")
        return renamed_gca
    else:
        ref_seq_gca_dir = os.path.join(refseq_dir, f"ncbi_dataset/data/{ref_seq_gca}/")
        renamed_gca = os.path.join(refseq_dir, f"{species_id}_{ref_seq_gca}_RefSeq.fasta")
        if not pd.isna(ref_seq_gca):
            os.chdir(refseq_dir)
            print(f"Downloading Reference Sequence Assembly for: {ref_seq_gca}")
            if not os.path.exists(renamed_gca):
                try:
                    ref_seq_gca = glob.glob(os.path.join(ref_seq_gca_dir, "*_genomic.fna"))[0]
                    if not os.path.exists(ref_seq_gca):
                        ref_seq_cmd = f"datasets download genome accession {ref_seq_gca} --include genome &&  unzip -o ncbi_dataset -d {refseq_dir}"
                        _ = run_subprocess_cmd(ref_seq_cmd, shell_check=True)
                except IndexError:
                    ref_seq_cmd = f"datasets download genome accession {ref_seq_gca} --include genome &&  unzip -o ncbi_dataset -d {refseq_dir}"
                    _ = run_subprocess_cmd(ref_seq_cmd, shell_check=True)
                pattern = os.path.join(ref_seq_gca_dir, "*_genomic.fna")
                ref_seq_gca = next(glob.iglob(pattern), None)                
                if ref_seq_gca is None:
                    raise FileNotFoundError(f"No ‘*_genomic.fna’ found in {ref_seq_gca_dir!r}")
                print(f"PASS:\tSuccessfully downloaded the GCA to: {ref_seq_gca}.")
                shutil.move(ref_seq_gca, renamed_gca)
                print(f"PASS:\tSuccessfully moved and renamed the GCA to: {renamed_gca}.")
            else:
                print(f"SKIP:\tREF_SEQ GCA already exists: {renamed_gca}")
    
    return renamed_gca


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
    
    renamed_gca = preprocess_refseq(sample_id, input_csv, output_dir, cpu_threads, ram_gb)
