# -*- coding: utf-8 -*-
"""
final_compress.py

PEP 8 Documentation

Created on Fri May 16 15:52:30 2025

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os, sys
import pandas as pd
from utilities import pigz_compress, get_current_row_data


def final_compress(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """
    PEP 8 Documentation
    """
    print(f"Compressing all FASTA and FASTQ files for {sample_id}...")

    # Read the CSV file and filter to the row corresponding to the sample of interest
    input_df = pd.read_csv(input_csv)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]  # Convert to Series (single row)

    # Identify read paths, reference, and BUSCO lineage info from CSV
    species_id = current_series["SPECIES_ID"]
    species_dir = os.path.join(output_dir, species_id)
    sample_dir = os.path.join(species_dir, sample_id)

    print(f"DEBUG - sample_dir - {sample_dir}")

    # Walk through directory and subdirectories and multi-thread compress ALL FASTA or FASTQ files
    for root, dirs, files in os.walk(sample_dir):
        for file in files:
            if file.endswith((".fasta", ".fastq")):
                full_path = os.path.join(root, file)
                print(f"Compressing: {full_path}")
                _ = pigz_compress(full_path, cpu_threads)

    print("PASS:\tAll FASTA and FASTQ successfully compressed!")
            
            
if __name__ == "__main__":
    # Handle command-line arguments
    if len(sys.argv) != 6:
        print("Usage: python3 final_compress.py <input_csv> "
              "<sample_id> <output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)

    final_compress(sys.argv[1],       # input_csv
                   sys.argv[2],       # sample_id
                   sys.argv[3],       # output_dir
                   str(sys.argv[4]),  # cpu_threads
                   str(sys.argv[5]))  # ram_gb