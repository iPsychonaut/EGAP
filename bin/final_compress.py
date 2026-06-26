# -*- coding: utf-8 -*-
"""
final_compress.py

Compress every remaining FASTA/FASTQ file under a sample's output directory.

Called at the very end of an EGAP run to reduce on-disk footprint once all
analyses that require uncompressed inputs have completed. Uses pigz for
parallel gzip compression; the sample directory is discovered from the
input TSV's SPECIES_ID / SAMPLE_ID columns rather than a user-supplied path.

Stage:
    Final Compression

Created on Fri May 16 15:52:30 2025

Updated on 2026-04-16

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os
import sys
import pandas as pd
from utilities import pigz_compress, get_current_row_data, read_sample_table


def final_compress(sample_id, input_tsv, output_dir, cpu_threads, ram_gb):
    """Compress every ``.fasta`` and ``.fastq`` file under the sample directory.

    The sample directory is resolved from the input TSV's SPECIES_ID and
    SAMPLE_ID columns as ``{output_dir}/{species_id}/{sample_id}``. Every
    uncompressed FASTA or FASTQ file under that tree is compressed in place
    with pigz using the supplied thread count.

    Args:
        sample_id (str): The SAMPLE_ID value identifying which TSV row (and
            therefore which sample directory) to compress.
        input_tsv (str): Path to the EGAP input TSV containing SPECIES_ID
            and SAMPLE_ID columns.
        output_dir (str): Root output directory for the EGAP run — species
            subdirectories live directly under this path.
        cpu_threads (int | str): Number of threads pigz may use when
            compressing each file.
        ram_gb (int | str): RAM budget passed for API compatibility with
            other EGAP steps. Not currently used by this step.

    Returns:
        None: Files are compressed in place; the function has no return value.
    """
    print(f"Compressing all FASTA and FASTQ files for {sample_id}...")

    # Read the sample table and filter to the row corresponding to the sample of interest
    input_df = read_sample_table(input_tsv)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]  # Convert to Series (single row)

    # Identify read paths, reference, and BUSCO lineage info from TSV
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
        print("Usage: python3 final_compress.py <input_tsv> "
              "<sample_id> <output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)

    final_compress(sys.argv[1],       # input_tsv
                   sys.argv[2],       # sample_id
                   sys.argv[3],       # output_dir
                   str(sys.argv[4]),  # cpu_threads
                   str(sys.argv[5]))  # ram_gb