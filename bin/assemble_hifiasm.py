#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
assemble_hifiasm.py

Updated on Sat Mar 29 2025

This script runs Hifiasm assembly with PacBio reads.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import os, sys
import pandas as pd
from utilities import run_subprocess_cmd, get_current_row_data
from qc_assessment import qc_assessment


# --------------------------------------------------------------
# Run Hifiasm assembly with PacBio reads
# --------------------------------------------------------------
def assemble_hifiasm(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """Assemble genomic data using Hifiasm with PacBio reads.

    Executes Hifiasm assembly with PacBio reads, converts the output GFA to FASTA,
    and performs quality control on the resulting assembly.

    Args:
        sample_id (str): Sample identifier.
        input_csv (str): Path to metadata CSV file.
        output_dir (str): Directory for output files.
        cpu_threads (int): Number of CPU threads to use.
        ram_gb (int): Available RAM in GB.

    Returns:
        str or None: Path to the gzipped Hifiasm assembly FASTA, or None if no PacBio reads are provided.
    """
    input_df = pd.read_csv(input_csv)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]  # Convert to Series (single row)

    # Identify read paths, reference, and BUSCO lineage info from CSV
    pacbio_sra = current_series["PACBIO_SRA"]
    pacbio_raw_reads = current_series["PACBIO_RAW_READS"]
    ref_seq_gca = current_series["REF_SEQ_GCA"]
    ref_seq = current_series["REF_SEQ"]
    species_id = current_series["SPECIES_ID"]

    species_dir = os.path.join(output_dir, species_id)
    
    if pd.isna(pacbio_sra) and pd.isna(pacbio_raw_reads):
        print("SKIP:\tNo PacBio reads files provided.")
        return None
    
    if pd.notna(pacbio_sra) and pd.isna(pacbio_raw_reads):
        pacbio_raw_reads = os.path.join(species_dir, "PacBio", f"{pacbio_sra}.fastq")

    print(f"DEBUG - pacbio_sra - {pacbio_sra}")
    print(f"DEBUG - pacbio_raw_reads - {pacbio_raw_reads}")
    print(f"DEBUG - ref_seq_gca - {ref_seq_gca}")
    print(f"DEBUG - ref_seq - {ref_seq}")   
    print(f"DEBUG - species_id - {species_id}")
    
    # Set long-read paths (ONT or PacBio), prefer prefiltered, fallback to raw
    highest_mean_qual_long_reads = None
    if pd.notna(pacbio_raw_reads):
        print("DEBUG - PACBIO RAW READS EXIST!")
        candidate = os.path.join(species_dir, "PacBio", f"{species_id}_PacBio_highest_mean_qual_long_reads.fastq")
        highest_mean_qual_long_reads = candidate if os.path.exists(candidate) else pacbio_raw_reads

    print(f"DEBUG - highest_mean_qual_long_reads    - {highest_mean_qual_long_reads}")

    sample_dir = os.path.join(species_dir, sample_id)
    hifiasm_out_dir = os.path.join(sample_dir, "hifiasm_assembly")
    os.makedirs(hifiasm_out_dir, exist_ok=True)
    
    egap_hifiasm_assembly_path = os.path.join(hifiasm_out_dir, f"{sample_id}_hifiasm.fasta")

    # Run Hifiasm
    pacbio_prefix = os.path.join(hifiasm_out_dir,
                                 os.path.basename(pacbio_raw_reads).split(".")[0])
    hifiasm_path = pacbio_prefix + ".asm"
    hifiasm_gfa  = pacbio_prefix + ".asm.bp.p_ctg.gfa"

    if os.path.exists(hifiasm_gfa):
        print(f"SKIP:\tHiFi Assembly already exists: {hifiasm_gfa}")
    else:
        # Ensure work directory output
        starting_work_dir = os.getcwd()
        if "work" not in starting_work_dir:
            current_work_dir = hifiasm_out_dir
        else:
            current_work_dir = starting_work_dir
        os.chdir(current_work_dir)
        
        hifiasm_cmd = ["hifiasm",
                       "-o", hifiasm_path,
                       "-t", str(cpu_threads),
                       highest_mean_qual_long_reads]
        _ = run_subprocess_cmd(hifiasm_cmd, shell_check=False)

    # Convert GFA to FASTA
    gfa_cmd = f"gfatools gfa2fa {hifiasm_gfa} > {egap_hifiasm_assembly_path}"
    _ = run_subprocess_cmd(gfa_cmd, shell_check=True)         
    
    egap_hifiasm_assembly_path, hifiasm_stats_list, _ = qc_assessment("hifiasm", input_csv, sample_id, output_dir, cpu_threads, ram_gb)        

    return egap_hifiasm_assembly_path


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python3 assemble_hifiasm.py <sample_id> <input_csv> "
            "<output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)
        
    egap_hifiasm_assembly_path = assemble_hifiasm(sys.argv[1],       # sample_id
                                                     sys.argv[2],       # input_csv
                                                     sys.argv[3],       # output_dir
                                                     str(sys.argv[4]),  # cpu_threads
                                                     str(sys.argv[5]))  # ram_gb