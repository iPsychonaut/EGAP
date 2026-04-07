#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
utilities.py

Module containing regularly used commands in various other EGAP scripts.

Created on Wed Aug 16 2023

Updated on Wed Sept 3 2025

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os
import subprocess
import platform
import shutil
import math
import hashlib
import tempfile
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from typing import List
from datetime import datetime

# --------------------------------------------------------------
# Catches and unzips compressed files for FASTQ
# --------------------------------------------------------------
def sum_fastq_bases_with_pigz_safe(fq_path: Path, cpu_threads: int) -> int:
    """
    Safely sum read lengths from a FASTQ that may be .gz by:
    - copying the .gz to a temp dir,
    - pigz_decompress() on the *copy*,
    - parsing the decompressed temp file,
    - removing the temp dir.
    """
    total = 0
    if str(fq_path).endswith(".gz"):
        with tempfile.TemporaryDirectory(prefix="egap_pigz_") as tdir:
            tmp_gz = Path(tdir) / fq_path.name
            shutil.copy2(str(fq_path), str(tmp_gz))                 # safe copy
            tmp_fastq = pigz_decompress(str(tmp_gz), cpu_threads)   # your function returns str path
            with open(tmp_fastq, "rt") as handle:
                for rec in SeqIO.parse(handle, "fastq"):
                    total += len(rec.seq)
            # temp dir cleanup happens automatically
    else:
        with open(fq_path, "rt") as handle:
            for rec in SeqIO.parse(handle, "fastq"):
                total += len(rec.seq)
    return total


# --------------------------------------------------------------
# Catches and unzips compressed files for FASTA
# --------------------------------------------------------------
def sum_fasta_bases_with_pigz_safe(fa_path: Path, cpu_threads: int) -> int:
    """
    Same safety pattern for FASTA/FA files that may be .gz.
    """
    total = 0
    if str(fa_path).endswith(".gz"):
        with tempfile.TemporaryDirectory(prefix="egap_pigz_") as tdir:
            tmp_gz = Path(tdir) / fa_path.name
            shutil.copy2(str(fa_path), str(tmp_gz))
            tmp_fa = pigz_decompress(str(tmp_gz), cpu_threads)      # returns str
            with open(tmp_fa, "rt") as handle:
                for rec in SeqIO.parse(handle, "fasta"):
                    total += len(rec.seq)
    else:
        with open(fa_path, "rt") as handle:
            for rec in SeqIO.parse(handle, "fasta"):
                total += len(rec.seq)
    return total


# --------------------------------------------------------------
# Calculates Coverage for a given Genome
# --------------------------------------------------------------
def calculate_genome_coverage(read_fastqs: List[str], assembly_fasta: str, cpu_threads: int) -> float:
    """
    Compute genome coverage = (total bases in all reads) / (total bases in the assembly).
    Uses your pigz_{compress,decompress} functions safely (no mutation of inputs).
    Supports .fastq/.fq(.gz) for reads and .fa/.fasta(.gz) for assembly.
    """
    total_bases = 0
    for fq in read_fastqs:
        fq_path = Path(fq)
        total_bases += sum_fastq_bases_with_pigz_safe(fq_path, cpu_threads)

    assembly_bases = sum_fasta_bases_with_pigz_safe(Path(assembly_fasta), cpu_threads)
    if assembly_bases == 0:
        raise ValueError(f"No contigs found in {assembly_fasta}")

    return total_bases / assembly_bases


# --------------------------------------------------------------
# Verify MD5 checksums for Illumina files
# --------------------------------------------------------------
def md5_check(folder_name, illumina_df):
    """Verify MD5 checksums for Illumina files in the specified folder.

    Compares computed MD5 checksums against those listed in ``MD5.txt``.
    Prints a PASS or ERROR message for each file.

    Parameters
    ----------
    folder_name : str
        Directory containing Illumina FASTQ files and ``MD5.txt``.
    illumina_df : pandas.DataFrame
        DataFrame used to accumulate MD5 and filename entries.
    """
    md5_file = os.path.join(folder_name, "MD5.txt")
    if not os.path.exists(md5_file):
        print(f"WARNING: MD5.txt not found in {folder_name}. Skipping MD5 check.")
        return
    with open(md5_file, "r") as f:
        for line in f:
            md5, filename = line.strip().split()
            illumina_df = illumina_df.append({"MD5": md5, "Filename": filename}, ignore_index=True)
    
    for index, row in illumina_df.iterrows():
        file_path = os.path.join(folder_name, row["Filename"])
        if os.path.exists(file_path):
            with open(file_path, "rb") as f:
                file_hash = hashlib.md5(f.read()).hexdigest()
            if file_hash == row["MD5"]:
                print(f"PASS: MD5 check passed for {row['Filename']}")
            else:
                print(f"ERROR: MD5 check failed for {row['Filename']}. Expected {row['MD5']}, got {file_hash}")
        else:
            print(f"ERROR: File not found for MD5 check: {file_path}")
            

# --------------------------------------------------------------
# Calculate resource allocation
# --------------------------------------------------------------
def get_resource_values(percent_resources, total_cpu, total_ram):
    """Calculate CPU threads and RAM based on a percentage of total resources.

    Parameters
    ----------
    percent_resources : float
        Fraction of total resources to allocate (0.0–1.0).
    total_cpu : int
        Total available CPU threads.
    total_ram : int
        Total available RAM in GB.

    Returns
    -------
    tuple of (int, int)
        ``(cpu_threads, ram_gb)`` — number of threads and RAM in GB.
    """
    cpu_threads = int(math.floor(total_cpu * percent_resources))
    ram_gb = int(total_ram * percent_resources)
    return cpu_threads, ram_gb


# --------------------------------------------------------------
# Execute a subprocess command and log its output
# --------------------------------------------------------------
def run_subprocess_cmd(cmd_list, shell_check):
    """Execute a subprocess command and stream its output to stdout.

    Runs the command using ``subprocess.Popen``, captures combined
    stdout/stderr, and streams each line in real time.  Returns the
    process exit code.

    Parameters
    ----------
    cmd_list : str or list of str
        Command to execute, either as a shell string or as an argument list.
    shell_check : bool
        If ``True``, execute the command through the shell (required when
        *cmd_list* is a string with shell operators).

    Returns
    -------
    int
        The subprocess return code (0 on success, 127 if the executable
        was not found).
    """
    cmd_display = cmd_list if isinstance(cmd_list, str) else ' '.join(cmd_list)
    print(f"CMD:\t{cmd_display}")
    try:
        process = subprocess.Popen(cmd_list, shell=shell_check, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True)
    except (FileNotFoundError, PermissionError) as exc:
        # The executable was not found on PATH or is not executable.
        # This happens when the pipeline is launched from a conda env that
        # does not have the required tool installed.  Return 127 (the
        # conventional "command not found" exit code) instead of crashing.
        tool = cmd_list.split()[0] if isinstance(cmd_list, str) else cmd_list[0]
        print(f"ERROR:\tCould not launch '{tool}': {exc}. "
              f"Make sure the correct conda environment is active "
              f"(e.g. conda activate EGAP_env).")
        return 127
    for line in process.stdout:
        print(line, end="")
    process.wait()
    if process.returncode != 0:
        print(f"NOTE:\tCommand failed with return code {process.returncode}")
    else:
        print(f"PASS:\tSuccessfully processed command: {cmd_display}")
    return process.returncode


# --------------------------------------------------------------
# Create a sample statistics dictionary from metadata
# --------------------------------------------------------------
def gen_sample_stats_dict(row):
    """Generate a sample statistics dictionary from a metadata row.

    Extracts key fields from a pandas Series and initializes placeholder
    values for all downstream quality metrics.

    Parameters
    ----------
    row : pandas.Series
        Single metadata row containing sample information columns.

    Returns
    -------
    dict
        Dictionary with all pipeline statistics fields initialized to
        ``None`` (except those extracted from *row*).
    """
    sample_stats_dict = {"SAMPLE_ID": row["SAMPLE_ID"],
                         "SPECIES_ID": row["SPECIES_ID"],
                         "ONT_SRA": row["ONT_SRA"] if isinstance(row["ONT_SRA"], str) else None,
                         "ONT": os.path.basename(row["ONT_RAW_READS"]) if isinstance(row["ONT_RAW_READS"], str) else None,
                         "ILLU_SRA": row["ILLUMINA_SRA"] if isinstance(row["ILLUMINA_SRA"], str) else None,
                         "ILLU_F": os.path.basename(row["ILLUMINA_RAW_F_READS"]) if isinstance(row["ILLUMINA_RAW_F_READS"], str) else None,
                         "ILLU_R": os.path.basename(row["ILLUMINA_RAW_R_READS"]) if isinstance(row["ILLUMINA_RAW_R_READS"], str) else None,
                         "PACBIO_SRA": row["PACBIO_SRA"] if isinstance(row["PACBIO_SRA"], str) else None,
                         "PACBIO": os.path.basename(row["PACBIO_RAW_READS"]) if isinstance(row["PACBIO_RAW_READS"], str) else None,
                         "REF_SEQ_GCA": row["REF_SEQ_GCA"] if isinstance(row["REF_SEQ_GCA"], str) else None,
                         "REF_SEQ": os.path.basename(row["REF_SEQ"]) if isinstance(row["REF_SEQ"], str) else None,
                         "RAW_ILLU_TOTAL_BASES": None,
                         "RAW_ILLU_COVERAGE": None,
                         "TRIMMED_ILLU_TOTAL_BASES": None,
                         "TRIMMED_ILLU_COVERAGE": None,
                         "DEDUPED_ILLU_TOTAL_BASES": None,
                         "DEDUPED_ILLU_COVERAGE": None,
                         "RAW_ONT_READS": None,
                         "RAW_ONT_MEAN_LENGTH": None,
                         "RAW_ONT_MEAN_QUAL": None,
                         "RAW_ONT_TOTAL_BASES": None,
                         "RAW_ONT_COVERAGE": None,
                         "FILT_ONT_READS": None,
                         "FILT_ONT_MEAN_LENGTH": None,
                         "FILT_ONT_MEAN_QUAL": None,
                         "FILT_ONT_TOTAL_BASES": None,
                         "FILT_ONT_COVERAGE": None,
                         "CORRECT_ONT_READS": None,
                         "CORRECT_ONT_MEAN_LENGTH": None,
                         "CORRECT_ONT_MEAN_QUAL": None,
                         "CORRECT_ONT_TOTAL_BASES": None,
                         "CORRECT_ONT_COVERAGE": None,
                         "KMER_COMPLETENESS": None,
                         "QUAL_VAL": None,
                         
                         "RAW_PACBIO_READS": None,
                         "RAW_PACBIO_MEAN_LENGTH": None,
                         "RAW_PACBIO_MEAN_QUAL": None,
                         "RAW_PACBIO_TOTAL_BASES": None,
                         "RAW_PACBIO_COVERAGE": None,
                         "HIFI_PACBIO_READS": None,
                         "HIFI_PACBIO_MEAN_LENGTH": None,
                         "HIFI_PACBIO_MEAN_QUAL": None,
                         "HIFI_PACBIO_TOTAL_BASES": None,
                         "HIFI_PACBIO_COVERAGE": None,
                         "FILT_PACBIO_READS": None,
                         "FILT_PACBIO_MEAN_LENGTH": None,
                         "FILT_PACBIO_MEAN_QUAL": None,
                         "FILT_PACBIO_TOTAL_BASES": None,
                         "FILT_PACBIO_COVERAGE": None,                        
                         
                         "FIRST_COMPLEASM_S": None,
                         "FIRST_COMPLEASM_D": None,
                         "FIRST_COMPLEASM_F": None,
                         "FIRST_COMPLEASM_M": None,
                         "FIRST_COMPLEASM_C": None,
                         "SECOND_COMPLEASM_S": None,
                         "SECOND_COMPLEASM_D": None,
                         "SECOND_COMPLEASM_F": None,
                         "SECOND_COMPLEASM_M": None,
                         "SECOND_COMPLEASM_C": None,
                         "GENOME_SIZE": None,
                         "ASSEMBLY_READS": None,
                         "ASSEMBLY_CONTIGS": None,
                         "ASSEMBLY_N50": None,
                         "ASSEMBLY_L50": None,
                         "ASSEMBLY_GC": None,
                         "MISASSEMBLIES": None,
                         "N_PER_100KBP": None,
                         "MIS_PER_100KBP": None,
                         "INDELS_PER_100KPB": None,
                         "FINAL_ASSEMBLY": None}
    return sample_stats_dict


# --------------------------------------------------------------
# Compress a file using pigz
# --------------------------------------------------------------
def pigz_compress(input_file, cpu_threads):
    """Compress a file using pigz with multiple threads.

    Parameters
    ----------
    input_file : str
        Path to the file to compress.
    cpu_threads : int
        Number of threads to use for compression.

    Returns
    -------
    str
        Path to the compressed ``.gz`` file.
    """
    pigz_cmd = f"pigz -p {cpu_threads} {input_file}"
    _ = run_subprocess_cmd(pigz_cmd, shell_check = True)
    gzip_file = input_file + ".gz"
    return gzip_file


# --------------------------------------------------------------
# Decompress a file using pigz
# --------------------------------------------------------------
def pigz_decompress(input_file, cpu_threads):
    """Decompress a file using pigz with multiple threads.

    Parameters
    ----------
    input_file : str
        Path to the ``.gz`` file to decompress.
    cpu_threads : int
        Number of threads to use for decompression.

    Returns
    -------
    str
        Path to the decompressed file (with ``.gz`` extension removed).
    """
    pigz_cmd = f"pigz -p {cpu_threads} -d -f {input_file}"
    _ = run_subprocess_cmd(pigz_cmd, shell_check = True)
    unzip_file = input_file.replace(".gz","")
    return unzip_file


# --------------------------------------------------------------
# Extract and process sample metadata
# --------------------------------------------------------------
def get_current_row_data(input_df, sample_id):
    """Extract row data for a sample ID and generate a statistics dictionary.

    Filters *input_df* to the row matching *sample_id*, replaces any
    literal ``"None"`` strings with ``pd.NA`` so that downstream
    ``pd.isna()`` / ``pd.notna()`` guards work correctly, and builds
    the initial sample statistics dictionary.

    Parameters
    ----------
    input_df : pandas.DataFrame
        Full metadata DataFrame loaded from the input CSV.
    sample_id : str
        Sample identifier to filter on the ``SAMPLE_ID`` column.

    Returns
    -------
    tuple of (pandas.DataFrame, list, dict)
        ``(current_row, current_index, sample_stats_dict)`` — the filtered
        single-row DataFrame, its integer index list, and the initialized
        statistics dictionary.
    """
    # Filter the DataFrame for rows where the "SAMPLE_ID" column equals the provided sample_id
    current_row = input_df[input_df["SAMPLE_ID"] == sample_id].copy()

    # Replace literal string "None" with actual NaN so downstream pd.isna()
    # checks work correctly (some CSV editors write "None" instead of leaving
    # cells empty).
    current_row = current_row.replace(to_replace="None", value=pd.NA)

    sample_stats_dict = gen_sample_stats_dict(current_row)
    current_index = current_row.index.tolist()

    return current_row, current_index, sample_stats_dict


# --------------------------------------------------------------
# Parse NanoPlot statistics for long reads
# --------------------------------------------------------------
def analyze_nanostats(READS_ORIGIN, nanoplot_out_file, sample_stats_dict):
    """Parse NanoPlot statistics and update the sample statistics dictionary.

    Reads NanoPlot output and extracts metrics based on read origin.

    Parameters
    ----------
    READS_ORIGIN : str
        Type and stage of reads (e.g., ``'Raw_ONT'``, ``'Filt_PacBio'``).
        Used to determine which keys in *sample_stats_dict* to populate.
    nanoplot_out_file : str
        Path to the NanoPlot ``NanoStats.txt`` output file.
    sample_stats_dict : dict
        Statistics dictionary to update in-place.

    Returns
    -------
    dict
        The updated *sample_stats_dict* with NanoPlot metrics filled in.
    """
    with open(nanoplot_out_file, "r") as nanostats:
        if "raw" in READS_ORIGIN.lower() and "ont" in READS_ORIGIN.lower():
            for line in nanostats:
                if "Number of reads:" in line:
                    sample_stats_dict["RAW_ONT_READS"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Mean read length:" in line:
                    sample_stats_dict["RAW_ONT_MEAN_LENGTH"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Mean read quality:" in line:
                    sample_stats_dict["RAW_ONT_MEAN_QUAL"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Total bases:" in line:
                    sample_stats_dict["RAW_ONT_TOTAL_BASES"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
        elif "filt" in READS_ORIGIN.lower() and "ont" in READS_ORIGIN.lower():    
            for line in nanostats:
                if "Number of reads:" in line:
                    sample_stats_dict["FILT_ONT_READS"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Mean read length:" in line:
                    sample_stats_dict["FILT_ONT_MEAN_LENGTH"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Mean read quality:" in line:
                    sample_stats_dict["FILT_ONT_MEAN_QUAL"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Total bases:" in line:
                    sample_stats_dict["FILT_ONT_TOTAL_BASES"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
        elif "cor" in READS_ORIGIN.lower() and "ont" in READS_ORIGIN.lower():  
            for line in nanostats:
                if "Number of reads:" in line:
                    sample_stats_dict["CORRECT_ONT_READS"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Mean read length:" in line:
                    sample_stats_dict["CORRECT_ONT_MEAN_LENGTH"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Mean read quality:" in line:
                    sample_stats_dict["CORRECT_ONT_MEAN_QUAL"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Total bases:" in line:
                    sample_stats_dict["CORRECT_ONT_TOTAL_BASES"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
        elif "raw" in READS_ORIGIN.lower() and "pacbio" in READS_ORIGIN.lower():
            for line in nanostats:
                if "Number of reads:" in line:
                    sample_stats_dict["RAW_PACBIO_READS"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Mean read length:" in line:
                    sample_stats_dict["RAW_PACBIO_MEAN_LENGTH"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Mean read quality:" in line:
                    sample_stats_dict["RAW_PACBIO_MEAN_QUAL"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Total bases:" in line:
                    sample_stats_dict["RAW_PACBIO_TOTAL_BASES"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
        elif "filt" in READS_ORIGIN.lower() and "pacbio" in READS_ORIGIN.lower():
            for line in nanostats:
                if "Number of reads:" in line:
                    sample_stats_dict["FILT_PACBIO_READS"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Mean read length:" in line:
                    sample_stats_dict["FILT_PACBIO_MEAN_LENGTH"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Mean read quality:" in line:
                    sample_stats_dict["FILT_PACBIO_MEAN_QUAL"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Total bases:" in line:
                    sample_stats_dict["FILT_PACBIO_TOTAL_BASES"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
    return sample_stats_dict


# --------------------------------------------------------------
# Move a file up the directory tree
# --------------------------------------------------------------
def move_file_up(input_file, up_count):
    """Move a file up the directory hierarchy by a specified number of levels.

    Parameters
    ----------
    input_file : str
        Path to the file to move.
    up_count : int
        Number of directory levels to ascend.

    Returns
    -------
    str
        New path to the moved file, or the original *input_file* path if
        the file does not exist.
    """
    if os.path.exists(input_file):
        move_dir = "/".join(os.path.dirname(input_file).split("/")[:-int(up_count)])
        os.makedirs(move_dir, exist_ok=True)
        new_path = os.path.join(move_dir, os.path.basename(input_file))
        shutil.move(input_file, new_path)
        return new_path
    else:
        print(f"ERROR:\tCannot move a non-existing file: {input_file}")
        return input_file


# --------------------------------------------------------------
# Select highest quality long reads
# --------------------------------------------------------------
def select_long_reads(output_dir, input_csv, sample_id, cpu_threads):
    """Select the highest-mean-quality long reads from ONT or PacBio data.

    Parses NanoPlot statistics to compare quality across raw, filtered, and
    corrected read sets, then copies the best set to a canonical
    ``*_highest_mean_qual_long_reads.fastq`` path.

    Parameters
    ----------
    output_dir : str
        Root output directory containing per-species subdirectories.
    input_csv : str
        Path to the metadata CSV file.
    sample_id : str
        Sample identifier used to look up the row in *input_csv*.
    cpu_threads : int
        Number of threads available for compression tasks.

    Returns
    -------
    str or None
        Absolute path to the selected highest-quality reads file, or
        ``None`` if no suitable file can be found.
    """
    input_df = pd.read_csv(input_csv)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]  # Convert to Series (single row)

    # Identify read paths, reference, and BUSCO lineage info from CSV
    ont_sra = current_series["ONT_SRA"]
    ont_raw_reads = current_series["ONT_RAW_READS"]
    pacbio_raw_reads = current_series["PACBIO_RAW_READS"]
    pacbio_sra = current_series["PACBIO_SRA"]
    species_id = current_series["SPECIES_ID"]

    print(f"DEBUG - species_id - {species_id}")
    print(f"DEBUG - sample_id - {sample_id}")
    
    if pd.notna(ont_sra) and pd.isna(ont_raw_reads):
        ont_raw_reads = os.path.join(output_dir, species_id, "ONT", f"{ont_sra}.fastq")
    if pd.notna(pacbio_sra) and pd.isna(pacbio_raw_reads):
        pacbio_raw_reads = os.path.join(output_dir, species_id, "PacBio", f"{pacbio_sra}.fastq")
 
    print(f"DEBUG - ont_raw_reads - {ont_raw_reads}")
    print(f"DEBUG - pacbio_raw_reads - {pacbio_raw_reads}")
    
    if pd.notna(ont_raw_reads):
        print("DEBUG - PROCESSING ONT HIGHEST MEAN QUAL")
        reads_type = "ONT"
        reads_dir = os.path.dirname(ont_raw_reads)
        filtered_reads = os.path.join(reads_dir, f"{species_id}_{reads_type}_filtered.fastq")
        corrected_reads = os.path.join(reads_dir, f"{species_id}_{reads_type}_corrected.fastq")
        reads_origin_list = ["Raw_ONT_", "Filt_ONT_", "Corr_ONT_"]
    elif pd.notna(pacbio_raw_reads):
        print("DEBUG - PROCESSING PACBIO HIGHEST MEAN QUAL")
        reads_type = "PacBio"
        reads_dir = os.path.dirname(pacbio_raw_reads)
        filtered_reads = os.path.join(reads_dir, f"{species_id}_{reads_type}_filtered.fastq")
        corrected_reads = os.path.join(reads_dir, f"{species_id}_{reads_type}_corrected.fastq")
        reads_origin_list = ["Raw_PacBio_", "Filt_PacBio_"]
    else:
        print(f"ERROR:\tUNABLE TO PARSE LONG READS AS BOTH ONT AND PACBIO RAW READS ARE NONE: {ont_raw_reads} & {pacbio_raw_reads}")
    
    for reads_origin in reads_origin_list:
        sample_stats_dict = analyze_nanostats(reads_origin, os.path.join(reads_dir, f"{reads_origin}nanoplot_analysis", f"{reads_origin}NanoStats.txt"), sample_stats_dict)    

    if pd.notna(ont_raw_reads):
        print("Selecting Highest Mean Quality Long reads...")
        highest_mean_qual_long_reads = corrected_reads
        if pd.notna(ont_raw_reads) and sample_stats_dict["CORRECT_ONT_MEAN_QUAL"] < sample_stats_dict["FILT_ONT_MEAN_QUAL"]:
            highest_mean_qual_long_reads = filtered_reads
            highest_mean_qual = sample_stats_dict["FILT_ONT_MEAN_QUAL"]
        else:
            highest_mean_qual = sample_stats_dict["CORRECT_ONT_MEAN_QUAL"]
    if pd.notna(pacbio_raw_reads):
        print("Selecting Highest Mean Quality Long reads...")
        highest_mean_qual_long_reads = pacbio_raw_reads
        if pd.notna(pacbio_raw_reads) and sample_stats_dict["RAW_PACBIO_MEAN_QUAL"] < sample_stats_dict["FILT_PACBIO_MEAN_QUAL"]:
            highest_mean_qual_long_reads = filtered_reads
            highest_mean_qual = sample_stats_dict["FILT_PACBIO_MEAN_QUAL"]
        else:
            highest_mean_qual = sample_stats_dict["RAW_ONT_MEAN_QUAL"]
    print(f"Highest Mean Quality Long reads: {highest_mean_qual_long_reads}")
    print(f"Mean Quality: {highest_mean_qual}")

    renamed_highest_mean_qual_long_reads = f"{species_id}_{reads_type}_highest_mean_qual_long_reads.fastq"
    if not os.path.exists(highest_mean_qual_long_reads):
        # try fallback: see if it's named like "Escherichia_coli_filtered.fastq"
        fallback_file = os.path.join(reads_dir, f"{species_id}_filtered.fastq")
        if os.path.exists(fallback_file):
            print(f"FALLBACK:\tFound fallback filtered file: {fallback_file}")
            highest_mean_qual_long_reads = fallback_file
        else:
            print("ERROR:\tNo usable highest-mean-quality long read file found.")
            return None
        
    renamed_highest_mean_qual_long_reads = os.path.join(reads_dir, f"{species_id}_{reads_type}_highest_mean_qual_long_reads.fastq")
    shutil.copy(highest_mean_qual_long_reads, renamed_highest_mean_qual_long_reads)

    print(f"NOTE:\tSelected highest quality long reads: {renamed_highest_mean_qual_long_reads} with mean quality {highest_mean_qual}")

    return renamed_highest_mean_qual_long_reads


# --------------------------------------------------------------
# Create or manage a log file
# --------------------------------------------------------------
def generate_log_file(log_file_path, use_numerical_suffix=False):
    """Generate a log file, optionally with a numerical suffix if it exists.

    Parameters
    ----------
    log_file_path : str
        Desired path for the log file.
    use_numerical_suffix : bool
        If ``True``, append an incrementing integer suffix rather than
        overwriting an existing file.

    Returns
    -------
    str
        Path to the created or selected log file.
    """
    if os.path.exists(log_file_path) and use_numerical_suffix:
        counter = 1
        base, ext = os.path.splitext(log_file_path)
        new_log_file_path = f"{base}_{counter}{ext}"
        while os.path.exists(new_log_file_path):
            counter += 1
            new_log_file_path = f"{base}_{counter}{ext}"
        log_file_path = new_log_file_path
    else:
        open(log_file_path, "w").close()
    return log_file_path


# --------------------------------------------------------------
# Log and print messages with color
# --------------------------------------------------------------
def log_print(input_message, log_file=None):
    """Log a message to a file and print it with ANSI-colored output.

    Prepends a timestamp, writes the message to the log file, and prints
    it in a color determined by the message prefix (e.g. ``ERROR`` →
    red, ``PASS`` → green).

    Parameters
    ----------
    input_message : str
        Message text to log and print.
    log_file : str, optional
        Path to the log file.  Defaults to the module-level
        ``DEFAULT_LOG_FILE`` set by ``initialize_logging_environment``.
    """
    global DEFAULT_LOG_FILE
    COLORS = {"grey": "\033[90m",
              "red": "\033[91m",
              "green": "\033[92m",
              "orange": "\033[38;5;208m",
              "yellow": "\033[93m",
              "blue": "\033[94m",
              "magenta": "\033[95m",
              "cyan": "\033[96m",
              "white": "\033[97m",
              "reset": "\033[0m"}
    if log_file is None:
        log_file = DEFAULT_LOG_FILE
    now = datetime.now()
    message = f"[{now:%Y-%m-%d %H:%M:%S}]\t{input_message}"
    message_type_dict = {"NOTE": "blue",
                         "CMD": "cyan",
                         "ERROR": "red",
                         "WARN": "yellow",
                         "PASS": "green",
                         "SKIP": "magenta",
                         "FAIL": "red"}
    print_color = "white"
    for key, value in message_type_dict.items():
        if key.lower() in input_message.lower():
            print_color = value
            break
    try:
        with open(log_file, "a") as file:
            print(message, file=file)
    except TypeError:
        print(f"UNLOGGED ERROR:\tUnable to load the log file provided: {log_file}")
    color_code = COLORS.get(print_color, COLORS["white"])
    print(f"{color_code}{message}{COLORS['reset']}")


# --------------------------------------------------------------
# Set up logging environment
# --------------------------------------------------------------
def initialize_logging_environment(INPUT_FOLDER, sample_id=None):
    """Initialize the logging environment based on the input folder.

    Sets the module-level ``DEFAULT_LOG_FILE`` and ``ENVIRONMENT_TYPE``
    globals and creates the log file.  The log file path is adjusted for
    WSL/Linux (``/mnt/<drive>/...``) when running on a non-Windows OS.

    Parameters
    ----------
    INPUT_FOLDER : str
        Output folder path used to determine the log file location.
        When *sample_id* is ``None``, the log file is written as
        ``<INPUT_FOLDER>/<basename>_log.txt``.
    sample_id : str, optional
        When provided, the log file is written per-sample as
        ``<INPUT_FOLDER>/<sample_id>_log.txt`` so that each sample in
        a multi-sample CSV run gets its own log file.
    """
    global DEFAULT_LOG_FILE, ENVIRONMENT_TYPE
    print(INPUT_FOLDER)
    if sample_id:
        input_file_path = f"{INPUT_FOLDER}/{sample_id}_log.txt"
    else:
        input_file_path = f"{INPUT_FOLDER}/{INPUT_FOLDER.split('/')[-1]}_log.txt"
    os_name = platform.system()
    if os_name == "Windows":
        print("UNLOGGED:\tWINDOWS ENVIRONMENT")
        ENVIRONMENT_TYPE = "WIN"
    elif os_name in ["Linux", "Darwin"]:
        drive, path_without_drive = os.path.splitdrive(input_file_path)
        if drive:
            drive_letter = drive.strip(":\\/")
            path_without_drive_mod = path_without_drive.replace("\\", "/")
            input_file_path = f"/mnt/{drive_letter.lower()}{path_without_drive_mod}"
        print("UNLOGGED:\tLINUX/WSL/MAC ENVIRONMENT")
        ENVIRONMENT_TYPE = "LINUX/WSL/MAC"
    else:
        print(f"UNLOGGED ERROR:\tUnsupported OS: {os_name}")
        return
    print(input_file_path)
    run_log = generate_log_file(input_file_path, use_numerical_suffix=False)
    DEFAULT_LOG_FILE = run_log
