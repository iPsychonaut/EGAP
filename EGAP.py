#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 19:57:46 2024

@author: ian.bollinger@entheome.org

EGAP (Entheome Genome Assembly Pipeline) is a versatile bioinformatics pipeline
for hybrid genome assembly from Oxford Nanopore, Illumina, and PacBio data.
It supports multiple input modes and assembly methods and determines the best 
based on multiple metrics: BUSCO Completeness (Single + Duplicated), Assembly 
Contig Count, Assembly N50, Assembly L50, and Assembly GC-content.
"""
# Library Imports
import math, platform, os, subprocess, multiprocessing, argparse, psutil, shutil, hashlib, re, gzip, glob
import numpy as np
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from bs4 import BeautifulSoup


# Global Variables
CPU_THREADS = 1  # Number of CPU threads to use.
RAM_GB = 1  # Amount of RAM in GB to use.
DEFAULT_LOG_FILE = None  # Default log file path.
ENVIRONMENT_TYPE = None  # Type of environment (e.g., "WIN", "LINUX/WSL/MAC").


def generate_log_file(log_file_path, use_numerical_suffix=False):
    """
    Generate a log file for the pipeline.

    If a file already exists and use_numerical_suffix is True, a new file with a numerical
    suffix is created; otherwise, the existing file is overwritten.

    Parameters:
        log_file_path (str): Path to the log file.
        use_numerical_suffix (bool): Whether to create a new file with a numerical suffix.

    Returns:
        log_file_path (str): The path to the log file.
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


def log_print(input_message, log_file=None):
    """
    Log a message to a file and print it with colored output.

    The message is timestamped and written to the specified log file (or a default one)
    and then printed in a color chosen based on the message type.

    Parameters:
        input_message (str): The message to log and print.
        log_file (str, optional): Path to the log file; if None, a default log file is used.
    """
    global DEFAULT_LOG_FILE
    COLORS = {"grey": "\033[90m",
              "red": "\033[91m",
              "green": "\033[92m",
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


def initialize_logging_environment(INPUT_FOLDER):
    """
    Initialize the logging environment based on the input folder.

    Sets the global DEFAULT_LOG_FILE and ENVIRONMENT_TYPE by creating a log file
    in a location derived from INPUT_FOLDER and the operating system.

    Parameters:
        INPUT_FOLDER (str): The folder used to influence log file generation.
    """
    global DEFAULT_LOG_FILE, ENVIRONMENT_TYPE
    print(INPUT_FOLDER)
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


def run_subprocess_cmd(cmd_list, shell_check):
    """
    Run a subprocess command and log its execution.

    Executes the command (as a string or list) using subprocess.Popen, streams its output,
    and logs whether it succeeded or failed.

    Parameters:
        cmd_list (str or list): The command to execute.
        shell_check (bool): Whether to execute the command through the shell.

    Returns:
        process.returncode (int): The return code of the subprocess.
    """
    if isinstance(cmd_list, str):
        log_print(f"CMD:\t{cmd_list}")
        process = subprocess.Popen(cmd_list, shell=shell_check, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True)
    else:
        log_print(f"CMD:\t{' '.join(cmd_list)}")
        process = subprocess.Popen(cmd_list, shell=shell_check, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True)
    for line in process.stdout:
        print(line, end="")
    process.wait()
    if process.returncode != 0:
        log_print(f"NOTE:\tCommand failed with return code {process.returncode}")
    else:
        log_print(f"PASS:\tSuccessfully processed command: {' '.join(cmd_list) if isinstance(cmd_list, list) else cmd_list}")
    return process.returncode


def get_resource_values(PERCENT_RESOURCES):
    """
    Calculate the number of CPU threads and GB of RAM to use based on a percentage.

    Parameters:
        PERCENT_RESOURCES (float): Fraction of total resources to use.

    Returns:
        cpu_threads (int): Number of threads avaialable for use.
        ram_gb (int): Amount of RAM avaialable for use.
    """
    num_cpus = multiprocessing.cpu_count()
    mem_info = psutil.virtual_memory()
    cpu_threads = int(math.floor(num_cpus * PERCENT_RESOURCES))
    ram_gb = int(mem_info.total / (1024.0 ** 3) * PERCENT_RESOURCES)
    return cpu_threads, ram_gb 


def find_file(filename):
    """
    Search the filesystem for a file with the given name.

    Walks the filesystem from the root (based on the OS) and returns the absolute path
    of the first occurrence of the file.

    Parameters:
        filename (str): The file name to search for.

    Returns:
        str or None: The absolute path if found; otherwise, None.
    """
    global ENVIRONMENT_TYPE
    log_print(f"Looking for {filename}")
    if ENVIRONMENT_TYPE == "WIN":
        root_directory = "C:\\"
    elif ENVIRONMENT_TYPE in ["LINUX/WSL/MAC"]:
        root_directory = "/"
    else:
        raise ValueError("Unknown ENVIRONMENT_TYPE")
    for root, dirs, files in os.walk(root_directory):
        if filename in files:
            return os.path.join(root, filename)
    return None


def find_ca_folder(input_folder):
    """
    Find the MaSuRCA CA folder within the given directory.

    Scans the subdirectories of input_folder for one whose name starts with "CA".
    If none is found, defaults to input_folder/CA.

    Parameters:
        input_folder (str): The directory to search.

    Returns:
        ca_folder (str): The path to the CA folder.
    """
    subfolders = [f.path for f in os.scandir(input_folder) if f.is_dir()]
    ca_folder = f"{input_folder}/CA"
    for folder in subfolders:
        if os.path.basename(folder).startswith("CA"):
            ca_folder = folder
            break
    return ca_folder


def md5_check(illumina_raw_data_dir, illumina_df):
    """
    Verify MD5 checksums for .fq.gz files against MD5.txt.

    Reads the MD5.txt file from illumina_raw_data_dir, updates the DataFrame with the
    expected checksums, computes the MD5 of each .gz file, and logs any mismatches.

    Parameters:
        illumina_raw_data_dir (str): Directory with Illumina raw reads and MD5.txt.
        illumina_df (DataFrame): DataFrame to update with checksum data.

    """
    md5_file_path = os.path.join(illumina_raw_data_dir, "MD5.txt")
    if os.path.exists(md5_file_path):
        with open(md5_file_path, "r") as f:
            md5_data = f.read().splitlines()
        md5_dict = dict(line.split() for line in md5_data)
        md5_df = pd.DataFrame(md5_dict.items(), columns=["MD5", "Filename"])
        illumina_df = pd.concat([illumina_df, md5_df], ignore_index=True)
        for file_name in os.listdir(illumina_raw_data_dir):
            if file_name.endswith(".gz"):
                gz_file_path = os.path.join(illumina_raw_data_dir, file_name)
                matching_row = illumina_df[illumina_df["Filename"].str.contains(file_name)].head(1)
                original_md5 = matching_row["MD5"].iloc[0] if not matching_row.empty else "None"
                if original_md5 != "None":
                    hash_md5 = hashlib.md5()
                    with open(gz_file_path, "rb") as f:
                        for chunk in iter(lambda: f.read(4096), b""):
                            hash_md5.update(chunk)
                    new_md5 = hash_md5.hexdigest()
                    if original_md5 != new_md5:
                        log_print(f"ERROR:\tMD5 checksum mismatch for {gz_file_path}: original {original_md5}, new {new_md5}")
                        break
                    else:
                        log_print(f"PASS:\tOriginal MD5 MATCHES for {os.path.basename(file_name)}")
                else:
                    log_print(f"ERROR:\tOriginal MD5 checksum not found for {os.path.basename(file_name)}")
                    break


def illumina_extract_and_check(folder_name, SPECIES_ID):
    """
    Combine paired-end Illumina reads after MD5 verification.

    Calls md5_check on folder_name, concatenates forward and reverse FASTQ files
    if not already present, gzips the combined files, and returns their paths.

    Parameters:
        folder_name (str): Directory with Illumina .fq.gz files and MD5.txt.
        SPECIES_ID (str): Identifier used to name the combined files.

    Returns:
        list: [combined_forward_file.gz, combined_reverse_file.gz]
    """
    log_print(f"Running MD5 Checksum Analysis on Raw Illumina FASTQ.GZ files in {folder_name}...")
    illumina_df = pd.DataFrame(columns=["MD5", "Filename"])
    base_data_dir = os.path.dirname(folder_name)
    combined_1_file = os.path.join(base_data_dir, f"{SPECIES_ID}_combined_1.fq")
    combined_2_file = os.path.join(base_data_dir, f"{SPECIES_ID}_combined_2.fq")
    combined_list = [combined_1_file, combined_2_file]
    combined_gz_list = [f"{combined_fq}.gz" for combined_fq in combined_list]
    if not os.path.isfile(combined_gz_list[0]) or not os.path.isfile(combined_gz_list[0]):
        if not os.path.isfile(combined_1_file) or not os.path.isfile(combined_2_file):
            md5_check(folder_name, illumina_df)
            raw_1_list = []
            raw_2_list = []
            for filename in os.listdir(folder_name):
                if "_1.fq.gz" in filename:
                    raw_1_list.append(os.path.join(folder_name, filename))
                elif "_2.fq.gz" in filename:
                    raw_2_list.append(os.path.join(folder_name, filename))
            fwd_cat_cmd = f"cat {raw_1_list[0]} {raw_1_list[1]} > {combined_1_file}"
            _ = run_subprocess_cmd(fwd_cat_cmd, shell_check=True)
            rev_cat_cmd = f"cat {raw_2_list[0]} {raw_2_list[1]} > {combined_2_file}"
            _ = run_subprocess_cmd(rev_cat_cmd, shell_check=True)
        else:
            log_print(f"SKIP:\tCombined FASTQ files already exist: {combined_list[0]}; {combined_list[1]}.")
        for combined_fq in combined_list:
            gzip_cmd = ["gzip", combined_fq]
            _ = run_subprocess_cmd(gzip_cmd, shell_check=False)
    else:
        log_print(f"SKIP:\tGzipped Combined FASTQ files already exist: {combined_gz_list[0]}; {combined_gz_list[1]}.")
    return combined_gz_list


def classify_metric(value, thresholds):
    """
    Classify a metric value based on preset thresholds.

    Parameters:
        value (float): The metric value.
        thresholds (dict): Thresholds with keys 'AMAZING', 'GREAT', 'OK', 'POOR'.

    Returns:
        str: One of "AMAZING", "GREAT", "OK", or "POOR".
    """
    if value >= thresholds["AMAZING"]:
        return "AMAZING"
    elif value >= thresholds["GREAT"]:
        return "GREAT"
    elif value >= thresholds["OK"]:
        return "OK"
    elif value >= thresholds["POOR"]:
        return "POOR"
    else:
        if value <= thresholds["AMAZING"]:
            return "AMAZING"
        elif value <= thresholds["GREAT"]:
            return "GREAT"
        elif value <= thresholds["OK"]:
            return "OK"
        elif value <= thresholds["POOR"]:
            return "POOR"


def classify_assembly(sample_stats):
    """
    Evaluate assembly quality based on multiple metrics.

    Classifies completeness, N50, and contig count from sample_stats and returns both
    individual and overall quality ratings.

    Parameters:
        sample_stats (dict): Assembly metrics (e.g., BUSCO scores, N50, contig count).

    Returns:
        results (dict): Quality classifications per metric and an overall rating.
    """
    results = {}
    if sample_stats["FIRST_COMPLEASM_C"] < 80.0:
        results["FIRST_COMPLEASM_C"] = "POOR"
    elif sample_stats["FIRST_COMPLEASM_C"] < 95.0:
        results["FIRST_COMPLEASM_C"] = "OK"
    elif sample_stats["FIRST_COMPLEASM_C"] < 98.5:
        results["FIRST_COMPLEASM_C"] = "GREAT"
    elif sample_stats["FIRST_COMPLEASM_C"] >= 98.5:
        results["FIRST_COMPLEASM_C"] = "AMAZING"
    if sample_stats["SECOND_COMPLEASM_C"] < 80.0:
        results["SECOND_COMPLEASM_C"] = "POOR"
    elif sample_stats["SECOND_COMPLEASM_C"] < 95.0:
        results["SECOND_COMPLEASM_C"] = "OK"
    elif sample_stats["SECOND_COMPLEASM_C"] < 98.5:
        results["SECOND_COMPLEASM_C"] = "GREAT"
    elif sample_stats["SECOND_COMPLEASM_C"] >= 98.5:
        results["SECOND_COMPLEASM_C"] = "AMAZING"
    if sample_stats["ASSEMBLY_CONTIGS"] <= 100:
        results["ASSEMBLY_CONTIGS"] = "AMAZING"
    elif sample_stats["ASSEMBLY_CONTIGS"] <= 1000:
        results["ASSEMBLY_CONTIGS"] = "GREAT"
    elif sample_stats["ASSEMBLY_CONTIGS"] <= 10000:
        results["ASSEMBLY_CONTIGS"] = "OK"
    elif sample_stats["ASSEMBLY_CONTIGS"] > 10000:
        results["ASSEMBLY_CONTIGS"] = "POOR"
    if sample_stats["ASSEMBLY_N50"] <= 100:
        results["ASSEMBLY_N50"] = "POOR"
    elif sample_stats["ASSEMBLY_N50"] <= 1000:
        results["ASSEMBLY_N50"] = "OK"
    elif sample_stats["ASSEMBLY_N50"] <= 10000:
        results["ASSEMBLY_N50"] = "GREAT"
    elif sample_stats["ASSEMBLY_N50"] > 10000:
        results["ASSEMBLY_N50"] = "AMAZING"
    if all(value == "AMAZING" for value in results.values()):
        results["OVERALL"] = "AMAZING"
    elif any(value == "POOR" for value in results.values()):
        results["OVERALL"] = "POOR"
    elif any(value == "OK" for value in results.values()):
        results["OVERALL"] = "OK"
    elif any(value == "GREAT" for value in results.values()):
        results["OVERALL"] = "GREAT"
    return results


def parse_bbmerge_output(insert_size_histogram_txt):
    """
    Parse BBMerge output to extract insert size statistics.

    Reads the histogram file to obtain the average insert size and its standard deviation.

    Parameters:
        insert_size_histogram_txt (str): Path to the BBMerge insert size histogram file.

    Returns:
        avg_insert (float): Average insert size.
        std_dev (float): Standard Deviation.

    Raises:
        ValueError: If the expected lines for mean and standard deviation are not found.
    """
    log_print(f"Processing insert size histogram: {insert_size_histogram_txt}...")
    avg_insert = None
    std_dev = None
    with open(insert_size_histogram_txt, "r") as file:
        for line in file:
            if "#Mean\t" in line:
                avg_insert = round(float(line.replace("#Mean\t", "").replace("\n", "")), 0)
            if "#STDev\t" in line:
                std_dev = round(float(line.replace("#STDev\t", "").replace("\n", "")), 0)
    if avg_insert is None or std_dev is None:
        raise ValueError("Could not find average insert size and/or standard deviation in the output.")
    return avg_insert, std_dev


def bbmap_stats(input_folder, reads_list):
    """
    Compute insert size statistics using BBMerge.

    If an existing histogram is found in input_folder, it is parsed; otherwise, BBMerge
    is run on the provided reads to generate the histogram and then parse it.

    Parameters:
        input_folder (str): Directory for BBMerge outputs.
        reads_list (list): List of FASTQ file paths (either paired-end or with an extra third read).

    Returns:
        avg_insert (float): Average insert size.
        std_dev (float): Standard Deviation.
    """
    bbmap_out_path = f"{input_folder}/bbmap_data.fq.gz"
    insert_size_histogram_txt = f"{input_folder}/insert_size_histogram.txt"
    log_print(f"NOTE:\tCurrent bbmap out path: {bbmap_out_path}")
    default_bbmerge_path = "bbmerge.sh"
    avg_insert = 251
    std_dev = 30
    if os.path.isfile(insert_size_histogram_txt):
        log_print(f"SKIP:\tbbmap insert size histogram output already exists: {insert_size_histogram_txt}")
        avg_insert, std_dev = parse_bbmerge_output(insert_size_histogram_txt)
    else:
        log_print("Processing fastq files for bbmap stats...")
        bbmerge_path = find_file(default_bbmerge_path)
        print(reads_list)
        bbmerge_cmd = [bbmerge_path,
                       f"in1={reads_list[1]}",
                       f"in2={reads_list[2]}",
                       f"out={bbmap_out_path}",
                       f"ihist={input_folder}/insert_size_histogram.txt"]
        _ = run_subprocess_cmd(bbmerge_cmd, False)
        try:
            avg_insert, std_dev = parse_bbmerge_output(insert_size_histogram_txt)
        except:
            log_print("NOTE:\tUnable to parse avg_insert or std_dev, using default values: 251 and 30 respectively.")
            avg_insert = 251
            std_dev = 30
    return avg_insert, std_dev


def skip_gap_closing_section(assembly_sh_path):
    """
    Modify the MaSuRCA assemble.sh script to skip gap closing.

    Replaces the gap-closing block in the script with code that logs a skip message and forces
    the terminator to "9-terminator". Writes the modified script as 'assemble_skip_gap.sh'.

    Parameters:
        assembly_sh_path (str): Path to the original assemble.sh script.
    """
    with open(assembly_sh_path, "r") as f_in:
        original_script = f_in.read()
    pattern = re.compile(r"(if \[ -s \$CA_DIR/9-terminator/genome\.scf\.fasta \];then)"
                         r"(.*?)"
                         r"(else\s+fail 'Assembly stopped or failed, see \$CA_DIR\.log'\nfi)",
                         re.DOTALL)
    replacement_snippet = (r"\1\n"
                           r"  # Force the final terminator to remain '9-terminator'\n"
                           r"  log \"Skipping gap closing step; using 9-terminator as final.\"\n"
                           r"  TERMINATOR='9-terminator'\n"
                           r"\3")
    modified_text = re.sub(pattern, replacement_snippet, original_script)
    with open("assemble_skip_gap.sh", "w") as f_out:
        f_out.write(modified_text)

def masurca_config_gen(input_folder, output_folder, 
                       ONT_RAW_READS, ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS, 
                       PACBIO_RAW_READS, 
                       clump_f_dedup_path, clump_r_dedup_path, 
                       EST_SIZE, CPU_THREADS, ram_gb, ref_seq=None):
    """
    Generate and run a MaSuRCA configuration for genome assembly using only the CABOG assembler.

    This version supports both Illumina-only and hybrid (Illumina + long-read) assemblies.
    It completely ignores the use of SOAP (and Flye) by always forcing SOAP_ASSEMBLY and
    FLYE_ASSEMBLY to 0. The output directory is determined via find_ca_folder(), ensuring
    that only CABOG is used.

    Parameters:
        input_folder (str): Directory containing input data.
        output_folder (str): Directory for MaSuRCA outputs.
        ONT_RAW_READS (str or NaN): Path to Nanopore reads file.
        ILLUMINA_RAW_F_READS (str or NaN): Path to Illumina forward reads.
        ILLUMINA_RAW_R_READS (str or NaN): Path to Illumina reverse reads.
        PACBIO_RAW_READS (str or NaN): Path to PacBio reads file.
        clump_f_dedup_path (str): Path to deduplicated Illumina forward reads.
        clump_r_dedup_path (str): Path to deduplicated Illumina reverse reads.
        EST_SIZE (str): Estimated genome size (e.g. "25M", "3.2G").
        CPU_THREADS (int): Number of CPU threads.
        ram_gb (int): Amount of available RAM in GB.
        ref_seq (str, optional): Reference genome FASTA path.

    Returns:
        tuple: (default_assembly_path, assembly_path, ref_assembly_path)
               where ref_assembly_path is None if ref_seq is not provided.
    """
    import os, shutil
    import pandas as pd  # for pd.notna and pd.isna

    # Unzip deduplicated files if they are gzipped.
    dedup_unzip_list = [clump_f_dedup_path.replace(".gz", ""),
                          clump_r_dedup_path.replace(".gz", "")]
    for input_file in [clump_f_dedup_path, clump_r_dedup_path]:
        if ".gz" in input_file:
            gunzip_file(input_file, input_file.replace(".gz", ""))

    # Calculate average insert size and standard deviation from the input read data.
    avg_insert, std_dev = bbmap_stats(input_folder, 
                                      [ONT_RAW_READS, ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS, PACBIO_RAW_READS])
    
    # Change to the output directory.
    os.chdir(output_folder)
    
    # Set the Jellyfish hash size (adjust if available RAM is lower than 62GB)
    jf_size = EST_SIZE * 20
    max_ram_tested = 62
    if ram_gb < max_ram_tested:
        adjustment_ratio = ram_gb / max_ram_tested
        jf_size = int(round(jf_size * adjustment_ratio, 0))
    
    # Define the expected assembly file paths.
    assembly_path = os.path.join(input_folder, "primary.genome.scf.fasta")
    ref_assembly_path = os.path.join(input_folder, "primary.genome.ref.fasta")
    
    # Determine if this is a hybrid assembly (Illumina paired-end + at least one long-read file).
    hybrid_assembly = False
    if pd.notna(ILLUMINA_RAW_F_READS) and pd.notna(ILLUMINA_RAW_R_READS) and (pd.notna(ONT_RAW_READS) or pd.notna(PACBIO_RAW_READS)):
        hybrid_assembly = True

    # Set assembler-specific flags.
    if pd.notna(ILLUMINA_RAW_F_READS) and pd.notna(ILLUMINA_RAW_R_READS):
        if hybrid_assembly:
            use_linking_mates = 0
            close_gaps = 0
        else:
            use_linking_mates = 1
            close_gaps = 1
        # Force use of CABOG only: disable SOAP and Flye.
        soap_assembly = 0
        flye_assembly = 0
        # Always use the CABOG output folder.
        data_output_folder = find_ca_folder(input_folder)
        default_assembly_path = os.path.join(data_output_folder, "primary.genome.scf.fasta")
    else:
        data_output_folder = output_folder
        default_assembly_path = os.path.join(data_output_folder, "unknown_assembly.fasta")

    # If the default assembly already exists, skip running MaSuRCA.
    if os.path.exists(default_assembly_path):
        log_print("PASS:\tSkipping MaSuRCA, moving output files")
        if ref_seq:
            shutil.move(default_assembly_path, ref_assembly_path)
        else:
            shutil.move(default_assembly_path, assembly_path)
    elif os.path.exists(assembly_path):
        log_print(f"SKIP:\tMaSuRCA Assembly, scaffolded assembly already exists: {assembly_path}.")
    else:
        # Build the configuration file.
        config_content = ["DATA\n"]
        
        # Add the Illumina paired-end reads line.
        if pd.notna(ILLUMINA_RAW_F_READS) and pd.notna(ILLUMINA_RAW_R_READS):
            config_content.append(
                f"PE= pe {int(avg_insert)} {int(std_dev)} "
                f"{clump_f_dedup_path.replace('.gz','')} {clump_r_dedup_path.replace('.gz','')}\n"
            )
        # For hybrid assemblies, add the long-read data line.
        if pd.notna(ONT_RAW_READS) and pd.isna(PACBIO_RAW_READS):
            config_content.append(f"NANOPORE={ONT_RAW_READS.replace('.gz','')}\n")
        elif pd.notna(PACBIO_RAW_READS) and pd.isna(ONT_RAW_READS):
            config_content.append(f"PACBIO={PACBIO_RAW_READS.replace('.gz','')}\n")
        elif pd.notna(ONT_RAW_READS) and pd.notna(PACBIO_RAW_READS):
            # If both long-read types are provided, prefer the PACBIO line.
            config_content.append(f"PACBIO={PACBIO_RAW_READS.replace('.gz','')}\n")
        
        # Optionally add a REFERENCE line.
        if ref_seq:
            config_content.append(f"REFERENCE={ref_seq.replace('.gbff','.fna.gz')}\n")
        config_content.append("END\n")
        
        # Add the PARAMETERS section.
        config_content.append("PARAMETERS\n")
        config_content.append("GRAPH_KMER_SIZE=auto\n")
        config_content.append(f"USE_LINKING_MATES={use_linking_mates}\n")
        config_content.append(f"CLOSE_GAPS={close_gaps}\n")
        config_content.append("MEGA_READS_ONE_PASS=0\n")
        config_content.append("LIMIT_JUMP_COVERAGE=300\n")
        config_content.append("CA_PARAMETERS=cgwErrorRate=0.15\n")
        config_content.append(f"NUM_THREADS={CPU_THREADS}\n")
        config_content.append(f"JF_SIZE={jf_size}\n")
        # Force CABOG: always disable SOAP and Flye.
        config_content.append(f"SOAP_ASSEMBLY={soap_assembly}\n")
        config_content.append(f"FLYE_ASSEMBLY={flye_assembly}\n")
        config_content.append("END\n")
        
        # Write the configuration file.
        config_path = os.path.join(output_folder, "masurca_config_file.txt")
        with open(config_path, "w") as file:
            for entry in config_content:
                file.write(entry)
        
        # Run the MaSuRCA configuration command.
        masurca_config_cmd = ["masurca", "masurca_config_file.txt"]
        _ = run_subprocess_cmd(masurca_config_cmd, False)
    
    # Modify the assemble.sh to skip gap closing if needed.
    assemble_sh_path = os.path.join(output_folder, "assemble.sh")
    skip_gap_closing_section(assemble_sh_path)
    
    # Run the assembly.
    masurca_assemble_cmd = ["bash", assemble_sh_path]
    return_code = run_subprocess_cmd(masurca_assemble_cmd, False)
    
    # Clean up any unzipped deduplication files.
    for output_file in dedup_unzip_list:
        try:
            os.remove(output_file)
            log_print(f"Removed: {output_file}")
        except FileNotFoundError:
            log_print(f"File not found: {output_file}")
        except PermissionError:
            log_print(f"Permission denied: {output_file}")
        except Exception as e:
            log_print(f"Error removing {output_file}: {e}")
    
    # If MaSuRCA assembly was interrupted intentionally.
    if return_code == 1:
        log_print("NOTE:\tMaSuRCA assembly interrupted intentionally; Gap Closing will be performed later.")
        default_assembly_path = os.path.join(data_output_folder, os.path.basename(default_assembly_path))
        return default_assembly_path, assembly_path, ref_assembly_path if ref_seq else None

    return default_assembly_path, assembly_path, ref_assembly_path if ref_seq else None


def bbduk_map(trimmo_f_pair_path, trimmo_r_pair_path):
    """
    Run BBDuk to trim and remove adapters from paired-end FASTQ files.

    Produces two new FASTQ files with "_mapped" added to the file names.

    Parameters:
        trimmo_f_pair_path (str): Path to the forward paired FASTQ file.
        trimmo_r_pair_path (str): Path to the reverse paired FASTQ file.

    Returns:
        forward_mapped_path (str): Path to the forward mapped FASTQ file.
        reverse_mapped_path (str): Path to the reverse mapped FASTQ file.
    """
    file_extension = trimmo_f_pair_path.split("_1_paired.")[1]
    bbduk_f_map_path = trimmo_f_pair_path.replace(f"_1_paired.{file_extension}", f"_forward_mapped.{file_extension}")
    bbduk_r_map_path = trimmo_r_pair_path.replace(f"_2_paired.{file_extension}", f"_reverse_mapped.{file_extension}")
    if os.path.exists(bbduk_f_map_path) and os.path.exists(bbduk_r_map_path):
        log_print(f"SKIP:\tbbduk Mapped outputs already exist: {bbduk_f_map_path}; {bbduk_r_map_path}.")
    else:
        default_adapters_path = "adapters.fa"
        default_bbduk_path = "bbduk.sh"
        adapters_path = find_file(default_adapters_path)
        bbduk_path = find_file(default_bbduk_path)
        bbduk_cmd = [bbduk_path,
                     f"in1={trimmo_f_pair_path}", f"in2={trimmo_r_pair_path}",
                     f"out1={bbduk_f_map_path}", f"out2={bbduk_r_map_path}",
                     f"ref={adapters_path}", "ktrim=r", "k=23", "mink=11", "hdist=1",
                     "tpe", "tbo", "qtrim=rl", "trimq=20"]
        _ = run_subprocess_cmd(bbduk_cmd, False)
    return bbduk_f_map_path, bbduk_r_map_path


def clumpify_dedup(bbduk_f_map_path, bbduk_r_map_path):
    """
    De-duplicate paired-end FASTQ files using Clumpify.

    If the deduplicated outputs already exist, the command is skipped; otherwise, Clumpify is run.

    Parameters:
        bbduk_f_map_path (str): Path to the forward mapped FASTQ file.
        bbduk_r_map_path (str): Path to the reverse mapped FASTQ file.

    Returns:
        dedup_forward_path (str): Path to the forward deduplicated and clumped FASTQ file.
        dedup_reverse_path (str): Path to the reverse deduplicated and clumped FASTQ file.
    """
    file_extension = bbduk_f_map_path.split("_forward_mapped.")[1]
    clump_f_dedup_path = bbduk_f_map_path.replace(f"_forward_mapped.{file_extension}", f"_forward_dedup.{file_extension}")
    clump_r_dedup_path = bbduk_r_map_path.replace(f"_reverse_mapped.{file_extension}", f"_reverse_dedup.{file_extension}")
    if os.path.exists(clump_f_dedup_path) and os.path.exists(clump_r_dedup_path):
        log_print(f"SKIP:\tclumpify deduplicated outputs already exist: {clump_f_dedup_path}.")
    else:
        default_clumpify_path = "clumpify.sh"
        clumpify_path = find_file(default_clumpify_path)
        clumpify_cmd = [clumpify_path,
                        f"in={bbduk_f_map_path}",
                        f"in2={bbduk_r_map_path}",
                        f"out={clump_f_dedup_path}",
                        f"out2={clump_r_dedup_path}",
                        "dedupe"]
        _ = run_subprocess_cmd(clumpify_cmd, False)
    return clump_f_dedup_path, clump_r_dedup_path


def get_total_bases(html_file):
    """
    Extract the 'Total Bases' value from an HTML report.

    Parses the HTML file to locate the table cell with "Total Bases" and returns
    the text from the adjacent cell.

    Parameters:
        html_file (str): Path to the HTML file.

    Returns:
        str or None: The total bases value (e.g., "4.8 Gbp") or None if not found.
    """
    with open(html_file, "r", encoding="utf-8") as f:
        contents = f.read()
    soup = BeautifulSoup(contents, "html.parser")
    total_bases_td = soup.find("td", string="Total Bases")
    if total_bases_td:
        value_td = total_bases_td.find_next_sibling("td")
        if value_td:
            return value_td.get_text(strip=True)
    return None


def process_read_file(read_path):
    """
    Ensure an Illumina read file uses the '.fastq.gz' extension.

    Renames or compresses the file as needed:
      - Renames '.fq.gz' to '.fastq.gz'
      - Gzips '.fq' or '.fastq' files

    Parameters:
        read_path (str): Path to the read file.

    Returns:
        read_path (str): The updated file path.
    """
    if not isinstance(read_path, str):
        log_print(f"Read path is not a string: {read_path}")
        return read_path

    dir_path, filename = os.path.split(read_path)
    basename_list = filename.split(".")
    basename = basename_list[0]
    new_read_path = read_path

    if read_path.endswith(".fastq.gz"):
        log_print(f"Read file already in .fastq.gz format: {read_path}")
        return read_path
    elif read_path.endswith(".fq"):
        new_read_path = os.path.join(dir_path, basename + ".fastq.gz")
        gzip_file(read_path, new_read_path)
        os.remove(read_path)
        log_print(f"Gzipped .fq to {new_read_path}")
        return new_read_path
    elif read_path.endswith(".fq.gz"):
        new_read_path = os.path.join(dir_path, basename + ".fastq.gz")
        os.rename(read_path, new_read_path)
        log_print(f"Renamed .fq.gz to .fastq.gz: {new_read_path}")
        return new_read_path
    elif read_path.endswith(".fastq"):
        new_read_path = os.path.join(dir_path, basename + ".fastq.gz")
        gzip_file(read_path, new_read_path)
        os.remove(read_path)
        log_print(f"Renamed .fastq to .fq and gzipped to {new_read_path}")
        return new_read_path
    else:
        log_print(f"Unrecognized file extension for read file: {read_path}")
        return read_path


def gzip_file(input_file, output_file):
    """
    Compress a file using gzip.

    Parameters:
        input_file (str): Path to the file to compress.
        output_file (str): Destination gzip file path.
    """
    with open(input_file, "rb") as f_in:
        with gzip.open(output_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def gunzip_file(input_file, output_file):
    """
    Decompress a gzip file.

    Parameters:
        input_file (str): Path to the gzip file.
        output_file (str): Destination file path for decompressed data.
    """
    with gzip.open(input_file, "rb") as f_in:
        with open(output_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def ont_combine_fastq_gz(ONT_FOLDER):
    """
    Combine multiple ONT FASTQ.GZ files into a single FASTQ file.

    Searches within ONT_FOLDER (or its subdirectories) for FASTQ.GZ files,
    concatenates the records, and writes them to a combined file.

    Parameters:
        ONT_FOLDER (str): Directory containing ONT FASTQ.GZ files.

    Returns:
        combined_ont_fastq_path (str): Path to the combined FASTQ file.
    """
    log_print("Combining ONT FASTQ.GZ files...")
    ont_raw_data_dir = next((subdir for subdir in glob.glob(os.path.join(ONT_FOLDER, "*"))
                            if os.path.isdir(subdir) and any(file.endswith(".gz") for file in os.listdir(subdir))), None)
    if ont_raw_data_dir is None:
        log_print(f"NOTE:\tNo directory containing '.gz' files found within '{ONT_FOLDER}'...Attempting with ONT_FOLDER '{ONT_FOLDER}'")
        ont_raw_data_dir = ONT_FOLDER
    base_name = ONT_FOLDER.split("/")[-3]
    if base_name == "ENTHEOME":
        base_name = ONT_FOLDER.split("/")[-2]
    combined_ont_fastq_path = os.path.join(ONT_FOLDER, f"{base_name}_ont_combined.fastq")
    raw_file_list = glob.glob(os.path.join(ont_raw_data_dir, "*.fastq.gz"))
    if os.path.isfile(combined_ont_fastq_path):
        log_print(f"NOTE:\tSkipping extraction & combination: Combined fastq file already exists: {combined_ont_fastq_path}")
        return combined_ont_fastq_path
    with open(combined_ont_fastq_path, "w") as combined_file:
        for filename in raw_file_list:
            with gzip.open(filename, "rt") as gz_file:
                for record in SeqIO.parse(gz_file, "fastq"):
                    try:
                        SeqIO.write(record, combined_file, "fastq")
                    except Exception as e:
                        log_print(f"ERROR:\tFound in FASTQ record: {e}")
                        raise e
    log_print(f"PASS: Successfully created combined fastq file: {combined_ont_fastq_path}")
    return combined_ont_fastq_path


def plot_busco(compleasm_odb, input_busco_tsv, input_fasta):
    """
    Generate BUSCO status plots from Compleasm output.

    Reads the BUSCO TSV file, prepares a stacked bar plot of BUSCO statuses per sequence,
    and saves the plot in both SVG and PNG formats.

    Parameters:
        compleasm_odb (str): BUSCO lineage identifier.
        input_busco_tsv (str): Path to the BUSCO TSV file.
        input_fasta (str): Path to the input FASTA file (used to name the outputs).
    """
    log_print(f"Generating BUSCO plot for {input_busco_tsv}...")
    busco_df = pd.read_csv(input_busco_tsv, sep="\t", header=0, dtype=str)
    if busco_df.empty:
        log_print("WARNING: BUSCO input file is empty. Skipping plot generation.")
        plt.figure(figsize=(12, 8))
        plt.text(0.5, 0.5, "No valid BUSCO data to plot", fontsize=14, ha='center', va='center')
        plt.xticks([])
        plt.yticks([])
        plt.title("BUSCO Status Plot - No Data Available")
        output_busco_svg = input_fasta.replace(".fasta", f"_{compleasm_odb}_busco.svg")
        plt.savefig(output_busco_svg, format="svg")
        output_busco_png = input_fasta.replace(".fasta", f"_{compleasm_odb}_busco.png")
        plt.savefig(output_busco_png, format="png")
        plt.close()
        return
    busco_genes = len(busco_df)
    busco_df['Status'] = busco_df['Status'].replace("Complete", "Single")
    status_counts = busco_df.pivot_table(index='Sequence', columns='Status', aggfunc='size', fill_value=0)
    desired_order = ['Single', 'Duplicated', 'Incomplete', 'Fragmented']
    status_counts = status_counts.reindex(columns=desired_order, fill_value=0)
    total_sequences = len(status_counts)
    excluded_sequences = len(status_counts.loc[status_counts.drop(columns='Duplicated', errors='ignore').sum(axis=1) == 0])
    included_sequences = total_sequences - excluded_sequences
    filtered_status_counts = status_counts.loc[status_counts.drop(columns='Duplicated', errors='ignore').sum(axis=1) > 0]
    if filtered_status_counts.empty:
        log_print("WARNING: No valid BUSCO data available for plotting. Using default autorange.")
        plt.figure(figsize=(12, 8))
        plt.text(0.5, 0.5, "No valid BUSCO data to plot", fontsize=14, ha='center', va='center')
        plt.xticks([])
        plt.yticks([])
        plt.title("BUSCO Status Plot - No Data Available")
        output_busco_svg = input_fasta.replace(".fasta", f"_{compleasm_odb}_busco.svg")
        plt.savefig(output_busco_svg, format="svg")
        output_busco_png = input_fasta.replace(".fasta", f"_{compleasm_odb}_busco.png")
        plt.savefig(output_busco_png, format="png")
        plt.close()
        return
    filtered_status_counts = filtered_status_counts.loc[filtered_status_counts.sum(axis=1).sort_values(ascending=False).index]
    status_totals = busco_df['Status'].value_counts()
    colors = {'Single': '#619B8AFF', 'Duplicated': '#A1C181FF', 'Incomplete': '#FE7F2DFF', 'Fragmented': '#FCCA46FF'}
    ax = filtered_status_counts.plot(kind='bar', stacked=True, figsize=(12, 8), color=[colors[col] for col in filtered_status_counts.columns])
    legend_labels = [f"{status} ({round((status_totals.get(status, 0)/busco_genes)*100, 2)}%)" for status in filtered_status_counts.columns]
    completeness_values = [round((status_totals.get(status, 0)/busco_genes)*100, 2) for status in filtered_status_counts.columns]
    completeness_calc = round(completeness_values[0] + completeness_values[1], 2)
    plt.title(f"Distribution of {compleasm_odb} BUSCO Status per Sequence\nCompleteness: {completeness_calc}%")
    plt.xlabel(f"Sequences (Contig/Scaffold/Chromosome)\nIncluded={included_sequences}, Excluded={excluded_sequences}")
    plt.ylabel(f"Number of BUSCO Matches (out of {busco_genes})")
    plt.xticks(rotation=45, ha='right')
    ax.legend(legend_labels, title="BUSCO Status", loc='upper right')
    plt.tight_layout()
    output_busco_svg = input_fasta.replace(".fasta", f"_{compleasm_odb}_busco.svg")
    plt.savefig(output_busco_svg, format="svg")
    output_busco_png = input_fasta.replace(".fasta", f"_{compleasm_odb}_busco.png")
    plt.savefig(output_busco_png, format="png")
    log_print(f"PASS:\tBUSCO {compleasm_odb} plot saved: {output_busco_svg} & {output_busco_png}")
    plt.close()


def comparative_plots(comparative_plot_dict, shared_root):
    """
    Generate a 2x2 comparative bar chart for assembly statistics.

    Creates bar charts for various assembly metrics using data in comparative_plot_dict,
    and saves the plots as SVG and PNG files in the shared_root directory.

    Parameters:
        comparative_plot_dict (dict): Statistics for each assembly method.
        shared_root (str): Directory in which to save the plots.
    """
    log_print("Generating comparative plots...")
    if not comparative_plot_dict:
        log_print("WARNING: Comparative plot data is missing. Skipping plot generation.")
        return
    masurca_stats = comparative_plot_dict.get("MaSuRCA", [0, 0, 0, 0])
    flye_stats = comparative_plot_dict.get("Flye", [0, 0, 0, 0])
    spades_stats = comparative_plot_dict.get("SPAdes", [0, 0, 0, 0])
    hifiasm_stats = comparative_plot_dict.get("hifiasm", [0, 0, 0, 0])
    shared_root = os.path.commonpath([p for p in [shared_root] if p is not None])
    masurca_color = (35/255, 61/255, 77/255)
    flye_color = (161/255, 193/255, 129/255)
    spades_color = (254/255, 127/255, 45/255)
    hifiasm_color = (97/255, 155/255, 138/255)
    asterisk_color_for = [(1.0, 1.0, 1.0),
                          (42/255, 32/255, 53/255),
                          (42/255, 32/255, 53/255),
                          (42/255, 32/255, 53/255)]
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), facecolor="white")
    def plot_stat(ax, title, stat_index, y_max=None, highlight_smallest=False):
        values = [masurca_stats[stat_index], flye_stats[stat_index], spades_stats[stat_index], hifiasm_stats[stat_index]]
        if all(v == 0 for v in values):
            ax.set_title(f"{title} (No Data Available)")
            ax.set_xticks([])
            ax.set_yticks([])
            return
        x = np.arange(4)
        colors = [masurca_color, flye_color, spades_color, hifiasm_color]
        ax.bar(x, values, color=colors, width=0.6)
        ax.set_title(title)
        ax.set_xticks(x)
        ax.set_xticklabels(["MaSuRCA", "Flye", "SPAdes", "hifiasm"])
        if y_max is not None:
            ax.set_ylim([0, y_max])
        highlighted_value = min(values) if highlight_smallest else max(values)
        idx = values.index(highlighted_value)
        offset = 0.02 * highlighted_value
        asterisk_y = highlighted_value - offset
        asterisk_color = asterisk_color_for[idx]
        ax.text(x[idx], asterisk_y, "*", color=asterisk_color, ha="center", va="top",
                fontsize=16, fontweight="bold")
    plot_stat(axes[0, 0], "First Compleasm C", stat_index=0, y_max=100)
    plot_stat(axes[0, 1], "Second Compleasm C", stat_index=1, y_max=100)
    plot_stat(axes[1, 0], "N50", stat_index=2)
    plot_stat(axes[1, 1], "Contig Count", stat_index=3, highlight_smallest=True)
    plt.tight_layout()
    output_svg = os.path.join(shared_root, "comparative_stats.svg")
    output_png = os.path.join(shared_root, "comparative_stats.png")
    plt.savefig(output_svg, format="svg")
    plt.savefig(output_png, format="png")
    plt.close()
    log_print(f"PASS: Comparative plot saved at {output_svg} & {output_png}")


def qc_assembly(final_assembly_path, shared_root, cwd, ONT_RAW_READS, ILLUMINA_RAW_F_READS,
                ILLUMINA_RAW_R_READS, PACBIO_RAW_READS, SPECIES_ID, first_compleasm_odb,
                second_compleasm_odb, REF_SEQ, karyote_id, kingdom_id, sample_stats_dict, results_df):
    """
    Perform quality control on an assembly using Compleasm, Quast, and coverage estimates.

    Runs BUSCO/Compleasm on two lineages, executes Quast to gather assembly statistics, and
    calculates coverage. Updates sample_stats_dict with the metrics and returns a list of key metrics.

    Parameters:
        final_assembly_path (str): Path to the assembly FASTA file.
        shared_root (str): Directory used for storing intermediate outputs.
        cwd (str): Current working directory.
        ONT_RAW_READS (str): Path to raw ONT reads.
        ILLUMINA_RAW_F_READS (str): Path to raw forward Illumina reads.
        ILLUMINA_RAW_R_READS (str): Path to raw reverse Illumina reads.
        PACBIO_RAW_READS (str): Path to raw PacBio reads.
        SPECIES_ID (str): Sample/species identifier.
        first_compleasm_odb (str): First Compleasm/BUSCO lineage identifier.
        second_compleasm_odb (str): Second Compleasm/BUSCO lineage identifier.
        REF_SEQ (str): Reference genome FASTA path.
        karyote_id (str): Organism karyote type.
        kingdom_id (str): Organism kingdom.
        sample_stats_dict (dict): Dictionary to update with QC metrics.
        results_df (DataFrame): DataFrame collecting results.

    Returns:
        final_assembly_path (str):
        sample_stats_list (str):
        sample_stats_dict (str):
    """
    first_compleasm_dir = final_assembly_path.replace(".fasta", f"_{first_compleasm_odb}_compleasm")
    if not os.path.exists(first_compleasm_dir):
        os.makedirs(first_compleasm_dir)
    first_compleasm_summary = os.path.join(first_compleasm_dir, "summary.txt")
    first_compleasm_tsv = os.path.join(first_compleasm_dir, first_compleasm_odb, "full_table_busco_format.tsv")
    if os.path.exists(first_compleasm_summary):
        log_print(f"SKIP\tFirst Compleasm Summary already exists: {first_compleasm_summary}.")
    else:
        first_compleasm_cmd = ["compleasm", "run", "-a", final_assembly_path,
                               "-o", first_compleasm_dir,
                               "--lineage", first_compleasm_odb, "-t", str(CPU_THREADS)]
        _ = run_subprocess_cmd(first_compleasm_cmd, shell_check=False)
    comp_1_busco_svg = final_assembly_path.replace(".fasta", f"_{first_compleasm_odb}_busco.svg")
    if not os.path.exists(comp_1_busco_svg):
        plot_busco(first_compleasm_odb, first_compleasm_tsv, final_assembly_path)

    second_compleasm_dir = final_assembly_path.replace(".fasta", f"_{second_compleasm_odb}_compleasm")
    if not os.path.exists(second_compleasm_dir):
        os.makedirs(second_compleasm_dir)
    second_compleasm_summary = os.path.join(second_compleasm_dir, "summary.txt")
    second_compleasm_tsv = os.path.join(second_compleasm_dir, second_compleasm_odb, "full_table_busco_format.tsv")
    if os.path.exists(second_compleasm_summary):
        log_print(f"SKIP\tSecond Compleasm Summary already exists: {second_compleasm_summary}.")
    else:
        second_compleasm_cmd = ["compleasm", "run", "-a", final_assembly_path,
                                "-o", second_compleasm_dir,
                                "--lineage", second_compleasm_odb, "-t", str(CPU_THREADS)]
        _ = run_subprocess_cmd(second_compleasm_cmd, shell_check=False)
    comp_2_busco_svg = final_assembly_path.replace(".fasta", f"_{second_compleasm_odb}_busco.svg")
    if not os.path.exists(comp_2_busco_svg):
        plot_busco(second_compleasm_odb, second_compleasm_tsv, final_assembly_path)

    quast_dir = final_assembly_path.replace(".fasta", "_quast")
    if not os.path.exists(quast_dir):
        os.makedirs(quast_dir)
    quast_report_tsv = os.path.join(quast_dir, "report.tsv")
    if os.path.exists(quast_report_tsv):
        log_print(f"SKIP\tQUAST Report already exists: {quast_report_tsv}.")
    else:
        if pd.isna(REF_SEQ):
            quast_cmd = ["quast", "--threads", str(CPU_THREADS),
                         f"--{karyote_id}", "-o", quast_dir, final_assembly_path]
        else:
            quast_cmd = ["quast", "--threads", str(CPU_THREADS),
                         "-r", REF_SEQ, f"--{karyote_id}",
                         "-o", quast_dir, final_assembly_path]
        if kingdom_id == "Funga":
            quast_cmd.append("--fungus")
        _ = run_subprocess_cmd(quast_cmd, shell_check=False)

    with open(first_compleasm_summary, "r") as first_compleasm_file:
        for line in first_compleasm_file:
            if "S:" in line:
                sample_stats_dict["FIRST_COMPLEASM_S"] = float(line.split("S:")[-1].split(", ")[0].replace("\n", "").replace("%", ""))
            elif "D:" in line:
                sample_stats_dict["FIRST_COMPLEASM_D"] = float(line.split("D:")[-1].split(", ")[0].replace("\n", "").replace("%", ""))
            elif "F:" in line:
                sample_stats_dict["FIRST_COMPLEASM_F"] = float(line.split("F:")[-1].split(", ")[0].replace("\n", "").replace("%", ""))
            elif "M:" in line:
                sample_stats_dict["FIRST_COMPLEASM_M"] = float(line.split("M:")[-1].split(", ")[0].replace("\n", "").replace("%", ""))
    sample_stats_dict["FIRST_COMPLEASM_C"] = sample_stats_dict["FIRST_COMPLEASM_S"] + sample_stats_dict["FIRST_COMPLEASM_D"]

    with open(second_compleasm_summary, "r") as second_compleasm_file:
        for line in second_compleasm_file:
            if "S:" in line:
                sample_stats_dict["SECOND_COMPLEASM_S"] = float(line.split("S:")[-1].split(", ")[0].replace("\n", "").replace("%", ""))
            elif "D:" in line:
                sample_stats_dict["SECOND_COMPLEASM_D"] = float(line.split("D:")[-1].split(", ")[0].replace("\n", "").replace("%", ""))
            elif "F:" in line:
                sample_stats_dict["SECOND_COMPLEASM_F"] = float(line.split("F:")[-1].split(", ")[0].replace("\n", "").replace("%", ""))
            elif "M:" in line:
                sample_stats_dict["SECOND_COMPLEASM_M"] = float(line.split("M:")[-1].split(", ")[0].replace("\n", "").replace("%", ""))
    sample_stats_dict["SECOND_COMPLEASM_C"] = sample_stats_dict["SECOND_COMPLEASM_S"] + sample_stats_dict["SECOND_COMPLEASM_D"]

    with open(quast_report_tsv, "r") as quast_file:
        for line in quast_file:
            if "Total length (>= 0 bp)" in line:
                sample_stats_dict["GENOME_SIZE"] = float(line.split("\t")[-1].replace("\n", ""))
            elif "# contigs" in line:
                sample_stats_dict["ASSEMBLY_CONTIGS"] = float(line.split("\t")[-1].replace("\n", ""))
            elif "N50" in line:
                sample_stats_dict["ASSEMBLY_N50"] = float(line.split("\t")[-1].replace("\n", ""))
            elif "L50" in line:
                sample_stats_dict["ASSEMBLY_L50"] = float(line.split("\t")[-1].replace("\n", ""))
            elif "GC (%)" in line:
                sample_stats_dict["ASSEMBLY_GC"] = float(line.split("\t")[-1].replace("\n", ""))
            if REF_SEQ is not None:
                if "# misassemblies" in line:
                    sample_stats_dict["MISASSEMBLIES"] = float(line.split("\t")[-1].replace("\n", ""))
                elif "# N's per 100 kbp" in line:
                    sample_stats_dict["N_PER_100KBP"] = float(line.split("\t")[-1].replace("\n", ""))
                elif "# mismatches per 100 kbp" in line:
                    sample_stats_dict["MIS_PER_100KBP"] = float(line.split("\t")[-1].replace("\n", ""))
                elif "# indels per 100 kbp" in line:
                    sample_stats_dict["INDELS_PER_100KPB"] = float(line.split("\t")[-1].replace("\n", ""))
    if not pd.isna(REF_SEQ):
        ref_total_bases = 0
        for record in SeqIO.parse(REF_SEQ, "fasta"):
            ref_total_bases += len(record.seq)
        if not pd.isna(ILLUMINA_RAW_F_READS) and not pd.isna(ILLUMINA_RAW_R_READS):
            sample_stats_dict["RAW_ILLU_COVERAGE"] = round(sample_stats_dict["RAW_ILLU_TOTAL_BASES"] / ref_total_bases, 2)
            sample_stats_dict["TRIMMED_ILLU_COVERAGE"] = round(sample_stats_dict["TRIMMED_ILLU_TOTAL_BASES"] / ref_total_bases, 2)
            sample_stats_dict["DEDUPED_ILLU_COVERAGE"] = round(sample_stats_dict["DEDUPED_ILLU_TOTAL_BASES"] / ref_total_bases, 2)
        if not pd.isna(ONT_RAW_READS):
            sample_stats_dict["RAW_ONT_COVERAGE"] = round(sample_stats_dict["RAW_ONT_TOTAL_BASES"] / ref_total_bases, 2)
            sample_stats_dict["FILT_ONT_COVERAGE"] = round(sample_stats_dict["FILT_ONT_TOTAL_BASES"] / ref_total_bases, 2)
            sample_stats_dict["CORRECT_ONT_COVERAGE"] = round(sample_stats_dict["CORRECT_ONT_TOTAL_BASES"] / ref_total_bases, 2)
        if not pd.isna(PACBIO_RAW_READS):
            sample_stats_dict["RAW_PACBIO_COVERAGE"] = round(sample_stats_dict["RAW_PACBIO_TOTAL_BASES"] / ref_total_bases, 2)
            sample_stats_dict["HIFI_PACBIO_COVERAGE"] = round(sample_stats_dict["HIFI_PACBIO_TOTAL_BASES"] / ref_total_bases, 2)
            sample_stats_dict["FILT_PACBIO_COVERAGE"] = round(sample_stats_dict["FILT_PACBIO_TOTAL_BASES"] / ref_total_bases, 2)

    else:
        if not pd.isna(ILLUMINA_RAW_F_READS) and not pd.isna(ILLUMINA_RAW_R_READS):
            sample_stats_dict["RAW_ILLU_COVERAGE"] = round(sample_stats_dict["RAW_ILLU_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
            sample_stats_dict["TRIMMED_ILLU_COVERAGE"] = round(sample_stats_dict["TRIMMED_ILLU_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
            sample_stats_dict["DEDUPED_ILLU_COVERAGE"] = round(sample_stats_dict["DEDUPED_ILLU_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
        if not pd.isna(ONT_RAW_READS):
            sample_stats_dict["RAW_ONT_COVERAGE"] = round(sample_stats_dict["RAW_ONT_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
            sample_stats_dict["FILT_ONT_COVERAGE"] = round(sample_stats_dict["FILT_ONT_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
            sample_stats_dict["CORRECT_ONT_COVERAGE"] = round(sample_stats_dict["CORRECT_ONT_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
        if not pd.isna(PACBIO_RAW_READS):
            sample_stats_dict["RAW_PACBIO_COVERAGE"] = round(sample_stats_dict["RAW_PACBIO_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
            sample_stats_dict["HIFI_PACBIO_COVERAGE"] = round(sample_stats_dict["HIFI_PACBIO_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
            sample_stats_dict["FILT_PACBIO_COVERAGE"] = round(sample_stats_dict["FILT_PACBIO_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
    
    first_compleasm_c = sample_stats_dict["FIRST_COMPLEASM_S"] + sample_stats_dict["FIRST_COMPLEASM_D"]
    second_compleasm_c = sample_stats_dict["SECOND_COMPLEASM_S"] + sample_stats_dict["SECOND_COMPLEASM_D"]
    n50 = sample_stats_dict["ASSEMBLY_N50"]
    contig_count = sample_stats_dict["ASSEMBLY_CONTIGS"]
    sample_stats_list = [first_compleasm_c, second_compleasm_c, n50, contig_count]
    return final_assembly_path, sample_stats_list, sample_stats_dict


def gen_sample_stats_dict(row):
    """
    Generate an initial sample statistics dictionary from a metadata row.

    Extracts key fields from the row and initializes placeholders for various metrics.

    Parameters:
        row (pandas.Series): A sample metadata row.

    Returns:
        sample_stats_dict (dict): A dictionary with keys for sample stats.
    """
    sample_stats_dict = {"SPECIES_ID": row["SPECIES_ID"],
                         "ONT_SRA": row["ONT_SRA"] if isinstance(row["ONT_SRA"], str) else None,
                         "ONT": os.path.basename(row["ONT_RAW_READS"]) if isinstance(row["ONT_RAW_READS"], str) else None,
                         "ILLU_SRA": row["ILLUMINA_SRA"] if isinstance(row["ILLUMINA_SRA"], str) else None,
                         "ILLU_F": os.path.basename(row["ILLUMINA_RAW_F_READS"]) if isinstance(row["ILLUMINA_RAW_F_READS"], str) else None,
                         "ILLU_R": os.path.basename(row["ILLUMINA_RAW_R_READS"]) if isinstance(row["ILLUMINA_RAW_R_READS"], str) else None,
                         "PACBIO_SRA": row["PACBIO_SRA"] if isinstance(row["PACBIO_SRA"], str) else None,
                         "PACBIO": os.path.basename(row["PACBIO_RAW_READS"]) if isinstance(row["PACBIO_RAW_READS"], str) else None,
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


def nanoplot_qc_reads(INPUT_READS, READS_ORIGIN, CPU_THREADS, sample_stats_dict):
    """
    Perform quality control on reads using NanoPlot.

    Parameters:
        INPUT_READS (str): Path to the input reads file.
        READS_ORIGIN (str): Origin of the reads (e.g., "Raw_ONT_", "Filt_ONT_").
        CPU_THREADS (int): Number of CPU threads to use.
        sample_stats_dict (dict): Dictionary to update with QC metrics.

    Returns:
        sample_stats_dict (dict): Updated sample_stats_dict with QC metrics.
    """
    output_dir = "/".join(INPUT_READS.split("/")[:-1]) + f"/{READS_ORIGIN}nanoplot_analysis"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    out_file = os.path.join(output_dir, "{READS_ORIGIN}NanoStats.txt")
    if os.path.exists(out_file):
        log_print(f"SKIP:\tNanoPlot output already exists: {out_file}.")
    else:    
        raw_nanoplot_cmd = [ "NanoPlot", "--fastq", INPUT_READS, "-t", str(CPU_THREADS),
                            "-o", output_dir, "--plots", "kde", "dot", "--loglength",
                            "--N50", "--title", f"{READS_ORIGIN} Reads: Preliminary Data",
                            "--prefix", READS_ORIGIN, "--verbose"]
        _ = run_subprocess_cmd(raw_nanoplot_cmd, shell_check = False)
    with open(out_file, "r") as nanostats:
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


def fastqc_reads(INPUT_F_READS, INPUT_R_READS, READS_ORIGIN, CPU_THREADS, sample_stats_dict):
    """
    DESCRIPTION

    Parameters:
        INPUT_F_READS (str): DESCRIPTION.
        INPUT_R_READS (str): DESCRIPTION.
        READS_ORIGIN (str): DESCRIPTION.
        CPU_THREADS (str): DESCRIPTION.
        sample_stats_dict (str): DESCRIPTION.

    Returns:
        sample_stats_dict (dict): DESCRIPTION.
    """
    fastqc_dir = "/".join(INPUT_F_READS.split("/")[:-1]) + f"/{READS_ORIGIN}fastqc_analysis"    
    if not os.path.exists(fastqc_dir):
        os.makedirs(fastqc_dir)
    fastqc_F_out_file = os.path.join(fastqc_dir, INPUT_F_READS.split("/")[-1].replace(".fastq.gz","_fastqc.html"))
    fastqc_R_out_file = os.path.join(fastqc_dir, INPUT_R_READS.split("/")[-1].replace(".fastq.gz","_fastqc.html"))
    if os.path.exists(fastqc_F_out_file) and os.path.exists(fastqc_R_out_file):
        log_print(f"SKIP:\tFastQC outputs already exist: {fastqc_F_out_file}; {fastqc_R_out_file}.")
    else:
        fastqc_cmd = ["fastqc", "-o", fastqc_dir, "-t", str(CPU_THREADS), INPUT_F_READS, INPUT_R_READS]
        _ = run_subprocess_cmd(fastqc_cmd, shell_check = False)
    if "raw" in READS_ORIGIN.lower():
        sample_stats_dict["RAW_ILLU_TOTAL_BASES"] =  round((float(get_total_bases(fastqc_F_out_file).split(" ")[0]) + float(get_total_bases(fastqc_R_out_file).split(" ")[0])) / 2, 2)
    if "trim" in READS_ORIGIN.lower():
        sample_stats_dict["TRIMMED_ILLU_TOTAL_BASES"] =  round((float(get_total_bases(fastqc_F_out_file).split(" ")[0]) + float(get_total_bases(fastqc_R_out_file).split(" ")[0])) / 2, 2)
    if "dedup" in READS_ORIGIN.lower():
        sample_stats_dict["DEDUPED_ILLU_TOTAL_BASES"] =  round((float(get_total_bases(fastqc_F_out_file).split(" ")[0]) + float(get_total_bases(fastqc_R_out_file).split(" ")[0])) / 2, 2)
    return sample_stats_dict


def racon_polish_assembly(input_assembly, highest_mean_qual_ont_reads, assembly_out_dir, iteration_count):
    """
    DESCRIPTION

    Parameters:
        input_assembly (str): DESCRIPTION.
        highest_mean_qual_ont_reads (str): DESCRIPTION.
        assembly_out_dir (str): DESCRIPTION.
        iteration_count (str): DESCRIPTION.

    Returns:
        racon_assembly (str): DESCRIPTION.
    """
    racon_paf = os.path.join(assembly_out_dir, f"racon_round{iteration_count}.paf")
    if os.path.exists(racon_paf):
        log_print(f"SKIP:\tRacon PAF {iteration_count} already exists: {racon_paf}.")
    else:
        minimap2_cmd = f"minimap2 -t {CPU_THREADS} -x map-ont {input_assembly} {highest_mean_qual_ont_reads} > {racon_paf}"
        _ = run_subprocess_cmd(minimap2_cmd, shell_check = True)
    racon_assembly = os.path.join(assembly_out_dir, f"assembly_racon{iteration_count}.fasta")
    if os.path.exists(racon_assembly):
        log_print(f"SKIP:\tRacon Assembly {iteration_count} already exists: {racon_assembly}.")
    else:
        racon_cmd = f"racon -t {CPU_THREADS} {highest_mean_qual_ont_reads} {racon_paf} {input_assembly} > {racon_assembly}"
        _ = run_subprocess_cmd(racon_cmd, shell_check = True)        
    return racon_assembly


def download_test_data(SPECIES_ID, ILLUMINA_SRA, ONT_SRA, PACBIO_SRA, REF_SEQ_GCA):
    """
    Download test data (Illumina, ONT, and reference sequences) from SRA/NCBI.

    Creates an "EGAP_Test_Data" folder (with a subfolder named after SPECIES_ID) and downloads
    the requested data using prefetch and fastq-dump for SRA data and datasets download for reference genomes.

    Parameters:
        SPECIES_ID (str): Sample/species identifier.
        ILLUMINA_SRA (str or NaN): Illumina SRA accession.
        ONT_SRA (str or NaN): ONT SRA accession.
        REF_SEQ_GCA (str or NaN): Reference genome GCA accession.

    Returns:
        illu_sra_f (str): Path to the Illumina forward SRA or None.
        illu_sra_r (str): Path to the Illumina reverse SRA or None.
        ont_sra (str): Path to the ONT SRA or None.
        pacbio_sra (str): Path to the PacBio SRA or None.
        ref_seq_gca (str): Path to the Reference Sequence Assembly GCA or None.
        EGAP_test_data_dir (str): Path to the directory data is being saved to.
    """
    cwd = os.getcwd()
    EGAP_test_dir = os.path.join(cwd, "EGAP_Test_Data")
    if not os.path.exists(EGAP_test_dir):
        os.mkdir(EGAP_test_dir)
    if len(SPECIES_ID.split("_")) == 2:
        EGAP_test_data_dir = os.path.join(cwd, "EGAP_Test_Data", SPECIES_ID)
    elif len(SPECIES_ID.split("_")) > 2:
        EGAP_test_data_dir = os.path.join(cwd, "EGAP_Test_Data", f"{SPECIES_ID.split('_')[0]}_{SPECIES_ID.split('_')[1]}", SPECIES_ID)
    if not os.path.exists(EGAP_test_data_dir):
        os.mkdir(EGAP_test_data_dir)
    illu_sra_f = None
    illu_sra_r = None
    ont_sra = None
    pacbio_sra = None
    ref_seq_gca = None
    if not pd.isna(PACBIO_SRA):
        pacbio_test_dir = os.path.join(EGAP_test_data_dir, "PacBio")
        pacbio_sra = os.path.join(pacbio_test_dir, f"{PACBIO_SRA}.fastq.gz")
        if not os.path.exists(pacbio_sra):
            if not os.path.exists(pacbio_test_dir):
                os.mkdir(pacbio_test_dir)
            os.chdir(pacbio_test_dir)
            pacbio_cmd = f"prefetch {PACBIO_SRA} && fastq-dump --gzip --split-files {PACBIO_SRA} && rm -rf {PACBIO_SRA}"
            _ = run_subprocess_cmd(pacbio_cmd, shell_check=True)
        else:
            log_print(f"SKIP:\tPacBio SRA already exists: {pacbio_sra}")
    os.chdir(EGAP_test_data_dir)
    if not pd.isna(ILLUMINA_SRA):
        illumina_test_dir = os.path.join(EGAP_test_data_dir, "Illumina")
        illu_sra_f = os.path.join(illumina_test_dir, f"{ILLUMINA_SRA}_1.fastq.gz")
        illu_sra_r = os.path.join(illumina_test_dir, f"{ILLUMINA_SRA}_2.fastq.gz")
        if not os.path.exists(illu_sra_f) and not os.path.exists(illu_sra_r):
            if not os.path.exists(illumina_test_dir):
                os.mkdir(illumina_test_dir)
            os.chdir(illumina_test_dir)
            illu_cmd = f"prefetch {ILLUMINA_SRA} && fastq-dump --gzip --split-files {ILLUMINA_SRA} && rm -rf {ILLUMINA_SRA}"
            _ = run_subprocess_cmd(illu_cmd, shell_check=True)
        else:
            log_print(f"SKIP:\tIllumina SRAs already exist: {illu_sra_f}; {illu_sra_r}")
    os.chdir(EGAP_test_data_dir)
    if not pd.isna(ONT_SRA):
        ont_test_dir = os.path.join(EGAP_test_data_dir, "ONT")
        ont_sra = os.path.join(ont_test_dir, f"{ONT_SRA}.fastq.gz")
        if not os.path.exists(ont_sra):
            if not os.path.exists(ont_test_dir):
                os.mkdir(ont_test_dir)
            os.chdir(ont_test_dir)
            ont_cmd = f"prefetch {ONT_SRA} && fastq-dump --gzip {ONT_SRA} && rm -rf {ONT_SRA}"
            _ = run_subprocess_cmd(ont_cmd, shell_check=True)
        else:
            log_print(f"SKIP:\tONT SRAs already exists: {ont_sra}")
    os.chdir(EGAP_test_data_dir)
    ref_seq_gca_dir = os.path.join(EGAP_test_data_dir, f"ncbi_dataset/data/{REF_SEQ_GCA}/")
    renamed_gca = os.path.join(EGAP_test_data_dir, f"{REF_SEQ_GCA}.fasta")
    if not pd.isna(REF_SEQ_GCA):
        if not os.path.exists(renamed_gca):
            try:
                ref_seq_gca = glob.glob(os.path.join(ref_seq_gca_dir, "*_genomic.fna"))[0]
                if not os.path.exists(ref_seq_gca):
                    ref_seq_cmd = f"datasets download genome accession {REF_SEQ_GCA} --include genome &&  unzip -o ncbi_dataset -d {EGAP_test_data_dir}"
                    _ = run_subprocess_cmd(ref_seq_cmd, shell_check=True)
            except IndexError:
                ref_seq_cmd = f"datasets download genome accession {REF_SEQ_GCA} --include genome &&  unzip -o ncbi_dataset -d {EGAP_test_data_dir}"
                _ = run_subprocess_cmd(ref_seq_cmd, shell_check=True)
            ref_seq_gca = glob.glob(os.path.join(ref_seq_gca_dir, "*_genomic.fna"))[0]
            shutil.move(ref_seq_gca, renamed_gca)
        else:
            log_print(f"SKIP:\tREF_SEQ GCA already exists: {renamed_gca}")
    return illu_sra_f, illu_sra_r, ont_sra, pacbio_sra, ref_seq_gca, EGAP_test_data_dir


def process_final_assembly(row, results_df, CPU_THREADS, RAM_GB, final_assembly_path=None, sample_stats_dict=None):
    """
    Finalize and assess the assembly, updating the results.

    Runs quality control (via qc_assembly), logs metric classifications, updates the results DataFrame,
    and returns the final assembly path along with updated results.

    Parameters:
        row (pandas.Series): Sample metadata.
        results_df (DataFrame): DataFrame for cumulative results.
        CPU_THREADS (int): Number of CPU threads.
        RAM_GB (int): RAM in GB.
        final_assembly_path (str, optional): Path to the assembly; if None, a default is used.
        sample_stats_dict (dict, optional): Dictionary of sample metrics.

    Returns:
        final_assembly_path (str): 
        results_df (dataframe): 
    """
    SPECIES_ID = row["SPECIES_ID"]
    REF_SEQ_GCA = row["REF_SEQ_GCA"]
    cwd = os.getcwd()
    EGAP_test_dir = os.path.join(cwd, "EGAP_Test_Data")
    if not os.path.exists(EGAP_test_dir):
        os.mkdir(EGAP_test_dir)
    EGAP_test_data_dir = os.path.join(cwd, "EGAP_Test_Data", SPECIES_ID)
    if not os.path.exists(EGAP_test_data_dir):
        os.mkdir(EGAP_test_data_dir)
    ref_seq_gca_dir = os.path.join(EGAP_test_data_dir, f"ncbi_dataset/data/{REF_SEQ_GCA}/")
    ref_seq_gca_fasta = glob.glob(os.path.join(ref_seq_gca_dir, "*_genomic.fna"))
    renamed_gca = os.path.join(EGAP_test_data_dir, f"{REF_SEQ_GCA}.fasta")
    if not pd.isna(REF_SEQ_GCA):
        if not os.path.exists(renamed_gca):
            if len(ref_seq_gca_fasta) == 0:
                ref_seq_cmd = f"datasets download genome accession {REF_SEQ_GCA} --include genome && unzip -o ncbi_dataset -d {EGAP_test_data_dir}"
                _ = run_subprocess_cmd(ref_seq_cmd, shell_check=True)
            ref_seq_gca_fasta = glob.glob(os.path.join(ref_seq_gca_dir, "*_genomic.fna"))[0]            
            shutil.move(ref_seq_gca_fasta, renamed_gca)
        else:
            log_print(f"SKIP:\tREF_SEQ GCA already exists: {renamed_gca}")
    if pd.isna(final_assembly_path):
        final_assembly_path = renamed_gca
    shared_root = EGAP_test_data_dir
    ONT_RAW_READS = row["ONT_RAW_READS"]
    ILLUMINA_RAW_F_READS = row["ILLUMINA_RAW_F_READS"]
    ILLUMINA_RAW_R_READS = row["ILLUMINA_RAW_R_READS"]
    PACBIO_RAW_READS = row["PACBIO_RAW_READS"]
    first_compleasm_odb = row["COMPLEASM_1"] + "_odb10"
    second_compleasm_odb = row["COMPLEASM_2"] + "_odb10"
    REF_SEQ = None
    karyote_id = row["ORGANISM_KARYOTE"]
    kingdom_id = row["ORGANISM_KINGDOM"]
    if pd.isna(sample_stats_dict):
        sample_stats_dict = gen_sample_stats_dict(row)
    final_assembly_path, sample_stats_list, sample_stats_dict = qc_assembly(final_assembly_path, shared_root, cwd,
                                                                            ONT_RAW_READS,
                                                                            ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS, 
                                                                            PACBIO_RAW_READS,
                                                                            SPECIES_ID,
                                                                            first_compleasm_odb, second_compleasm_odb,
                                                                            REF_SEQ, karyote_id, kingdom_id,
                                                                            sample_stats_dict, results_df)
    quality_classifications = classify_assembly(sample_stats_dict)
    for metric, classification in quality_classifications.items():
        log_print(f"{metric}: {classification}")
    # (Additional result processing omitted for brevity.)
    log_print(f"Assessment of Final Assembly: {final_assembly_path}")
    log_print(f"PASS:\tEGAP Final Assembly Complete: {final_assembly_path}")
    log_print("This was produced with the help of the Entheogen Genome (Entheome) Foundation\n")
    log_print("If this was useful, please support us at https://entheome.org/\n")
    print("\n\n\n")
    return final_assembly_path, results_df


def egap_sample(row, results_df, INPUT_CSV, CPU_THREADS, RAM_GB):
    """
    Run the EGAP pipeline on a single sample.

    Performs read processing, assembly, polishing, and curation for the sample,
    updating the results DataFrame.

    Parameters:
        row (pandas.Series): Sample metadata.
        results_df (DataFrame): Cumulative results DataFrame.
        INPUT_CSV (str): Path to the input CSV (if any).
        CPU_THREADS (int): Number of CPU threads to use.
        RAM_GB (int): Amount of RAM in GB.

    Returns:
        tuple: (final_assembly_path, results_df)
    """
    cwd = os.getcwd()
    ONT_SRA = row["ONT_SRA"]
    ONT_RAW_DIR = row["ONT_RAW_DIR"]
    ONT_RAW_READS = row["ONT_RAW_READS"]
    ILLU_SRA = row["ILLUMINA_SRA"]
    ILLU_RAW_DIR = row["ILLUMINA_RAW_DIR"]
    ILLUMINA_RAW_F_READS = row["ILLUMINA_RAW_F_READS"]
    ILLUMINA_RAW_R_READS = row["ILLUMINA_RAW_R_READS"]
    PACBIO_SRA = row["PACBIO_SRA"]
    PACBIO_RAW_DIR = row["PACBIO_RAW_DIR"]
    PACBIO_RAW_READS = row["PACBIO_RAW_READS"]
    SPECIES_ID = row["SPECIES_ID"]
    EST_SIZE = row["EST_SIZE"]
    REF_SEQ_GCA = row["REF_SEQ_GCA"]
    REF_SEQ = row["REF_SEQ"]
    
    if pd.notna(REF_SEQ_GCA) or pd.notna(ILLU_SRA) or pd.notna(ONT_SRA) or pd.notna(PACBIO_SRA):
        if pd.isna(ONT_SRA):
            # NO ONT NOR PACBIO AKA ILLU ONLY
            if pd.isna(PACBIO_SRA):
                ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS, _, _, REF_SEQ, EGAP_test_data_dir = download_test_data(SPECIES_ID, ILLU_SRA, ONT_SRA, PACBIO_SRA, REF_SEQ_GCA)
            # NO ONT
            else:
                ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS, _, PACBIO_RAW_READS, REF_SEQ, EGAP_test_data_dir = download_test_data(SPECIES_ID, ILLU_SRA, ONT_SRA, PACBIO_SRA, REF_SEQ_GCA)                                
        elif pd.isna(PACBIO_SRA):
            # NO PACBIO
            ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS, ONT_RAW_READS, _, REF_SEQ, EGAP_test_data_dir = download_test_data(SPECIES_ID, ILLU_SRA, ONT_SRA, PACBIO_SRA, REF_SEQ_GCA)                                
        else:
            # ALL
            ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS, ONT_RAW_READS, PACBIO_RAW_READS, REF_SEQ, EGAP_test_data_dir = download_test_data(SPECIES_ID, ILLU_SRA, ONT_SRA, PACBIO_SRA, REF_SEQ_GCA)

    if type(PACBIO_RAW_DIR) == str:
        shared_root = PACBIO_RAW_DIR
        initialize_logging_environment(shared_root)
        log_print(f"Running Entheome Genome Assembly Pipeline on: {shared_root}")
        log_print("UNFINISHED PACBIO RAW READS DIRECTORY PROCESSING")
        # TODO: run ccs on Raw PacBio reads to turn them into HiFi reads using the bam file

    if type(ONT_RAW_DIR) == str:
        if type(ILLU_RAW_DIR) != str:
            shared_root = os.path.join(ONT_RAW_DIR, SPECIES_ID)
        else:
            shared_root = os.path.commonpath([ONT_RAW_DIR, ILLU_RAW_DIR])
        initialize_logging_environment(shared_root)
        log_print(f"Running Entheome Genome Assembly Pipeline on: {shared_root}")
        ONT_RAW_READS = ont_combine_fastq_gz(ONT_RAW_DIR)
    if type(ILLU_RAW_DIR) == str:
        if type(ONT_RAW_DIR) == str:
            shared_root = os.path.commonpath([ONT_RAW_DIR, ILLU_RAW_DIR])
        if type(PACBIO_RAW_DIR) == str:
            shared_root = os.path.commonpath([PACBIO_RAW_DIR, ILLU_RAW_DIR])
        else:
            shared_root = ILLU_RAW_DIR
        initialize_logging_environment(shared_root)
        log_print(f"Running Entheome Genome Assembly Pipeline on: {shared_root}")
        illu_combined_files = illumina_extract_and_check(ILLU_RAW_DIR, SPECIES_ID)
        ILLUMINA_RAW_F_READS = illu_combined_files[0]
        ILLUMINA_RAW_R_READS = illu_combined_files[1]
    
    if type(PACBIO_RAW_READS) == str:
        if type(ILLUMINA_RAW_F_READS) != str and type(ILLUMINA_RAW_R_READS) != str:
            shared_root = "/".join(PACBIO_RAW_READS.split("/")[:-1])
        else:
            shared_root = os.path.commonpath([PACBIO_RAW_READS, ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS])
        shared_root = shared_root.replace(".","_") if "." in shared_root else shared_root       
        initialize_logging_environment(shared_root)
        log_print(f"Running Entheome Genome Assembly Pipeline on: {shared_root}")
    elif type(ONT_RAW_READS) == str:
        if type(ILLUMINA_RAW_F_READS) != str and type(ILLUMINA_RAW_R_READS) != str:
            shared_root = os.path.join(os.path.dirname(ONT_RAW_READS), SPECIES_ID)
        else:
            shared_root = os.path.commonpath([ONT_RAW_READS, ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS])
        initialize_logging_environment(shared_root)
        log_print(f"Running Entheome Genome Assembly Pipeline on: {shared_root}")
    else:
        if pd.notna(INPUT_CSV):
            shared_root = EGAP_test_data_dir
        elif type(ILLUMINA_RAW_F_READS) == str and type(ILLUMINA_RAW_R_READS) == str:
            shared_root = os.path.commonpath([ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS])
        initialize_logging_environment(shared_root)
        log_print(f"Running Entheome Genome Assembly Pipeline on: {shared_root}")

    if pd.notna(REF_SEQ):
        ref_seq_path = REF_SEQ
        ref_dir, ref_filename = os.path.split(ref_seq_path)
        ref_basename, ref_ext = os.path.splitext(ref_filename)
        if ref_basename.endswith('.fna'):
            ref_basename, _ = os.path.splitext(ref_basename)
            ref_ext = os.path.splitext(ref_ext)[1] + ref_ext
        new_ref_path = None
        if ref_seq_path.endswith('.fna.gz'):
            new_ref_filename = ref_basename + '.fa.gz'
            new_ref_path = os.path.join(ref_dir, new_ref_filename)
            if os.path.exists(ref_seq_path):
                try:
                    os.rename(ref_seq_path, new_ref_path)
                    REF_SEQ = new_ref_path
                    log_print(f"NOTE:\tRenamed REF_SEQ from {ref_seq_path} to {new_ref_path}")
                except Exception as e:
                    log_print(f"ERROR:\tUnable to rename {ref_seq_path} to {new_ref_path}: {e}")
            else:
                fa_gz_path = os.path.join(ref_dir, ref_basename + '.fa.gz')
                if os.path.exists(fa_gz_path):
                    REF_SEQ = fa_gz_path
                    log_print(f"NOTE:\tOriginal file {ref_seq_path} not found. Using existing .fa.gz file: {fa_gz_path}")
                else:
                    log_print(f"NOTE:\tNeither {ref_seq_path} nor {fa_gz_path} exists. REF_SEQ cannot be processed.")
        elif ref_seq_path.endswith('.fna'):
            new_ref_filename = ref_basename + '.fa'
            new_ref_path = os.path.join(ref_dir, new_ref_filename)
            if os.path.exists(ref_seq_path):
                try:
                    os.rename(ref_seq_path, new_ref_path)
                    REF_SEQ = new_ref_path
                    log_print(f"NOTE:\tRenamed REF_SEQ from {ref_seq_path} to {new_ref_path}")
                except Exception as e:
                    log_print(f"ERROR:\tUnable to rename {ref_seq_path} to {new_ref_path}: {e}")
            else:
                fa_path = os.path.join(ref_dir, ref_basename + '.fa')
                if os.path.exists(fa_path):
                    REF_SEQ = fa_path
                    log_print(f"NOTE:\tOriginal file {ref_seq_path} not found. Using existing .fa file: {fa_path}")
                else:
                    log_print(f"NOTE:\tNeither {ref_seq_path} nor {fa_path} exists. REF_SEQ cannot be processed.")
        else:
            log_print(f"NOTE:\tREF_SEQ file {ref_seq_path} does not match expected extensions (.fna or .fna.gz). Skipping renaming.")
    else:
        log_print("NOTE:\tNo REF_SEQ provided; skipping REF_SEQ processing.")

    if pd.notna(ONT_RAW_READS):
        ONT_RAW_READS = process_read_file(ONT_RAW_READS)
    if pd.notna(ILLUMINA_RAW_F_READS):
        ILLUMINA_RAW_F_READS = process_read_file(ILLUMINA_RAW_F_READS)
    if pd.notna(ILLUMINA_RAW_R_READS):
        ILLUMINA_RAW_R_READS = process_read_file(ILLUMINA_RAW_R_READS)
    if pd.notna(PACBIO_RAW_READS):
        PACBIO_RAW_READS = process_read_file(PACBIO_RAW_READS)
    sample_stats_dict = gen_sample_stats_dict(row)
    first_compleasm_odb = f"{row['COMPLEASM_1']}_odb10"
    second_compleasm_odb = f"{row['COMPLEASM_2']}_odb10"
    kingdom_id = row["ORGANISM_KINGDOM"].lower()
    karyote_id = row["ORGANISM_KARYOTE"].lower()
    
    # Convert the estimated genome size (e.g., "25M") to an integer (base pairs)
    size_str = EST_SIZE.lower().strip()
    multipliers = {'m': 10**6, 'g': 10**9}
    if size_str[-1] in multipliers:
        estimated_genome_size = int(float(size_str[:-1]) * multipliers[size_str[-1]])
    else:
        log_print(f"NOTE:\tUnable to parse input estimated size {EST_SIZE}, using default: 25000000")
        estimated_genome_size = 25000000
    EST_SIZE = estimated_genome_size

###############################################################################
    # Reads Pre-Processing
###############################################################################

    if pd.notna(ILLUMINA_RAW_F_READS) and pd.notna(ILLUMINA_RAW_R_READS):
        # FastQC Illumina Raw Reads
        sample_stats_dict = fastqc_reads(ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS, "Raw_Illumina_", CPU_THREADS, sample_stats_dict)
    
        # Trimmomatic of Illumina Raw Reads
        trimmo_f_pair_path = ILLUMINA_RAW_F_READS.replace(".fastq.gz", "_paired.fastq.gz")
        fwd_unpaired_out = trimmo_f_pair_path.replace("paired","unpaired")
        trimmo_r_pair_path = ILLUMINA_RAW_R_READS.replace(".fastq.gz", "_paired.fastq.gz")
        rev_unpaired_out = trimmo_r_pair_path.replace("paired","unpaired")
        HEADCROP = "HEADCROP:10"
        CROP = "CROP:145"
        SLIDINGWINDOW = "SLIDINGWINDOW:50:25"
        MINLEN = "MINLEN:125"
        if os.path.exists(trimmo_f_pair_path) and os.path.exists(trimmo_r_pair_path):
            log_print("SKIP:\tTrimmomatic Paired outputs already exist: {trimmo_f_pair_path}; {trimmo_r_pair_path}")
        else:
            trimmomatic_cmd = ["trimmomatic", "PE", "-threads", str(CPU_THREADS), "-phred33",
                                ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS, 
                                trimmo_f_pair_path, fwd_unpaired_out,
                                trimmo_r_pair_path, rev_unpaired_out,
                                f"ILLUMINACLIP:{find_file('TruSeq3-PE.fa')}:2:30:10:11",
                                HEADCROP, CROP, SLIDINGWINDOW, MINLEN]
            _ = run_subprocess_cmd(trimmomatic_cmd, shell_check = False)
        if trimmo_f_pair_path == None and trimmo_r_pair_path == None:
            trimmo_f_pair_path = ILLUMINA_RAW_F_READS
            trimmo_r_pair_path = ILLUMINA_RAW_R_READS
    
        # FastQC Illumina Trimmed Reads
        sample_stats_dict = fastqc_reads(trimmo_f_pair_path, trimmo_r_pair_path, "Trimmed_Illumina_", CPU_THREADS, sample_stats_dict)
                
        # Run bbduk on the trimmed files then clumpify
        bbduk_f_map_path, bbduk_r_map_path = bbduk_map(trimmo_f_pair_path, trimmo_r_pair_path)
        clump_f_dedup_path, clump_r_dedup_path = clumpify_dedup(bbduk_f_map_path, bbduk_r_map_path)
        
        # FastQC Illumina Deduplicated Reads
        sample_stats_dict = fastqc_reads(clump_f_dedup_path, clump_r_dedup_path, "Deduplicated_Illumina_", CPU_THREADS, sample_stats_dict)

    # NanoPlot ONT Raw Reads    
    if not pd.isna(ONT_RAW_READS):
        # NanoPlot ONT Raw Reads
        sample_stats_dict = nanoplot_qc_reads(ONT_RAW_READS, "Raw_ONT_", CPU_THREADS, sample_stats_dict)
    
        # Filtlong ONT Raw Reads
        filtered_ONT_reads = ONT_RAW_READS.replace(".fastq.gz","_filtered.fastq")
        gzipped_filtered_ONT_reads = filtered_ONT_reads.replace(".fasta",".fasta.gz")
        if os.path.exists(gzipped_filtered_ONT_reads):
            log_print(f"SKIP:\tGzipped FiltLong output already exists: {gzipped_filtered_ONT_reads}.")
        else:
            if os.path.exits(filtered_ONT_reads):
                log_print(f"SKIP:\tFiltLong output already exists: {filtered_ONT_reads}.")
            else:
                filtlong_cmd = f"filtlong --min_length 1000 --min_mean_q 8 --target_bases 500000000 --trim -1 {clump_f_dedup_path} -2 {clump_r_dedup_path} {ONT_RAW_READS} > {filtered_ONT_reads}"
                _ = run_subprocess_cmd(filtlong_cmd, shell_check = True)
                
        # NanoPlot ONT Filtered Reads
        sample_stats_dict = nanoplot_qc_reads(filtered_ONT_reads, "Filt_ONT_", CPU_THREADS, sample_stats_dict)        
        gzip_file(filtered_ONT_reads, gzipped_filtered_ONT_reads)
    
        # Error Correct ONT Filtered Reads With Illumina Reads
        corrected_ONT_out = gzipped_filtered_ONT_reads.replace("filtered.fastq.gz","corrected")
        corrected_ONT_Reads = corrected_ONT_out + ".fastq.gz"
        if os.path.exists(corrected_ONT_Reads):
            log_print(f"SKIP:\tRatatosk Corrected Reads already exists: {corrected_ONT_Reads}.")
        else:
            ratatosk_cmd = ["Ratatosk", "correct",
                            "-s", clump_f_dedup_path, "-s", clump_r_dedup_path,
                            "-l", gzipped_filtered_ONT_reads,
                            "-o", corrected_ONT_out,
                            "-c", str(CPU_THREADS), "-G", "-v"]
            _ = run_subprocess_cmd(ratatosk_cmd, shell_check = False)
        
        # NanoPlot ONT Corrected Reads
        sample_stats_dict = nanoplot_qc_reads(corrected_ONT_Reads, "Corrected_ONT", CPU_THREADS, sample_stats_dict)
    
        log_print("Selecting Highest Mean Quality ONT reads...")
        highest_mean_qual_ont_reads = corrected_ONT_Reads
        if sample_stats_dict["CORRECT_ONT_MEAN_QUAL"] < sample_stats_dict["FILT_ONT_MEAN_QUAL"]:
            highest_mean_qual_ont_reads = gzipped_filtered_ONT_reads
        log_print(f"Highest Mean Quality ONT reads: {highest_mean_qual_ont_reads}")
        
    if pd.notna(PACBIO_RAW_READS):
        # NanoPlot QC Raw PacBio Reads
        sample_stats_dict = nanoplot_qc_reads(PACBIO_RAW_READS, "RawPacBio", CPU_THREADS, sample_stats_dict)

        # TODO: Determine the files necessary to use hifiasm and get hifi reads from the PACBIO_RAW_READS
        hifi_reads = PACBIO_RAW_READS
       
        # Filtlong PacBio Reads      
        filtered_PACBIO_reads = PACBIO_RAW_READS.replace(".fastq.gz","_filtered.fastq")
        gzipped_filtered_PACBIO_reads = filtered_PACBIO_reads.replace(".fastq",".fastq.gz")
        if os.path.exists(gzipped_filtered_PACBIO_reads):
            log_print(f"SKIP:\tGzipped FiltLong output already exists: {gzipped_filtered_PACBIO_reads}.")
        else:   
            if os.path.exists(filtered_PACBIO_reads):
                log_print(f"SKIP:\t FiltLong output already exists: {filtered_PACBIO_reads}.")
            else:
                coverage = 300 # can be changed later if needed
                target_bases = EST_SIZE * coverage
                filtlong_cmd = f"filtlong --min_length 10000 --keep_percent 90 --target_bases {target_bases} {hifi_reads} > {filtered_PACBIO_reads}"
                _ = run_subprocess_cmd(filtlong_cmd, shell_check=True)
                
        # NanoPlot Flitered Reads
        sample_stats_dict = nanoplot_qc_reads(filtered_PACBIO_reads, "FiltPacBio", CPU_THREADS, sample_stats_dict)
        gzip_file(filtered_PACBIO_reads, gzipped_filtered_PACBIO_reads)
        
###############################################################################
    # Assembly
###############################################################################

    # MaSuRCA de Novo Assembly based on input Reads     
    masurca_out_dir = os.path.join(shared_root, "masurca_assembly")
    if not os.path.exists(masurca_out_dir):
        os.makedirs(masurca_out_dir)
    final_masurca_path = os.path.join(shared_root, f"{SPECIES_ID}_masurca.fasta")
    if os.path.exists(final_masurca_path):
        log_print(f"SKIP:\tFinal Masurca Assembly already exists: {final_masurca_path}.")
    else:
        if pd.notna(ONT_RAW_READS) and pd.isna(PACBIO_RAW_READS):
            default_assembly_path, assembly_path, ref_assembly_path = masurca_config_gen(shared_root, masurca_out_dir,
                                                                                      highest_mean_qual_ont_reads,
                                                                                      ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS,
                                                                                      PACBIO_RAW_READS,
                                                                                      clump_f_dedup_path, clump_r_dedup_path,
                                                                                      EST_SIZE, CPU_THREADS, RAM_GB, REF_SEQ)
        elif pd.notna(PACBIO_RAW_READS) and pd.isna(ONT_RAW_READS):
            default_assembly_path, assembly_path, ref_assembly_path = masurca_config_gen(shared_root, masurca_out_dir,
                                                                                      ONT_RAW_READS,
                                                                                      ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS,
                                                                                      gzipped_filtered_PACBIO_reads,
                                                                                      clump_f_dedup_path, clump_r_dedup_path,
                                                                                      EST_SIZE, CPU_THREADS, RAM_GB, REF_SEQ)

        else:            
            default_assembly_path, assembly_path, ref_assembly_path = masurca_config_gen(shared_root, masurca_out_dir,
                                                                                      ONT_RAW_READS, ILLUMINA_RAW_F_READS,
                                                                                      ILLUMINA_RAW_R_READS, PACBIO_RAW_READS,
                                                                                      clump_f_dedup_path, clump_r_dedup_path,
                                                                                      EST_SIZE, CPU_THREADS, RAM_GB, REF_SEQ)       
        gap_assembly = os.path.join(find_ca_folder(masurca_out_dir), "9-terminator", "genome.scf.fasta")
        shutil.move(gap_assembly, final_masurca_path)
    
    final_masurca_path, masurca_stats_list, _ = qc_assembly(final_masurca_path, shared_root, cwd,
                                                            ONT_RAW_READS,
                                                            ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS,
                                                            PACBIO_RAW_READS, SPECIES_ID,
                                                            first_compleasm_odb, second_compleasm_odb,
                                                            REF_SEQ, karyote_id, kingdom_id,
                                                            sample_stats_dict, results_df)
    os.chdir(shared_root)
    kmer_list = ["21", "33", "55", "77", "99", "127"]        
    if pd.notna(ONT_RAW_READS) or pd.notna(PACBIO_RAW_READS):
        flye_out_dir = os.path.join(shared_root, "flye_assembly")
        if not os.path.exists(flye_out_dir):
            os.makedirs(flye_out_dir)
        flye_path = os.path.join(flye_out_dir, "assembly.fasta")
        final_flye_path = os.path.join("/".join(flye_path.split("/")[:-1]), f"{SPECIES_ID}_flye.fasta")
        if os.path.exists(final_flye_path):
            log_print(f"SKIP:\tFinal Flye Assembly already exists: {final_flye_path}.")
        else:
            if not os.path.exists(flye_path):
                if pd.notna(ONT_RAW_READS):
                    flye_cmd = ["flye",
                                "--nano-corr", highest_mean_qual_ont_reads,
                                "--out-dir", flye_out_dir, 
                                "--genome-size", EST_SIZE, "--threads", str(CPU_THREADS),
                                "--iterations", "3", "--keep-haplotypes"]
                elif pd.notna(PACBIO_RAW_READS):
                    flye_cmd = ["flye",
                                "--pacbio-corr", gzipped_filtered_PACBIO_reads,
                                "--out-dir", flye_out_dir, 
                                "--genome-size", EST_SIZE, "--threads", str(CPU_THREADS),
                                "--iterations", "3", "--keep-haplotypes"]
                _ = run_subprocess_cmd(flye_cmd, shell_check = False)
            shutil.move(flye_path, final_flye_path)
        final_flye_path, flye_stats_list, _ = qc_assembly(final_flye_path, shared_root, cwd,
                                                          ONT_RAW_READS,
                                                          ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS,
                                                          PACBIO_RAW_READS, SPECIES_ID,
                                                          first_compleasm_odb, second_compleasm_odb,
                                                          REF_SEQ, karyote_id, kingdom_id,
                                                          sample_stats_dict, results_df)
    else:
        final_flye_path, flye_stats_list = None, [0,0,0,100000]

    spades_out_dir = os.path.join(shared_root, "spades_assembly")
    if not os.path.exists(spades_out_dir):
        os.makedirs(spades_out_dir)
    spades_path = os.path.join(spades_out_dir, "scaffolds.fasta")
    final_spades_path = os.path.join("/".join(spades_path.split("/")[:-1]), f"{SPECIES_ID}_spades.fasta")
    if os.path.exists(final_spades_path):
        log_print(f"SKIP:\tFinal SPAdes Assembly already exists: {final_spades_path}.")
    else:   
        if not os.path.exists(spades_path):
            if pd.notna(ONT_RAW_READS):
                spades_cmd = ["spades.py",
                              "-1", clump_f_dedup_path,
                              "-2", clump_r_dedup_path,
                              "--nanopore", highest_mean_qual_ont_reads,
                              "-o", spades_out_dir,
                              "-t", str(CPU_THREADS),
                              "-m", str(RAM_GB),
                              "--careful", "--cov-cutoff", "auto",
                              "-k", ",".join(kmer_list)]
            elif pd.notna(PACBIO_RAW_READS): # TODO ADJUST FOR PACBIO
                spades_cmd = ["spades.py",
                              "--careful",
                              "-1", clump_f_dedup_path,
                              "-2", clump_r_dedup_path,
                              "-o", spades_out_dir,
                              "-t", str(CPU_THREADS),
                              "-m", str(RAM_GB),
                              "--cov-cutoff", "auto",
                              "-k", ",".join(kmer_list)]
            else:
                spades_cmd = ["spades.py",
                              "--careful",
                              "-1", clump_f_dedup_path,
                              "-2", clump_r_dedup_path,
                              "-o", spades_out_dir,
                              "-t", str(CPU_THREADS),
                              "-m", str(RAM_GB),
                              "--cov-cutoff", "auto",
                              "-k", ",".join(kmer_list)]
            if not pd.isna(REF_SEQ):
                spades_cmd.append(f"--trusted-contigs {REF_SEQ}")
            _ = run_subprocess_cmd(spades_cmd, shell_check = False)
        shutil.move(spades_path, final_spades_path)
    final_spades_path, spades_stats_list, _ = qc_assembly(final_spades_path, shared_root, cwd,
                                                          ONT_RAW_READS,
                                                          ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS,
                                                          PACBIO_RAW_READS, SPECIES_ID,
                                                          first_compleasm_odb, second_compleasm_odb,
                                                          REF_SEQ, karyote_id, kingdom_id,
                                                          sample_stats_dict, results_df)
    
    if pd.notna(PACBIO_RAW_READS):
        hifiasm_out_dir = os.path.join(shared_root, "hifiasm_assembly")
        if not os.path.exists(hifiasm_out_dir):
            os.mkdirs(hifiasm_out_dir)
        pacbio_prefix = os.path.join(hifiasm_out_dir, os.path.basename(PACBIO_RAW_READS.split(".")[0]))
        hifiasm_path = pacbio_prefix + ".asm"
        hifiasm_gfa_path = pacbio_prefix + ".p_ctg.gfa"
        final_hifiasm_path = os.path.join("/".join(hifiasm_path.split("/")[:-1]), f"{SPECIES_ID}_hifiasm.fasta")
        if os.path.exists(final_hifiasm_path):
            log_print(f"SKIP:\tFinal HiFi Assembly already exists: {final_hifiasm_path}.")            
        else:
            if os.path.exists(hifiasm_path):
                log_print(f"SKIP:\tHiFi Assembly already exists: {hifiasm_path}.")
            else:
                hifhiasm_cmd = ["hifiasm", "-o", hifiasm_path, "-t", str(CPU_THREADS), gzipped_filtered_PACBIO_reads]
                _ = run_subprocess_cmd(hifhiasm_cmd, shell_check=False)
            gfa_cmd = f"gfatools gfa2a -s {hifiasm_gfa_path} > {final_hifiasm_path}"
            _ = run_subprocess_cmd(gfa_cmd, shell_check=False)            
        final_hifiasm_path, hifiasm_stats_list, _ = qc_assembly(final_hifiasm_path, shared_root, cwd,
                                                                ONT_RAW_READS,
                                                                ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS,
                                                                PACBIO_RAW_READS, SPECIES_ID,
                                                                first_compleasm_odb, second_compleasm_odb,
                                                                REF_SEQ, karyote_id, kingdom_id,
                                                                sample_stats_dict, results_df)
    else:
        final_hifiasm_path, hifiasm_stats_list = None, [0,0,0,100000]

    stats_combined = zip(masurca_stats_list, flye_stats_list, spades_stats_list, hifiasm_stats_list)
    custom_stats = []
    for index, values in enumerate(stats_combined):
        if index == len(masurca_stats_list) - 1:
            min_value = min(values)
            min_index = values.index(min_value)
            methods = ["MaSuRCA", "Flye", "SPAdes", "hifiasm"]
            custom_stats.append(methods[min_index])
        else:  # For all other items, find the largest
            max_value = max(values)
            max_index = values.index(max_value)
            methods = ["MaSuRCA", "Flye", "SPAdes", "hifiasm"]
            custom_stats.append(methods[max_index])
    method_counts = {}
    for method in custom_stats:
        method_counts[method] = method_counts.get(method, 0) + 1
    shared_root
    most_represented_method = max(method_counts, key=method_counts.get)
    comparative_plot_dict = {"MaSuRCA": masurca_stats_list,
                             "Flye": flye_stats_list,
                             "SPAdes": spades_stats_list,
                             "hifiasm": hifiasm_stats_list}
    log_print(f"Best Initial Assembly: {most_represented_method}.")
    if pd.isna(final_flye_path):
        common_path = os.path.commonpath([final_masurca_path, final_spades_path])
    else:
        common_path = os.path.commonpath([final_masurca_path, final_flye_path, final_spades_path])
    comparative_plots(comparative_plot_dict, common_path)
    if most_represented_method == "Flye":
        assembly_out_dir = flye_out_dir
        de_novo_assembly = final_flye_path
    elif most_represented_method == "SPAdes":
        assembly_out_dir = spades_out_dir
        de_novo_assembly = final_spades_path
    elif most_represented_method == "hifiasm":
        assembly_out_dir = hifiasm_out_dir
        de_novo_assembly = final_hifiasm_path
    else:
        assembly_out_dir = masurca_out_dir
        de_novo_assembly = final_masurca_path
    os.chdir(shared_root)
    
###############################################################################
    # Assembly Polishing
###############################################################################

    # 2x Polish with Racon
    if not pd.isna(ONT_RAW_READS):
        first_racon_assembly = racon_polish_assembly(de_novo_assembly, highest_mean_qual_ont_reads, assembly_out_dir, 1)
        pilon_polish_target = racon_polish_assembly(first_racon_assembly, highest_mean_qual_ont_reads, assembly_out_dir, 2)
    else:
        log_print("SKIP:\t2x-Racon Polish as no ONT reads were provided.")
        pilon_polish_target = de_novo_assembly
    
    # Pilon Polish Preparation
    polish_bam = pilon_polish_target.replace(".fasta",".bam")# os.path.join(assembly_out_dir, "second_polish_assembly_illumina_sorted.bam")
    output_sam = polish_bam.replace("bam","sam")
    sorted_bam = polish_bam.replace(".bam","_sorted.bam")
    if not os.path.exists(output_sam):
        bwa_index_cmd = ["bwa", "index", pilon_polish_target]
        _ = run_subprocess_cmd(bwa_index_cmd, shell_check = False)
        bwa_cmd = f"bwa mem {pilon_polish_target} {clump_f_dedup_path} {clump_r_dedup_path} > {output_sam}"
        _ = run_subprocess_cmd(bwa_cmd, shell_check = True)
    if not os.path.exists(sorted_bam):
        samview_cmd = f"samtools view -@ {CPU_THREADS} -S -b {output_sam} > {sorted_bam}"
        _ = run_subprocess_cmd(samview_cmd, shell_check = True)
    if not os.path.exists(polish_bam):    
        bamsort_cmd = ["bamtools", "sort", "-in", sorted_bam, "-out", polish_bam]
        _ = run_subprocess_cmd(bamsort_cmd, shell_check = False)
        samtools_index_cmd = ["samtools", "index", polish_bam]
        _ = run_subprocess_cmd(samtools_index_cmd, shell_check = False)   

    # Pilon Polish with Illumina Reads
    pilon_out_prefix = "pilon_assembly"
    pilon_out_dir = os.path.join(shared_root, "pilon_polished_assembly")
    pilon_out_fasta = os.path.join(pilon_out_dir, "pilon_assembly.fasta")
    pilon_renamed_fasta = de_novo_assembly.replace(".fasta","_pilon.fasta")
    if os.path.exists(pilon_out_fasta):
        shutil.move(pilon_out_fasta, pilon_renamed_fasta)
        log_print(f"SKIP:\tPilon Polished Assembly already exists: {pilon_out_fasta}.")
    elif os.path.exists(pilon_renamed_fasta):
        log_print(f"SKIP:\tPilon Polished Assembly already exists: {pilon_renamed_fasta}.")
    else:
        pilon_cmd = ["pilon", f"-Xmx{RAM_GB}g",
                     "--genome", pilon_polish_target,
                     "--frags", polish_bam,
                     "--output", pilon_out_prefix,
                     "--outdir", pilon_out_dir,
                     "--changes", "--vcf",
                     "--chunksize", str(5000000)]
        _ = run_subprocess_cmd(pilon_cmd, shell_check = False)
        shutil.move(pilon_out_fasta, pilon_renamed_fasta)
    pilon_out_fasta = pilon_renamed_fasta

###############################################################################
    # Assembly Curation
###############################################################################    
    
    # Purge haplotigs and overlaps with purge_dups (if ONT Reads exist)
    if not pd.isna(ONT_RAW_READS):
        # Define paths and variables
        pd_work_dir = os.path.join(shared_root, "purge_dups_work")
        if not os.path.exists(pd_work_dir):
            os.mkdir(pd_work_dir)
        os.chdir(pd_work_dir)
        
        # Create the FOFN file for PacBio reads
        pd_fofn = os.path.join(pd_work_dir, "ont_reads.fofn")
        with open(pd_fofn, "w") as f:
            f.write(highest_mean_qual_ont_reads)
        
        # Define paths for output files
        pd_json = os.path.join(pd_work_dir, "purge_dups_config.json")
        dup_purged_assembly = os.path.join(pd_work_dir, f"{SPECIES_ID}_{most_represented_method.lower()}_pilon/seqs/{SPECIES_ID}_{most_represented_method.lower()}_pilon.purged.fa")
        updated_dup_purged_assembly = dup_purged_assembly.replace(".fa", ".fasta")
        pd_config_path = find_file("pd_config.py")
        pd_config_dir = os.path.dirname(os.path.dirname(pd_config_path))
        pd_path = os.path.join(pd_config_dir, "bin")
      
        # Step 1: Generate the purge_dups configuration file
        if os.path.exists(pd_json):
            log_print(f"SKIP:\tPurge Dupes JSON already exists: {pd_json}.")
            log_print(f"NOTE:\tCurrent pd_path: {pd_path}.")            
        else:
            purge_dupes_config_cmd = ["python", pd_config_path,
                                      pilon_out_fasta,
                                      pd_fofn,
                                      "-l", pd_work_dir,
                                      "-n", pd_json]
            _ = run_subprocess_cmd(purge_dupes_config_cmd, shell_check=False)
        
        # Step 2: Run the purge_dups pipeline        
        if os.path.exists(dup_purged_assembly):
            shutil.move(dup_purged_assembly, updated_dup_purged_assembly)
            log_print(f"SKIP:\tDupe Purged Assembly already exists: {updated_dup_purged_assembly}.")
        elif os.path.exists(updated_dup_purged_assembly):
            log_print(f"SKIP:\tDupe Purged Assembly already exists: {updated_dup_purged_assembly}.")
        else:
            # Run the purge_dups pipeline
            purge_dupes_cmd = ["python", find_file("run_purge_dups.py"),
                               pd_json,
                               pd_path,
                               SPECIES_ID,
                               "-p", "bash"]
            _ = run_subprocess_cmd(purge_dupes_cmd, shell_check=False)
        
            # Move the output file to the updated path
            if os.path.exists(dup_purged_assembly):
                shutil.move(dup_purged_assembly, updated_dup_purged_assembly)
        
        # Restore the original working directory
        os.chdir(cwd)
        dup_purged_assembly = updated_dup_purged_assembly
    else:
        dup_purged_assembly = pilon_out_fasta

    # RagTag Correct Reads (if REF_SEQ exists)
    ragtag_ref_assembly = os.path.join("/".join(os.path.dirname(dup_purged_assembly).split("/")[:-1]), os.path.basename(dup_purged_assembly).replace(".fasta","_ragtag_final.fasta"))
    ragtag_patched = ragtag_ref_assembly.replace("_pilon_ragtag_final.fasta","_patched")
    patch_fasta = os.path.join(ragtag_patched, "ragtag.patch.fasta")
    ragtag_scaff = ragtag_patched.replace("patched","scaffold")
    scaff_fasta = os.path.join(ragtag_scaff, "ragtag.scaffold.fasta")
    ragtag_corr = ragtag_patched.replace("patched","corrected")
    corr_fasta = os.path.join(ragtag_corr, "ragtag.correct.fasta")
    if not pd.isna(REF_SEQ):
        if not os.path.exists(ragtag_patched):
            os.mkdir(ragtag_patched)
        if not os.path.exists(ragtag_scaff):
            os.mkdir(ragtag_scaff)
        if not os.path.exists(ragtag_corr):
            os.mkdir(ragtag_corr)
        if os.path.exists(ragtag_ref_assembly):
            log_print(f"SKIP:\tRagTag Reference Corrected Assembly already exists: {ragtag_ref_assembly}")
        else:
            if os.path.exists(scaff_fasta):
                log_print(f"SKIP:\tRagTag Scaffolded Assembly already exists: {scaff_fasta}")
            else:
                ragscaf_cmd = ["ragtag.py", "scaffold", REF_SEQ, dup_purged_assembly, "-o", ragtag_scaff, "-t", str(CPU_THREADS), "-C", "-u"]
                _ = run_subprocess_cmd(ragscaf_cmd, shell_check = False)
            if os.path.exists(corr_fasta):
                log_print(f"SKIP:\tRagTag Corrected Assembly already exists: {corr_fasta}")
            else:
                ragcorr_cmd = ["ragtag.py", "correct", REF_SEQ, scaff_fasta, "-o", ragtag_corr, "-t", str(CPU_THREADS), "-u"]
                _ = run_subprocess_cmd(ragcorr_cmd, shell_check = False)
            ragpatch_cmd = ["ragtag.py", "patch", REF_SEQ, corr_fasta, "-o", ragtag_patched, "-t", str(CPU_THREADS), "-u"]
            _ = run_subprocess_cmd(ragpatch_cmd, shell_check = False)         
        
            shutil.copyfile(patch_fasta, ragtag_ref_assembly)
    else:
        ragtag_ref_assembly = dup_purged_assembly

    # Close Gaps with TGS-GapCloser (if ONT Reads exist) or Abyss-Sealer
    if not pd.isna(ONT_RAW_READS):
        tgs_gapcloser_dir = os.path.join(shared_root, "tgs_gapcloser")
        if not os.path.exists(tgs_gapcloser_dir):
            os.makedirs(tgs_gapcloser_dir)
        os.chdir(tgs_gapcloser_dir)
        tgs_gapcloser_prefix =  os.path.join(tgs_gapcloser_dir, "purged_gapclosed")
        tgs_gapcloser_output = os.path.join(tgs_gapcloser_dir, "purged_gapclosed.scaff_seqs")
        final_assembly_path = os.path.join(shared_root, f"{SPECIES_ID}_EGAP_assembly.fasta")
        closed_gaps_path = os.path.join(tgs_gapcloser_dir, "purged_gapclosed.scaff_seqs")
        if os.path.exists(tgs_gapcloser_output):
            log_print(f"SKIP:\tTGS Gapcloser Assembly already exists: {tgs_gapcloser_output}.")
        else:
            tgs_gapcloser_cmd = ["tgsgapcloser", "--scaff", ragtag_ref_assembly,
                                 "--reads", highest_mean_qual_ont_reads,
                                 "--output", tgs_gapcloser_prefix,
                                 "--tgstype", "ont", "--thread", str(CPU_THREADS), "--ne"]
            _ = run_subprocess_cmd(tgs_gapcloser_cmd, shell_check = False)
        os.chdir(cwd)
        if os.path.exists(final_assembly_path):
            log_print(f"SKIP:\tFinal Assembly already exists: {final_assembly_path}.")
        else:        
            shutil.copyfile(closed_gaps_path, final_assembly_path)
    else:
        output_prefix = "/".join(shared_root.split("/")[:-1]) + f"/{ragtag_ref_assembly.split('/')[-1].replace('_pilon_ragtag_final.fasta','_sealed').replace('_pilon.fasta','_sealed')}"
        sealer_output_file = f"{output_prefix}_scaffold.fa"        
        renamed_sealer_output_file = sealer_output_file.replace(".fa",".fasta")
        if os.path.isfile(sealer_output_file):        
            log_print(f"SKIP:\tABySS sealer output file already exists: {sealer_output_file}.")
        else:
            if not os.path.exists(renamed_sealer_output_file):
                kmer_sizes = [55,75,95]
                abyss_sealer_cmd = ["abyss-sealer", "-o", output_prefix, "-S", ragtag_ref_assembly,
                                    "-L", str(400), "-G", str(1000),
                                    "-j", str(CPU_THREADS), "-b", "500M"]            
                for k in kmer_sizes:
                    abyss_sealer_cmd.extend(["-k", str(k)])
                abyss_sealer_cmd.extend([clump_f_dedup_path, clump_r_dedup_path])
                _ = run_subprocess_cmd(abyss_sealer_cmd, shell_check = False)
                shutil.move(sealer_output_file, renamed_sealer_output_file)                
        final_assembly_path = renamed_sealer_output_file
    final_gz_assembly_path = final_assembly_path + ".gz"
    gzip_file(final_assembly_path, final_gz_assembly_path)

###############################################################################
    # Final Assembly Assessment
###############################################################################
    final_assembly_path, results_df = process_final_assembly(row, results_df,
                                                             CPU_THREADS, RAM_GB,
                                                             final_assembly_path,
                                                             sample_stats_dict)    
    return final_assembly_path, results_df


if __name__ == "__main__":
    # Argument Parsing & Test Data
    parser = argparse.ArgumentParser(description="Run Entheome Genome Assembly Pipeline (EGAP)")
    default_input_csv = None
    default_ont_sra = None
    default_raw_ont_dir = None
    default_ont_reads = None
    default_illu_sra = None
    default_raw_illu_dir = None
    default_raw_illu_reads_1 = None
    default_raw_illu_reads_2 = None
    default_pacbio_sra = None
    default_raw_pacbio_dir = None
    default_raw_pacbio_reads = None
    default_species_id = None
    default_organism_kingdom = None
    default_organism_karyote = None
    default_compleasm_1 = None
    default_compleasm_2 = None
    default_estimated_genome_size = None
    default_reference_sequence = None
    default_reference_sequence_gca = None
    default_percent_resources = 0.75

    parser.add_argument("--input_csv", "-csv", type=str, default=default_input_csv,
                        help=f"Path to a csv containing multiple sample data. (default: {default_input_csv})")
    parser.add_argument("--ont_sra", "-osra", type=str, default=default_ont_sra,
                        help=f"Oxford Nanopore SRA Accession number. (default: {default_ont_sra})")
    parser.add_argument("--raw_ont_dir", "-odir", type=str, default=default_raw_ont_dir,
                        help=f"Directory containing Raw ONT Reads. (default: {default_raw_ont_dir})")
    parser.add_argument("--raw_ont_reads", "-i0", type=str, default=default_ont_reads,
                        help=f"Path to combined Raw ONT fastq reads. (default: {default_ont_reads})")
    parser.add_argument("--illu_sra", "-isra", type=str, default=default_illu_sra,
                        help=f"Illumina SRA Accession number (paired-end). (default: {default_illu_sra})")
    parser.add_argument("--raw_illu_dir", "-idir", type=str, default=default_raw_illu_dir,
                        help=f"Directory containing Raw Illumina Reads. (default: {default_raw_illu_dir})")
    parser.add_argument("--raw_illu_reads_1", "-i1", type=str, default=default_raw_illu_reads_1,
                        help=f"Path to Raw Forward Illumina Reads. (default: {default_raw_illu_reads_1})")
    parser.add_argument("--raw_illu_reads_2", "-i2", type=str, default=default_raw_illu_reads_2,
                        help=f"Path to Raw Reverse Illumina Reads. (default: {default_raw_illu_reads_2})")
    parser.add_argument("--pacbio_sra", "-pbsra", type=str, default=default_pacbio_sra,
                        help=f"PacBio SRA Accession number. (default: {default_pacbio_sra})")
    parser.add_argument("--raw_pacbio_dir", "-pbdir", type=str, default=default_raw_pacbio_dir,
                        help=f"Directory containing Raw PacBio Reads. (default: {default_raw_pacbio_dir})")
    parser.add_argument("--raw_pacbio_reads", "-pb", type=str, default=default_raw_pacbio_reads,
                        help=f"Path to Raw PacBio Reads. (default: {default_raw_pacbio_reads})")
    parser.add_argument("--species_id", "-ID", type=str, default=default_species_id,
                        help=f"Species ID (format: <2-letters of Genus>_<full species name>). (default: {default_species_id})")
    parser.add_argument("--organism_kingdom", "-Kg", type=str, default=default_organism_kingdom,
                        help=f"Organism Kingdom. (default: {default_organism_kingdom})")
    parser.add_argument("--organism_karyote", "-Ka", type=str, default=default_organism_karyote,
                        help=f"Organism Karyote type. (default: {default_organism_karyote})")
    parser.add_argument("--compleasm_1", "-c1", type=str, default=default_compleasm_1,
                        help=f"First Compleasm/BUSCO database. (default: {default_compleasm_1})")
    parser.add_argument("--compleasm_2", "-c2", type=str, default=default_compleasm_2,
                        help=f"Second Compleasm/BUSCO database. (default: {default_compleasm_2})")
    parser.add_argument("--est_size", "-es", type=str, default=default_estimated_genome_size,
                        help=f"Estimated genome size in Mbp. (default: {default_estimated_genome_size})")
    parser.add_argument("--ref_seq_gca", "-rgca", type=str, default=default_reference_sequence_gca,
                        help=f"Reference Genome Assembly (GCA) Accession number. (default: {default_reference_sequence_gca})")
    parser.add_argument("--ref_seq", "-rf", type=str, default=default_reference_sequence,
                        help=f"Path to the reference genome. (default: {default_reference_sequence})")
    parser.add_argument("--percent_resources", "-R", type=float, default=default_percent_resources,
                        help=f"Fraction of system resources to use. (default: {default_percent_resources})")
    
    args = parser.parse_args()
    INPUT_CSV = args.input_csv
    if INPUT_CSV is not None:
        input_csv_df = pd.read_csv(INPUT_CSV)
    else:
        sample_dict = {"ONT_SRA": [args.ont_sra],
                       "ONT_RAW_DIR": [args.raw_ont_dir],
                       "ONT_RAW_READS": [args.raw_ont_reads],
                       "ILLUMINA_SRA": [args.illu_sra],
                       "ILLUMINA_RAW_DIR": [args.raw_illu_dir],
                       "ILLUMINA_RAW_F_READS": [args.raw_illu_reads_1],
                       "ILLUMINA_RAW_R_READS": [args.raw_illu_reads_2],
                       "PACBIO_SRA": [args.ont_sra],
                       "PACBIO_RAW_DIR": [args.raw_ont_dir],
                       "PACBIO_RAW_READS": [args.raw_ont_reads],
                       "SPECIES_ID": [args.species_id],
                       "ORGANISM_KINGDOM": [args.organism_kingdom],
                       "ORGANISM_KARYOTE": [args.organism_karyote],
                       "COMPLEASM_1": [args.compleasm_1],
                       "COMPLEASM_2": [args.compleasm_2],
                       "EST_SIZE": [args.est_size],
                       "REF_SEQ_GCA": [args.ref_seq_gca],
                       "REF_SEQ": [args.ref_seq]}
        input_csv_df = pd.DataFrame.from_dict(sample_dict)
    PERCENT_RESOURCES = args.percent_resources
    CPU_THREADS, RAM_GB = get_resource_values(PERCENT_RESOURCES)

    results_df = pd.DataFrame()
    for index, row in input_csv_df.iterrows():
        if pd.isna(row["ONT_SRA"]) and pd.isna(row["ILLUMINA_SRA"]) and pd.isna(row["EST_SIZE"]) and not pd.isna(row["REF_SEQ_GCA"]):
            final_assembly_path, results_df = process_final_assembly(row, results_df, CPU_THREADS, RAM_GB)
        else:
            final_assembly_path, results_df = egap_sample(row, results_df, INPUT_CSV, CPU_THREADS, RAM_GB)
    
    if INPUT_CSV is not None:
        final_csv_filename = INPUT_CSV.replace(".csv", "_final_assembly_stats.csv")
    else:
        final_csv_filename = input_csv_df.iloc[1]["FINAL_ASSEMBLY"].replace(".fasta", "_stats.csv")
    input_csv_df.to_csv(final_csv_filename, index=False)
