#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
utilities.py

Module containing regularly used commands in various other EGAP scripts.

Created on Tue Apr  8 22:13:08 2025

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os, subprocess, datetime, platform, shutil
import pandas as pd


# --------------------------------------------------------------
# Execute a subprocess command and log its output
# --------------------------------------------------------------
def run_subprocess_cmd(cmd_list, shell_check):
    """Execute a subprocess command and log its execution.

    Runs the command (as a string or list) using subprocess.Popen, captures and streams
    its output in real-time, and logs success or failure.

    Args:
        cmd_list (str or list): Command to execute, as a string or list of arguments.
        shell_check (bool): If True, execute the command through the shell.

    Returns:
        int: The subprocess return code.
    """
    if isinstance(cmd_list, str):
        print(f"CMD:\t{cmd_list}")
        process = subprocess.Popen(cmd_list, shell=shell_check, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True)
    else:
        print(f"CMD:\t{' '.join(cmd_list)}")
        process = subprocess.Popen(cmd_list, shell=shell_check, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True)
    for line in process.stdout:
        print(line, end="")
    process.wait()
    if process.returncode != 0:
        print(f"NOTE:\tCommand failed with return code {process.returncode}")
    else:
        print(f"PASS:\tSuccessfully processed command: {' '.join(cmd_list) if isinstance(cmd_list, list) else cmd_list}")
    return process.returncode


# --------------------------------------------------------------
# Create a sample statistics dictionary from metadata
# --------------------------------------------------------------
def gen_sample_stats_dict(row):
    """Generate a sample statistics dictionary from a metadata row.

    Extracts key fields from a pandas Series and initializes placeholders for metrics.

    Args:
        row (pandas.Series): Metadata row containing sample information.

    Returns:
        dict: Dictionary with initialized sample statistics.
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

    Args:
        input_file (str): Path to the file to compress.
        cpu_threads (int): Number of threads for compression.

    Returns:
        str: Path to the compressed .gz file.
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

    Args:
        input_file (str): Path to the .gz file to decompress.
        cpu_threads (int): Number of threads for decompression.

    Returns:
        str: Path to the decompressed file.
    """
    pigz_cmd = f"pigz -p {cpu_threads} -d -f {input_file}"
    _ = run_subprocess_cmd(pigz_cmd, shell_check = True)
    unzip_file = input_file.replace(".gz","")
    return unzip_file


# --------------------------------------------------------------
# Extract and process sample metadata
# --------------------------------------------------------------
def get_current_row_data(input_df, sample_id):
    """Extract row data for a sample ID and generate a stats dictionary.

    Filters a DataFrame for a specific sample ID and creates a statistics dictionary.

    Args:
        input_df (pandas.DataFrame): DataFrame with sample metadata.
        sample_id (str): Sample identifier to filter.

    Returns:
        tuple: (filtered DataFrame row, row index list, sample statistics dictionary).
    """
    # Filter the DataFrame for rows where the "SAMPLE_ID" column equals the provided sample_id
    current_row = input_df[input_df["SAMPLE_ID"] == sample_id]
    sample_stats_dict = gen_sample_stats_dict(current_row)
    current_index = current_row.index.tolist()
    
    return current_row, current_index, sample_stats_dict

# --------------------------------------------------------------
# Parse NanoPlot statistics for long reads
# --------------------------------------------------------------
def analyze_nanostats(READS_ORIGIN, nanoplot_out_file, sample_stats_dict):
    """Parse NanoPlot statistics and update the sample statistics dictionary.

    Reads NanoPlot output and extracts metrics based on read origin (e.g., raw ONT).

    Args:
        READS_ORIGIN (str): Type and stage of reads (e.g., 'Raw_ONT', 'Filt_PacBio').
        nanoplot_out_file (str): Path to NanoPlot statistics file.
        sample_stats_dict (dict): Dictionary to update with statistics.

    Returns:
        dict: Updated sample statistics dictionary.
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
# Search filesystem for a file
# --------------------------------------------------------------
def find_file(filename, folder=None):
    """Search the filesystem for a file by name.

    Walks the filesystem from a root directory, skipping certain paths, to find the
    first occurrence of the file.

    Args:
        filename (str): Name of the file to search for.
        folder (str, optional): Starting directory for search. Defaults to OS-specific root.

    Returns:
        str or None: Absolute path to the file if found, else None.
    """
    global ENVIRONMENT_TYPE
    print(f"Looking for {filename}")
    root_directory = "/"

    for root, dirs, files in os.walk(root_directory):
        # Skip paths containing "$RECYCLE.BIN"
        if "$RECYCLE.BIN" in root:
            continue
        if "ncbi" in dirs:
            continue

        if filename in files:
            return os.path.join(root, filename)

    return None


# --------------------------------------------------------------
# Move a file up the directory tree
# --------------------------------------------------------------
def move_file_up(input_file, up_count):
    """Move a file up the directory hierarchy by a specified number of levels.

    Args:
        input_file (str): Path to the file to move.
        up_count (int): Number of directory levels to move up.

    Returns:
        str: New path to the moved file, or original path if file doesn't exist.
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
    """Select the highest quality long reads from ONT or PacBio data.

    Processes metadata to identify, filter, and select the best long reads based on
    mean quality.

    Args:
        output_dir (str): Directory for output files.
        input_csv (str): Path to metadata CSV file.
        sample_id (str): Sample identifier.
        cpu_threads (int): Number of threads for compression tasks.

    Returns:
        str or None: Path to the selected high-quality reads file, or None if not found.
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
        ont_raw_reads = os.path.join(output_dir, species_id, "ONT", f"{ont_sra}.fastq.gz")
    if pd.notna(pacbio_sra) and pd.isna(pacbio_raw_reads):
        pacbio_raw_reads = os.path.join(output_dir, species_id, "PacBio", f"{pacbio_sra}.fastq.gz")
 
    print(f"DEBUG - ont_raw_reads - {ont_raw_reads}")
    print(f"DEBUG - pacbio_raw_reads - {pacbio_raw_reads}")
    
    if pd.notna(ont_raw_reads):
        print("DEBUG - PROCESSING ONT HIGHEST MEAN QUAL")
        reads_type = "ONT"
        reads_dir = os.path.dirname(ont_raw_reads)
        filtered_reads = os.path.join(reads_dir, f"{species_id}_{reads_type}_filtered.fastq.gz")
        corrected_reads = os.path.join(reads_dir, f"{species_id}_{reads_type}_corrected.fastq.gz")
        reads_origin_list = ["Raw_ONT_", "Filt_ONT_", "Corr_ONT_"]
    elif pd.notna(pacbio_raw_reads):
        print("DEBUG - PROCESSING PACBIO HIGHEST MEAN QUAL")
        reads_type = "PacBio"
        reads_dir = os.path.dirname(pacbio_raw_reads)
        filtered_reads = os.path.join(reads_dir, f"{species_id}_{reads_type}_filtered.fastq.gz")
        corrected_reads = os.path.join(reads_dir, f"{species_id}_{reads_type}_corrected.fastq.gz")
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

    if ".gz" not in highest_mean_qual_long_reads:
        highest_mean_qual_long_reads_gz = pigz_compress(highest_mean_qual_long_reads, cpu_threads)
    else:
        highest_mean_qual_long_reads_gz = highest_mean_qual_long_reads

    renamed_highest_mean_qual_long_reads_gz = f"{species_id}_{reads_type}_highest_mean_qual_long_reads.fastq.gz"
    if not os.path.exists(highest_mean_qual_long_reads):
        # try fallback: see if it's named like "Escherichia_coli_filtered.fastq.gz"
        fallback_file = os.path.join(reads_dir, f"{species_id}_filtered.fastq.gz")
        if os.path.exists(fallback_file):
            print(f"FALLBACK:\tFound fallback filtered file: {fallback_file}")
            highest_mean_qual_long_reads = fallback_file
        else:
            print("ERROR:\tNo usable highest-mean-quality long read file found.")
            return None
    
    if ".gz" not in highest_mean_qual_long_reads:
        highest_mean_qual_long_reads_gz = pigz_compress(highest_mean_qual_long_reads, cpu_threads)
    else:
        highest_mean_qual_long_reads_gz = highest_mean_qual_long_reads
    
    renamed_highest_mean_qual_long_reads_gz = os.path.join(reads_dir, f"{species_id}_{reads_type}_highest_mean_qual_long_reads.fastq.gz")
    shutil.copy(highest_mean_qual_long_reads_gz, renamed_highest_mean_qual_long_reads_gz)
    print(f"NOTE:\tSelected highest quality long reads: {renamed_highest_mean_qual_long_reads_gz} with mean quality {highest_mean_qual}")
    return renamed_highest_mean_qual_long_reads_gz


# --------------------------------------------------------------
# Create or manage a log file
# --------------------------------------------------------------
def generate_log_file(log_file_path, use_numerical_suffix=False):
    """Generate a log file, optionally with a numerical suffix if it exists.

    Args:
        log_file_path (str): Desired path for the log file.
        use_numerical_suffix (bool): If True, append a numerical suffix to avoid overwriting.

    Returns:
        str: Path to the created or selected log file.
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
    """Log a message to a file and print it with colored output.

    Timestamps the message, writes it to a log file, and prints it in a color based
    on message type.

    Args:
        input_message (str): Message to log and print.
        log_file (str, optional): Path to the log file. Defaults to DEFAULT_LOG_FILE.
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
def initialize_logging_environment(INPUT_FOLDER):
    """Initialize the logging environment based on the input folder.

    Sets global logging variables and creates a log file based on the OS and input folder.

    Args:
        INPUT_FOLDER (str): Folder used to determine log file location.
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