# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 19:57:46 2024

@author: ian.bollinger@entheome.org

EGAP (Entheome Genome Assembly Pipeline) is a versatile bioinformatics pipeline
developed for assembling high-quality hybrid genomes using Oxford Nanopore
Technologies (ONT) and Illumina sequencing data. It also supports de novo and
reference-based assemblies using Illumina data alone. The pipeline encompasses
comprehensive steps for read quality control, trimming, genome assembly, polishing,
and scaffolding. While optimized for fungal genomes, EGAP can be customized to
work with other types of organisms.

Modified From:
    Bollinger IM, Singer H, Jacobs J, Tyler M, Scott K, Pauli CS, Miller DR,
    Barlow C, Rockefeller A, Slot JC, Angel-Mosti V. High-quality draft genomes
    of ecologically and geographically diverse Psilocybe species. Microbiol Resour
    Announc 0:e00250-24. https://doi.org/10.1128/mra.00250-24
    
    Muñoz-Barrera A, Rubio-Rodríguez LA, Jáspez D, Corrales A , Marcelino-Rodriguez I,
    Lorenzo-Salazar JM, González-Montelongo R, Flores C. Benchmarking of bioinformatics
    tools for the hybrid de novo assembly of human whole-genome sequencing data.
    bioRxiv 2024.05.28.595812; doi: https://doi.org/10.1101/2024.05.28.595812 

CLI EXAMPLE: python EGAP.py ...

    Parameters: 
        --input_csv, -csv (str): Path to a csv containing multiple sample data. (default = None)
        --raw_ont_dir, -odir (str): Path to a directory containing all Raw ONT Reads. (if -csv = None; else REQUIRED)
        --raw_ont_reads, -i0 (str): Path to the combined Raw ONT fastq reads. (if -csv = None; else REQUIRED)
        --raw_illu_dir, -idir (str): Path to a directory containing all Raw Illumina Reads. (if -csv = None; else REQUIRED)
        --raw_illu_reads_1, -i1 (str): Path to the Raw Forward Illumina Reads. (if -csv = None; else REQUIRED)
        --raw_illu_reads_2, -i2 (str): Path to the Raw Reverse Illumina Reads. (if -csv = None; else REQUIRED)
        --species_id, -ID (str): Species ID formatted: <2-letters of Genus>_<full species name>. (if -csv = None; else REQUIRED)
        --organism_kingdom, -Kg (str): Phylogenetic Kingdom the current organism data belongs to. (default: Funga)
        --organism_karyote, -Ka (str): Karyote type of the organism. (default: Eukaryote)
        --compleasm_1, -c1 (str): Name of the first organism compleasm/BUSCO database to compare to. (default: basidiomycota)
        --compleasm_2, -c2 (str): Name of the second organism compleasm/BUSCO database to compare to. (default: agaricales)
        --est_size, -es (str): Estimaged size of the genome in Mbp (aka million-base-pairs). (default: 60m)
        --ref_seq, -rf (str): Path to the reference genome for assembly. (default: None)
        --percent_resources, -R (float): Percentage of resources for processing. (default: 1.00)

CSV EXAMPLE:
    
    ONT_RAW_DIR,ONT_RAW_READS,ILLUMINA_RAW_DIR,ILLUMINA_RAW_F_READS,ILLUMINA_RAW_R_READS,SPECIES_ID,ORGANISM_KINGDOM,ORGANISM_KARYOTE,COMPLEASM_1,COMPLEASM_2,EST_SIZE,REF_SEQ
    None,/mnt/d/TESTING_SPACE/Ps_zapotecorum/ONT_MinION/SRR25932369.fq.gz,None,/mnt/d/TESTING_SPACE/Ps_zapotecorum/Illumina_PE150/SRR25932370_1.fq.gz,/mnt/d/TESTING_SPACE/Ps_zapotecorum/Illumina_PE150/SRR25932370_2.fq.gz,Ps_zapotecorum,Funga,Eukaryote,basidiomycota,agaricales,60m,None
    None,/mnt/d/TESTING_SPACE/Ps_gandalfiana/ONT_MinION/SRR27945396.fq.gz,/mnt/d/TESTING_SPACE/Ps_gandalfiana/Illumina_PE150/B1_3,None,None,Ps_gandalfiana,Funga,Eukaryote,basidiomycota,agaricales,60m,/mnt/d/TESTING_SPACE/Ps_cubensis/GCF_017499595_1_MGC_Penvy_REF_SEQ/GCF_017499595_1_MGC_Penvy_1_genomic.fna

"""
# Python Imports
import math, platform, os, subprocess, multiprocessing, argparse, psutil, shutil, hashlib, re, gzip, glob
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from bs4 import BeautifulSoup

# Global Variables
CPU_THREADS = 1
RAM_GB = 1
DEFAULT_LOG_FILE = None
ENVIRONMENT_TYPE = None


# Functions Section
def generate_log_file(log_file_path, use_numerical_suffix=False):
    """
    Generates or clears a log file based on the given parameters.
    
    This function either creates a new log file or clears an existing one, depending on the specified parameters. 
    If the 'use_numerical_suffix' parameter is True, and the file already exists, a new file with a numerical suffix 
    will be created. Otherwise, the existing file will be cleared.
    
    Parameters:
        log_file_path (str): Path to the log file.
        use_numerical_suffix (bool): If True, creates new files with numerical suffixes if the file exists; otherwise, clears the existing file.
    
    Returns:
        str: Path to the log file.
    
    Notes:
        - This function is essential for managing log file versions, especially in long-running applications or in situations where log file preservation is crucial.
        - The numerical suffix increments for each new file created in the same location.
        - When 'use_numerical_suffix' is False, it refreshes the log file by clearing existing content.
    
    Considerations:
        - Ensure the directory for the log file exists, or handle directory creation within the function or externally.
    
    Examples:
        log_file_path = "logs/my_log.txt"
        generate_log_file(log_file_path, use_numerical_suffix=True)
    """
    if os.path.exists(log_file_path) and use_numerical_suffix:
        counter = 1
        new_log_file_path = f"{log_file_path.rsplit('.', 1)[0]}_{counter}.txt"
        while os.path.exists(new_log_file_path):
            counter += 1
            new_log_file_path = f"{log_file_path.rsplit('.', 1)[0]}_{counter}.txt"
        log_file_path = new_log_file_path
    else:
        open(log_file_path, "w").close()
    return log_file_path


def log_print(input_message, log_file=None):
    """
    Logs a message to a file and prints it to the console with appropriate coloring.
    
    This function takes a message and logs it to the specified file. Additionally, the message is printed to the 
    console, potentially with specific coloring depending on the context.
    
    Parameters:
        input_message (str): Message to be logged and printed.
        log_file (str): Path to the log file.

    Notes:
        - This function serves as a centralized way to manage both logging and console output, ensuring consistency across the application.
        - The function uses a global default log file if none is specified.
        - Timestamps each log entry for easy tracking.
        - Utilizes color coding in the console to distinguish between different types of messages (e.g., errors, warnings).
        - Supports color coding for specific message types: NOTE, CMD, ERROR, WARN, and PASS.
        - Falls back to default (white) color if the message type is unrecognized.
    
    Considerations:
        - Consider the security implications of logging sensitive information.
        
    Examples:
        log_print("NOTE: Starting process")
        log_print("ERROR: An error occurred", log_file="error_log.txt")
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
    print_color = "white"  # Default color
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
    Initializes the logging environment based on the given input file path.

    This function sets up the logging environment by adjusting file paths according to the operating system in use, 
    ensuring file existence, and then generating a log file. It sets the global DEFAULT_LOG_FILE variable to the path 
    of the generated log file.

    Parameters:
        INPUT_FOLDER (str): Path to the folder which influences the log file generation.

    Global Variables:
        DEFAULT_LOG_FILE (str): The default path for the log file used throughout the logging process.

    Notes:
        - Adapts to different operating systems, making the logging system more flexible.
        - Prints unlogged messages to the console regarding environment detection and file existence.
        - Modifies the global DEFAULT_LOG_FILE variable.
        
    Considerations:
        - Verify the input folder's existence and accessibility before attempting to create a log file.

    Examples:
        input_folder = "data/input_data"
        initialize_logging_environment(input_folder)
    """
    global DEFAULT_LOG_FILE, ENVIRONMENT_TYPE
    input_file_path = f"{INPUT_FOLDER}/{INPUT_FOLDER.split('/')[-1]}_log.txt"
    os_name = platform.system()
    if os_name == "Windows":
        print("UNLOGGED:\tWINDOWS ENVIRONMENT")
        ENVIRONMENT_TYPE = "WIN"
    elif os_name in ["Linux", "Darwin"]:  # Darwin is the system name for macOS
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
    run_log = generate_log_file(input_file_path, use_numerical_suffix=False)
    DEFAULT_LOG_FILE = run_log


def run_subprocess_cmd(cmd_list, shell_check):
    """
    Executes a command using the subprocess.Popen and displays its output in real-time.

    This function is designed to execute shell commands from within a Python script. It uses subprocess.Popen to
    provide real-time output of the command to the command line window. It also logs the command execution details.

    Parameters:
        cmd_list (str or list of str): The command to be executed. Can be a single string or a list of strings
                                       representing the command and its arguments.
        shell_check (bool): If True, the command is executed through the shell. This is necessary for some 
                            commands, especially those that are built into the shell or involve shell features
                            like wildcard expansion.

    Features:
        - Uses subprocess.Popen for more control over command execution.
        - Captures the command's standard output and errors in real-time and displays them as they occur.
        - Waits for the command to complete, checks the return code to determine success or failure, and returns it.
        - Logs the command, its real-time output, and any execution errors.

    Usage and Considerations:
        - Useful for executing commands where live feedback is important, especially for long-running commands.
        - Requires careful use of 'shell_check' due to potential security risks with shell commands.

    Example:
        result_code = run_subprocess_cmd(["ls", "-l"], shell_check=False)
        print("Return code:", result_code)
    """
    if isinstance(cmd_list, str):
        log_print(f"CMD:\t{cmd_list}")    
        process = subprocess.Popen(cmd_list, shell=shell_check, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    else:
        log_print(f"CMD:\t{' '.join(cmd_list)}")    
        process = subprocess.Popen(cmd_list, shell=shell_check, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    for line in process.stdout:
        print(line, end="")
    process.wait()
    if process.returncode != 0:
        log_print(f"NOTE:\tCommand failed with return code {process.returncode}")
    else:
        log_print(f"PASS:\tSuccessfully processed command: {' '.join(cmd_list)}" if isinstance(cmd_list, list) else cmd_list)
    return process.returncode
        

def get_resource_values(PERCENT_RESOURCES):
    """
    Converts user input PERCENT_RESOURCES into usuable cpu_threads and ram_gb values.

    Parameters:
        PERCENT_RESOURCES (float): Percentage of resources to use.

    Returns:
        cpu_threads (str): A count of the CPUs available for processing.
        ram_gb (str): A count of the RAM (in GB) available for processing.

    Notes:
        - Allows dynamic allocation of resources based on the system's current state.

    Considerations:
        - Ensure that the PERCENT_RESOURCES value is within an acceptable range to avoid over-utilization of system resources.

    Examples:
        cpu_threads, ram_gb = get_resource_values(0.5)  # Use 50% of system resources
    """
    num_cpus = multiprocessing.cpu_count()
    mem_info = psutil.virtual_memory()
    cpu_threads = int(math.floor(num_cpus * PERCENT_RESOURCES))
    ram_gb = int(mem_info.total / (1024.0 ** 3) * PERCENT_RESOURCES)
    return cpu_threads, ram_gb 


def find_file(filename):
    """
    Searches for a file within the current working directory.

    Parameters:
        filename (str): The name of the file to search for.

    Returns:
        str: The path to the first instance of the file, if found.
        None: If the file is not found.
        
    Notes:
        - Useful for setups where the exact location of a file might vary.
    
    Considerations:
        - This function might take a long time in large file systems.
    
    Examples:
        file_path = find_file("config.json")
        if file_path:
            print(f"Found config at {file_path}")
        else:
            print("Config file not found")
    """
    global ENVIRONMENT_TYPE
    log_print(f"Looking for {filename}")    
    if ENVIRONMENT_TYPE == "WIN":
        root_directory = "C:\\"  # Adjust if necessary for different drives
    elif ENVIRONMENT_TYPE in ["LINUX/WSL/MAC"]:
        root_directory = "/"
    else:
        raise ValueError("Unknown ENVIRONMENT_TYPE")
    for root, dirs, files in os.walk(root_directory):
        if filename in files:
            return os.path.join(root, filename)
    return None


def md5_check(illumina_raw_data_dir, illumina_df):
    """
    Run MD5 checksums on all `.FQ.GZ` files in the provided directory, referencing expected checksums from an `MD5.txt` file.

    Parameters:
        illumina_raw_data_dir (str): 
            Path to the main Illumina Raw Reads Data Directory, which should contain a `MD5.txt` file 
            and the `.FQ.GZ` files to be checked.
        illumina_df (DataFrame): 
            A pandas DataFrame intended to be populated with `MD5` checksums and filenames.

    Returns:
        DataFrame:
            Updated `illumina_df` containing the original MD5 checksums (from `MD5.txt`) and any newly discovered `.FQ.GZ` files 
            that match those checksums.

    Notes:
        - The function reads a `MD5.txt` file, splits its contents into a dictionary, then populates the provided DataFrame.
        - It iterates over all `.gz` files, computes the new MD5 and compares it to the original MD5 to detect mismatches.

    Considerations:
        - This function will raise an error only through log messages (`log_print`) if the MD5 is missing or mismatched.
        - If a mismatch is found, the function breaks out of the loop but does not raise an exception. 
          Adjust logic if you need strict error handling.

    Examples:
        >>> # Example usage:
        >>> illumina_df = pd.DataFrame(columns=["MD5", "Filename"])
        >>> updated_df = md5_check("/path/to/illumina/raw_data", illumina_df)
        >>> print(updated_df.head())
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
                        log_print(f"PASS:\tOriginal MD5 MATCHES for {file_name.split('/')[-1]}")
                else:
                    log_print(f"ERROR:\tOriginal MD5 checksum not found for {file_name.split('/')[-1]}")
                    break
    

def illumina_extract_and_check(folder_name, SPECIES_ID):
    """
    Extract and verify Illumina MD5 checksums, then combine paired-end reads into single forward and reverse files.

    Parameters:
        folder_name (str): 
            Path to the folder containing Illumina data (`.fq.gz` files and `MD5.txt`).
        SPECIES_ID (str): 
            A unique identifier for the species. Used to name the combined output files.

    Returns:
        list:
            A list of two gzipped combined FASTQ files. The first (`_combined_1.fq.gz`) corresponds to forward reads, 
            and the second (`_combined_2.fq.gz`) to reverse reads.

    Notes:
        - This function calls `md5_check` to verify checksums and populate a DataFrame with MD5 data.
        - After MD5 verification, if combined files do not already exist, it concatenates raw forward and reverse 
          FASTQ files, and then gzips them.

    Considerations:
        - If the combined or gzipped files already exist, the function will skip regeneration steps and simply return 
          the existing file paths.
        - Adjust logging or error-handling logic as needed.

    Examples:
        >>> # Example usage:
        >>> combined_files = illumina_extract_and_check("/path/to/illumina/folder", "MySpecies")
        >>> print(combined_files)
        ['/path/to/illumina/../MySpecies_combined_1.fq.gz',
         '/path/to/illumina/../MySpecies_combined_2.fq.gz']
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
                    raw_1_list.append(os.path.join(folder_name,filename))
                elif "_2.fq.gz" in filename:
                    raw_2_list.append(os.path.join(folder_name,filename))
            fwd_cat_cmd = f"cat {raw_1_list[0]} {raw_1_list[1]} > {combined_1_file}"            
            _ = run_subprocess_cmd(fwd_cat_cmd, shell_check = True)
            rev_cat_cmd = f"cat {raw_2_list[0]} {raw_2_list[1]} > {combined_2_file}"            
            _ = run_subprocess_cmd(rev_cat_cmd, shell_check = True)
        else:
            log_print(f"SKIP:\tCombined FASTQ files already exist: {combined_list[0]}; {combined_list[1]}.")
        for combined_fq in combined_list:
            gzip_cmd = ["gzip", combined_fq]
            _ = run_subprocess_cmd(gzip_cmd, shell_check = False)
    else:
        log_print(f"SKIP:\tGzipped Combined FASTQ files already exist: {combined_gz_list[0]}; {combined_gz_list[1]}.")
    return combined_gz_list 


def classify_metric(value, thresholds):
    """
    Classifies a given value according to specified thresholds for genome assembly metrics.

    Parameters:
        value (float): The metric value to classify.
        thresholds (dict): A dictionary specifying the classification thresholds for AMAZING, GREAT, OK, and POOR.

    Returns:
        str: A string classification of the value ("AMAZING", "GREAT", "OK", "POOR").

    Notes:
        - This function simplifies the process of metric evaluation across different genome assembly parameters.

    Considerations:
        - Ensure that the 'thresholds' dictionary is correctly set up with appropriate limits for each category.
    
    Examples:
        n50_classification = classify_metric(1500000, {"AMAZING": 1000000, "GREAT": 100000, "OK": 1000, "POOR": 100})
        print(n50_classification)  # Output should be "AMAZING"
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
    Classifies the genome assembly quality based on multiple metrics, including BUSCO scores, N50, and number of contigs.

    Parameters:
        sample_stats (dict): A dictionary containing the values for various assembly metrics such as "FIRST_COMPLEASM_C",
                             "SECOND_COMPLEASM_C", "ASSEMBLY_N50", and "ASSEMBLY_CONTIGS".

    Returns:
        dict: A dictionary with the classification for each metric and an overall classification based on combined criteria.

    Notes:
        - This function utilizes `classify_metric` for individual assessments and determines the overall quality.
        - Provides a granular and comprehensive evaluation of genome assembly quality.

    Considerations:
        - The overall assembly quality heavily depends on the individual classifications. A single "POOR" rating can impact
          the overall rating.
        - Adjusting the individual thresholds will affect the sensitivity of the classification.

    Examples:
        sample_stats_dict = {"FIRST_COMPLEASM_C": 99, "SECOND_COMPLEASM_C": 97, "ASSEMBLY_N50": 1500000, "ASSEMBLY_CONTIGS": 50}
        assembly_quality = classify_assembly(sample_stats_dict)
        print(assembly_quality)
        # Output: {"FIRST_COMPLEASM_C": "AMAZING", "SECOND_COMPLEASM_C": "AMAZING", "ASSEMBLY_N50": "AMAZING", "ASSEMBLY_CONTIGS": "AMAZING", "OVERALL": "AMAZING"}
    """
    results = {}
    if sample_stats["FIRST_COMPLEASM_C"] <= 80.0:
        results["FIRST_COMPLEASM_C"] = "POOR"
    elif sample_stats["FIRST_COMPLEASM_C"] <= 90.0:
        results["FIRST_COMPLEASM_C"] = "OK"
    elif sample_stats["FIRST_COMPLEASM_C"] <= 95.0:
        results["FIRST_COMPLEASM_C"] = "GREAT"
    elif sample_stats["FIRST_COMPLEASM_C"] <= 98.5:
        results["FIRST_COMPLEASM_C"] = "AMAZING"
    if sample_stats["SECOND_COMPLEASM_C"] <= 80.0:
        results["SECOND_COMPLEASM_C"] = "POOR"
    elif sample_stats["SECOND_COMPLEASM_C"] <= 90.0:
        results["SECOND_COMPLEASM_C"] = "OK"
    elif sample_stats["SECOND_COMPLEASM_C"] <= 95.0:
        results["SECOND_COMPLEASM_C"] = "GREAT"
    elif sample_stats["SECOND_COMPLEASM_C"] <= 98.5:
        results["SECOND_COMPLEASM_C"] = "AMAZING"
    if sample_stats["ASSEMBLY_CONTIGS"] <= 100:
        results["ASSEMBLY_CONTIGS"] = "AMAZING"
    elif sample_stats["ASSEMBLY_CONTIGS"] <= 1000:
        results["ASSEMBLY_CONTIGS"] = "GREAT"
    elif sample_stats["ASSEMBLY_CONTIGS"] <= 10000:
        results["ASSEMBLY_CONTIGS"] = "OK"
    elif sample_stats["ASSEMBLY_CONTIGS"] <= 100000:
        results["ASSEMBLY_CONTIGS"] = "POOR"
    if sample_stats["ASSEMBLY_N50"] <= 100:
        results["ASSEMBLY_N50"] = "POOR"
    elif sample_stats["ASSEMBLY_N50"] <= 1000:
        results["ASSEMBLY_N50"] = "OK"
    elif sample_stats["ASSEMBLY_N50"] <= 10000:
        results["ASSEMBLY_N50"] = "GREAT"
    elif sample_stats["ASSEMBLY_N50"] <= 100000:
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


def plot_classification_table(df):
    """
    Generate a color-coded table from a DataFrame where each cell’s color corresponds to a classification label.

    Parameters:
        df (DataFrame): 
            A pandas DataFrame where indices are sample identifiers and columns are metrics. 
            The cell values should be classification labels (e.g., "AMAZING", "GREAT", "OK", "POOR").

    Returns:
        None:
            This function displays a matplotlib figure showing the table, but does not return a value.

    Notes:
        - This function uses matplotlib’s table functionality to create a color-mapped table.
        - The `color_map` dictionary defines the relationship between classification labels and cell colors.

    Considerations:
        - Ensure that the cell values in `df` match the keys in `color_map` if you want custom colors.
        - Adjust the figure size or layout for optimal display of large tables.

    Examples:
        >>> # Example usage:
        >>> data = {
        ...     "Metric1": ["AMAZING", "GREAT", "POOR"],
        ...     "Metric2": ["OK", "GREAT", "AMAZING"]
        ... }
        >>> sample_df = pd.DataFrame(data, index=["SampleA", "SampleB", "SampleC"])
        >>> plot_classification_table(sample_df)
        # A matplotlib figure should pop up showing a color-coded table.
    """
    color_map = {"AMAZING": "lightorchid",
                 "GREAT": "green",
                 "OK": "cornflowerblue",
                 "POOR": "orange",}
    colors = df.applymap(lambda x: color_map.get(x, "white"))
    fig, ax = plt.subplots(figsize=(10, 2))  # Adjust the figure size as necessary
    ax.axis("tight")
    ax.axis("off")
    ax.table(cellText = df.values, colLabels = df.columns,
             cellColours = colors.values, cellLoc = "center", loc = "center")
    plt.show()


def find_ca_folder(input_folder):
    """
    Determine the specific location of the MaSuRCA output folder named `CA`.

    Parameters:
        input_folder (str):
            Path to the folder containing MaSuRCA outputs. It may have subdirectories, 
            one of which might begin with "CA" (e.g., "CA", "CA_12345", etc.).

    Returns:
        str:
            The path to the first subfolder starting with "CA". If none is found, 
            returns a default path `<input_folder>/CA`.

    Notes:
        - This function scans the immediate subfolders of `input_folder` and returns 
          the one that starts with "CA".
        - If no subfolder starts with "CA", it defaults to `<input_folder>/CA`.

    Considerations:
        - Ensure `input_folder` is a valid path containing directories. 
        - If multiple folders begin with "CA", it returns the first one found.

    Examples:
        >>> # Example usage:
        >>> ca_path = find_ca_folder("/path/to/masurca_outputs")
        >>> print(ca_path)
        "/path/to/masurca_outputs/CA_20230101"
    """
    subfolders = [f.path for f in os.scandir(input_folder) if f.is_dir()]
    ca_folder = f"{input_folder}/CA"
    for folder in subfolders:
        if os.path.basename(folder).startswith("CA"):
            ca_folder = folder
            break
    return ca_folder

def find_soap_folder(input_folder):
    """
    Determine the specific location of the MaSuRCA output folder named `CA`.

    Parameters:
        input_folder (str):
            Path to the folder containing MaSuRCA outputs. It may have subdirectories, 
            one of which might begin with "CA" (e.g., "CA", "CA_12345", etc.).

    Returns:
        str:
            The path to the first subfolder starting with "CA". If none is found, 
            returns a default path `<input_folder>/CA`.

    Notes:
        - This function scans the immediate subfolders of `input_folder` and returns 
          the one that starts with "CA".
        - If no subfolder starts with "CA", it defaults to `<input_folder>/CA`.

    Considerations:
        - Ensure `input_folder` is a valid path containing directories. 
        - If multiple folders begin with "CA", it returns the first one found.

    Examples:
        >>> # Example usage:
        >>> ca_path = find_soap_folder("/path/to/masurca_outputs")
        >>> print(ca_path)
        "/path/to/masurca_outputs/CA_20230101"
    """
    subfolders = [f.path for f in os.scandir(input_folder) if f.is_dir()]
    soap_folder = f"{input_folder}/SOAP_assembly"
    for folder in subfolders:
        if os.path.basename(folder).startswith("CA"):
            soap_folder = folder
            break
    return soap_folder


def parse_bbmerge_output(insert_size_histogram_txt):
    """
    Extract average insert size and standard deviation from a BBMerge insert size histogram.

    Parameters:
        insert_size_histogram_txt (str):
            Path to the `insert_size_histogram.txt` file produced by BBMerge.

    Returns:
        tuple:
            A 2-tuple containing:
                - avg_insert (float): The average insert size.
                - std_dev (float): The standard deviation of the insert size.

    Notes:
        - The file is expected to contain lines starting with `#Mean` and `#STDev`.
        - If those lines are not found, a ValueError is raised.

    Considerations:
        - If the insert size file format changes in newer BBMerge versions, 
          this function may need updating.

    Examples:
        >>> # Suppose the insert_size_histogram.txt has lines:
        >>> #Mean    250.7
        >>> #STDev   30.2
        >>> avg, std = parse_bbmerge_output("/path/to/insert_size_histogram.txt")
        >>> print(avg, std)
        (251.0, 30.0)
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
    Use BBMap (specifically BBMerge) to calculate statistics like average insert size and standard deviation.

    Parameters:
        input_folder (str):
            Path to the folder where intermediate outputs (e.g., `bbmap_data.fq.gz`, `insert_size_histogram.txt`) 
            will be stored or read from.
        reads_list (list):
            A list of paths to FASTQ files. 
            - If there are two items, they are treated as forward and reverse reads. 
            - If there is one item, it is treated as a single-end read.

    Returns:
        tuple:
            A 2-tuple of (avg_insert, std_dev) representing the average insert size and standard deviation.

    Notes:
        - If an existing `insert_size_histogram.txt` is found, the function reuses it and skips rerunning BBMerge.
        - If the file is missing, BBMerge is invoked to generate it.

    Considerations:
        - By default, if parsing fails, the function logs a note and uses the default `(251, 30)` values.
        - Ensure `BBMerge` or `bbmerge.sh` is installed and accessible in `PATH` or specify the full path as needed.

    Examples:
        >>> # Example usage with paired-end reads:
        >>> stats = bbmap_stats("/path/to/folder", ["/path/to/read1.fq.gz", "/path/to/read2.fq.gz"])
        >>> print(stats)
        (250, 30)
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
        if len(reads_list) == 3:
            bbmerge_cmd = [bbmerge_path,
                           f"in1={reads_list[1]}",
                           f"in2={reads_list[2]}",
                           f"out={bbmap_out_path}",
                           f"ihist={input_folder}/insert_size_histogram.txt"]
        elif len(reads_list) == 2:
            bbmerge_cmd = [bbmerge_path,
                           f"in1={reads_list[0]}",
                           f"in2={reads_list[1]}",
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
    Modify a MaSuRCA `assemble.sh` script to skip the gap closing step, forcing the final output to remain at 
    stage "9-terminator".

    Parameters:
        assembly_sh_path (str):
            Path to the `assemble.sh` script generated by MaSuRCA.

    Returns:
        None:
            This function writes out a modified script named `assemble_skip_gap.sh` in the current working directory 
            but does not return a value.

    Notes:
        - Uses a regular expression to locate the `if [ -s $CA_DIR/9-terminator/genome.scf.fasta ];then ... else ... fi` 
          block and replace its content with a snippet that logs a skip message and sets `TERMINATOR="9-terminator"`.
        - The original script is read in full, and only that specific section is replaced. 
          The rest of the script is left unmodified.

    Considerations:
        - Ensure the original `assemble.sh` contains the expected pattern. If the pattern changes in newer versions, 
          this function needs updating.
        - The new script is always named `assemble_skip_gap.sh`; be mindful that running this multiple times 
          might overwrite existing files.

    Examples:
        >>> # Example usage:
        >>> skip_gap_closing_section("/path/to/assemble.sh")
        >>> # A file named "assemble_skip_gap.sh" is created with the gap closing step skipped.
    """
    with open(assembly_sh_path, "r") as f_in:
        original_script = f_in.read()
    # Regex Explanation:
    #
    #   (if \[ -s \$CA_DIR/9-terminator/genome\.scf\.fasta \];then)
    #      Captures the exact `if [ -s $CA_DIR/9-terminator/genome.scf.fasta ];then` line as group 1
    #
    #   ( .*? )
    #      Lazily matches everything in between, as group 2
    #
    #   (else\s+fail "Assembly stopped or failed, see \$CA_DIR\.log"\nfi)
    #      Captures the part from `else` down to the `fi` as group 3, preserving that code
    #
    # We use DOTALL so that `.` also matches newlines.
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


def masurca_config_gen(input_folder, output_folder, input_fq_list, clump_f_dedup_path, clump_r_dedup_path, CPU_THREADS, ram_gb, ref_seq=None):
    """
    DESCRIPTION:
        Generates a MaSuRCA configuration file and runs a genome assembly workflow. 
        This includes creating a config file based on input parameters, optionally skipping
        gap closing, and running the assembly. The assembled output is moved or renamed 
        depending on whether a reference sequence is used.

    Parameters:
        input_folder (str):
            The directory where input data (and possibly the CA sub-folder) is located.
        output_folder (str):
            The directory where the MaSuRCA config file and assembly scripts will be generated.
        input_fq_list (list of str):
            A list of input FASTQ files (often 2 for paired-end reads).
        clump_f_dedup_path (str):
            The path to the forward de-duplicated FASTQ file (e.g., from Clumpify).
        clump_r_dedup_path (str):
            The path to the reverse de-duplicated FASTQ file (e.g., from Clumpify).
        CPU_THREADS (int):
            The number of threads available for the assembly process.
        ram_gb (int):
            The amount of RAM (in GB) available to the assembly process.
        ref_seq (str, optional):
            An optional path to a reference genome. If provided, the assembled output 
            will be renamed to "primary.genome.ref.fasta".

    Returns:
        tuple:
            A tuple of strings indicating:
            (1) `default_assembly_path`: The default path to the primary assembled genome 
                (e.g., "CAfolder/primary.genome.scf.fasta").
            (2) `assembly_path`: The path to the primary assembled genome if reference 
                is not used (e.g., "input_folder/primary.genome.scf.fasta").
            (3) `ref_assembly_path` (or None if `ref_seq` is not provided).

    Notes:
        - The function relies on a correct folder structure where a "CA" sub-folder 
          (created by MaSuRCA or another pipeline step) is located in `input_folder`.
        - Adjustments to the JF_SIZE parameter are made based on available RAM 
          to potentially improve efficiency.

    Considerations:
        - The function skips re-running MaSuRCA if it detects that the assembled output 
          files already exist.
        - It modifies the `assemble.sh` script via `skip_gap_closing_section` to skip 
          certain gap-closing steps, which can shorten runtime but may reduce assembly contiguity.
        - Error handling for MaSuRCA is performed by checking the return code of 
          the assembly command; if it’s 1, an error is logged, but no further exception is raised.

    Examples:
        >>> # Example usage:
        >>> input_folder = "/path/to/input_data"
        >>> output_folder = "/path/to/output_data"
        >>> input_fq_list = ["reads_1.fastq", "reads_2.fastq"]
        >>> clump_f_dedup_path = "/path/to/reads_forward_dedup.fastq"
        >>> clump_r_dedup_path = "/path/to/reads_reverse_dedup.fastq"
        >>> CPU_THREADS = 16
        >>> ram_gb = 64
        >>> ref_seq = "/path/to/reference.fna"
        >>> result_paths = masurca_config_gen(
        ...     input_folder, output_folder, input_fq_list,
        ...     clump_f_dedup_path, clump_r_dedup_path,
        ...     CPU_THREADS, ram_gb, ref_seq
        ... )
        >>> print(result_paths)
        ("/path/to/CAfolder/primary.genome.scf.fasta",
         "/path/to/input_data/primary.genome.scf.fasta",
         "/path/to/input_data/primary.genome.ref.fasta")
    """
    # If clump_f_dedup_path and/or clump_r_dedup_path contain ".gz" in their name, unzip them with gzip
    dedup_unzip_list = [input_file.replace(".gz","") for input_file in [clump_f_dedup_path, clump_r_dedup_path]]
    for input_file in [clump_f_dedup_path, clump_r_dedup_path]:
        if ".gz" in input_file:
            gunzip_file(input_file, input_file.replace(".gz",""))
    
    avg_insert, std_dev = bbmap_stats(input_folder, input_fq_list)
    os.chdir(output_folder)
    jf_size = 2500000000 # BASED ON estimated_genome_size*20
    max_ram_tested = 62
    if ram_gb < max_ram_tested:
        adjustment_ratio = ram_gb / max_ram_tested
        jf_size = int(round(jf_size * adjustment_ratio, 0))
    ca_folder = find_ca_folder(input_folder)
    soap_folder = find_soap_folder(input_folder)
    data_output_folder = ca_folder
    default_assembly_path = os.path.join(ca_folder, "primary.genome.scf.fasta")
    assembly_path = os.path.join(input_folder, "primary.genome.scf.fasta")
    ref_assembly_path = os.path.join(input_folder, "primary.genome.ref.fasta")
    if len(input_fq_list) == 2:
        illu_only_mates = 1
        illu_only_gaps = 1
        illu_only_soap = 0
    elif len(input_fq_list) >= 2:
        illu_only_mates = 0
        illu_only_gaps = 0
        illu_only_soap = 0
    if os.path.exists(default_assembly_path):
        log_print("PASS:\tSkipping MaSuRCA, moving output files")
        if ref_seq:
            shutil.move(default_assembly_path, ref_assembly_path)
        else:
            shutil.move(default_assembly_path, assembly_path)
    elif os.path.exists(assembly_path):
        log_print(f"SKIP:\tMaSuRCA Assembly, scaffolded assembly already exists: {assembly_path}.")
    else:
        config_content = ["DATA\n",
                          f"PE= pe {avg_insert} {std_dev} {clump_f_dedup_path.replace('.gz','')} {clump_r_dedup_path.replace('.gz','')}\n"]
        if ref_seq:
            config_content.append(f"REFERENCE={ref_seq.replace('.gbff','.fna.gz')}\n")
        config_content.append("END\n")
        config_content.append("PARAMETERS\n")
        config_content.append("GRAPH_KMER_SIZE=auto\n")
        config_content.append(f"USE_LINKING_MATES={illu_only_mates}\n") # IF ILLUMINA ONLY SET TO 1, ELSE 0
        config_content.append(f"CLOSE_GAPS={illu_only_gaps}\n") # IF ILLUMINA ONLY SET TO 1, ELSE 0
        config_content.append("MEGA_READS_ONE_PASS=0\n")
        config_content.append("LIMIT_JUMP_COVERAGE=300\n")
        config_content.append("CA_PARAMETERS=cgwErrorRate=0.15\n")
        config_content.append(f"NUM_THREADS={CPU_THREADS}\n")
        config_content.append(f"JF_SIZE={jf_size}\n")
        config_content.append(f"SOAP_ASSEMBLY={illu_only_soap}\n") # IF ILLUMINA ONLY SET TO 0, ELSE 1
        config_content.append("END\n")
        config_path = f"{output_folder}/masurca_config_file.txt"
        with open(config_path, "w") as file:
            for entry in config_content:
                file.write(entry)
        masurca_config_cmd = ["masurca", "masurca_config_file.txt"]
        _ = run_subprocess_cmd(masurca_config_cmd, False)       
    assemble_sh_path = f"{output_folder}/assemble.sh"
    skip_gap_closing_section(assemble_sh_path)
    masurca_assemble_cmd = ["bash", assemble_sh_path]
    return_code = run_subprocess_cmd(masurca_assemble_cmd, False)    
    for output_file in dedup_unzip_list:
        try:
            os.remove(output_file)  # Remove the file
            print(f"Removed: {output_file}")
        except FileNotFoundError:
            print(f"File not found: {output_file}")
        except PermissionError:
            print(f"Permission denied: {output_file}")
        except Exception as e:
            print(f"Error removing {output_file}: {e}")
    if return_code == 1:
        log_print("NOTE:\tMaSuRCA assembly interrupted intentionally; Gap Closing will be performed later.")
        default_assembly_path = os.path.join(data_output_folder, os.path.basename(default_assembly_path))
        return default_assembly_path, assembly_path, ref_assembly_path if ref_seq else None
    return default_assembly_path, assembly_path, ref_assembly_path if ref_seq else None


def bbduk_map(trimmo_f_pair_path, trimmo_r_pair_path):
    """
    DESCRIPTION:
        Runs BBDuk to perform quality trimming and adapter removal on paired-end FASTQ files.
        The function will generate two new FASTQ files (forward and reverse) with "_mapped" 
        in their file names to indicate they have been processed by BBDuk.

    Parameters:
        trimmo_f_pair_path (str):
            The file path to the forward paired FASTQ file (e.g., "sample_1_paired.fastq").
        trimmo_r_pair_path (str):
            The file path to the reverse paired FASTQ file (e.g., "sample_2_paired.fastq").

    Returns:
        tuple of str:
            A tuple containing the paths to the forward and reverse mapped FASTQ files 
            (e.g., "sample_forward_mapped.fastq", "sample_reverse_mapped.fastq").

    Notes:
        - This function assumes your forward reads end in "_1_paired.{extension}" 
          and your reverse reads end in "_2_paired.{extension}".
        - It also expects that BBDuk ("bbduk.sh") and a file of adapters ("adapters.fa")
          are accessible via the `find_file` utility function.
        - The BBDuk command used here includes parameters for trimming based on quality score (qtrim=rl, trimq=20)
          and adapter removal (ktrim=r, ref=adapters.fa, etc.). Adjust them as needed.

    Considerations:
        - If BBDuk or the adapters file are not found, the function will fail unless handled in `find_file`.
        - If the output mapped files already exist, the function skips re-processing and logs a message.
        - Ensure the user has read/write access to the specified paths.

    Examples:
        >>> # Example usage:
        >>> forward_in = "sample_1_paired.fastq"
        >>> reverse_in = "sample_2_paired.fastq"
        >>> f_mapped, r_mapped = bbduk_map(forward_in, reverse_in)
        >>> print(f_mapped, r_mapped)
        sample_forward_mapped.fastq sample_reverse_mapped.fastq
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
    DESCRIPTION:
        De-duplicates paired-end FASTQ files using Clumpify. If the output files already exist,
        the command is skipped, otherwise Clumpify is invoked to generate de-duplicated versions.

    Parameters:
        bbduk_f_map_path (str):
            The path to the forward mapped FASTQ file (e.g., "reads_forward_mapped.fastq").
        bbduk_r_map_path (str):
            The path to the reverse mapped FASTQ file (e.g., "reads_reverse_mapped.fastq").

    Returns:
        tuple of str:
            A tuple containing the paths to the forward and reverse de-duplicated FASTQ files.

    Notes:
        - This function assumes that the forward file name contains "_forward_mapped."
          and the reverse file name contains "_reverse_mapped." within their paths.
        - Relies on the external Clumpify program ("clumpify.sh") being accessible in the 
          system's PATH or in a known location that `find_file` can discover.

    Considerations:
        - Ensure that the user has permission to read and write to the specified paths.
        - The deduplicated output files will replace "mapped" with "dedup" in their file names.
        - If Clumpify is not found, the process will likely fail unless handled by `find_file`.

    Examples:
        >>> # Example usage:
        >>> # Suppose we have two mapped FASTQ files:
        >>> forward_path = "sample_forward_mapped.fastq"
        >>> reverse_path = "sample_reverse_mapped.fastq"
        >>> f_dedup, r_dedup = clumpify_dedup(forward_path, reverse_path)
        >>> print(f_dedup, r_dedup)
        sample_forward_dedup.fastq sample_reverse_dedup.fastq
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
    DESCRIPTION:
        Opens an HTML file, parses it to find the row that contains "Total Bases",
        then returns the string value in the subsequent cell (e.g., "4.8 Gbp").
    
    Parameters:
        html_file (str):
            The path to the HTML file that contains the "Total Bases" table row.
    
    Returns:
        total_bases_value (str or None):
            The value in the adjacent cell to "Total Bases" (e.g., "4.8 Gbp").
            Returns None if "Total Bases" is not found or if the adjacent cell is missing.
    
    Notes:
        - This function uses BeautifulSoup to parse the HTML, which must be installed 
          (e.g., `pip install beautifulsoup4`).
        - Expects well-formed HTML with a <td> containing exactly the text "Total Bases".
    
    Considerations:
        - If the HTML structure changes significantly, the function may fail to locate
          the required text.
        - The file must be encoded in UTF-8 or compatible with the specified encoding.
    
    Examples:
        >>> # Example usage:
        >>> value = get_total_bases('example.html')
        >>> if value:
        ...     print(f"Total Bases: {value}")
        ... else:
        ...     print("Total Bases entry not found.")
    """
    # Open and read the contents of the HTML file
    with open(html_file, "r", encoding="utf-8") as f:
        contents = f.read()
    
    # Parse the HTML contents with BeautifulSoup
    soup = BeautifulSoup(contents, "html.parser")
    
    # Find the table cell that contains the text "Total Bases"
    # Then get its sibling (the next <td>), which contains the value (e.g., "4.8 Gbp")
    total_bases_td = soup.find("td", string="Total Bases")
    
    if total_bases_td:
        value_td = total_bases_td.find_next_sibling("td")
        if value_td:
            return value_td.get_text(strip=True)
    
    return None


def process_read_file(read_path):
    """
    Processes an Illumina read file to ensure it has a .fq.gz extension.
    
    Parameters:
    -----------
    read_path : str
        Path to the read file.
    
    Returns:
    --------
    new_read_path : str
        Path to the processed read file with the correct extension.
    """
    if not isinstance(read_path, str):
        log_print(f"Read path is not a string: {read_path}")
        return read_path  # Return as is if not a string

    dir_path, filename = os.path.split(read_path)
    basename, ext = os.path.splitext(filename)
    new_read_path = read_path  # Initialize with original path

    # Handle double extensions like .fastq.gz
    if ext == ".gz":
        basename, ext1 = os.path.splitext(basename)
        ext = ext1 + ext  # e.g., ".fq.gz" or ".fastq.gz"
    if read_path.endswith(".fq.gz"):
        # Correct extension; do nothing
        log_print(f"Read file already in .fq.gz format: {read_path}")
        return read_path
    
    elif read_path.endswith(".fq"):
        # Gzip the file
        new_read_path = os.path.join(dir_path, basename + ".fq.gz")
        with open(read_path, "rb") as f_in, gzip.open(new_read_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(read_path)  # Remove the original .fq file
        log_print(f"Gzipped .fq to {new_read_path}")
        return new_read_path
    
    elif read_path.endswith(".fastq.gz"):
        # Rename to .fq.gz
        new_read_path = os.path.join(dir_path, basename + ".fq.gz")
        os.rename(read_path, new_read_path)
        log_print(f"Renamed .fastq.gz to .fq.gz: {new_read_path}")
        return new_read_path
    
    elif read_path.endswith(".fastq"):
        # Rename to .fq and gzip
        temp_path = os.path.join(dir_path, basename + ".fq")
        os.rename(read_path, temp_path)
        with open(temp_path, "rb") as f_in, gzip.open(new_read_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(temp_path)  # Remove the intermediate .fq file
        log_print(f"Renamed .fastq to .fq and gzipped to {new_read_path}")
        return new_read_path
    
    else:
        log_print(f"Unrecognized file extension for read file: {read_path}")
        return read_path  # Return as is if extension is unrecognized


def gzip_file(input_file, output_file):
    """
    Compresses a file using gzip compression.

    Parameters:
    - input_file (str): Path to the input file to be compressed.
    - output_file (str): Path to the output gzip file.

    Returns:
    - None

    Considerations:
    - The function overwrites the output file if it already exists.
    - Ensure the input file exists and the program has permission to read it.
    - Ensure the program has permission to write to the output file's directory.
    - Gzipping is best suited for text or binary files and may not significantly reduce the size of already compressed files.

    Example:
    >>> input_file = "example.txt"
    >>> output_file = "example.txt.gz"
    >>> gzip_file(input_file, output_file)
    >>> print(f"{input_file} has been compressed to {output_file}")
    """
    with open(input_file, "rb") as f_in:  # Open the original file in binary read mode
        with gzip.open(output_file, "wb") as f_out:  # Open the gzip file in binary write mode
            shutil.copyfileobj(f_in, f_out)  # Copy the content to the gzip file

def gunzip_file(input_file, output_file):
    """
    Decompresses a gzip file back to its original format.

    Parameters:
    - input_file (str): Path to the gzip file to be decompressed.
    - output_file (str): Path to the output decompressed file.

    Returns:
    - None

    Considerations:
    - The function overwrites the output file if it already exists.
    - Ensure the input file exists and the program has permission to read it.
    - Ensure the program has permission to write to the output file's directory.
    - The function assumes the input file is a valid gzip file.

    Example:
    >>> input_file = "example.txt.gz"
    >>> output_file = "example.txt"
    >>> gunzip_file(input_file, output_file)
    >>> print(f"{input_file} has been decompressed to {output_file}")
    """    
    with gzip.open(input_file, "rb") as f_in:  # Open the gzip file in binary read mode
        with open(output_file, "wb") as f_out:  # Open the output file in binary write mode
            shutil.copyfileobj(f_in, f_out)  # Copy the decompressed content to the output file


def ont_combine_fastq_gz(ONT_FOLDER):
    """
    Combine multiple ONT FASTQ.GZ files into a single FASTQ file.

    Args:
        ONT_FOLDER (str): Path to the folder containing ONT FASTQ.GZ files.

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
        base_name  = ONT_FOLDER.split("/")[-2]
    combined_ont_fastq_path = os.path.join(ONT_FOLDER, f"{base_name}_ont_combined.fastq")
    raw_file_list = glob.glob(os.path.join(ont_raw_data_dir, "*.fastq.gz"))
    if os.path.isfile(combined_ont_fastq_path):
        log_print(f"NOTE:\tSkipping extraction & combination: Combined fastq file: {combined_ont_fastq_path} already exists")
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
    log_print(f"Generating BUSCO plot for {input_busco_tsv}...")
    # Load the BUSCO TSV into a pandas DataFrame
    busco_df = pd.read_csv(input_busco_tsv, sep="\t", header=0, dtype=str)
    
    busco_genes = len(busco_df)
    
    # Replace "Complete" with "Single" in the "Status" column
    busco_df['Status'] = busco_df['Status'].replace("Complete", "Single")
    
    # Create a pivot table to count Status values for each Sequence
    status_counts = busco_df.pivot_table(index='Sequence', columns='Status', aggfunc='size', fill_value=0)
    
    # Reorder the columns to maintain the original order
    desired_order = ['Single', 'Duplicated', 'Incomplete', 'Fragmented']
    status_counts = status_counts.reindex(columns=desired_order, fill_value=0)
    
    # Count total sequences and sequences to exclude
    total_sequences = len(status_counts)
    excluded_sequences = len(status_counts.loc[status_counts.drop(columns='Duplicated', errors='ignore').sum(axis=1) == 0])
    included_sequences = total_sequences - excluded_sequences
    
    # Remove sequences containing only Duplicated statuses
    filtered_status_counts = status_counts.loc[status_counts.drop(columns='Duplicated', errors='ignore').sum(axis=1) > 0]
    
    # Sort sequences by the total count of all statuses (optional)
    filtered_status_counts = filtered_status_counts.loc[filtered_status_counts.sum(axis=1).sort_values(ascending=False).index]
    
    # Calculate the total counts for each status for the legend
    status_totals = busco_df['Status'].value_counts()
    
    # Plot the data as a stacked bar chart with specified colors
    colors = {'Single': '#619B8AFF', 'Duplicated': '#A1C181FF', 'Incomplete': '#FE7F2DFF', 'Fragmented': '#FCCA46FF'}
    ax = filtered_status_counts.plot(
        kind='bar',
        stacked=True,
        figsize=(12, 8),
        color=[colors[col] for col in filtered_status_counts.columns])
    
    # Update the legend with totals and percentages
    legend_labels = [f"{status} ({round((status_totals.get(status, 0)/busco_genes)*100, 2)}%)" 
                     for status in filtered_status_counts.columns]
    
    completeness_values = [round((status_totals.get(status, 0)/busco_genes)*100, 2) 
                           for status in filtered_status_counts.columns]
    completeness_calc = round(completeness_values[0] + completeness_values[1],2)
    
    # Add chart labels and title with n and x counts
    plt.title(f"Distribution of {compleasm_odb} BUSCO Status per Sequence\nCompleteness: {completeness_calc}%")
    plt.xlabel(f"Sequences (Contig/Scaffold/Chromosome)\nIncluded={included_sequences}, Excluded={excluded_sequences}")
    plt.ylabel(f"Number of BUSCO Matches (out of {busco_genes})")
    plt.xticks(rotation=45, ha='right')  # Tilt x-axis labels for readability
    
    ax.legend(legend_labels, title="BUSCO Status", loc='upper right')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the file as an svg based on input_fasta
    output_busco_svg = input_fasta.replace(".fasta", "").replace(".fa", f"_{compleasm_odb}_busco.svg")
    plt.savefig(output_busco_svg, format="svg")
    output_busco_png = input_fasta.replace(".fasta", "").replace(".fa", f"_{compleasm_odb}_busco.png")
    plt.savefig(output_busco_png, format="png")
    
    log_print(f"PASS:\tBUSCO {compleasm_odb} plot saved: {output_busco_svg} & output_busco_png")



def egap_sample(row, results_df, CPU_THREADS, RAM_GB):
    """
    Run the Entheome Genome Hybrid Assembly Pipeline on a single sample (row of metadata),
    performing quality checks, trimming, assembly, polishing, and final curation steps.

    This function orchestrates a series of steps to:
    1.  Validate and combine Illumina reads.
    2.  Perform quality checks on raw, trimmed, and deduplicated Illumina reads.
    3.  Process and filter ONT reads (if available).
    4.  Correct ONT reads (if available) with Illumina data.
    5.  Assemble the genome using a hybrid approach (if ONT is provided) or Illumina-only
        approach (if ONT is not provided).
    6.  Polish the resulting assembly with Racon twice (if ONT is provided).
    7.  Perform a final polish with Pilon.
    8.  Remove haplotigs with purge_dedups (if ONT is provided).
    9.  RagTag scaffold and patch the assembly with a Reference Sequence (if provided)
    10. Close gaps using TGS-GapCloser (if ONT is provided) or Abyss-Sealer (Illumina-only).
    11. Evaluate the assembly using Merqury, BUSCO/Compleasm, and Quast.
    12. Update a shared `results_df` with computed metrics and quality classifications.

    Parameters:
    -----------
    row : pandas.Series or dict
        A metadata row describing the sample. Must contain keys such as:
        - "ONT_RAW_DIR": path to a directory with raw ONT files (string or NaN)       
        - "ONT_RAW_READS": path to ONT raw reads (string or NaN if none)
        - "ILLUMINA_RAW_DIR": path to a directory with raw Illumina files (string or NaN)
        - "ILLUMINA_RAW_F_READS"/"ILLUMINA_RAW_R_READS": forward/reverse Illumina reads
        - "SPECIES_ID": sample/species identifier
        - "ORGANISM_KINGDOM": e.g., "Funga", "Flora", etc.
        - "ORGANISM_KARYOTE": e.g., "Prokaryote" or "Eukaryote"
        - "COMPLEASM_1": e.g., "basidiomycota", "agaricales", etc.
        - "COMPLEASM_2": e.g., "basidiomycota", "agaricales", etc.
        - "EST_SIZE": estimated genome size (e.g., "50m" for 50 Mb)
        - "REF_SEQ": reference sequence FASTA (if any) or NaN
        (Other keys are also extracted/used within this function.)
    
    results_df : pandas.DataFrame
        A cumulative results table to which this function will append or update
        assembly and QC metrics for each sample row.

    CPU_THREADS : str or int
        Number of CPU threads to use in various external tools (e.g., FastQC, NanoPlot,
        Flye, MaSuRCA, etc.).

    RAM_GB : str or int
        Amount of RAM in gigabytes that certain tools (e.g., Trimmomatic, Pilon)
        can allocate.

    Returns:
    --------
    final_assembly_path : str
        Path to the final (optionally gzipped) assembled genome FASTA file.

    results_df : pandas.DataFrame
        The updated results DataFrame with new metrics for this specific sample (row).

    Notes:
    ------
    - This function calls various external tools (FastQC, NanoPlot, Trimmomatic, Flye,
      MaSuRCA, Racon, Pilon, RagTag, TGS-GapCloser, Abyss-Sealer, Merqury, etc.). 
      Make sure they are installed and accessible via the system $PATH.
    - BUSCO-like checks are performed through `compleasm.py`. Ensure `compleasm.py`
      is installed and that the lineage databases (e.g., `basidiomycota_odb10`,
      `agaricales_odb10`) are downloaded.
    - The function assumes certain naming conventions for FASTQ/FASTA files and
      automatically creates or skips directories and analysis outputs if they already
      exist.
    - Coverage is calculated either using a provided reference sequence size or 
      the final assembly size, depending on availability.

    Considerations:
    ---------------
    - If ONT reads are not provided, the pipeline skips steps involving long reads
      (ONT read correction, Racon polishing, TGS-GapCloser, etc.).
    - If a reference genome is not provided (`REF_SEQ` is NaN), scaffold/patch steps
      (RagTag) and some QUAST steps are skipped.
    - Make sure the input metadata row contains correct paths—missing or malformed
      paths will cause errors.
    - The function logs progress and skips steps if outputs from a prior run already
      exist to improve reproducibility and speed for subsequent runs.

    Examples:
    ---------
    # Example 1: Hybrid assembly with ONT + Illumina
    >>> import pandas as pd
    >>> data = {
    ...     "ONT_RAW_READS": "/path/to/ONT_reads.fq.gz",
    ...     "ILLUMINA_RAW_DIR": "/path/to/illumina_reads",
    ...     "SPECIES_ID": "SampleX",
    ...     "ORGANISM_KINGDOM": "Funga",
    ...     "ORGANISM_KARYOTE": "Eukaryote",
    ...     "ILLUMINA_RAW_F_READS": "/path/to/illumina_F.fq.gz",
    ...     "ILLUMINA_RAW_R_READS": "/path/to/illumina_R.fq.gz",
    ...     "COMPLEASM_1": "basidiomycota",
    ...     "COMPLEASM_2": "agaricales",
    ...     "EST_SIZE": "50m",
    ...     "REF_SEQ": "/path/to/reference.fna"
    ... }
    >>> sample_series = pd.Series(data)
    >>> results_df = pd.DataFrame()
    >>> final_path, updated_results = egap_sample(sample_series, results_df, 8, 16)
    >>> print(final_path)
    /shared/root/SampleX_EGAP_assembly.fna.gz

    # Example 2: Illumina-only assembly (no ONT, no reference)
    >>> data_illumina_only = {
    ...     "ONT_RAW_READS": float('nan'),
    ...     "ILLUMINA_RAW_DIR": float('nan'),  # or a valid path if you have multiple Illumina files
    ...     "SPECIES_ID": "SampleY",
    ...     "ORGANISM_KINGDOM": "Funga",
    ...     "ORGANISM_KARYOTE": "Eukaryote",
    ...     "ILLUMINA_RAW_F_READS": "/path/to/illuminaY_F.fq.gz",
    ...     "ILLUMINA_RAW_R_READS": "/path/to/illuminaY_R.fq.gz",
    ...     "COMPLEASM_1": "basidiomycota",
    ...     "COMPLEASM_2": "agaricales",
    ...     "EST_SIZE": "60m",
    ...     "REF_SEQ": float('nan')
    ... }
    >>> sample_series_illumina_only = pd.Series(data_illumina_only)
    >>> results_df_illumina_only = pd.DataFrame()
    >>> final_path, updated_results = egap_sample(sample_series_illumina_only, results_df_illumina_only, '4', '8')
    >>> print(final_path)
    /shared/root/SampleY_EGAP_assembly.fna.gz

    """
    cwd = os.getcwd()
    ONT_RAW_DIR = row["ONT_RAW_DIR"]
    ONT_RAW_READS = row["ONT_RAW_READS"]
    ILLU_RAW_DIR = row["ILLUMINA_RAW_DIR"]
    ILLUMINA_RAW_F_READS = row["ILLUMINA_RAW_F_READS"]
    ILLUMINA_RAW_R_READS = row["ILLUMINA_RAW_R_READS"]
    SPECIES_ID = row["SPECIES_ID"]

    if type(ONT_RAW_DIR) == str:
        shared_root = os.path.commonpath([ONT_RAW_DIR, ILLU_RAW_DIR])
        initialize_logging_environment(shared_root)
        log_print(f"Running Entheome Genome Assembly Pipeline on: {shared_root}")
        ONT_RAW_READS = ont_combine_fastq_gz(ONT_RAW_DIR)
    if type(ILLU_RAW_DIR) == str:
        shared_root = os.path.commonpath([ONT_RAW_READS, ILLU_RAW_DIR])
        initialize_logging_environment(shared_root)
        log_print(f"Running Entheome Genome Assembly Pipeline on: {shared_root}")
        illu_combined_files = illumina_extract_and_check(ILLU_RAW_DIR, SPECIES_ID)
        ILLUMINA_RAW_F_READS = illu_combined_files[0]
        ILLUMINA_RAW_R_READS = illu_combined_files[1]
    elif type(ONT_RAW_READS) == str:    
        shared_root = os.path.commonpath([ONT_RAW_READS, ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS])
        initialize_logging_environment(shared_root)
        log_print(f"Running Entheome Genome Assembly Pipeline on: {shared_root}")
    else:
        shared_root = os.path.commonpath([ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS])
        initialize_logging_environment(shared_root)
        log_print(f"Running Entheome Genome Assembly Pipeline on: {shared_root}")    

    EST_SIZE = row["EST_SIZE"]
    REF_SEQ = row["REF_SEQ"]

    if pd.notna(REF_SEQ):
        ref_seq_path = REF_SEQ
        ref_dir, ref_filename = os.path.split(ref_seq_path)
        ref_basename, ref_ext = os.path.splitext(ref_filename)
        
        # Handle double extensions like .fna.gz
        if ref_basename.endswith('.fna'):
            ref_basename, _ = os.path.splitext(ref_basename)
            ref_ext = os.path.splitext(ref_ext)[1] + ref_ext  # e.g., '.fa.gz'

        # Initialize new_ref_path
        new_ref_path = None

        # Case 1: .fna.gz
        if ref_seq_path.endswith('.fna.gz'):
            new_ref_filename = ref_basename + '.fa.gz'
            new_ref_path = os.path.join(ref_dir, new_ref_filename)
            
            if os.path.exists(ref_seq_path):
                try:
                    os.rename(ref_seq_path, new_ref_path)
                    REF_SEQ = new_ref_path  # Update REF_SEQ in the row
                    log_print(f"NOTE:\tRenamed REF_SEQ from {ref_seq_path} to {new_ref_path}")
                except Exception as e:
                    log_print(f"ERROR:\tUnable to rename {ref_seq_path} to {new_ref_path}: {e}")
            else:
                # If .fna.gz does not exist, try using existing .fa.gz
                fa_gz_path = os.path.join(ref_dir, ref_basename + '.fa.gz')
                if os.path.exists(fa_gz_path):
                    REF_SEQ = fa_gz_path
                    log_print(f"NOTE:\tOriginal file {ref_seq_path} not found. Using existing .fa.gz file: {fa_gz_path}")
                else:
                    log_print(f"NOTE:\tNeither {ref_seq_path} nor {fa_gz_path} exists. REF_SEQ cannot be processed.")
        
        # Case 2: .fna
        elif ref_seq_path.endswith('.fna'):
            new_ref_filename = ref_basename + '.fa'
            new_ref_path = os.path.join(ref_dir, new_ref_filename)
            
            if os.path.exists(ref_seq_path):
                try:
                    os.rename(ref_seq_path, new_ref_path)
                    REF_SEQ = new_ref_path  # Update REF_SEQ in the row
                    log_print(f"NOTE:\tRenamed REF_SEQ from {ref_seq_path} to {new_ref_path}")
                except Exception as e:
                    log_print(f"ERROR:\tUnable to rename {ref_seq_path} to {new_ref_path}: {e}")
            else:
                # If .fna does not exist, try using existing .fa
                fa_path = os.path.join(ref_dir, ref_basename + '.fa')
                if os.path.exists(fa_path):
                    REF_SEQ = fa_path
                    log_print(f"NOTE:\tOriginal file {ref_seq_path} not found. Using existing .fa file: {fa_path}")
                else:
                    log_print(f"NOTE:\tNeither {ref_seq_path} nor {fa_path} exists. REF_SEQ cannot be processed.")
        
        # Handle other extensions if necessary
        else:
            log_print(f"NOTE:\tREF_SEQ file {ref_seq_path} does not match expected extensions (.fna or .fna.gz). Skipping renaming.")
        
    
    else:
        log_print("NOTE:\tNo REF_SEQ provided; skipping REF_SEQ processing.")


    # Process input read files
    if pd.notna(ONT_RAW_READS):
        ONT_RAW_READS = process_read_file(ONT_RAW_READS)
    if pd.notna(ILLUMINA_RAW_F_READS):
        ILLUMINA_RAW_F_READS = process_read_file(ILLUMINA_RAW_F_READS)    
    if pd.notna(ILLUMINA_RAW_R_READS):
        ILLUMINA_RAW_R_READS = process_read_file(ILLUMINA_RAW_R_READS)

    sample_stats_dict = {"SPECIES_ID": SPECIES_ID,
                         "ONT": os.path.basename(ONT_RAW_READS) if isinstance(ONT_RAW_READS, str) else None,
                         "ILLU_F": os.path.basename(ILLUMINA_RAW_F_READS),
                         "ILLU_R": os.path.basename(ILLUMINA_RAW_R_READS),
    
                         "RAW_ILLU_TOTAL_BASES": None, # FROM FastQC
                         "RAW_ILLU_COVERAGE": None, # FROM Calculated later based on REF_SEQ or final_assembly
                         "TRIMMED_ILLU_TOTAL_BASES": None, # FROM FastQC
                         "TRIMMED_ILLU_COVERAGE": None, # FROM Calculated later based on REF_SEQ or final_assembly
                         "DEDUPED_ILLU_TOTAL_BASES": None, # FROM FastQC
                         "DEDUPED_ILLU_COVERAGE": None, # FROM Calculated later based on REF_SEQ or final_assembly

                         "RAW_ONT_READS": None, # FROM Raw NanoPlot
                         "RAW_ONT_MEAN_LENGTH": None, # FROM Raw NanoPlot
                         "RAW_ONT_MEAN_QUAL": None, # FROM Raw NanoPlot
                         "RAW_ONT_TOTAL_BASES": None, # FROM Raw NanoPlot
                         "RAW_ONT_COVERAGE": None, # FROM Calculated later based on REF_SEQ or final_assembly
        
                         "FILT_ONT_READS": None, # FROM Filtered NanoPlot
                         "FILT_ONT_MEAN_LENGTH": None, # FROM Filtered NanoPlot
                         "FILT_ONT_MEAN_QUAL": None, # FROM Filtered NanoPlot
                         "FILT_ONT_TOTAL_BASES": None, # FROM Filtered NanoPlot
                         "FILT_ONT_COVERAGE": None, # FROM Calculated later based on REF_SEQ or final_assembly
        
                         "CORRECT_ONT_READS": None, # FROM Corrected NanoPlot
                         "CORRECT_ONT_MEAN_LENGTH": None, # FROM CORRECT NanoPlot
                         "CORRECT_ONT_MEAN_QUAL": None, # FROM CORRECT NanoPlot
                         "CORRECT_ONT_TOTAL_BASES": None, # FROM CORRECT NanoPlot
                         "CORRECT_ONT_COVERAGE": None, # FROM Calculated later based on REF_SEQ or final_assembly
        
                         "KMER_COMPLETENESS": None, # FROM Merqury best k-mer completeness >97%
                         "QUAL_VAL": None, # FROM Merqury Quality Value >43
        
                         "FIRST_COMPLEASM_S": None, # FROM First Compleasm
                         "FIRST_COMPLEASM_D": None, # FROM First Compleasm
                         "FIRST_COMPLEASM_F": None, # FROM First Compleasm
                         "FIRST_COMPLEASM_M": None, # FROM First Compleasm
                         "FIRST_COMPLEASM_C": None, # FROM First Compleasm
                            
                         "SECOND_COMPLEASM_S": None, # FROM Second Compleasm
                         "SECOND_COMPLEASM_D": None, # FROM Second Compleasm
                         "SECOND_COMPLEASM_F": None, # FROM Second Compleasm
                         "SECOND_COMPLEASM_M": None, # FROM Second Compleasm
                         "SECOND_COMPLEASM_C": None, # FROM Second Compleasm
        
                         "GENOME_SIZE": None, # FROM QUAST
                         "ASSEMBLY_READS": None, # FROM QUAST
                         "ASSEMBLY_CONTIGS": None, # FROM QUAST
                         "ASSEMBLY_N50": None, # FROM QUAST
                         "ASSEMBLY_L50": None, # FROM QUAST
                         "ASSEMBLY_GC": None, # FROM QUAST
                         "MISASSEMBLIES": None, # FROM QUAST if REF_SEQ != None
                         "N_PER_100KBP": None, # FROM QUAST if REF_SEQ != None
                         "MIS_PER_100KBP": None, # FROM QUAST if REF_SEQ != None
                         "INDELS_PER_100KPB": None, # FROM QUAST if REF_SEQ != None
                            
                         "FINAL_ASSEMBLY": None}
    first_compleasm_odb = f"{row['COMPLEASM_1']}_odb10"
    second_compleasm_odb = f"{row['COMPLEASM_2']}_odb10"
    kingdom_id = row["ORGANISM_KINGDOM"].lower()
    karyote_id = row["ORGANISM_KARYOTE"].lower()
    
###############################################################################
    # Reads Pre-Processing
###############################################################################

    # FastQC Illumina Raw Reads
    raw_fastqc_dir = "/".join(ILLUMINA_RAW_F_READS.split("/")[:-1]) + "/raw_fastqc_analysis"
    if not os.path.exists(raw_fastqc_dir):
        os.makedirs(raw_fastqc_dir)
    raw_fastqc_F_out_file = os.path.join(raw_fastqc_dir, ILLUMINA_RAW_F_READS.split("/")[-1].replace(".fq.gz","_fastqc.html"))
    raw_fastqc_R_out_file = os.path.join(raw_fastqc_dir, ILLUMINA_RAW_R_READS.split("/")[-1].replace(".fq.gz","_fastqc.html"))
    if os.path.exists(raw_fastqc_F_out_file) and os.path.exists(raw_fastqc_R_out_file):
        log_print(f"SKIP:\tTrimmed FastQC outputs already exist: {raw_fastqc_F_out_file}; {raw_fastqc_R_out_file}.")
    else:
        fastqc_cmd = ["fastqc", "-o", raw_fastqc_dir, "-t", str(CPU_THREADS), ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS]
        _ = run_subprocess_cmd(fastqc_cmd, shell_check = False)
    sample_stats_dict["RAW_ILLU_TOTAL_BASES"] =  round((float(get_total_bases(raw_fastqc_F_out_file).split(" ")[0]) + float(get_total_bases(raw_fastqc_R_out_file).split(" ")[0])) / 2, 2)

    # Trimmomatic of Illumina Raw Reads
    trimmo_f_pair_path = ILLUMINA_RAW_F_READS.replace(".fq.gz", "_paired.fq.gz")
    fwd_unpaired_out = trimmo_f_pair_path.replace("paired","unpaired")
    trimmo_r_pair_path = ILLUMINA_RAW_R_READS.replace(".fq.gz", "_paired.fq.gz")
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
    trimmed_fastqc_dir = "/".join(trimmo_f_pair_path.split("/")[:-1]) + "/trimmed_fastqc_analysis"
    if not os.path.exists(trimmed_fastqc_dir):
        os.makedirs(trimmed_fastqc_dir)
    fastqc_F_out_file = os.path.join(trimmed_fastqc_dir, trimmo_f_pair_path.split("/")[-1].replace(".fq.gz","_fastqc.html"))
    fastqc_R_out_file = os.path.join(trimmed_fastqc_dir, trimmo_r_pair_path.split("/")[-1].replace(".fq.gz","_fastqc.html"))
    if os.path.exists(fastqc_F_out_file) and os.path.exists(fastqc_R_out_file):
        log_print(f"SKIP:\tTrimmed FastQC outputs already exist: {fastqc_F_out_file}; {fastqc_R_out_file}.")
    else:
        fastqc_cmd = ["fastqc", "-o", trimmed_fastqc_dir, "-t", str(CPU_THREADS), trimmo_f_pair_path, trimmo_r_pair_path]
        _ = run_subprocess_cmd(fastqc_cmd, shell_check = False)
    sample_stats_dict["TRIMMED_ILLU_TOTAL_BASES"] =  round((float(get_total_bases(fastqc_F_out_file).split(" ")[0]) + float(get_total_bases(fastqc_R_out_file).split(" ")[0])) / 2, 2)
    
    # Run bbduk on the trimmed files
    bbduk_f_map_path, bbduk_r_map_path = bbduk_map(trimmo_f_pair_path, trimmo_r_pair_path)

    # Run clumpify on the mapped files
    clump_f_dedup_path, clump_r_dedup_path = clumpify_dedup(bbduk_f_map_path, bbduk_r_map_path)
    
    # FastQC Illumina Deduplicated Reads
    dedup_fastqc_dir = "/".join(clump_f_dedup_path.split("/")[:-1]) + "/dedup_fastqc_analysis"
    if not os.path.exists(dedup_fastqc_dir):
        os.makedirs(dedup_fastqc_dir)
    fastqc_dedup_F_out_file = os.path.join(dedup_fastqc_dir, clump_f_dedup_path.split("/")[-1].replace(".fq.gz","_fastqc.html")).replace("trimmed","dedup")
    fastqc_dedup_R_out_file = os.path.join(dedup_fastqc_dir, clump_r_dedup_path.split("/")[-1].replace(".fq.gz","_fastqc.html")).replace("trimmed","dedup")
    if os.path.exists(fastqc_dedup_F_out_file) and os.path.exists(fastqc_dedup_R_out_file):
        log_print(f"SKIP:\tTrimmed FastQC outputs already exist: {fastqc_dedup_F_out_file}; {fastqc_dedup_R_out_file}.")
    else:
        fastqc_cmd = ["fastqc", "-o", dedup_fastqc_dir, "-t", str(CPU_THREADS), clump_f_dedup_path, clump_r_dedup_path]
        _ = run_subprocess_cmd(fastqc_cmd, shell_check = False)
    sample_stats_dict["DEDUPED_ILLU_TOTAL_BASES"] =  round((float(get_total_bases(fastqc_dedup_F_out_file).split(" ")[0]) + float(get_total_bases(fastqc_dedup_R_out_file).split(" ")[0])) / 2, 2)

    # NanoPlot ONT Raw Reads    
    if not pd.isna(ONT_RAW_READS):
        raw_nanoplot_dir = "/".join(ONT_RAW_READS.split("/")[:-1]) + "/raw_nanoplot_analysis"
        if not os.path.exists(raw_nanoplot_dir):
            os.makedirs(raw_nanoplot_dir)
        raw_nanoplot_out_file = os.path.join(raw_nanoplot_dir,"RawONTNanoStats.txt")
        if os.path.exists(raw_nanoplot_out_file):
            log_print(f"SKIP:\tRaw NanoPlot output already exists: {raw_nanoplot_out_file}.")
        else:    
            raw_nanoplot_cmd = [ "NanoPlot", "--fastq", ONT_RAW_READS, "-t", str(CPU_THREADS),
                                "-o", raw_nanoplot_dir, "--plots", "kde", "dot", "--loglength",
                                "--N50", "--title", "Raw ONT Reads: Preliminary Data",
                                "--prefix", "RawONT", "--verbose"]
            _ = run_subprocess_cmd(raw_nanoplot_cmd, shell_check = False)
        with open(raw_nanoplot_out_file, "r") as raw_nanostats:
            for line in raw_nanostats:
                if "Number of reads:" in line:
                    sample_stats_dict["RAW_ONT_READS"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Mean read length:" in line:
                    sample_stats_dict["RAW_ONT_MEAN_LENGTH"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Mean read quality:" in line:
                    sample_stats_dict["RAW_ONT_MEAN_QUAL"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Total bases:" in line:
                    sample_stats_dict["RAW_ONT_TOTAL_BASES"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
    
        # Filtlong ONT Raw Reads
        filtered_ONT_reads = ONT_RAW_READS.replace(".fq.gz","_filtered.fq")
        gzipped_filtered_ONT_reads = filtered_ONT_reads + ".gz"
        if os.path.exists(gzipped_filtered_ONT_reads):
            log_print(f"SKIP:\tGzipped FiltLong output already exists: {gzipped_filtered_ONT_reads}.")
        else:    
            filtlong_cmd = f"filtlong --min_length 1000 --min_mean_q 8 --target_bases 500000000 --trim -1 {clump_f_dedup_path} -2 {clump_r_dedup_path} {ONT_RAW_READS} > {filtered_ONT_reads}"
            _ = run_subprocess_cmd(filtlong_cmd, shell_check = True)
            gzip_file(filtered_ONT_reads, gzipped_filtered_ONT_reads)
        
        # NanoPlot ONT Filtered Reads
        filtered_nanoplot_dir = raw_nanoplot_dir.replace("raw","filtered")
        if not os.path.exists(filtered_nanoplot_dir):
            os.makedirs(filtered_nanoplot_dir)
        filtered_nanoplot_out_file = os.path.join(filtered_nanoplot_dir,"FilteredONTNanoStats.txt")
        if os.path.exists(filtered_nanoplot_out_file):
            log_print(f"SKIP:\tFiltered NanoPlot output already exists: {filtered_nanoplot_out_file}.")    
        else:
            filtered_nanoplot_cmd = [ "NanoPlot", "--fastq", gzipped_filtered_ONT_reads, "-t", str(CPU_THREADS),
                                "-o", filtered_nanoplot_dir, "--plots", "kde", "dot", "--loglength",
                                "--N50", "--title", "Filtered ONT Reads: Secondary Data",
                                "--prefix", "FilteredONT", "--verbose"]
            _ = run_subprocess_cmd(filtered_nanoplot_cmd, shell_check = False)
        with open(filtered_nanoplot_out_file, "r") as filt_nanostats:
            for line in filt_nanostats:
                if "Number of reads:" in line:
                    sample_stats_dict["FILT_ONT_READS"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Mean read length:" in line:
                    sample_stats_dict["FILT_ONT_MEAN_LENGTH"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Mean read quality:" in line:
                    sample_stats_dict["FILT_ONT_MEAN_QUAL"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Total bases:" in line:
                    sample_stats_dict["FILT_ONT_TOTAL_BASES"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
    
        # Error Correct ONT Filtered Reads With Illumina Reads
        corrected_ONT_out = gzipped_filtered_ONT_reads.replace("filtered.fq.gz","corrected")
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
        corrected_nanoplot_dir = raw_nanoplot_dir.replace("raw","corrected")
        if not os.path.exists(corrected_nanoplot_dir):
            os.makedirs(corrected_nanoplot_dir)
        corrected_nanoplot_out_file = os.path.join(corrected_nanoplot_dir,"CorrectedONTNanoStats.txt")
        if os.path.exists(corrected_nanoplot_out_file):
            log_print(f"SKIP:\tCorrected NanoPlot output already exists: {corrected_nanoplot_out_file}.")    
        else:    
            corrected_nanoplot_cmd = [ "NanoPlot", "--fastq", corrected_ONT_Reads, "-t", str(CPU_THREADS),
                                "-o", corrected_nanoplot_dir, "--plots", "kde", "dot", "--loglength",
                                "--N50", "--title", "Corrected ONT Reads: Tertiary Data",
                                "--prefix", "CorrectedONT", "--verbose"]
            _ = run_subprocess_cmd(corrected_nanoplot_cmd, shell_check = False)
        with open(corrected_nanoplot_out_file, "r") as corr_nanostats:
            for line in corr_nanostats:
                if "Number of reads:" in line:
                    sample_stats_dict["CORRECT_ONT_READS"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Mean read length:" in line:
                    sample_stats_dict["CORRECT_ONT_MEAN_LENGTH"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Mean read quality:" in line:
                    sample_stats_dict["CORRECT_ONT_MEAN_QUAL"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
                elif "Total bases:" in line:
                    sample_stats_dict["CORRECT_ONT_TOTAL_BASES"] = float(line.split(":")[-1].replace(" ","").replace(",","").replace("\n",""))
    
        log_print("Selecting Highest Mean Quality ONT reads...")
        highest_mean_qual_ont_reads = corrected_ONT_Reads
        if sample_stats_dict["CORRECT_ONT_MEAN_QUAL"] < sample_stats_dict["FILT_ONT_MEAN_QUAL"]:
            highest_mean_qual_ont_reads = gzipped_filtered_ONT_reads
        log_print(f"Highest Mean Quality ONT reads: {highest_mean_qual_ont_reads}")
        
###############################################################################
    # Assembly
###############################################################################   
    
    # MaSuRCA de Novo Assembly based on input Reads     
    masurca_out_dir = os.path.join(shared_root, "masurca_assembly")
    if not os.path.exists(masurca_out_dir):
        os.makedirs(masurca_out_dir)
    final_masurca_path = os.path.join(shared_root, f"{SPECIES_ID}_masurca.fasta")
    if not pd.isna(ONT_RAW_READS):
        default_assembly_path, assembly_path, ref_assembly_path = masurca_config_gen(shared_root, masurca_out_dir,
                                                                                 [highest_mean_qual_ont_reads,
                                                                                  ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS],
                                                                                 clump_f_dedup_path, clump_r_dedup_path,
                                                                                 CPU_THREADS, RAM_GB, REF_SEQ)
    else:            
        default_assembly_path, assembly_path, ref_assembly_path = masurca_config_gen(shared_root, masurca_out_dir,
                                                                                 [ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS],
                                                                                 clump_f_dedup_path, clump_r_dedup_path,
                                                                                 CPU_THREADS, RAM_GB, REF_SEQ)
    gap_assembly = os.path.join(find_ca_folder(masurca_out_dir), "9-terminator", "genome.scf.fasta")
    shutil.copyfile(gap_assembly, final_masurca_path)
    de_novo_assembly = final_masurca_path
    assembly_out_dir = masurca_out_dir

###############################################################################
    # Assembly Polishing
###############################################################################

    # 2x Polish with Racon
    if not pd.isna(ONT_RAW_READS):
        first_racon_paf = os.path.join(assembly_out_dir,"racon_round1.paf")
        if os.path.exists(first_racon_paf):
            log_print(f"SKIP:\tFirst Racon PAF already exists: {first_racon_paf}.")
        else:
            first_minimap2_cmd = f"minimap2 -t {CPU_THREADS} -x map-ont {de_novo_assembly} {highest_mean_qual_ont_reads} > {first_racon_paf}"
            _ = run_subprocess_cmd(first_minimap2_cmd, shell_check = True)
        first_racon_assembly = os.path.join(assembly_out_dir, "assembly_racon1.fasta")
        if os.path.exists(first_racon_assembly):
            log_print(f"SKIP:\tFirst Racon Assembly already exists: {first_racon_assembly}.")
        else:
            first_racon_cmd = f"racon -t {CPU_THREADS} {highest_mean_qual_ont_reads} {first_racon_paf} {de_novo_assembly} > {first_racon_assembly}"
            _ = run_subprocess_cmd(first_racon_cmd, shell_check = True)        
        second_racon_paf = os.path.join(assembly_out_dir,"racon_round2.paf")
        if os.path.exists(second_racon_paf):
            log_print(f"SKIP:\tSecond Racon PAF already exists: {second_racon_paf}.")
        else:
            second_minimap2_cmd = f"minimap2 -t {CPU_THREADS} -x map-ont {de_novo_assembly} {highest_mean_qual_ont_reads} > {second_racon_paf}"
            _ = run_subprocess_cmd(second_minimap2_cmd, shell_check = True)
        second_racon_assembly = os.path.join(assembly_out_dir, "assembly_racon2.fasta")
        if os.path.exists(second_racon_assembly):
            log_print(f"SKIP:\tSecond Racon Assembly already exists: {second_racon_assembly}.")
        else:
            second_racon_cmd = f"racon -t {CPU_THREADS} {highest_mean_qual_ont_reads} {second_racon_paf} {de_novo_assembly} > {second_racon_assembly}"
            _ = run_subprocess_cmd(second_racon_cmd, shell_check = True)
        pilon_polish_target = second_racon_assembly
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
    pilon_renamed_fasta = de_novo_assembly.replace("_masurca.fasta","_pilon.fasta")
    if os.path.exists(pilon_out_fasta):
        log_print(f"SKIP:\tPilon Polished Assembly already exists: {pilon_out_fasta}.")
    else:
        pilon_cmd = ["pilon", f"-Xmx{RAM_GB}g",
                     "--genome", pilon_polish_target,
                     "--frags", polish_bam,
                     "--output", pilon_out_prefix,
                     "--outdir", pilon_out_dir,
                     "--changes", "--vcf",
                     "--chunksize", str(5000000)]
        _ = run_subprocess_cmd(pilon_cmd, shell_check = False)
    if os.path.exists(pilon_renamed_fasta):
        log_print(f"SKIP:\tPolished Assembly already exists: {pilon_renamed_fasta}.")
    else:        
        shutil.copyfile(pilon_out_fasta, pilon_renamed_fasta)
    pilon_out_fasta = pilon_renamed_fasta
    
###############################################################################
    # Assembly Curation
###############################################################################    
    
# Purge haplotigs and overlaps with purge_dups (if ONT Reads exist)
    if not pd.isna(ONT_RAW_READS):
        pd_work_dir = os.path.join(shared_root,"purge_dups_work")
        if not os.path.exists(pd_work_dir):
            os.mkdir(pd_work_dir)
        os.chdir(pd_work_dir)
        pd_fofn = os.path.join(pd_work_dir,"ont_reads.fofn")
        with open(pd_fofn, "w") as f:
            f.write(highest_mean_qual_ont_reads)
        pd_json = os.path.join(pd_work_dir, "purge_dups_config.json")
        dup_purged_assembly = os.path.join(pd_work_dir, f"{SPECIES_ID}_pilon/seqs/{SPECIES_ID}_pilon.purged.fa")
        pd_config_path = find_file("pd_config.py")
        if os.path.exists(pd_json):
            log_print(f"SKIP:\tPurge Dupes JSON already exists: {pd_json}.")
        else:
            purge_dupes_config_cmd = ["python", pd_config_path, pilon_out_fasta,
                                      pd_fofn, "-l", pd_work_dir, "-n", pd_json]
            _ = run_subprocess_cmd(purge_dupes_config_cmd, shell_check = False)
        pd_path = pd_config_path.replace("scripts/pd_config.py","bin") # "./purge_dups/bin"
        if os.path.exists(dup_purged_assembly):
            log_print(f"SKIP:\tDupe Purged Assembly already exists: {dup_purged_assembly}.")
        else:
            purge_dupes_cmd = ["python", find_file("run_purge_dups.py"),
                               pd_json, pd_path, SPECIES_ID, "-p", "bash"]
            print(purge_dupes_cmd)
            _ = run_subprocess_cmd(purge_dupes_cmd, shell_check = False)       
        os.chdir(cwd)
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
        final_assembly_path = os.path.join(shared_root, f"{SPECIES_ID}_EGAP_assembly.fa")
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
        output_prefix = "/".join(shared_root.split("/")[:-1]) + f"/{ragtag_ref_assembly.split('/')[-1].replace('_pilon_ragtag_final.fasta','_sealed')}"
        sealer_output_file = f"{output_prefix}_scaffold.fa"        
        if os.path.isfile(sealer_output_file):        
            log_print(f"SKIP:\tABySS sealer output file already exists: {sealer_output_file}.")
        else:
            kmer_sizes = [55,75,95]
            abyss_sealer_cmd = ["abyss-sealer", "-o", output_prefix, "-S", ragtag_ref_assembly,
                                "-L", str(400), "-G", str(1000),
                                "-j", str(CPU_THREADS), "-b", "500M"]            
            for k in kmer_sizes:
                abyss_sealer_cmd.extend(["-k", str(k)])
            abyss_sealer_cmd.extend([clump_f_dedup_path, clump_r_dedup_path])
            _ = run_subprocess_cmd(abyss_sealer_cmd, shell_check = False)
        final_assembly_path = sealer_output_file

###############################################################################
    # Final Assembly Assessment
###############################################################################

    # Use Illumina Reads with Merqury
    os.chdir(shared_root)
    meryl_db_out = os.path.join(os.path.dirname(shared_root), f"{SPECIES_ID}.meryl")
    if os.path.exists(meryl_db_out):
        log_print(f"SKIP:\tMeryl Database already exists: {meryl_db_out}.")
    else:
        meryl_cmd = ["meryl", "k=21", "count", "output", meryl_db_out,
                      clump_f_dedup_path, clump_r_dedup_path]
        _ = run_subprocess_cmd(meryl_cmd, shell_check = False)
    
    merqury_out = meryl_db_out.replace(".meryl","_merqury")
    if not os.path.exists(merqury_out):
        os.makedirs(merqury_out)
    merqury_cmd = ["merqury.sh", meryl_db_out, final_assembly_path, merqury_out]
    _ = run_subprocess_cmd(merqury_cmd, shell_check = False)
    os.chdir(cwd)
    
    # Run Compelasm using first_odb10
    first_compleasm_dir = os.path.join(os.path.dirname(shared_root), f"{SPECIES_ID}_{first_compleasm_odb}_compleasm")
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
        _ = run_subprocess_cmd(first_compleasm_cmd, shell_check = False)
    plot_busco(first_compleasm_odb, first_compleasm_tsv, final_assembly_path)

    # Run Compelasm using second_odb10
    second_compleasm_dir = os.path.join(os.path.dirname(shared_root), f"{SPECIES_ID}_{second_compleasm_odb}_compleasm")
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
        _ = run_subprocess_cmd(second_compleasm_cmd, shell_check = False)
    plot_busco(second_compleasm_odb, second_compleasm_tsv, final_assembly_path)

    # Run Quast
    quast_dir = os.path.join(os.path.dirname(shared_root), f"{SPECIES_ID}_quast")
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
        _ = run_subprocess_cmd(quast_cmd, shell_check = False)
    final_gz_assembly_path = final_assembly_path + ".gz"
    gzip_file(final_assembly_path, final_gz_assembly_path)

###############################################################################
# Assessment of Input Reads, Process, & Final Assembly
###############################################################################

    # Parse First Compleasm BUSCO Report
    with open(first_compleasm_summary, "r") as first_compleasm_file:
        for line in first_compleasm_file:
            if "S:" in line:
                sample_stats_dict["FIRST_COMPLEASM_S"] = float(line.split("S:")[-1].split(", ")[0].replace("\n","").replace("%",""))
            elif "D:" in line:
                sample_stats_dict["FIRST_COMPLEASM_D"] = float(line.split("D:")[-1].split(", ")[0].replace("\n","").replace("%",""))
            elif "F:" in line:
                sample_stats_dict["FIRST_COMPLEASM_F"] = float(line.split("F:")[-1].split(", ")[0].replace("\n","").replace("%",""))
            elif "M:" in line:
                sample_stats_dict["FIRST_COMPLEASM_M"] = float(line.split("M:")[-1].split(", ")[0].replace("\n","").replace("%",""))
    sample_stats_dict["FIRST_COMPLEASM_C"] = sample_stats_dict["FIRST_COMPLEASM_S"] + sample_stats_dict["FIRST_COMPLEASM_D"]

    # Parse Second Compleasm BUSCO Report
    with open(second_compleasm_summary, "r") as second_compleasm_file:
        for line in second_compleasm_file:
            if "S:" in line:
                sample_stats_dict["SECOND_COMPLEASM_S"] = float(line.split("S:")[-1].split(", ")[0].replace("\n","").replace("%",""))
            elif "D:" in line:
                sample_stats_dict["SECOND_COMPLEASM_D"] = float(line.split("D:")[-1].split(", ")[0].replace("\n","").replace("%",""))
            elif "F:" in line:
                sample_stats_dict["SECOND_COMPLEASM_F"] = float(line.split("F:")[-1].split(", ")[0].replace("\n","").replace("%",""))
            elif "M:" in line:
                sample_stats_dict["SECOND_COMPLEASM_M"] = float(line.split("M:")[-1].split(", ")[0].replace("\n","").replace("%",""))
    sample_stats_dict["SECOND_COMPLEASM_C"] = sample_stats_dict["SECOND_COMPLEASM_S"] + sample_stats_dict["SECOND_COMPLEASM_D"]

    # Parse QUAST Report
    with open(quast_report_tsv, "r") as quast_file:
        for line in quast_file:
            if "Total length (>= 0 bp)" in line:
                sample_stats_dict["GENOME_SIZE"] = float(line.split("\t")[-1].replace("\n",""))
            elif "# contigs" in line:
                sample_stats_dict["ASSEMBLY_CONTIGS"] = float(line.split("\t")[-1].replace("\n",""))
            elif "N50" in line:
                sample_stats_dict["ASSEMBLY_N50"] = float(line.split("\t")[-1].replace("\n",""))
            elif "L50" in line:
                sample_stats_dict["ASSEMBLY_L50"] = float(line.split("\t")[-1].replace("\n",""))
            elif "GC (%)" in line:
                sample_stats_dict["ASSEMBLY_GC"] = float(line.split("\t")[-1].replace("\n",""))
            if REF_SEQ != None:
                if "# misassemblies" in line:
                    sample_stats_dict["MISASSEMBLIES"] = float(line.split("\t")[-1].replace("\n",""))
                elif "# N's per 100 kbp" in line:
                    sample_stats_dict["N_PER_100KBP"] = float(line.split("\t")[-1].replace("\n",""))
                elif "# mismatches per 100 kbp" in line:
                    sample_stats_dict["MIS_PER_100KBP"] = float(line.split("\t")[-1].replace("\n",""))
                elif "# indels per 100 kbp" in line:
                    sample_stats_dict["INDELS_PER_100KPB"] = float(line.split("\t")[-1].replace("\n",""))

    # Calculate Coverage Based on if REF_SEQ != None then use final_assembly
    if not pd.isna(REF_SEQ):
        ref_total_bases = 0
        for record in SeqIO.parse(REF_SEQ, "fasta"):
            ref_total_bases += len(record.seq)
        sample_stats_dict["RAW_ILLU_COVERAGE"] = round(sample_stats_dict["RAW_ILLU_TOTAL_BASES"] / ref_total_bases, 2) # Calculated based on REF_SEQ or final_assembly
        sample_stats_dict["TRIMMED_ILLU_COVERAGE"] = round(sample_stats_dict["TRIMMED_ILLU_TOTAL_BASES"] / ref_total_bases, 2) # Calculated based on REF_SEQ or final_assembly
        sample_stats_dict["DEDUPED_ILLU_COVERAGE"] = round(sample_stats_dict["DEDUPED_ILLU_TOTAL_BASES"] / ref_total_bases, 2) # Calculated based on REF_SEQ or final_assembly
        if not pd.isna(ONT_RAW_READS):
            sample_stats_dict["RAW_ONT_COVERAGE"] = round(sample_stats_dict["RAW_ONT_TOTAL_BASES"] / ref_total_bases, 2) # Calculated based on REF_SEQ or final_assembly
            sample_stats_dict["FILT_ONT_COVERAGE"] = round(sample_stats_dict["FILT_ONT_TOTAL_BASES"] / ref_total_bases, 2) # Calculated based on REF_SEQ or final_assembly
            sample_stats_dict["CORRECT_ONT_COVERAGE"] = round(sample_stats_dict["CORRECT_ONT_TOTAL_BASES"] / ref_total_bases, 2) # Calculated based on REF_SEQ or final_assembly
    else:
        sample_stats_dict["RAW_ILLU_COVERAGE"] = round(sample_stats_dict["RAW_ILLU_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2) # Calculated based on REF_SEQ or final_assembly
        sample_stats_dict["TRIMMED_ILLU_COVERAGE"] =  round(sample_stats_dict["TRIMMED_ILLU_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2) # Calculated based on REF_SEQ or final_assembly
        sample_stats_dict["DEDUPED_ILLU_COVERAGE"] = round(sample_stats_dict["DEDUPED_ILLU_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2) # Calculated based on REF_SEQ or final_assembly
        if not pd.isna(ONT_RAW_READS):
            sample_stats_dict["RAW_ONT_COVERAGE"] =  round(sample_stats_dict["RAW_ONT_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2) # Calculated based on REF_SEQ or final_assembly
            sample_stats_dict["FILT_ONT_COVERAGE"] = round(sample_stats_dict["FILT_ONT_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2) # Calculated based on REF_SEQ or final_assembly
            sample_stats_dict["CORRECT_ONT_COVERAGE"] = round(sample_stats_dict["CORRECT_ONT_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2) # Calculated based on REF_SEQ or final_assembly
    sample_stats_dict["FINAL_ASSEMBLY"] = final_assembly_path
    quality_classifications = classify_assembly(sample_stats_dict)
    for metric, classification in quality_classifications.items():
        log_print(f"{metric}: {classification}")
    result_row = pd.DataFrame([quality_classifications], index=[index])
    results_df = pd.concat([results_df, result_row])
    for key, value in sample_stats_dict.items():
        input_csv_df.loc[index, key] = value
    log_print(f"Assessment of Final Assembly: {final_assembly_path}")
    log_print(f"PASS:\tEGAP Final Assembly Complete: {final_assembly_path}")
    log_print("This was produced with the help of the Entheogen Genome (Entheome) Foundation\n")
    log_print("If this was useful, please support us at https://entheome.org/\n")
    print("\n\n\n")    
    return final_assembly_path, results_df


if __name__ == "__main__":
    # Argument Parsing
    parser = argparse.ArgumentParser(description="Run Entheome Genome Assembly Pipeline (EGAP)")

    # Default values
    default_input_csv = None
    default_raw_ont_dir = None
    default_ont_reads = None
    default_raw_illu_dir = None
    default_raw_illu_reads_1 = "/mnt/d/ENTHEOME/Ps_mexicana/IlluminaPE150/SRR22434202_1.fq.gz"
    default_raw_illu_reads_2 = "/mnt/d/ENTHEOME/Ps_mexicana/IlluminaPE150/SRR22434202_2.fq.gz"
    default_species_id = "Ps_mexicana" # Format: <2-letters of Genus>_<full species name>
    default_organism_kingdom = "Funga"
    default_organism_karyote = "Eukaryote"
    default_compleasm_1 = "basidiomycota"
    default_compleasm_2 = "agaricales"
    default_estimated_genome_size = "60m"
    default_reference_sequence= None
    default_percent_resources = 0.75
    
    # Add arguments with default values
    parser.add_argument("--input_csv", "-csv",
                        type = str, default = default_input_csv,
                        help = f"Path to a csv containing multiple sample data. (default: {default_input_csv})")
    parser.add_argument("--raw_ont_dir", "-odir",
                        type = str, default = default_raw_ont_dir,
                        help = f"Path to a directory containing all Raw ONT Reads. (default: {default_raw_ont_dir})")
    parser.add_argument("--raw_ont_reads", "-i0",
                        type = str, default = default_ont_reads,
                        help = f"Path to the combined Raw ONT fastq reads. (default: {default_ont_reads})")
    parser.add_argument("--raw_illu_dir", "-idir",
                        type = str, default = default_raw_illu_dir,
                        help = f"Path to a directory containing all Raw Illumina Reads. (default: {default_raw_illu_dir})")
    parser.add_argument("--raw_illu_reads_1", "-i1",
                        type = str, default = default_raw_illu_reads_1,
                        help = f"Path to the Raw Forward Illumina Reads. (default: {default_raw_illu_reads_1})")
    parser.add_argument("--raw_illu_reads_2", "-i2",
                        type = str, default = default_raw_illu_reads_2,
                        help = f"Path to the Raw Reverse Illumina Reads. (default: {default_raw_illu_reads_2})")
    parser.add_argument("--species_id", "-ID",
                        type = str, default = default_species_id,
                        help = f"Species ID formatted: <2-letters of Genus>_<full species name>. (default: {default_species_id})")
    parser.add_argument("--organism_kingdom", "-Kg",
                        type = str, default = default_organism_kingdom,
                        help = f"Phylogenetic Kingdom the organism belongs to. (default: {default_organism_kingdom})")
    parser.add_argument("--organism_karyote", "-Ka",
                        type = str, default = default_organism_karyote,
                        help = f"Karyote type of the organism. (default: {default_organism_karyote})")
    parser.add_argument("--compleasm_1", "-c1",
                        type = str, default = default_compleasm_1,
                        help = f"Name of the first organism compleasm/BUSCO database to compare to. (default: {default_compleasm_1})")
    parser.add_argument("--compleasm_2", "-c2",
                        type = str, default = default_compleasm_2,
                        help = f"Name of the second organism compleasm/BUSCO database to compare to. (default: {default_compleasm_2})")
    parser.add_argument("--est_size", "-es",
                        type = str, default = default_estimated_genome_size,
                        help="Estimaged size of the genome in Mbp (aka million-base-pairs). (default: {default_estimated_genome_size})")
    parser.add_argument("--ref_seq", "-rf",
                        type = str, default = default_reference_sequence,
                        help = "Path to the reference genome for assembly. (default: {default_reference_sequence})")
    parser.add_argument("--percent_resources", "-R",
                        type = float, default = default_percent_resources,
                        help = f"Percentage of resources for processing. (default: {default_percent_resources})")
    
    # Parse the arguments into Data Frame for iterative processing
    args = parser.parse_args()
    INPUT_CSV = args.input_csv
    if INPUT_CSV != None:
        input_csv_df = pd.read_csv(INPUT_CSV)
    else:
        sample_dict = {"ONT_RAW_DIR": [args.raw_ont_dir],
                       "ONT_RAW_READS": [args.raw_ont_reads],
                       "ILLUMINA_RAW_DIR": [args.raw_illu_dir],
                       "ILLUMINA_RAW_F_READS": [args.raw_illu_reads_1],
                       "ILLUMINA_RAW_R_READS": [args.raw_illu_reads_2],
                       "SPECIES_ID": [args.species_id],
                       "ORGANISM_KINGDOM": [args.organism_kingdom],
                       "ORGANISM_KARYOTE": [args.organism_karyote],
                       "COMPLEASM_1": [args.compleasm_1],
                       "COMPLEASM_2": [args.compleasm_2],
                       "EST_SIZE": [args.est_size],
                       "REF_SEQ": [args.ref_seq]}
        input_csv_df = pd.DataFrame.from_dict(sample_dict)
    PERCENT_RESOURCES = args.percent_resources
    CPU_THREADS, RAM_GB = get_resource_values(PERCENT_RESOURCES)

    # Parse samples in sample_dict & input_csv_df
    results_df = pd.DataFrame()
    print(input_csv_df)
    for index, row in input_csv_df.iterrows():
        final_assembly_path, results_df = egap_sample(row, results_df, CPU_THREADS, RAM_GB)
    
    # plot_classification_table(results_df)

    if not pd.isna(INPUT_CSV):
        final_csv_filename = INPUT_CSV.replace(".csv", "_final_assembly_stats.csv")
    else:
        final_csv_filename = os.path.join("/".join(os.path.dirname(default_raw_illu_reads_1).split("/")[:-1]), f"{args.species_id}_final_assembly_stats.csv")
    input_csv_df.to_csv(final_csv_filename, index=False)
