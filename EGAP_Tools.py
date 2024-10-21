# -*- coding: utf-8 -*-
"""
Created on Wed Aug 2 11:39:38 2023

@author: ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting with the help of ChatGPT 4.0

"""
# Base Python Imports
import math, platform, os, subprocess, shutil, multiprocessing
from datetime import datetime


# Required Python Imports
import psutil


# Global variables
CPU_THREADS = 1
RAM_GB = 1
DEFAULT_LOG_FILE = None
ENVIRONMENT_TYPE = None


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
        # If using numerical suffixes, increment until a new filename is found
        counter = 1
        new_log_file_path = f"{log_file_path.rsplit('.', 1)[0]}_{counter}.txt"
        while os.path.exists(new_log_file_path):
            counter += 1
            new_log_file_path = f"{log_file_path.rsplit('.', 1)[0]}_{counter}.txt"
        log_file_path = new_log_file_path
    else:
        # Clear the existing log file or create a new one
        open(log_file_path, 'w').close()
    
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
    # Access the global variable
    global DEFAULT_LOG_FILE
    # ANSI escape sequences for colors
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
    message = f'[{now:%Y-%m-%d %H:%M:%S}]\t{input_message}'

    # Determine the print color
    message_type_dict = {'NOTE': 'blue',
                         'CMD': 'cyan',
                         'ERROR': 'red',
                         'WARN': 'yellow',
                         'PASS': 'green',}
    print_color = 'white'  # Default color
    for key, value in message_type_dict.items():
        if key.lower() in input_message.lower():
            print_color = value
            break

    # Writing the message to the log file
    try:
        with open(log_file, 'a') as file:
            print(message, file=file)
    except TypeError:
        print(f"UNLOGGED ERROR:\tUnable to load the log file provided: {log_file}")

    # Print the message with color
    color_code = COLORS.get(print_color, COLORS['white'])
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

    # Change the file extension to .txt
    input_file_path = f"{INPUT_FOLDER}/{INPUT_FOLDER.split('/')[-1]}_log.txt"
 
    # Determine the operating system
    os_name = platform.system()
    
    # Depending on the operating system, handle the input file path differently
    if os_name == "Windows":
        print('UNLOGGED:\tWINDOWS ENVIRONMENT')
        ENVIRONMENT_TYPE = "WIN"
    elif os_name in ["Linux", "Darwin"]:  # Darwin is the system name for macOS
        drive, path_without_drive = os.path.splitdrive(input_file_path)
        if drive:
            drive_letter = drive.strip(":\\/")
            path_without_drive_mod = path_without_drive.replace("\\", "/")
            input_file_path = f'/mnt/{drive_letter.lower()}{path_without_drive_mod}'
        print('UNLOGGED:\tLINUX/WSL/MAC ENVIRONMENT')
        ENVIRONMENT_TYPE = "LINUX/WSL/MAC"
    else:
        print(f'UNLOGGED ERROR:\tUnsupported OS: {os_name}')
        return
    
    # Generate the log file based on the input file
    run_log = generate_log_file(input_file_path, use_numerical_suffix=False)
    
    # Set the default log file to the generated run_log
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
        - Waits for the command to complete and checks the return code to determine success or failure.
        - Logs the command, its real-time output, and any execution errors.

    Usage and Considerations:
        - Useful for executing commands where live feedback is important, especially for long-running commands.
        - Requires careful use of 'shell_check' due to potential security risks with shell commands.

    Example:
        run_subprocess_cmd(["ls", "-l"], shell_check=False)
        # This would execute 'ls -l' and display its output in real-time, while handling logging.
    """
    if isinstance(cmd_list, str):
        log_print(f"CMD:\t{cmd_list}")    
        process = subprocess.Popen(cmd_list, shell=shell_check, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    else:
        log_print(f"CMD:\t{' '.join(cmd_list)}")    
        process = subprocess.Popen(cmd_list, shell=shell_check, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    # Iterating over the output
    for line in process.stdout:
        print(line, end='')

    # Wait for the process to complete and fetch the return code
    process.wait()

    if process.returncode != 0:
        log_print(f"ERROR:\tCommand failed with return code {process.returncode}")
    else:
        log_print(f"PASS:\tSuccessfully processed command: {' '.join(cmd_list)}" if isinstance(cmd_list, list) else cmd_list)
        

def get_resource_values(PERCENT_RESOURCES):
    """
    Converts user input PERCENT_RESOURCES into usuable cpu_threads and ram_gb values.

    Arg:
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
    # Get the number of CPUs available on the system
    num_cpus = multiprocessing.cpu_count()
    
    # Get the amount of RAM (GB) currently available
    mem_info = psutil.virtual_memory()
    
    # Calculate the number of threads as 80% of available CPUs & RAM
    cpu_threads = int(math.floor(num_cpus * PERCENT_RESOURCES))
    ram_gb = int(mem_info.total / (1024.0 ** 3) * PERCENT_RESOURCES)
    
    return cpu_threads, ram_gb 


def generate_sam(assembly_file, input_reads, CPU_THREADS):
    """
    Generates a SAM Map of the ONT Cleaned Flye de novo Assembly with the Trimmed Illumina Paired Forward & Reverse Reads.

    Args:
        assembly_file (str): Path to an Assembly FASTA file.
        input_reads (list): List of paths to the FASTQ reads files.
        CPU_THREADS (int): Number of threads available.

    Returns:
        output_sam (str): A SAM Map File of the ONT Cleaned Flye de novo Assembly with the Trimmed Illumina Paired Forward & Reverse Reads.
    """
    log_print(f"Generating SAM from: {input_reads}")    
    log_print(f"Length of items in Input Reads: {len(input_reads)}")

    if len(input_reads) == 1:
        output_sam = assembly_file.replace(".fasta", "_ont_minimap_aligned.sam")
        if os.path.isfile(output_sam):
            log_print(f"PASS:\tSkipping SAM Generation, {output_sam} already exists")
        else:
            sam_cmd = f"minimap2 -ax sr {assembly_file} {input_reads[0]} > {output_sam}"
            
            _ = run_subprocess_cmd(sam_cmd, True)
        
    elif len(input_reads) == 2:
        output_sam = assembly_file.replace(".fasta", "_illu_minimap_aligned.sam")
        if os.path.isfile(output_sam):
            log_print(f"PASS:\tSkipping SAM Generation, {output_sam} already exists")
        else:
            sam_cmd = f"minimap2 -ax sr {assembly_file} {input_reads[0]} {input_reads[1]} > {output_sam}"
        
            _ = run_subprocess_cmd(sam_cmd, True)
            
    return output_sam


def convert_sam_to_bam(sam_file, CPU_THREADS, output_bam=None):
    """
    Generates Binary Alignment Map (BAM file) from the Illumina Sequence Alignment Map with BWA.

    Args:
        sam_file (str): A SAM Map File of the ONT Cleaned Flye de novo Assembly with the Trimmed Illumina Paired Forward & Reverse Reads.
        CPU_THREADS (int): Number of threads available.

    Returns:
        output_bam (str): A Binary Alignment Map (BAM file) from the Illumina Sequence Alignment Map with BWA. 
    """
    log_print(f'Converting {sam_file} into BAM file with Samtools...')

    # Check if the BAM file already exists
    if output_bam == None:
        output_bam = sam_file.replace(".sam", "_sorted.bam")    

    if os.path.isfile(output_bam):
        log_print(f"PASS:\tSkipping BAM Generation, {output_bam} already exists")
    else:        
        # First command: samtools view
        view_cmd = f"samtools view -@ {CPU_THREADS} -S {sam_file} -b > {output_bam}"
        log_print(f"CMD:\t{view_cmd}")
    
        view_result = subprocess.run(view_cmd, shell=True, stderr=subprocess.PIPE, executable="/bin/bash")
    
        # Check if SAM->BAM conversion ran successfully
        if view_result.returncode != 0:
            log_print(f"ERROR:\t{view_result.stderr.decode()}")
            return None
        else:
            log_print("PASS:\tSuccessfully generated BAM file")

        # Second command: samtools sort
        sort_cmd = f"bamtools sort -in {output_bam} -out {output_bam}"
        log_print(f"CMD:\t{sort_cmd}")
        sort_result = subprocess.run(sort_cmd, shell=True, stderr=subprocess.PIPE, executable="/bin/bash")
        
        # Check if BAM sorting ran successfully
        if sort_result.returncode != 0:
            log_print(f"ERROR:\t{sort_result.stderr.decode()}")
            return None
        else:
            log_print("PASS:\tSuccessfully sorted BAM file")
        
        # Third command: samtools index
        index_cmd = f"samtools index {output_bam}"
        log_print(f"CMD:\t{index_cmd}")
        index_result = subprocess.run(index_cmd, shell=True, stderr=subprocess.PIPE, executable="/bin/bash")
        
        # Check if BAM indexing ran successfully
        if index_result.returncode != 0:
            log_print(f"ERROR:\t{index_result.stderr.decode()}")
            return None
        else:
            log_print("PASS:\tSuccessfully indexed BAM file")

        log_print("PASS:\tSuccessfully generated BAM file")
    
    return output_bam


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
        root_directory = 'C:\\'  # Adjust if necessary for different drives
    elif ENVIRONMENT_TYPE in ["LINUX/WSL/MAC"]:
        root_directory = '/'
    else:
        raise ValueError("Unknown ENVIRONMENT_TYPE")

    for root, dirs, files in os.walk(root_directory):
        if filename in files:
            return os.path.join(root, filename)
    return None


def cleanup(keep_paths, input_folder, log_file):
    """
    Removes unwanted files and folders after successful QC checks.

    Args:
        keep_paths (list): List of Path strings of files to keep.
        input_folder (str): Path to the input folder containing original data and temporary files.
        log_file (str): Path to the log file to keep.
    """

    # Traverse the input folder and delete files/folders not in the keep list
    for root, dirs, files in os.walk(input_folder, topdown=False):
        # Check each file in the directory
        for name in files:
            file_path = os.path.join(root, name)
            if file_path not in keep_paths:
                os.remove(file_path)
                log_print(f"Deleted file: {file_path}")
        
        # Check each directory
        for name in dirs:
            dir_path = os.path.join(root, name)
            if dir_path not in keep_paths:
                shutil.rmtree(dir_path)
                log_print(f"Deleted directory: {dir_path}")

    log_print("Cleanup completed: Removed all unnecessary files and folders.")
    

# Debuging Main Space & Example
if __name__ == "__main__":
    print("EGAP_Tools.py DEBUG")

    # Example for initializing the logging environment
    print("\n-- Testing initialize_logging_environment --")
    test_input_folder = "test_data"
    initialize_logging_environment(test_input_folder)
    print(f"Default Log File: {DEFAULT_LOG_FILE}")

    # Example for running a subprocess command
    log_print("\n-- Testing run_subprocess_cmd --")
    test_command = ["echo", "Hello, world!"]
    run_subprocess_cmd(test_command, shell_check=False)

    # Example for converting resources to usable values
    log_print("\n-- Testing get_resource_values --")
    test_percent_resources = 0.8
    cpu_threads, ram_gb = get_resource_values(test_percent_resources)
    log_print(f"CPU Threads: {cpu_threads}, RAM (GB): {ram_gb}")

    # TODO: Example for generating a sam file from a reads_list and assembly_fasta

    # TODO: Example for converting a sam into a bam file

    # Example for finding a file
    log_print("\n-- Testing find_file --")
    test_filename = "EGAP.py"
    found_file = find_file(test_filename)
    log_print(f"Found File: {found_file}")

    # TODO: Example for file cleanup
