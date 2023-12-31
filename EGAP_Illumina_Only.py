# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 12:20:03 2023

@author: ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting

Command Line Example:
    python Illumina_Pipe.py -i /path/to/folder -r FLOAT

The -i, --input_folder must have input FASTQ.GZ file, matching primers.txt and Index.txt
The -r, --resources Percentage of resources to use (0.01-1.00; default: 0.2) 

"""
import os, subprocess, multiprocessing, math, platform, argparse, gzip, shutil
from datetime import datetime

# Global output_area variable
CPU_THREADS = 1
DEFAULT_LOG_FILE = None

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
    - Useful for managing log file versions without overwriting existing logs.
    - The numerical suffix increments for each new file created in the same location.
    - When 'use_numerical_suffix' is False, it refreshes the log file by clearing existing content.
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
        - The function uses a global default log file if none is specified.
        - Timestamps each log entry for easy tracking.
        - Utilizes color coding in the console to distinguish between different types of messages (e.g., errors, warnings).
        - Supports color coding for specific message types: NOTE, CMD, ERROR, WARN, and PASS.
        - Falls back to default (white) color if the message type is unrecognized.
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
    message_type_dict = {
        'NOTE': 'blue',
        'CMD': 'cyan',
        'ERROR': 'red',
        'WARN': 'yellow',
        'PASS': 'green',
    }
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
        input_file_path (str): Path to the input file which influences the log file generation.

    Global Variables:
        DEFAULT_LOG_FILE (str): The default path for the log file used throughout the logging process.

    Notes:
        - Supports Windows, Linux/WSL, and Darwin (macOS) environments.
        - Prints unlogged messages to the console regarding environment detection and file existence.
        - Modifies the global DEFAULT_LOG_FILE variable.
    """
    global DEFAULT_LOG_FILE

    # Change the file extension to .txt
    input_file_path = f"{INPUT_FOLDER}/{INPUT_FOLDER.split('/')[-1]}_log.txt"
 
    # Determine the operating system
    os_name = platform.system()
    
    # Depending on the operating system, handle the input file path differently
    if os_name == "Windows":
        print('UNLOGGED:\tWINDOWS ENVIRONMENT')
    elif os_name in ["Linux", "Darwin"]:  # Darwin is the system name for macOS
        drive, path_without_drive = os.path.splitdrive(input_file_path)
        if drive:
            drive_letter = drive.strip(":\\/")
            path_without_drive_mod = path_without_drive.replace("\\", "/")
            input_file_path = f'/mnt/{drive_letter.lower()}{path_without_drive_mod}'
        print('UNLOGGED:\tLINUX/WSL/MAC ENVIRONMENT')
    else:
        print(f'UNLOGGED ERROR:\tUnsupported OS: {os_name}')
        return
    
    # Generate the log file based on the input file
    run_log = generate_log_file(input_file_path, use_numerical_suffix=False)
    
    # Set the default log file to the generated run_log
    DEFAULT_LOG_FILE = run_log

def run_subprocess_cmd(cmd_list, shell_check):
    """
    Executes a command using the subprocess module and logs the output.
    
    This function runs a given command (either as a string or a list of strings) using the subprocess module. 
    It logs the command being executed and handles any output or errors generated by the command. The function 
    supports both single-string commands and commands split into arguments within a list.
    
    Parameters:
        cmd_list (str or list of str): The command to be executed. Can be a single string or a list of strings 
                                       representing the command and its arguments.
        shell_check (bool): If True, the command is executed through the shell. This is necessary for some 
                            commands, especially those that are built into the shell or involve shell features 
                            like wildcard expansion.
    
    Notes:
        - Uses the 'log_print' function to log and print the command and its output.
        - Captures and logs both standard output and standard error.
        - Logs an error message if the command execution fails (non-zero return code).
        - When 'cmd_list' is a list, it joins the list elements into a single string for logging purposes.
        - Requires careful use of 'shell_check' due to potential security risks when executing shell commands.
        - The function is designed to be versatile for various command execution scenarios in different environments.
    """
    if isinstance(cmd_list, str):
        log_print(f"CMD:\t{cmd_list}")    
        process = subprocess.run(cmd_list, text=True, shell=shell_check, capture_output=True)
        if process.returncode != 0:
            log_print(f"ERROR:\t{process.stderr}")
        else:
            log_print(f"PASS:\tSuccessfully processed command: {cmd_list}")
        return process.stdout
    else:        
        log_print(f"CMD:\t{' '.join(cmd_list)}")    
        process = subprocess.run(cmd_list, text=True, shell=shell_check, capture_output=True)
        if process.returncode != 0:
            log_print(f"ERROR:\t{process.stderr}")
        else:
            log_print(f"PASS:\tSuccessfully processed command: {' '.join(cmd_list)}")
        return process.stdout

def find_fq_files(input_folder):
    """
    Searches for specific files in a given folder and returns their paths.
    
    This function looks for two specific types of files within the specified input folder: files ending with "_1.fq" 
    and "_2.fq". It stores and returns the paths of these files in a list, ensuring that the file with "_1.fq" is at 
    position 0 and the one with "_2.fq" is at position 1 in the returned list.
    
    Parameters:
        input_folder (str): The path to the folder where the files will be searched.
    
    Returns:
        list of str: A list containing the paths of the "_1.fq" and "_2.fq" files. The list elements are ordered 
                     such that the "_1.fq" file path is at index 0 and the "_2.fq" file path is at index 1. If 
                     either file is not found, `None` is placed in the respective position.
    
    Notes:
        - The function specifically looks for files with the exact endings "_1.fq" and "_2.fq".
        - It is assumed that there are no more than one of each file type in the folder.
        - If either "_1.fq" or "_2.fq" files are not present in the folder, the corresponding list element will be `None`.
        - This function is typically used in contexts where paired-end sequence data files are expected, often in 
          bioinformatics workflows.
    """
    fq_files = {'_1.fq': None, '_2.fq': None}

    # Iterate through all files in the given folder
    for file in os.listdir(input_folder):
        file_path = os.path.join(input_folder, file)

        # Check if the file ends with '_1.fq.gz' or '_2.fq.gz' and unzip it
        for end in ['_1.fq', '_2.fq']:
            if file.endswith(end + '.gz'):
                with gzip.open(file_path, 'rb') as f_in:
                    with open(file_path[:-3], 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                fq_files[end] = file_path[:-3]
            elif file.endswith(end):
                fq_files[end] = file_path

    return [fq_files['_1.fq'], fq_files['_2.fq']]


def find_file(filename):
    """
    Searches for a file within the current working directory.

    Parameters:
        filename (str): The name of the file to search for.

    Returns:
        str: The path to the first instance of the file, if found.
        None: If the file is not found.
    """
    current_working_directory = os.getcwd()

    for root, dirs, files in os.walk(current_working_directory):
        if filename in files:
            return os.path.join(root, filename)
    return None

def find_drives():
    """
    Identifies and returns a list of root directories for all available drives on the system.
    
    This function is designed to be cross-platform, supporting both Unix-like and Windows systems. 
    It categorizes and returns the system's root directories differently based on the operating system.
    
    On Unix-like systems (identified by os.name == 'posix'), it returns a list containing only the root ('/'), 
    as these systems mount all drives under the root directory.
    
    On Windows systems (identified by os.name == 'nt'), it generates a list of drive letters (e.g., 'C:\\', 'D:\\', ...) 
    by checking each letter from A to Z to see if it corresponds to an existing drive. 
    
    Returns:
        list of str: A list of strings, each representing the root directory of an available drive on the system. 
        The list will contain just '/' for Unix-like systems or a list of drive letters like ['C:\\', 'D:\\', ...] for Windows systems.
         
        If the operating system is neither Unix-like nor Windows, or if an error occurs, it returns an empty list.
    """
    # For Unix and Unix-like systems, return '/'
    if os.name == 'posix':
        return ['/']
    # For Windows, return a list of drive letters
    elif os.name == 'nt':
        drives = []
        for drive_letter in range(65, 91): # ASCII values for A-Z
            drive = chr(drive_letter) + ':\\'
            if os.path.exists(drive):
                drives.append(drive)
        return drives
    else:
        return []

def trimmomatic_prep(input_fq_list, CPU_THREADS):
    """
    Prepares and runs the Trimmomatic command for paired-end sequence trimming.

    This function constructs the command for running Trimmomatic on a pair of FASTQ files. It sets up the paths for 
    output files (both paired and unpaired) and includes parameters for Trimmomatic such as adapter clipping, 
    sliding window trimming, and minimum length filtering. The paths to the adapter sequences used by Trimmomatic 
    are also determined.

    Parameters:
        input_fq_list (list of str): A list containing the paths to the input FASTQ files. The first element is 
                                     expected to end with "_1.fq" and the second with "_2.fq".

    Returns:
        tuple of str: A tuple containing the paths to the forward and reverse paired output FASTQ files.

    Raises:
        FileNotFoundError: If the Trimmomatic adapters file is not found at the default location.

    Notes:
        - The function uses a global 'CPU_THREADS' variable to set the number of threads for Trimmomatic.
        - The Trimmomatic command includes fixed parameters for adapter clipping, sliding window, and minimum length. 
          These could be modified or made configurable if different settings are required.
        - It's important to have the correct path to the adapter file used by Trimmomatic; this function uses a 
          default path and raises an error if the file is not found there.
    """    
    # Build output file path names based on input_fq_list
    trimmo_f_pair_path = input_fq_list[0].replace('_1.fq','_forward_paired.fq')
    trimmo_f_unpair_path = input_fq_list[0].replace('_1.fq','_forward_unpaired.fq')
    trimmo_r_pair_path = input_fq_list[1].replace('_2.fq','_reverse_paired.fq')
    trimmo_r_unpair_path = input_fq_list[1].replace('_2.fq','_reverse_unpaired.fq')
    
    # Default path to Trimmomatic paths
    default_trimmo_path = "./Trimmomatic-0.39/trimmomatic-0.39.jar"
    default_trimmo_adapters_path = "./Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
    
    if not os.path.exists(default_trimmo_adapters_path):
        # Perform a search for the file TruSeq3-PE.fa
        trimmo_adapters_path = find_file("TruSeq3-PE.fa", "/")
        trimmo_path = find_file("trimmomatic-0.39.jar","/")
        if trimmo_adapters_path is None or trimmo_path is None:
            # try:
            drives = find_drives()
            for drive in drives:
                for root, dirs, _ in os.walk(drive):
                    if 'Trimmomatic-0.39' in dirs:
                        trimmo_path = os.path.join(root, 'Trimmomatic-0.39/trimmomatic-0.39.jar')
                        trimmo_adapters_path =  os.path.join(root, 'Trimmomatic-0.39/adapters/TruSeq3-PE.fa')
            # except:
            #     raise FileNotFoundError("Trimmomatic adapters file not found")
    else:
        trimmo_path = default_trimmo_path
        trimmo_adapters_path = default_trimmo_adapters_path
    
    
    # Make variables for other commands as needed: {trimmo_adapters_path}:{}:{}:{}:{}, SLIDINGWINDOW, MINLEN
    # Joining command list into a single string
    trimmo_cmd = f"java -jar {trimmo_path} PE -phred33 -threads {CPU_THREADS} " \
                f"{input_fq_list[0]} {input_fq_list[1]} " \
                f"{trimmo_f_pair_path} {trimmo_f_unpair_path} " \
                f"{trimmo_r_pair_path} {trimmo_r_unpair_path} " \
                f"ILLUMINACLIP:{trimmo_adapters_path}:2:30:10:11 SLIDINGWINDOW:50:25 MINLEN:125"

    # Executing the command
    _ = run_subprocess_cmd(trimmo_cmd, True)
    
    return trimmo_f_pair_path, trimmo_r_pair_path

def bbduk_map(trimmo_f_pair_path, trimmo_r_pair_path):
    """
    Maps and trims Illumina reads using BBDuk and returns the paths to the mapped files.
    
    This function uses BBDuk to perform quality trimming and adapter removal on paired-end FASTQ files. It constructs 
    the command for BBDuk with specific parameters and paths to the input and output files. The function also specifies 
    the path to the adapter sequences file used by BBDuk.
    
    Parameters:
        trimmo_f_pair_path (str): The path to the forward paired FASTQ file, typically an output from Trimmomatic.
        trimmo_r_pair_path (str): The path to the reverse paired FASTQ file, typically an output from Trimmomatic.
    
    Returns:
        tuple of str: A tuple containing the paths to the forward and reverse mapped output FASTQ files.
    
    Raises:
        FileNotFoundError: If the BBDuk adapters file is not found at the default location.
    
    Notes:
        - The BBDuk command includes parameters for trimming quality and removing adapters, with specific settings 
          for 'ktrim', 'k', 'mink', 'hdist', 'qtrim', and 'trimq'. These settings are fixed in this function but could 
          be made configurable if different settings are required.
        - It's crucial to have the correct path to the adapter file used by BBDuk; this function uses a default path 
          and raises an error if the file is not found there.
    """
    bbduk_f_map_path = trimmo_f_pair_path.replace('_forward_paired.fq','_forward_mapped.fq')
    bbduk_r_map_path = trimmo_r_pair_path.replace('_reverse_paired.fq','_reverse_mapped.fq')
    
    # Default path to BBDuk adapters
    default_adapters_path = "/home/eye/miniconda3/envs/entheome_env/bbtools/lib/resources/adapters.fa"
    
    if not os.path.exists(default_adapters_path):
        # Perform a search for the file TruSeq3-PE.fa
        adapters_path = find_file("adapters.fa", "/")
        if adapters_path is None:
            # try:
            conda_env_path = os.environ['CONDA_PREFIX']
            # Walk through the directories in the Conda environment
            for root, dirs, files in os.walk(conda_env_path):
                if 'bbtools' in dirs:
                    # Return the full path if bbtools is found
                    adapters_path = os.path.join(root, 'bbtools/lib/resources/adapters.fa')
            # except:
            #     raise FileNotFoundError("bbduk adapters file not found")
    else:
        adapters_path = default_adapters_path
        
    # bbduk Map Illumina Reads with run_subprocess_cmd shell_check = True
    bbduk_cmd = ["bbduk.sh",
                 f"in1={trimmo_f_pair_path}", f"in2={trimmo_r_pair_path}",
                 f"out1={bbduk_f_map_path}", f"out2={bbduk_r_map_path}",
                 f"ref={adapters_path}",
                 "ktrim=r", "k=23", "mink=11", "hdist=1", "tpe", "tbo", "qtrim=rl", "trimq=20"]
    _ = run_subprocess_cmd(bbduk_cmd, True)
    return bbduk_f_map_path, bbduk_r_map_path

def clumpify_dedup(bbduk_f_map_path, bbduk_r_map_path):
    """
    Runs Clumpify for de-duplication of mapped FASTQ files.
    
    This function uses Clumpify to de-duplicate paired-end FASTQ files. It creates new file paths for the output 
    de-duplicated files based on the input mapped file paths. The function then constructs and runs the Clumpify 
    command with the specified input and output file paths and the 'dedupe' option.
    
    Parameters:
        bbduk_f_map_path (str): The path to the forward mapped FASTQ file.
        bbduk_r_map_path (str): The path to the reverse mapped FASTQ file.
    
    Returns:
        tuple of str: A tuple containing the paths to the forward and reverse de-duplicated output FASTQ files.
    
    Notes:
        - The 'dedupe' option in Clumpify is used to remove duplicate reads, which is a common step in many bioinformatics 
          workflows to reduce redundancy and improve downstream analysis accuracy.
    """
    clump_f_dedup_path = bbduk_f_map_path.replace('_forward_mapped.fq','_forward_dedup.fq')
    clump_r_dedup_path = bbduk_r_map_path.replace('_reverse_mapped.fq','_reverse_dedup.fq')
    
    clumpify_cmd = ["clumpify.sh",
                    f"in={bbduk_f_map_path}",
                    f"in2={bbduk_r_map_path}",
                    f"out={clump_f_dedup_path}",
                    f"out2={clump_r_dedup_path}",
                    "dedupe"]
    _ = run_subprocess_cmd(clumpify_cmd, True)
    
    return clump_f_dedup_path, clump_r_dedup_path

def parse_bbmerge_output(output):
    """
    Parses the output from the BBMerge command to extract average insert size and standard deviation.

    Parameters:
        output (str): The standard output from the BBMerge command.

    Returns:
        tuple: A tuple containing the average insert size and the standard deviation.
    """
    avg_insert = None
    std_dev = None

    for line in output.splitlines():
        if "Avg Insert:" in line:
            avg_insert = float(line.split()[-1])
        elif "Standard Deviation:" in line:
            std_dev = float(line.split()[-1])

    if avg_insert is None or std_dev is None:
        raise ValueError("Could not find average insert size and/or standard deviation in the output.")

    return avg_insert, std_dev


def masurca_config_gen(input_folder, clump_f_dedup_path, clump_r_dedup_path, CPU_THREADS):
    """
    Generates a configuration file for MaSuRCA and executes genome assembly.
    
    This function performs multiple steps in preparing and executing a de novo genome assembly using MaSuRCA. It 
    first runs BBMerge to determine the average insert size and standard deviation of the input paired-end FASTQ files.
    These metrics are then used to generate a configuration file for MaSuRCA. The function handles the assembly 
    process by running MaSuRCA with the generated configuration file and then executing the assembly script produced 
    by MaSuRCA.
    
    Parameters:
        input_folder (str): The path to the folder where the MaSuRCA configuration file will be created.
        clump_f_dedup_path (str): The path to the forward de-duplicated FASTQ file.
        clump_r_dedup_path (str): The path to the reverse de-duplicated FASTQ file.
    
    Returns:
        str: The path to the scaffolded assembly file generated by MaSuRCA.
    
    Raises:
        Exception: If the assembly process does not complete successfully.
    
    Notes:
        - BBMerge is used to calculate the average insert size and standard deviation, which are crucial parameters for 
          the assembly process.
        - If BBMerge fails to provide these metrics, default values are used.
        - The configuration file for MaSuRCA is dynamically generated based on the insert size metrics and other 
          predefined parameters. These parameters can be adjusted or made configurable for different assembly requirements.
        - The function initiates the MaSuRCA assembly process, which includes generating an assembly script ('assemble.sh') 
          and then executing it.
          - The function uses the global `CPU_THREADS` variable to set the number of threads for MaSurCA.
    
    TODO:
        - Implement error handling to catch and log any issues during the assembly process.
        - Make additional parameters (like USE_LINKING_MATES, LIMIT_JUMP_COVERAGE, etc.) configurable if needed.
    """
    bbmap_out_path = clump_f_dedup_path.replace('_forward_dedup.fq','_bbmap')
    
    # bbmerge to Determine average insert size and standard deviation
    bbmerge_cmd = ["bbmerge.sh",
                   f"in1={clump_f_dedup_path}",
                   f"in2={clump_r_dedup_path}",
                   f"out={bbmap_out_path}",
                   "ihist=insert_size_histogram.txt"]
    bbmerge_output = run_subprocess_cmd(bbmerge_cmd, True)
   
    # Parse standard out text for avg_insert and std_dev
    try:
        avg_insert, std_dev = parse_bbmerge_output(bbmerge_output)
    except:
        avg_insert = 251
        std_dev = 30
        
    # Define the content to be written to the file
    # TODO: Make variables for other commands as needed: USE_LINKING_MATES, LIMIT_JUMP_COVERAGE, CA_PARAMETERS, MEGA_READS_ONE_PASS, LHE_COVERAGE, KMER_COUNT_THRESHOLD, JF_SIZE, DO_HOMOPOLYMER_TRIM
    config_content = f"""DATA
                       PE= pe {avg_insert} {std_dev} {clump_f_dedup_path} {clump_r_dedup_path}
                       END
                       PARAMETERS
                       GRAPH_KMER_SIZE = auto
                       USE_LINKING_MATES = 1
                       LIMIT_JUMP_COVERAGE = 300
                       CA_PARAMETERS = cgwErrorRate=0.15
                       MEGA_READS_ONE_PASS=0
                       LHE_COVERAGE=35
                       KMER_COUNT_THRESHOLD = 1
                       NUM_THREADS = {CPU_THREADS}
                       JF_SIZE = 2500000000
                       DO_HOMOPOLYMER_TRIM = 0
                       END"""
    
    # Specify the file name
    config_path = f"{input_folder}/masurca_config_file.txt"
    
    # Open the file for writing
    with open(config_path, 'w') as file:
        file.write(config_content)

    # Masurca Assembly commands in succession
    masurca_config_cmd = ["masurca", config_path]
    _ = run_subprocess_cmd(masurca_config_cmd, True)
    
    masurca_assemble_cmd = [f"{input_folder}/assemble.sh"]
    _ = run_subprocess_cmd(masurca_assemble_cmd, True)

    # TODO: If completed successfully return scaffolded_assmebly_path, else error out and print path to log file for review (TODO: generate log file)
    scaffolded_assmebly_path = f"{os.path.dirname(input_fq_list[0])}/CA/primary.genome.scf.fasta"
    return scaffolded_assmebly_path

def qc_checks(scaffolded_assmebly_path, CPU_THREADS):
    """
    Performs quality control checks on the scaffolded assembly using QUAST and CompleAsm.
    
    This function runs two quality control tools: QUAST and CompleAsm. QUAST is used for assessing the quality of 
    the assembled sequences, while CompleAsm is run twice with different lineage specifications (Basidiomycota and 
    Agaricales) for completeness analysis of the assembly.
    
    Parameters:
        scaffolded_assmebly_path (str): The path to the scaffolded genome assembly file.
    
    Notes:
        - The output of QUAST is stored in a directory named after the scaffolded assembly file with '_quast' appended.
        - CompleAsm is used to analyze the completeness of the assembly relative to specific lineages, namely 
          Basidiomycota and Agaricales. The output for each lineage is stored in separate directories.
        - The function uses the global `CPU_THREADS` variable to set the number of threads for both QUAST and CompleAsm.
    """
    quast_out = scaffolded_assmebly_path.replace(".fasta","_quast")
    quast_cmd = ["quast",
                 "-o", quast_out,
                 "-t", str(CPU_THREADS),
                 scaffolded_assmebly_path]
    _ = run_subprocess_cmd(quast_cmd, True)
    
    # Default path to Compleasm
    default_compleasm_path = "./Compleasm/compleasm_kit/compleasm.py" 
    
    if not os.path.exists(default_compleasm_path):
        # Perform a search for the file compleasm.py
        compleasm_path = find_file("compleasm.py", "/")
        if compleasm_path is None:
            # try:
            drives = find_drives()
            for drive in drives:
                for root, dirs, _ in os.walk(drive):
                    if 'Compleasm' in dirs:
                        compleasm_path = os.path.join(root, 'compleasm_kit/compleasm.py')
            # except:
            #     raise FileNotFoundError("Compleasm python file not found")
    else:
        compleasm_path = default_compleasm_path
    
    basidio_out = scaffolded_assmebly_path.replace(".fasta","_compleasm_basidiomycota")
    compleasm_basidio_cmd = ["python",compleasm_path, "run", "-a", scaffolded_assmebly_path,
                             "-o", basidio_out,
                             "-l", "basidiomycota",
                             "-t", str(CPU_THREADS)]
    _ = run_subprocess_cmd(compleasm_basidio_cmd, True)
    
    agaricales_out = scaffolded_assmebly_path.replace(".fasta","_compleasm_agaricales")
    compleasm_agaricales_cmd = ["python",compleasm_path, "run", "-a", scaffolded_assmebly_path,
                                "-o", agaricales_out,
                                "-l", "agaricales",
                                "-t", str(CPU_THREADS)]
    _ = run_subprocess_cmd(compleasm_agaricales_cmd, True)

def barcode_extractor(scaffolded_assmebly_path, barcodes_path):
    """
    Extracts barcode sequences from the scaffolded assembly using Exonerate.

    This function uses Exonerate to align barcode sequences to the scaffolded genome assembly. The aim is to 
    identify and extract regions in the assembly that correspond to specific barcode sequences provided in the 
    barcodes file.

    Parameters:
        scaffolded_assmebly_path (str): The path to the scaffolded genome assembly file.
        barcodes_path (str): The path to the file containing barcode sequences.

    Notes:
        - The output of Exonerate, which includes alignments and potentially identified barcodes, is stored in a 
          file named after the scaffolded assembly file with the basename of the barcodes file appended.
        - The Exonerate command is configured to use the 'est2genome' model, which is suitable for aligning 
          expressed sequence tags (ESTs) to genomic sequences. This may be appropriate for barcode sequences 
          depending on their nature and origin.
    """
    exonerate_output_path = scaffolded_assmebly_path.replace(".fasta",f"_{os.path.basename(barcodes_path)}")

    exonerate_cmd = ["exonerate", "--model", "est2genome",
                     "--query", barcodes_path,
                     "--target", scaffolded_assmebly_path,
                     "--showvulgar", "true",
                     "--showalignment", "true",
                     ">", exonerate_output_path]
    
    _ = run_subprocess_cmd(exonerate_cmd, True)

def illumina_only_main(INPUT_FOLDER, PERCENT_RESOURCES):
    # CPU Threads count setup
    num_cpus = multiprocessing.cpu_count()
    CPU_THREADS = int(math.floor(num_cpus * PERCENT_RESOURCES))
    
    # Generate list of raw input fq files
    input_fq_list = find_fq_files(INPUT_FOLDER)
    initialize_logging_environment(INPUT_FOLDER)
    
    # Run Trimmomatic on the raw input files
    trimmo_f_pair_path, trimmo_r_pair_path = trimmomatic_prep(input_fq_list, CPU_THREADS)
    
    # Run bbduk on the trimmed files
    bbduk_f_map_path, bbduk_r_map_path = bbduk_map(trimmo_f_pair_path, trimmo_r_pair_path)
    
    # Run clumpify on the mapped files
    clump_f_dedup_path, clump_r_dedup_path = clumpify_dedup(bbduk_f_map_path, bbduk_r_map_path)
    
    # Run Masurca to generate a Scaffolded Assembly
    scaffolded_assmebly_path = masurca_config_gen(INPUT_FOLDER, clump_f_dedup_path, clump_r_dedup_path, CPU_THREADS)
    
    # Run QC Checks on Scaffolded Assembly
    qc_checks(scaffolded_assmebly_path, CPU_THREADS)
    
    return scaffolded_assmebly_path

if __name__ == "__main__":
    # Set Default Values
    default_input_folder = "/mnt/d/ENTHEOME/Ps_zapatecorum/RECREATION"
    default_barcodes_path = "/mnt/d/ENTHEOME/Assembled Genes/Fungal_Barcodes.txt"
    default_percent_resources = 0.2
    
    # Create the parser
    parser = argparse.ArgumentParser(description="Run bioinformatics pipeline")
    
    # Add arguments
    parser.add_argument("-i", "--input_folder", type=str, default=default_input_folder,
                        help="Path to the input folder containing FASTQ files")
    parser.add_argument("-r", "--resources", type=float, default=default_percent_resources,
                        help="Percentage of resources to use (0.01-1.00; default: 0.2)")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Set the input folder and barcodes path from the arguments
    INPUT_FOLDER = args.input_folder
    PERCENT_RESOURCES = args.resources
    
    # Main function
    scaffolded_assmebly_path = illumina_only_main(INPUT_FOLDER, PERCENT_RESOURCES)
