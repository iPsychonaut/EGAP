# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 11:39:38 2023

@author: ian.michael.bollinger@gmail.com with the help of ChatGPT 4.0
"""
# Base Python Imports
import math, platform, os, subprocess, shutil, pkg_resources, multiprocessing

# Required Python Imports
import psutil, hashlib

from datetime import datetime

# Global output_area variable
CPU_THREADS = 1
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

def find_file(filename):
    """
    Searches for a file within the current working directory.

    Parameters:
        filename (str): The name of the file to search for.

    Returns:
        str: The path to the first instance of the file, if found.
        None: If the file is not found.
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

# Function to convert the user input PERCENT_RESOURCES into usuable cpu_threads and ram_gb values
def get_resource_values(PERCENT_RESOURCES):
    """
    Converts user input PERCENT_RESOURCES into usuable cpu_threads and ram_gb values.

    Arg:
        PERCENT_RESOURCES (str): A string that is between 1-100. 

    Returns
    -------
    cpu_threads (str): A count of the CPUs available for processing.
    ram_gb (str): A count of the RAM (in GB) available for processing.

    """
    # Get the number of CPUs available on the system
    num_cpus = multiprocessing.cpu_count()
    
    # Get the amount of RAM (GB) currently available
    mem_info = psutil.virtual_memory()
    
    # Calculate the number of threads as 80% of available CPUs & RAM
    cpu_threads = int(math.floor(num_cpus * PERCENT_RESOURCES))
    ram_gb = int(mem_info.total / (1024.0 ** 3) * PERCENT_RESOURCES)
    
    return cpu_threads, ram_gb 

# Function to move a file up one directory
def move_file_up(file_path, move_bool=True):
    """
    Moves the given file up one directory.
    
    Args:
        file_path (str): Path to a file that needs to be moved up one directory.
        move_bool (bool): Boolean; when 'False' will just return the path of the file one directory above.

    Returns:
        new_path (str): Path to the new file location.        
    """
    # Get the directory containing the file and its parent directory
    current_directory = os.path.dirname(file_path)
    parent_directory = os.path.dirname(current_directory)
    
    # Construct the new path for the file in the parent directory
    new_path = f"{parent_directory}/{os.path.basename(file_path)}"
    
    if move_bool == True:
        # Move the file
        print(f"UNLOGGED:\tAttempting to move {file_path} one directory up...")
        shutil.move(file_path, new_path)
        print(f"UNLOGGED PASS:\tSuccessfully moved {file_path} to {new_path}")
    if move_bool == False:
        # Just return the new_path
        pass
    
    return new_path

# Function to generate an md5 checsum for a given file
def get_md5(file_path):
    """
    Calculate the MD5 checksum of a given file.

    Args:
        file_path (str): Path to the file.

    Returns:
        hash_md5.hexdigest (str): MD5 checksum of the file.
    """
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

# Function to search for a directory containing a specific file
def search_directory_for_file(start_dir, target_file):
    """
    Search recursively for a directory containing the target file.

    Args:
        start_dir (str): Directory to start the search.
        target_file (str): Name of the file to search for.

    Returns:
        root (str): Path to the directory containing the file, or None if not found.
    """
    for root, _, files in os.walk(start_dir):
        if target_file in files:
            return root
    return None

# # Function to find a file in a directory of its sub directories
# def find_file(root, filename):
#     """
#     Find a file in a directory or its subdirectories.
    
#     Args:
#         root (str): Description.
#         filename (str): Description.
    
#     Returns:
#         found_path (str): Confirmed path found.
#     """
#     for dirpath, dirnames, filenames in os.walk(root):
#         if filename in filenames:
#             found_path = os.path.join(dirpath, filename)
#             return found_path
#     return None

# Function to determine the current operating system and set the environment directory and command prefix
def get_env_dir(BASE_FOLDER):
    """
    Determine the operating system and set environment directory and command prefix.

    Returns:
        str: The environment directory for the detected operating system.

    Raises:
        Exception: If the OS is neither Windows nor Linux.
    """
    # Determine the operating system in use
    os_name = os.name
    platform_system = platform.system()

    # Mapping OS to environment directory and command prefix
    mounted_drive = f"{'/'.join(BASE_FOLDER.split('/')[:3]).replace('/mnt/','')}"
    os_mapping = {('nt', 'Windows'): (f"{mounted_drive}:", "wsl "),
                  ('posix', 'Linux'): (f"/mnt/{mounted_drive}", "")}

    # Check the current OS against the mapping and set the environment directory and command prefix
    for os_keys, (env_dir, cmd_prefix) in os_mapping.items():
        if os_name in os_keys or platform_system in os_keys:
            global environment_dir, environment_cmd_prefix
            environment_dir, environment_cmd_prefix = env_dir, cmd_prefix
            break
    else:
        # If the operating system is neither Windows nor Linux, raise an Exception
        raise Exception("ERROR: OS NOT TESTED WITH THIS CODE")

    # Print the detected operating system and determined environment directory
    print(f'UNLOGGED:\tOperating System: {platform_system}')
    print(f'UNLOGGED:\tEnvironment Directory: {environment_dir}')

    return environment_dir

# Function to check that all Python Libraries are installed, and if not installs them using conda
def libraries_check(libraries):
    """
    Check that all Python Libraries in the provided list are installed, and if not installs them using conda.

    Args:
        libraries (list): List of Python libraries to make sure are installed.
    """
    # Print the start of the library check
    print("UNLOGGED:\tChecking Python Library Prerequisites...")
    
    # Retrieve all installed Python libraries using pkg_resources, a runtime facilities for Python libraries
    installed = {pkg.key for pkg in pkg_resources.working_set}
    
    # Check which libraries from the list are not installed yet, by comparing with the installed libraries
    missing = [lib for lib in libraries if lib not in installed]

    # If there are missing libraries
    if missing:
        print(f"UNLOGGED WARN:\tMissing libraries: {missing}. Installing...")
        
        try:
            # Run a subprocess calling conda install for each missing library
            # subprocess.check_call will raise an error if the installation fails
            # stdout=subprocess.DEVNULL will prevent the output from appearing on your console
            missing_cmd = ['mamba', 'install', '-y',
                           '-c', 'bioconda',
                           '-c', 'agbiome',
                           '-c', 'prkrekel',
                           '-c', 'conda-forge',
                           *missing]
            print(f'UNLOGGED CMD:\t{" ".join(missing_cmd)}')
            subprocess.check_call(missing_cmd, stdout=subprocess.DEVNULL)
        except:
            # Print an error if something goes wrong during the installation
            print(f"UNLOGGED ERROR:\t Unable to Instal {missing} library")
    else:
        # If all libraries are installed, print a confirmation message
        print("UNLOGGED PASS:\tAll libraries are installed")

# Debuging Main Space & Example
if __name__ == "__main__":
    print("check_tools.py DEBUG")
