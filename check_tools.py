# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 11:39:38 2023

@author: ian.michael.bollinger@gmail.com with the help of ChatGPT 4.0
"""
# Base Python Imports
import math, platform, os, subprocess, shutil, pkg_resources, fnmatch, multiprocessing
from concurrent.futures import ThreadPoolExecutor

# Required Python Imports
import psutil, hashlib, requests, zipfile

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

# Function to find a file in a directory of its sub directories
def find_file(root, filename):
    """
    Find a file in a directory or its subdirectories.
    
    Args:
        root (str): Description.
        filename (str): Description.
    
    Returns:
        found_path (str): Confirmed path found.
    """
    for dirpath, dirnames, filenames in os.walk(root):
        if filename in filenames:
            found_path = os.path.join(dirpath, filename)
            return found_path
    return None

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
