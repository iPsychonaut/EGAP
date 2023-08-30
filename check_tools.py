# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 11:39:38 2023

@author: ian.michael.bollinger@gmail.com with the help of ChatGPT 4.0
"""
import platform, os, subprocess, hashlib, pkg_resources, fnmatch, shutil
from log_print import log_print, generate_log_file

# Function to move a file up one directory
def move_file_up(file_path, log_file, move_bool=True):
    """
    Moves the given file up one directory.
    
    Args:
        file_path (str): Path to a file that needs to be moved up one directory.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.
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
        log_print(f"Attempting to move {file_path} one directory up...", log_file)
        shutil.move(file_path, new_path)
        log_print(f"PASS:\tSuccessfully moved {file_path} to {new_path}", log_file)
    if move_bool == False:
        # Just return the new_path
        log_print(f"NOTE:\tSkipping move, returning {new_path}", log_file)
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
    mounted_drive = mounted_drive = f"{'/'.join(BASE_FOLDER.split('/')[:3])}/"
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
    print(f'Operating System: {platform_system}')
    print(f'Environment Directory: {environment_dir}')

    return environment_dir

# Function to check that all Python Libraries are installed, and if not installs them using conda
def libraries_check(libraries, log_file):
    """
    Check that all Python Libraries in the provided list are installed, and if not installs them using conda.

    Args:
        libraries (list): List of Python libraries to make sure are installed.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.
    """
    # Print the start of the library check
    log_print("Checking Python Library Prerequisites...", log_file)
    
    # Retrieve all installed Python libraries using pkg_resources, a runtime facilities for Python libraries
    installed = {pkg.key for pkg in pkg_resources.working_set}
    
    # Check which libraries from the list are not installed yet, by comparing with the installed libraries
    missing = [lib for lib in libraries if lib not in installed]

    # If there are missing libraries
    if missing:
        log_print(f"WARN:\tMissing libraries: {missing}. Installing...", log_file)
        
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
            log_print(f'CMD:\t{" ".join(missing_cmd)}',log_file)
            subprocess.check_call(missing_cmd, stdout=subprocess.DEVNULL)
        except:
            # Print an error if something goes wrong during the installation
            log_print(f"ERROR:\t Unable to Instal {missing} library", log_file)
    else:
        # If all libraries are installed, print a confirmation message
        log_print("PASS:\tAll libraries are installed", log_file)

# Function to check if third-party programs' specific JAR file is installed in the current working environment
def check_for_jars(program_dict, log_file):
    """
    Check if third-party programs' specific JAR files are installed within the immediate subdirectories of the base directory.
    
    Args:
        program_dict (dict): Dictionary containing program names as keys and tuples containing their git repos and jar file patterns as values.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.
    
    Returns:
        jar_paths_dict (dict): Dictionary containing program names as keys and paths to the respective JAR files as values.
    """
    # TODO: INSTALL INTO CURRENT WORKING DIRECTORY (SHOULD BE )
    base_dir = os.path.expanduser("~")
    jar_paths_dict = {}
    
    # List all immediate subdirectories of the base directory
    subdirs = [os.path.join(base_dir, d) for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    
    for program_name, (git_repo, jar_pattern) in program_dict.items():
    
        # Searching for the JAR file in immediate subdirectories
        file_path = None
        for subdir in subdirs:
            for root, _, files in os.walk(subdir):
                for filename in fnmatch.filter(files, jar_pattern):
                    file_path = os.path.join(root, filename)
                    break
                if file_path:  # break the loop as soon as we found the file
                    break
            if file_path:  # break the outer loop as well if we found the file
                break
    
        if file_path:
            message = f"PASS:\t{program_name} JAR file is installed at {file_path}"
            jar_paths_dict[program_name] = file_path
        else:
            log_print(f"WARN:\t{program_name} JAR file is not installed. Installing from {git_repo}...", log_file)
            
            # Install the program from git
            install_path = os.path.join(base_dir, program_name)
            jar_cmd = ['git', 'clone', git_repo, install_path]
            log_print(f"CMD:\t{' '.join(jar_cmd)}", log_file)
            subprocess.run(jar_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            # Now, search for the JAR file within the installed directory
            for root, _, files in os.walk(install_path):
                for filename in fnmatch.filter(files, jar_pattern):
                    file_path = os.path.join(root, filename)
                    break
                if file_path:  # break the loop as soon as we found the file
                    break
            
            if file_path:
                message = f"PASS:\t{program_name} JAR file is now installed at {file_path}"
                jar_paths_dict[program_name] = file_path
            else:
                message = f"ERROR:\tFailed to find {program_name} JAR file after installation"
                jar_paths_dict[program_name] = None
        
        log_print(message, log_file)

    return jar_paths_dict, base_dir 

def check_prereqs_installed(prerequisites, log_file):
    """
    Check if third-party programs are installed in the current working environment.

    Args:
        prerequisites (list): List of third-party programs to check for installation.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.
    """
    log_print("Checking Third-Party Prerequisites...", log_file)

    # Convert list into a single string for command
    prereqs_string = ' '.join(prerequisites)

    # Using f-string to correctly format the string with the prerequisites
    prereq_cmd = ['bash', '-c', 
                f'''for cmd in {prereqs_string}; \
                    do \
                        command -v $cmd >/dev/null 2>&1 && {{ echo "$cmd is installed"; }} || {{ echo >&2 "$cmd is not installed"; exit 1; }}; \
                    done;''']
    log_print(f"Checking for: {', '.join(prerequisites)}", log_file)
    prereq_result = subprocess.run(prereq_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if prereq_result.returncode != 0:
        log_print(f"ERROR:\tThere was a problem checking the installation of Third-Party Prerequisites. Return code: {prereq_result.returncode}", log_file)
        
        # Print detailed output
        output_messages = prereq_result.stdout.decode().splitlines()
        for message in output_messages:
            log_print(f"DETAIL:\t{message}", log_file)
        
        # If there's any error messages, log them as well
        error_messages = prereq_result.stderr.decode().splitlines()
        for message in error_messages:
            log_print(f"ERROR DETAIL:\t{message}", log_file)
            
    else:
        log_print("PASS:\tAll Third-Party Prerequisites successfully found", log_file)

# Debuging Main Space & Example
if __name__ == "__main__":
    print("check_tools.py DEBUG")

    # Generate Main Logfile
    debug_log = 'check_tools_log.tsv'
    log_file = generate_log_file(debug_log, use_numerical_suffix=False)
    
    # Ensure mamba is installed
    mamba_cmd = ['conda', 'install', '-y', 'conda-forge', 'mamba==1.5.0']
    print(f"UNLOGGED CMD:\t{' '.join(mamba_cmd)}")
    mamba_result = subprocess.check_call(mamba_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if mamba_result == 0:
        print(f'UNLOGGED PASS:\tSuccessfully installed: mamba')
    else:
        print(f'UNLOGGED ERROR:\tUnable to install: mamba')
    
    # Check Python Libraries
    libraries = ['gdown==4.7.1','busco==4.1.2','openjdk==20.0.0', 'nanoq==0.10.0', 'pandas==2.0.3',
                 'biopython==1.81', 'tqdm==4.38.0', 'psutil==5.9.5', 'termcolor==2.3.0',
                 'beautifulsoup4==4.12.2', 'fastqc==0.11.8', 'quast==5.2.0', 'nanostat==1.6.0',
                 'flye==2.9.2', 'bbtools==37.62', 'metaeuk==6.a5d39d9', 'blast==2.14.1',
                 'bwa==0.7.17', 'minimap2==2.26', 'pysam==0.21.0', 'samtools==1.17',
                 'arcs==1.2.5', 'tigmint==1.2.10', 'abyss==2.3.7', 'racon==1.5.0', 'spades==3.15.3']
    libraries_check(libraries, log_file)

    # Check for specific Third-Party JAR Files with adjusted program_dict
    program_dict = {"trimmomatic": ("https://github.com/usadellab/Trimmomatic", "trimmomatic-*.jar"),
                    "pilon": ("https://github.com/broadinstitute/pilon", "pilon*.jar")}
    jar_paths_dict, _ = check_for_jars(program_dict, log_file)

    # Check Pipeline Third-Party Prerequisites
    prerequisites = ["java", "fastqc", "quast.py", "busco", "samtools", "bwa",
                     "makeblastdb", "blastn", "racon", "nanoq", "NanoStat",
                     "spades.py", "tblastn", "flye", "minimap2"]
    check_prereqs_installed(prerequisites, log_file)
    ### TODO: ADD KELSEY'S PREREQUISITE: abyss-sealer
