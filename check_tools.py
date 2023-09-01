# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 11:39:38 2023

@author: ian.michael.bollinger@gmail.com with the help of ChatGPT 4.0
"""
import math, platform, os, subprocess, hashlib, shutil, pkg_resources, fnmatch, requests, zipfile, psutil, multiprocessing
from concurrent.futures import ThreadPoolExecutor

# Function to search for a folder named 'target_folder' within a specified maximum depth in the directory tree
def find_folder(target_folder, max_depth=5):
    """
    Searches for a folder named 'target_folder' within a specified maximum depth in the directory tree.

    Args:
        target_folder (str): Name of the folder to search for.
        max_depth (int, optional): Maximum depth of folders to search in; default is 5.

    Returns:
        result (str): Path to the folder of interest.
    """
    
    # Nested function to search a single drive
    def search_drive(drive):
        depth = 0  # Initialize depth counter
        # os.walk generates the file names in a directory tree by walking either top-down or bottom-up through the directory tree.
        for root, dirs, files in os.walk(drive):
            if depth > max_depth:
                break  # Break if maximum depth reached
            if target_folder in dirs:
                return os.path.join(root, target_folder)  # Return full path if target_folder is found
            depth += 1  # Increment depth for each new level in the directory tree

    # List all available drives in Windows
    all_drives = [f"{d}:\\" for d in "ABCDEFGHIJKLMNOPQRSTUVWXYZ" if os.path.exists(f"{d}:\\")]
    # List all available drives in Linux (WSL)
    all_drives += [f"/mnt/{d.lower()}/" for d in "ABCDEFGHIJKLMNOPQRSTUVWXYZ" if os.path.exists(f"/mnt/{d.lower()}/")]

    # Use multi-threading to search all drives in parallel
    with ThreadPoolExecutor() as executor:
        for result in executor.map(search_drive, all_drives):
            if result is not None:
                return result  # Return the first found result

    return None  # Return None if target_folder is not found

# Function to check for JAR files, downloads them if they are not found, and then moves them into the correct folder
def check_for_jars(program_dict, search_dir):
    """
    Checks for JAR files, downloads them if they are not found, and then moves them into the correct folder.

    Args:
        program_dict (dict): Dictionary of the JAR's to locate.
        search_dir (str): Name of the program install folder to store JARs in.

    Returns:
        jar_paths_dict (dict): 
        install_dir (str): Path to the installed prgram directory.
    """
    
    print(f"UNLOGGED:\tSearching for {search_dir} install folder...")
    # Find the base install directory
    install_dir = find_folder(search_dir)
    
    if install_dir is None:
        print(f'UNLOGGED ERROR: Unable to find base {search_dir} install folder')
        return None, None  # Return None values if install directory not found

    jar_paths_dict = {}
    # List subdirectories under the install directory
    subdirs = [os.path.join(install_dir, d) for d in os.listdir(install_dir) if os.path.isdir(os.path.join(install_dir, d))]

    # Loop through each program to check for JAR files
    for program_name, (git_repo, jar_pattern, zip_url) in program_dict.items():
        file_path = None
        # Search for JAR files in the subdirectories
        for subdir in subdirs:
            for root, _, files in os.walk(subdir):
                for filename in fnmatch.filter(files, jar_pattern):
                    file_path = os.path.join(root, filename)
                    break  # Break if JAR file is found
                if file_path:
                    break  # Break if JAR file is found
            if file_path:
                break  # Break if JAR file is found
                
        # If JAR file is found, update the message and dictionary
        if file_path:
            message = f"UNLOGGED PASS:\t{program_name} JAR file is installed at {file_path}"
        else:  # If JAR file is not found, download and install it
            install_path = os.path.join(install_dir, program_name)
            os.makedirs(install_path, exist_ok=True)

            # Clone the git repository if it exists
            if git_repo:
                git_cmd = ['git', 'clone', git_repo, install_path]
                print(f"UNLOGGED CMD:\t{' '.join(git_cmd)}")
                subprocess.run(git_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            # Download the JAR file
            download_path = os.path.join(install_dir, f"{program_name}.jar")
            print(f"UNLOGGED CMD:\tDownloading {zip_url} to {download_path}")
            response = requests.get(zip_url)
            with open(download_path, 'wb') as f:
                f.write(response.content)
            
            # Move the JAR file to the correct folder
            shutil.move(download_path, os.path.join(install_path, f"{program_name}.jar"))
            
            # Update the message and dictionary
            message = f"UNLOGGED PASS:\t{program_name} JAR file is now installed at {os.path.join(install_path, f'{program_name}.jar')}"
        
        jar_paths_dict[program_name] = file_path or os.path.join(install_path, f"{program_name}.jar")
        print(message)
        
    return jar_paths_dict, install_dir  # Return the dictionary of JAR paths and the install directory

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

def check_prereqs_installed(prerequisites):
    """
    Check if third-party programs are installed in the current working environment.

    Args:
        prerequisites (list): List of third-party programs to check for installation.
    """
    log_print("Checking Third-Party Prerequisites...")

    # Convert list into a single string for command
    prereqs_string = ' '.join(prerequisites)

    # Using f-string to correctly format the string with the prerequisites
    prereq_cmd = ['bash', '-c', 
                f'''for cmd in {prereqs_string}; \
                    do \
                        command -v $cmd >/dev/null 2>&1 && {{ echo "$cmd is installed"; }} || {{ echo >&2 "$cmd is not installed"; exit 1; }}; \
                    done;''']
    print(f"UNLOGGED:\tChecking for: {', '.join(prerequisites)}")
    prereq_result = subprocess.run(prereq_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if prereq_result.returncode != 0:
        print(f"UNLOGGED ERROR:\tThere was a problem checking the installation of Third-Party Prerequisites. Return code: {prereq_result.returncode}")
        
        # Print detailed output
        output_messages = prereq_result.stdout.decode().splitlines()
        for message in output_messages:
            print(f"UNLOGGED DETAIL:\t{message}")
        
        # If there's any error messages, log them as well
        error_messages = prereq_result.stderr.decode().splitlines()
        for message in error_messages:
            print(f"UNLOGGED ERROR DETAIL:\t{message}")
            
    else:
        print("UNLOGGED PASS:\tAll Third-Party Prerequisites successfully found")

# Debuging Main Space & Example
if __name__ == "__main__":
    print("check_tools.py DEBUG")
    
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
    libraries_check(libraries)

    # Check for specific Third-Party JAR Files with adjusted program_dict
    program_dict = {"trimmomatic": ("https://github.com/usadellab/Trimmomatic",
                                "trimmomatic-*.jar",
                                "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip"),
                    "pilon": ("https://github.com/broadinstitute/pilon",
                              "pilon-*.jar",
                              "https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar")}
    jar_paths_dict, _ = check_for_jars(program_dict, 'EGAP')

    # Check Pipeline Third-Party Prerequisites
    prerequisites = ["java", "fastqc", "quast.py", "busco", "samtools", "bwa",
                     "makeblastdb", "blastn", "racon", "nanoq", "NanoStat",
                     "spades.py", "tblastn", "flye", "minimap2"]
    check_prereqs_installed(prerequisites)
    ### TODO: ADD KELSEY'S PREREQUISITE: abyss-sealer
