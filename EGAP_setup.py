# -*- coding: utf-8 -*-
"""
Created on Sat Sep  2 11:46:49 2023

@author: theda
"""
# Base Python Imports
import os, re, tarfile, zipfile, subprocess
from concurrent.futures import ThreadPoolExecutor

if __name__ == "__main__":
    # Make sure mamba is installed
    module_name = 'mamba==1.5.0'
    print(f"UNLOGGED:\tAttempting to install: {module_name}...")
    install_cmd = ['conda', 'install', '-y', '-c', 'conda-forge', module_name]
    print(f"UNLOGGED CMD:\t{' '.join(install_cmd)}")
    intsall_result = subprocess.check_call(install_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if intsall_result == 0:
        print(f'UNLOGGED PASS:\tSuccessfully installed: {module_name}')
    else:
        print(f'UNLOGGED ERROR:\tUnable to install: {module_name}')
    
    # Make sure requests and gdown are installed
    module_list = ['requests==2.31.0', 'gdown==4.7.1' ]
    print(f"UNLOGGED:\tAttempting to install: {module_list}...")
    install_cmd = ['mamba', 'install', '-y', '-c', 'conda-forge', *module_list]
    print(f"UNLOGGED CMD:\t{' '.join(install_cmd)}")
    intsall_result = subprocess.check_call(install_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if intsall_result == 0:
        print(f'UNLOGGED PASS:\tSuccessfully installed: {module_list}')
    else:
        print(f'UNLOGGED ERROR:\tUnable to install: {module_list}')

# Required Python Imports
import requests, gdown

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
    print(f"UNLOGGED:\tSearching for the target_folder: {target_folder}...")
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
                print(f"UNLOGGED PASS:\tFound: {target_folder}")
                return result  # Return the first found result
    print(f"UNLOGGED ERROR:\tFailed to find the target_folder: {target_folder}")
    return None  # Return None if target_folder is not found

# Function to download and setup files based on the URL provided
def download_and_setup(install_dir, url):
    """
    Download and setup files based on the URL provided.
    
    Parameters:
    install_dir (str): Installation directory where the files will be saved.
    url (str): URL of the file to download.

    Returns:
    jar_file_path (str): The path to the JAR file.
    """   
    print(f'UNLOGGED:\tAttempting to ensure {install_dir} is installed...')
 
    # Initialize jar_file_path
    jar_file_path = None 

    # Determine the file type based on the URL
    if url.endswith('.jar'):
        file_type = 'jar'
    elif url.endswith('.zip'):
        file_type = 'zip'
    else:
        print("Unsupported file type. Only .jar and .zip are supported.")
        return jar_file_path

    # Check if the install directory exists
    if not os.path.exists(install_dir):
        os.makedirs(install_dir, exist_ok=True)
    else:
        # Check if JAR file already exists
        expected_jar_name = url.split('/')[-1].replace('.zip', '.jar')
        existing_jar_path = os.path.join(install_dir, expected_jar_name)
        if os.path.exists(existing_jar_path):
            print(f"UNLOGGED PASS:\tJAR file already exists at {existing_jar_path}. Skipping download.")
            return existing_jar_path

    # Create the install directory if it doesn't exist
    os.makedirs(install_dir, exist_ok=True)

    # Download the file
    response = requests.get(url)
    
    # Check for HTTP errors
    response.raise_for_status()

    # Define the file name from the URL
    file_name = url.split('/')[-1]

    # Define the full file path
    file_path = os.path.join(install_dir, file_name)
    
    # Save the downloaded file to the install directory
    with open(file_path, 'wb') as f:
        f.write(response.content)

    # If it's a zip file, extract it
    if file_type == 'zip':
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(install_dir)

        # Remove the ZIP file after extraction
        os.remove(file_path)

        # Assuming the JAR file has the same base name as the ZIP file
        jar_file_name = file_name.replace('.zip', '.jar')
        jar_file_path = os.path.join(install_dir, os.path.splitext(jar_file_name)[0], jar_file_name)

    # If it's a jar file, it's ready to use and you can move or rename it as needed
    elif file_type == 'jar':
        jar_file_path = file_path  # The JAR file path is the same as the downloaded file path

    # Return the path to the JAR file
    print(f"UNLOGGED PASS:\tPath to JAR file: {jar_file_path}")
    return jar_file_path

# Function to check if third-party programs are installed in the current working environment
def check_prereqs_installed(prerequisites):
    """
    Check if third-party programs are installed in the current working environment.

    Args:
        prerequisites (list): List of third-party programs to check for installation.
    """
    print("UNLOGGED:\tChecking Third-Party Prerequisites...")
    
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
            print(f"UNLOGGED:\t{message}")
        
        # If there's any error messages, log them as well
        error_messages = prereq_result.stderr.decode().splitlines()
        for message in error_messages:
            print(f"UNLOGGED ERROR DETAIL:\t{message}")
            
    else:
        print("UNLOGGED PASS:\tAll Third-Party Prerequisites successfully found")

# Function to install a specific python module
def install_module(module_name):
    """
    Installs a provided python module with mamba, accepts versions as well.

    Args:
    module_name (str): The name of the module that can be installed with mamba, versions are acceptable.
    """
    # Make sure module is installed
    install_cmd = ['mamba', 'install', '-y', '-c', 'bioconda', '-c', 'agbiome', '-c', 'prkrekel', '-c', 'conda-forge', module_name]
    print(f"UNLOGGED CMD:\t{' '.join(install_cmd)}")
    intsall_result = subprocess.check_call(install_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if intsall_result == 0:
        print(f'UNLOGGED PASS:\tSuccessfully installed: {module_name}')
    else:
        print(f'UNLOGGED ERROR:\tUnable to install: {module_name}')

# Function to download BUSCO databases if not currently in the base_install_dir
def download_busco_dbs(base_install_dir):
    """
    Download BUSCO databases if not currently in the base_install_dir.

    Args:
        base_install_dir (str). Path to the directory the module is installed in.

    Returns:
        busco_db_dict (dict): Dictionary of each Kingdom's ODB database paths.
    """
    print(f'UNLOGGED:\tAttempting to download BUSCO databases...')
    
    # Generate BUSCO database information
    busco_db_website = 'https://busco-data.ezlab.org/v5/data/lineages/'

    # Fetch the webpage content
    try:
        response = requests.get(busco_db_website)
        response.raise_for_status()
        webpage_content = response.text
    except requests.RequestException as e:
        print(f"UNLOGGED ERROR:\tFailed to fetch the list of available databases: {e}")
        return None
    
    # Use regular expressions to extract all hyperlinks (simplified)
    hyperlinks = re.findall('href=[\'"]?([^\'" >]+)', webpage_content)
    
    # Filter only the hyperlinks that are tar.gz (assuming all databases are tar.gz files)
    available_dbs = [link.split('/')[-1].replace('.tar.gz', '') for link in hyperlinks if link.endswith('.tar.gz')]
    
    busco_db_dict = {'Archaea':  ['archaea_odb10',
                                  'euryarchaeota_odb10'],
                     'Bacteria': ['actinobacteria_phylum_odb10',
                                  'proteobacteria_odb10',],
                     'Fauna':    ['vertebrata_odb10',
                                  'arthropoda_odb10',],
                     'Flora':    ['eudicots_odb10',
                                  'liliopsida_odb10'],
                     'Funga':    ['basidiomycota_odb10',
                                  'agaricales_odb10'],
                     'Protista': ['alveolata_odb10',
                                  'euglenozoa_odb10']}
    
    # Check if the folder "/EGAP_Databases/BUSCO_Databases" exists in the base_install_dir, if not create it
    busco_db_path = os.path.join(base_install_dir, "EGAP_Databases", "BUSCO_Databases")
    if not os.path.exists(busco_db_path):
        os.makedirs(busco_db_path)
    
        # Download and extract BUSCO databases if not already present
        for category, databases in busco_db_dict.items():
            for i, db in enumerate(databases):
                # Using regular expression to find the database with a date suffix
                pattern = re.compile(f"{re.escape(db)}\.\d{{4}}-\d{{2}}-\d{{2}}")
                matching_dbs = list(filter(pattern.match, available_dbs))
                
                if not matching_dbs:
                    print(f"UNLOGGED ERROR:\tDatabase {db} is not available for download.")
                    continue
               
                # Choose the first matching database (or you could choose the latest one)
                chosen_db = matching_dbs[0]
                
                # Check if the database folder already exists
                db_path = os.path.join(busco_db_path, chosen_db)  # Use 'chosen_db' instead of 'db'
                if not os.path.exists(db_path):
                    download_url = f"{busco_db_website}{chosen_db}.tar.gz"  # Use 'chosen_db' instead of 'db_with_date'
                    download_path = os.path.join(busco_db_path, f"{chosen_db}.tar.gz")  # Use 'chosen_db' instead of 'db_with_date'
                              
                    try:
                        # Download the database
                        response = requests.get(download_url)
                        response.raise_for_status()  # Raise an exception for HTTP errors
                        
                        # Save the downloaded content to a .tar.gz file
                        with open(download_path, 'wb') as f:
                            f.write(response.content)
                        
                        # Extract the .tar.gz file
                        with tarfile.open(download_path, 'r:gz') as tar:
                            tar.extractall(path=busco_db_path)
                        
                        # Remove the .tar.gz file after extraction
                        os.remove(download_path)
                        
                    except requests.RequestException as e:
                        print(f"UNLOGGED ERROR:\tFailed to download {db} due to network issue: {e}")
                        continue
                    except (tarfile.TarError, IOError):
                        print(f"UNLOGGED ERROR:\tFailed to extract {db}")
                        continue
        
                # Update the dictionary entry with the full path
                busco_db_dict[category][i] = db_path
                print(f"UNLOGGED PASS:\tUpdated dictionary with {db_path}")
    else:
        # Get the list of all folders (not sub-folders) in busco_db_path
        db_path_folder_list = [d for d in os.listdir(busco_db_path) if os.path.isdir(os.path.join(busco_db_path, d))]
        
        # Replace each entry with the appropriate path from db_path_folder_list
        for category, databases in busco_db_dict.items():
            for i, db in enumerate(databases):
                # Using regex to match database names in the folder list
                pattern = re.compile(f"{re.escape(db)}\.\d{{4}}-\d{{2}}-\d{{2}}")
                matching_dbs = list(filter(pattern.match, db_path_folder_list))
                
                if matching_dbs:
                    # If a matching database folder exists, update the dictionary with its path
                    chosen_db = matching_dbs[0]
                    busco_db_dict[category][i] = os.path.join(busco_db_path, chosen_db)
        print(f"UNLOGGED PASS:\tUpdated dictionary with {db_path_folder_list}")
    return busco_db_dict

# Function to...
def download_assembled_dbs(file_id, base_install_dir):
    """
    Description

    Args:
        file_id (): DESCRIPTION.
        base_install_dir (): DESCRIPTION.
    """
    print(f'UNLOGGED:\tAttempting to download Assembled databases for cleaning...')
    
    # File ID extracted from the Google Drive link for the Databases
    zip_output_path = f'{base_install_dir}/EGAP_Databases/file.zip'  # Output filename
    url = f'https://drive.google.com/uc?id={file_id}'

    # Check if the EGAP_Database directory already exists
    if os.path.exists('~/EGAP/EGAP_Databases/Assembled_Databases'):
        print(f"UNLOGGED PASS:\tSkipping download and extraction: EGAP_Database directory already exists")
    else:
        # Download Databases from Google Drive then Unzip the downloaded and Remove the zip file
        print(f"UNLOGGED:\tAttempting to download EGAP_Databases.zip from Google Drive...")
    
        # Download Zip with gdown
        gdown.download(url, zip_output_path, quiet=False)
        
        # Unzip contents and then remove zip file
        if os.path.exists(zip_output_path):
            try:
                with zipfile.ZipFile(zip_output_path, 'r') as zip_ref:
                    zip_ref.extractall(zip_output_path.replace('/file.zip',''))
            except zipfile.BadZipFile:
                print("The file is not a valid zip file.")
        
        # Remove zip file after downloading
        os.remove(zip_output_path)

def download_and_extract_compleasm(base_install_dir):
    """
    Download and extract compleasm and its dependencies.

    Args:
        base_install_dir (str): Path of the directory the module is installed in.
    """
    print(f'UNLOGGED:\tAttempting to download and extract compleasm...')
    
    # Specify the download URL and output path
    compleasm_url = "https://github.com/huangnengCSU/compleasm/releases/download/v0.2.2/compleasm-0.2.2_x64-linux.tar.bz2"
    tar_output_path = os.path.join(base_install_dir, "compleasm-0.2.2_x64-linux.tar.bz2")
    
    # Download compleasm
    response = requests.get(compleasm_url)
    response.raise_for_status()  # Raise an exception for HTTP errors
    
    # Save the downloaded content to a .tar.bz2 file
    with open(tar_output_path, 'wb') as f:
        f.write(response.content)
    
    # Extract the .tar.bz2 file
    extract_dir = os.path.join(base_install_dir, "compleasm")
    os.makedirs(extract_dir, exist_ok=True)
    with tarfile.open(tar_output_path, 'r:bz2') as tar:
        tar.extractall(path=extract_dir)
    
    # Remove the .tar.bz2 file after extraction
    os.remove(tar_output_path)
    print(f'UNLOGGED PASS:\tSuccessfully downloaded and extracted compleasm to {extract_dir}')

# Function to ensure Python libraries are installed for base_install_dir
def setup_module(base_install_dir, python_libraries, java_program_dict, prereq_list):
    """
    Ensures Python libraries are installed for base_install_dir.

    Args:
        base_install_dir (str): Path of the directory the module is installed in.
        python_libraries (list): List of python libraries to check/install with mamba.
        java_program_dict (dict): Dictionary of Java based programs that require jar files that need checking.
        prereq_list (list): List of third-party programs called that need checking.
    """        
    print(f'UNLOGGED:\tAttempting to setup {base_install_dir}...')
    
    # Check if BUSCO databases are located there, if not download them there
    busco_db_dict = download_busco_dbs(base_install_dir)
       
    # Ensure all other Python libraries are installed
    print(f'UNLOGGED:\tAttempting to install the following Python libraries: {python_libraries}...')
    try:
        for library in python_libraries:
            # Run a subprocess calling conda install for each missing library
            install_module(library)
    except:
        # Print an error if something goes wrong during the installation
        print(f"UNLOGGED ERROR:\t Unable to Install {library}")      
    
    # Check for specific Third-Party JAR Files with adjusted program_dict
    for key, entry in java_program_dict.items():
        # Generate install directories
        install_dir = os.path.join(base_install_dir, key)
        jar_path = download_and_setup(install_dir, entry)
    
    # Check Pipeline Third-Party Prerequisites
    check_prereqs_installed(prereq_list)
    
    # Download the assembled databases for cleaning
    download_assembled_dbs(assembled_db_zip_file_id, base_install_dir)
    
    # Download and install compleasm
    download_and_extract_compleasm(base_install_dir)
    
    print(f'UNLOGGED:\t Completed {base_install_dir} setup successfully')

if __name__ == "__main__":
    # Gather module information
    module_name = 'EGAP'
    python_libraries = ['busco==5.5.0','openjdk==20.0.0', 'nanoq==0.10.0', 'pandas==2.0.3',
                        'biopython==1.81', 'tqdm==4.38.0', 'termcolor==2.3.0', 'beautifulsoup4==4.12.2',
                        'fastqc==0.11.8', 'quast==5.2.0', 'nanostat==1.6.0', 'flye==2.9.2',
                        'bbtools==37.62', 'metaeuk==6.a5d39d9', 'blast==2.14.1', 'bwa==0.7.17',
                        'minimap2==2.26', 'pysam==0.21.0', 'samtools==1.17', 'arcs==1.2.5',
                        'tigmint==1.2.10', 'abyss==2.3.7', 'racon==1.5.0', 'spades==3.15.3',
                        'gdown==4.7.1', 'psutil==5.9.5']
    java_program_dict = {"trimmomatic": "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip",
                         "pilon": "https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar"}
    prereq_list = ["java", "fastqc", "quast.py", "busco", "samtools", "bwa",
                   "makeblastdb", "blastn", "racon", "nanoq", "NanoStat",
                   "spades.py", "tblastn", "flye", "minimap2", "metaeuk", "tqdm"]
    assembled_db_zip_file_id = '1Hj-8tFlJPiOoP_8_Sp4pyWipaQ3zr745'
    ### TODO: ADD KELSEY'S PREREQUISITE: abyss-sealer
    
    # Find the base install directory
    base_install_dir = find_folder(module_name)
    
    # Run main setup for the module
    setup_module(base_install_dir, python_libraries, java_program_dict, prereq_list)
