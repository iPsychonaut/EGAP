# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 11:27:53 2023

@author: ian.michael.bollinger@gmail.com with the help of ChatGPT 4.0

Command Line Example:
    python EGAP_main.py -i /path/to/folder -k STRING -g INTEGER -p STRING -d [same/different] -r INTEGER -a [0/1]

The -i, --input_dir must have sub-folders with names containing either 'illumina' (sub-folder with raw PE150 .fq.gz files and their matching MD5.txt file) AND a folder with 'ont' (sub-folder with raw pass .fastq.gz files)
The -k, --organism_kingdom must be from the following: Archaea, Bacteria, Fauna, Flora, Funga, or Protista
The -g, --genome_size is a number in Mega-Bytes/Bases that the expected genome is to be. 
The -p, --primer_type must be a string similar to 'TruSeq3-PE' to represent the Illumina primer type to use with trimmomatic.
"""
# Base Python Imports
import sys, subprocess, argparse, multiprocessing, math, os, platform, shutil
from threading import Thread

# Required Python Imports
import psutil

# Custom Python Imports
from EGAP_setup import setup_module, install_module

# Define and parse command line arguments
parser = argparse.ArgumentParser(description="Run Entheome Illumina+ONT Pipeline")

# Default values
default_folder = f'/mnt/e/Entheome/Ps_aff_hopii/MODULAR_TEST'
default_organism_kingdom = 'Funga'
default_genome_size = 60
default_primer = 'TruSeq3-PE'
default_organism_data = 'same' # same = Ian's same organism pipeline; different = Kelsey's Slot Lab different organism pipeline
default_percent_resources = 80
default_install = '0' # 0 = Skip install; 1 = Attempt to install

# Add organism data argument
parser.add_argument('--input_dir', '-i',
                    type = str, default = default_folder,
                    help='The root folder containing a folder with "illumina" in its name and one with "ont" in its name')
parser.add_argument('--organism_kingdom', '-k', 
                    type = str, default = default_organism_kingdom,
                    choices = ['Archaea', 'Bacteria', 'Fauna', 'Flora', 'Funga', 'Protista'],
                    help = f'Kingdom the current organism data belongs to. (default: {default_organism_kingdom})')
parser.add_argument('--genome_size', '-g',
                    type = int, default = default_genome_size,
                    help = f'Genome Size. (default: {default_genome_size})')
parser.add_argument('--primer_type', '-p',
                    type = str, default = default_primer,
                    help = f'Type of Illumina Primers used. (default: {default_primer})')
parser.add_argument('--organism_data', '-d',
                    type = str, default = default_organism_data,
                    choices = ["same", "different"],
                    help = 'Indicate if the provided data are generated from the same organism or different organisms')
parser.add_argument('--resource_use', '-r',
                    type = int, default = default_percent_resources,
                    help = f'Percent of Resources to use. (default: {default_percent_resources})')
parser.add_argument('--attempt_install', '-a',
                    type = str, default = default_install,
                    choices = ['0', '1'], 
                    help = f'Flag to indicate if the installation should be attempted 0 = No; 1 = Yes (default: {default_install}).')

# Parse the arguments
args = parser.parse_args()
BASE_FOLDER = args.input_dir
CURRENT_ORGANISM_KINGDOM = args.organism_kingdom
GENOME_SIZE = args.genome_size
ILLU_PRIMER_TYPE = args.primer_type
ORGANISM_DATA = args.organism_data
PERCENT_RESOURCES = (args.resource_use/100)
EGAP_ATTEMPT_INSTALL = args.attempt_install

# Check if already tried installing required Python libraries
if EGAP_ATTEMPT_INSTALL == '1':  
    # Gather module information
    module_name = 'EGAP'
    python_libraries = ['requests==2.31.0', 'busco==5.5.0','openjdk==20.0.0', 'nanoq==0.10.0', 'pandas==2.0.3',
                        'biopython==1.81', 'tqdm==4.38.0', 'termcolor==2.3.0', 'beautifulsoup4==4.12.2',
                        'fastqc==0.11.8', 'quast==5.2.0', 'nanostat==1.6.0', 'flye==2.9.2',
                        'bbtools==37.62', 'metaeuk==6.a5d39d9', 'blast==2.14.1', 'bwa==0.7.17',
                        'minimap2==2.26', 'pysam==0.21.0', 'samtools==1.17', 'arcs==1.2.5',
                        'tigmint==1.2.10', 'abyss==2.3.7', 'racon==1.5.0', 'spades==3.15.3',
                        'psutil==5.9.5']
    java_program_dict = {"trimmomatic": "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip",
                         "pilon": "https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar"}
    prereq_list = ["java", "fastqc", "quast.py", "busco", "samtools", "bwa",
                   "makeblastdb", "blastn", "racon", "nanoq", "NanoStat",
                   "spades.py", "tblastn", "flye", "minimap2", "metaeuk", "tqdm"]
    ### TODO: ADD KELSEY'S PREREQUISITE: abyss-sealer
    
    # Run main setup for the module
    setup_module(module_name, python_libraries, java_program_dict, prereq_list)
    
    # Set the flag to indicate we've attempted installation
    os.environ['EGAP_ATTEMPT_INSTALL'] = '0'
    
    # Reset the script with loaded libraries
    print('UNLOGGED PASS:\tRestarting EGAP')
    os.execv(sys.executable, ['python'] + sys.argv + ['--attempt_install', '0'])
elif EGAP_ATTEMPT_INSTALL == '0':
    print(f'UNLOGGED:\tSkipping Module setup and prerequisite installation')
    
# Custom Python Imports
from check_tools import get_env_dir
from log_print import log_print, generate_log_file
from EGAP_ONT import process_ONT
from EGAP_illumina import process_illumina
from EGAP_pilon_polish import final_pilon_polish
from EGAP_qc import assess_with_fastqc, assess_with_quast, assess_with_compleasm

# Main function to run the Entheome Pipeline developed by Ian
def PILON_POLISH_PIPELINE(BASE_FOLDER, CURRENT_ORGANISM_KINGDOM, GENOME_SIZE, ILLU_PRIMER_TYPE, PERCENT_RESOURCES, busco_db_dict, log_file):
    """
    Assembles a final fasta file based on data from a single lineage of an organism using Ian's Pipeline.

    Args:
        BASE_FOLDER (str): Path the main folder containaing all sub-folder data for processing.
        CURRENT_ORGANISM_KINGDOM (str): Single word Kingdom Description: Archaea, Bacteria, Fauna, Flora, Funga, or Protista.
        GENOME_SIZE (int): Expected size of the organims genome in Mega-Bytes/Bases.
        ILLU_PRIMER_TYPE (str): A string representing the Illumina primer type to use with trimmomatic.
        PERCENT_RESOURCES (int): Amount of resources availalbe for use during processing.
        busco_db_dict (dict): Dictionary of all Kingdoms and their respective BUSCO databases.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.

    Returns:
        pilon_output (str): Path to the final Polished Pilon Assembly.
    """
    # Get the number of CPUs available on the system
    num_cpus = multiprocessing.cpu_count()
        
    # Get the amount of RAM (GB) currently available
    mem_info = psutil.virtual_memory()
    
    # Calculate the number of threads as based on PERCENT_RESOURCES of available CPUs & RAM
    cpu_threads = round(int(num_cpus * PERCENT_RESOURCES),0)
    ont_cpus = round(int(cpu_threads/2),0)
    illu_cpus = cpu_threads - ont_cpus

    ram_gb = int(mem_info.total / (1024.0 ** 3) * PERCENT_RESOURCES)
    ont_ram_gb = round(int(ram_gb/2),2)
    illu_ram_gb = ram_gb - ont_cpus
    
    ont_percent_resources = float(round(PERCENT_RESOURCES/2, 1))
    illu_percent_resources = PERCENT_RESOURCES - ont_percent_resources
    
    # Set the boolean checks for each folder to False
    processed_same_illumina = False
    processed_same_ont = False

    # Walk through the root folder
    for folder_name, subfolders, file_names in os.walk(BASE_FOLDER):
        if folder_name == BASE_FOLDER:
            continue        
        elif 'illumina' in folder_name.lower() and not processed_same_illumina:
            log_print(f"Starting Illumina Raw Data Pipeline on {folder_name}...", log_file)
            # Run main Illumina Trimming function
            fq_paired_list, fastqc_output_dirs, data_type = process_illumina(folder_name, ILLU_PRIMER_TYPE, illu_percent_resources, log_file)
            
            # Quality Control Check A FastQC on Trimmed Illumina Reads
            r_fastqc_cpus = int(math.floor(illu_cpus/2))
            f_fastqc_cpus = cpu_threads - r_fastqc_cpus
            
            if data_type == 'PE':
                f_fastqc_thread = Thread(target = assess_with_fastqc, args = (fq_paired_list[0], fastqc_output_dirs[0], log_file, f_fastqc_cpus))
                f_fastqc_thread.start()            
                r_fastqc_thread = Thread(target = assess_with_fastqc, args = (fq_paired_list[1], fastqc_output_dirs[1], log_file, r_fastqc_cpus))
                r_fastqc_thread.start()
            elif data_type == 'SE':
                fastqc_thread = Thread(target = assess_with_fastqc, args = (fq_paired_list[0], fastqc_output_dirs[0], log_file, f_fastqc_cpus))
                fastqc_thread.start()     
            # Register that Illumina Folder has been processed
            processed_same_illumina = True

        # Check if 'ont' is in the base folder name, and if so, run the ont Raw Data Pipeline
        elif 'ont' in folder_name.lower() and not processed_same_ont:
            log_print(f"Starting ONT Raw Data Pipeline on {folder_name}...", log_file)
            # Run main ONT Cleaning function
            cleaned_ont_assembly = process_ONT(folder_name, CURRENT_ORGANISM_KINGDOM, GENOME_SIZE, ont_percent_resources, busco_db_dict, log_file)
            
            # Quality Control Check Cleaned ONT Assembly with QUAST
            ont_quast_thread = Thread(target = assess_with_quast, args = (cleaned_ont_assembly, log_file, ont_cpus))
            ont_quast_thread.start()
            
            # Quality Control Check Flye de Novo Assembly with BUSCO agasint first database
            first_busco_thread = Thread(target = assess_with_compleasm, args = (cleaned_ont_assembly, log_file, busco_db_dict[CURRENT_ORGANISM_KINGDOM][0]))
            first_busco_thread.start()
            
            # Quality Control Check Flye de Novo Assembly with BUSCO agasint second database
            second_busco_thread = Thread(target = assess_with_compleasm, args = (cleaned_ont_assembly, log_file, busco_db_dict[CURRENT_ORGANISM_KINGDOM][1]))
            second_busco_thread.start()

            # Register that ONT Folder has been processed
            processed_same_ont = True
    
    # Wait for all QC threads to finish
    if data_type == 'PE':
        f_fastqc_thread.join()
        r_fastqc_thread.join()
        ont_quast_thread.join()
        first_busco_thread.join()
        second_busco_thread.join()
    elif data_type == 'SE':
        fastqc_thread.join()
        ont_quast_thread.join()
        first_busco_thread.join()
        second_busco_thread.join()

    log_print(f"PASS:\tAll Initial QC Checks Complete for {BASE_FOLDER}...", log_file)
    
    # If the Trimmed Illumina Forward & Reverse Reads completed AND the ONT Reads were cleaned -> perform Pilon Polish
    if all(os.path.exists(file) for file in fq_paired_list) and os.path.exists(cleaned_ont_assembly):
        log_print(f'Starting Polishing Process on {fq_paired_list} and {cleaned_ont_assembly}', log_file)      
        
        if data_type == 'PE':
            # Run main Pilon Polish with Paired-End reads
            pilon_output = final_pilon_polish(cleaned_ont_assembly, fq_paired_list, CURRENT_ORGANISM_KINGDOM, PERCENT_RESOURCES, busco_db_dict, log_file)
        elif data_type == 'SE':
            # Run main Pilon Polish with Paired-End reads
            pilon_output = final_pilon_polish(cleaned_ont_assembly, fq_paired_list, CURRENT_ORGANISM_KINGDOM, PERCENT_RESOURCES, busco_db_dict, log_file)

        # Quality Control Check Pilon Polished Assembly with QUAST
        pilon_quast_thread = Thread(target = assess_with_quast, args = (pilon_output, log_file, cpu_threads))
        pilon_quast_thread.start()

        # Quality Control Check Pilon Polished Assembly with BUSCO agasint first database
        first_pilon_busco_thread = Thread(target = assess_with_compleasm, args = (pilon_output, log_file, busco_db_dict[CURRENT_ORGANISM_KINGDOM][0]))
        first_pilon_busco_thread.start()

        # Quality Control Check Pilon Polished Assembly with BUSCO agasint the second database
        second_pilon_busco_thread = Thread(target = assess_with_compleasm, args = (pilon_output, log_file, busco_db_dict[CURRENT_ORGANISM_KINGDOM][1]))
        second_pilon_busco_thread.start()   
    else:
        log_print('ERROR:\tcleaned_ont_assembly or fq_paired_list not found', log_file)

    # Wait for all QC threads to finish
    pilon_quast_thread.join()
    first_pilon_busco_thread.join()
    second_pilon_busco_thread.join()
    
    log_print(f"PASS:\tAll Final QC Checks Complete for {BASE_FOLDER}", log_file)
    log_print(f"EGAP PIPELINE COMPLETE", log_file)
    
    return pilon_output

# Main function to run the Slot Lab Hybrid Organism Pipeline developed by Kelsey
def SPADES_HYBRID_PIPELINE(BASE_FOLDER, CURRENT_ORGANISM_KINGDOM, GENOME_SIZE, ILLU_PRIMER_TYPE, PERCENT_RESOURCES, busco_db_dict, log_file):
    """
    Assembles a final fasta file based on data from two different lineages of the same organism using the Slot Lab Pipeline.

    Args:
        BASE_FOLDER (str): Path the main folder containaing all sub-folder data for processing.
        CURRENT_ORGANISM_KINGDOM (str): Single word Kingdom Description: Archaea, Bacteria, Fauna, Flora, Funga, or Protista.
        GENOME_SIZE (int): Expected size of the organims genome in Mega-Bytes/Bases.
        ILLU_PRIMER_TYPE (str): A string representing the Illumina primer type to use with trimmomatic.
        PERCENT_RESOURCES (int): Amount of resources availalbe for use during processing.
        busco_db_dict (dict): Dictionary of all Kingdoms and their respective BUSCO databases.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.

    Returns:
        spades_output (str): Path to the final SPAdes Assembly.
    """
    # Get the number of CPUs available on the system
    num_cpus = multiprocessing.cpu_count()
    
    # Get the amount of RAM (GB) currently available
    mem_info = psutil.virtual_memory()
    
    # Calculate the number of threads as 80% of available CPUs & RAM
    cpu_threads = int(math.floor(num_cpus * PERCENT_RESOURCES))
    ram_gb = int(mem_info.total / (1024.0 ** 3) * PERCENT_RESOURCES)
    
    # Set the boolean checks for each folder to False
    processed_different_illumina = False
    processed_different_ont = False

    # Walk through the root folder
    for folder_name, subfolders, file_names in os.walk(BASE_FOLDER):
        if folder_name == BASE_FOLDER:
            continue        

        # Check if 'ont' is in the base folder name, and if so, run the ont Raw Data Pipeline
        elif 'ont' in folder_name.lower() and not processed_different_ont:
            log_print(f"Starting ONT Raw Data Pipeline on {folder_name}...", log_file)

            # TODO: GENERATION MODULE AND FUNCTIONS FOR:
            """
            Nanoq filtering of ONT Basecalled Reads processed with to remove reads below Q10

            NanoStat analysis of Filtered ONT FASTQ

            Flye de novo Assembly of Filtered ONT FASTQ
            """

            # TODO: QUAST analysis of ONT Flye Assembly

            # TODO: BUSCO analysis of ONT Flye Assembly against first database

            # TODO: BUSCO analysis of ONT Flye Assembly against second database
            
            # Register that Illumina Folder has been processed
            processed_different_ont = True
            
        elif 'illumina' in folder_name.lower() and not processed_different_illumina:
            log_print(f"Starting Illumina Raw Data Pipeline on {folder_name}...", log_file)    
            
            # TODO: GENERATION MODULE AND FUNCTIONS FOR:
            """
            Raw Illumina PE150 Reads trimmed with Trimmomatic
            """
            
            # TODO: FastQC Analysis of Trimmed Illumina Forward Paired Reads
            
            # TODO: FastQC Analysis of Trimmed Illumina Reverse Paired Reads
            
            
            # TODO: GENERATION MODULE AND FUNCTIONS FOR:
            """
            SPAdes Assembly of Raw Illumina PE150 Reads
                --careful, -k 21,33,55,77,99,121
            """
            
            # TODO: QUAST analysis of Illumina SPAdes Assembly

            # TODO: BUSCO analysis of Illumina SPAdes Assembly against first database

            # TODO: BUSCO analysis of Illumina SPAdes Assembly against second database
                        
            # Register that ONT Folder has been processed
            processed_different_illumina = True

    if all(os.path.exists(file) for file in kelsey_ont_output) and os.path.exists(kelsey_illumina_output):
        log_print(f'Starting Scaffolding on {kelsey_ont_output} and {kelsey_illumina_output}', log_file)
        
        # TODO: GENERATION MODULE AND FUNCTIONS FOR:
        """
        SPAdes Iteratively Scaffolding of the Illumina SPAdes Assembly four (4x) times to the Filtered ONT FASTQ
            option = scaffold

        abyss-sealer gap filling of Final Scaffolded Hybrid Assembly
            bloom-size = 20G, kmer=90,80,70,60,50,40,30

        minimap2 mapping of Trimmed Illumina Forward & Reverse Paired Reads to the Scaffold Hybrid Assembly
        
        Racon Iteratively Polishing of the Assembly with the Illumina Map two (2x) times
        """
        
        # TODO: QUAST analysis of Racon Polished SPAdes Assembly

        # TODO: BUSCO analysis of Racon Polished SPAdes Assembly against Basidiomycota

        # TODO: BUSCO analysis of Racon Polished SPAdes Assembly against Agaricales

    else:
        log_print(f'ERROR:\t{kelsey_ont_output} or {kelsey_illumina_output} not found.\n', log_file)
    
    spades_output = '' # PLACEHOLDER
    return spades_output

## Debuging Main Space & Example
if __name__ == "__main__":    
    # Get working environment information
    environment_dir = get_env_dir(BASE_FOLDER)
    
    # Generate Main Logfile
    debug_log = f'{BASE_FOLDER}EGAP_log.tsv'
    log_file = generate_log_file(debug_log, use_numerical_suffix=False)
    log_print('RUNNING ENTHEOME GENOME EXTRACTION PIPELINE - EGAP', log_file)
    
    # Generate BUSCO Database Dictionary
    busco_db_dict = {'Archaea':  [f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/archaea_odb10',
                                  f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/euryarchaeota_odb10',],
                     'Bacteria': [f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/actinobacteria_phylum_odb10',
                                  f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/proteobacteria_odb10',],
                     'Fauna':    [f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/vertebrata_odb10',
                                  f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/arthropoda_odb10',],
                     'Flora':    [f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/eudicots_odb10',
                                  f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/liliopsida_odb10'],
                     'Funga':    [f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/basidiomycota_odb10',
                                  f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/agaricales_odb10'],
                     'Protista': [f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/alveolata_odb10',
                                  f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/euglenozoa_odb10']}
   
    # Modify your pipeline function to accept the new argument and adjust its behavior accordingly
    if ORGANISM_DATA == 'same':
        log_print("RUN IAN'S SAME ORGANISM DATA PIPELINE", log_file)
        pilon_output = PILON_POLISH_PIPELINE(BASE_FOLDER, CURRENT_ORGANISM_KINGDOM, GENOME_SIZE, ILLU_PRIMER_TYPE, PERCENT_RESOURCES, busco_db_dict, log_file)
    elif ORGANISM_DATA == 'different':
        log_print("RUN KELSEY'S DIFFERENT ORGANISM DATA PIPELINE", log_file)
        spades_output = SPADES_HYBRID_PIPELINE(BASE_FOLDER, CURRENT_ORGANISM_KINGDOM, GENOME_SIZE, ILLU_PRIMER_TYPE, PERCENT_RESOURCES, busco_db_dict, log_file)
