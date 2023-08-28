"""
Created on Mon Jul 17 11:49:17 2023

@author: ian.michael.bollinger@gmail.com with the help of ChatGPT 4.0

Command Line Example:
    python EGAP_ONT.py --ont_folder /path/to/folder --organism_kingdom STRING --genome_size INTEGER --resource_use INTEGER

The --organism_kingdom must be from the following: Archaea, Bacteria, Fauna, Flora, Funga, or Protista
"""
import os, subprocess, glob, shutil, gzip, random, tempfile, multiprocessing, math, psutil, csv, argparse
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from threading import Thread
from log_print import log_print, generate_log_file
from check_tools import get_env_dir, move_file_up
from EGAP_qc import assess_with_quast, assess_with_busco, assess_with_nanostat
from EGAP_cleaner import clean_dirty_fasta

# Function to extract and combine multiple FASTQ.GZ files into a single FASTQ file
def ont_combine_fastq_gz(ONT_FOLDER, log_file):
    """
    Combine multiple ONT FASTQ.GZ files into a single FASTQ file.

    Args:
        ONT_FOLDER (str): Path to the folder containing ONT FASTQ.GZ files.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.

    Returns:
        combined_ont_fastq_path (str): Path to the combined FASTQ file.
    """
    log_print('Combining ONT FASTQ.GZ files...', log_file)

    # Search for the folder containing multiple .gz files
    ont_raw_data_dir = next((subdir for subdir in glob.glob(os.path.join(ONT_FOLDER, "*"))
                            if os.path.isdir(subdir) and any(file.endswith(".gz") for file in os.listdir(subdir))), None)
    if ont_raw_data_dir is None:
        log_print(f"WARN: No directory containing '.gz' files found within '{ONT_FOLDER}'", log_file)
        log_print(f"Attempting with ONT_FOLDER '{ONT_FOLDER}'", log_file)
        ont_raw_data_dir = ONT_FOLDER

    # Get the base name from the path and append the '.fastq' extension
    base_name = os.path.basename(ont_raw_data_dir)
    combined_ont_fastq_path = os.path.join(ONT_FOLDER, f'{base_name}_ont_combined.fastq')
    return_raw_data_dir = os.path.join(ONT_FOLDER, base_name)
    
    # Check if the combined file already exists
    if os.path.exists(combined_ont_fastq_path):
        log_print(f"PASS: Skipping extraction & combination: Combined fastq file: {combined_ont_fastq_path} already exists", log_file)
        return return_raw_data_dir, combined_ont_fastq_path

    # Combine the FASTQ files
    file_list = glob.glob(os.path.join(ont_raw_data_dir, '*.fastq.gz'))
    with open(combined_ont_fastq_path, 'w') as combined_file:
        for filename in tqdm(file_list, desc='Combining files'):
            with gzip.open(filename, 'rt') as gz_file:
                for record in SeqIO.parse(gz_file, "fastq"):
                    try:
                        SeqIO.write(record, combined_file, "fastq")
                    except Exception as e:
                        log_print(f"ERROR: in FASTQ record: {e}", log_file)
                        raise e

    log_print(f"PASS: Successfully created combined fastq file: {combined_ont_fastq_path}", log_file)
    return return_raw_data_dir, combined_ont_fastq_path

# Function to generate a de novo assembly from ONT Reads with Flye
def assemble_ont_flye(input_fastq, cpu_threads, log_file, GENOME_SIZE):
    """
    Generate a de novo assembly from ONT Reads with Flye.

    Args:
        input_fastq (str): Path to the combined ONT fastq reads.
        cpu_threads (int): Number of threads available for processing.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.
        GENOME_SIZE (int): Expected Mega-Base/Byte (MB) size of genome.      

    Returns:
        assembly_file_path (str): Path to the assembled FASTA file.
        output_directory (str): Path to the Flye Output Directory.
    """
    log_print(f'Generating Flye de novo Assembly from {input_fastq}...', log_file)
    
    # Create the output directory if it does not exist
    output_directory = input_fastq.replace('_ont_combined.fastq','_ont_flye_output')
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    base_name = os.path.basename(input_fastq)
    
    # Path to the assembled file
    assembly_file_path = f"{output_directory}/{base_name.replace('combined.fastq','flye.fasta')}"
    filtered_file_path = assembly_file_path.replace('flye.fasta','flye_filtered.fasta')
    final_assembly_path = move_file_up(assembly_file_path, log_file, move_bool = False)
    final_filtered_path = move_file_up(filtered_file_path, log_file, move_bool = False)
    
    # Check if the assembly already exists
    if os.path.isfile(assembly_file_path):
        log_print(f"PASS:\tSkipping Flye Assembly: {assembly_file_path} already exists", log_file)    
        return assembly_file_path, output_directory
    elif os.path.isfile(final_assembly_path):
        log_print(f"PASS:\tSkipping Flye Assembly: {final_assembly_path} already exists", log_file)    
        return final_assembly_path, output_directory
    elif os.path.isfile(filtered_file_path):
        log_print(f"PASS:\tSkipping Flye Assembly: {filtered_file_path} already exists", log_file)    
        return filtered_file_path, output_directory
    elif os.path.isfile(final_filtered_path):
        log_print(f"PASS:\tSkipping Flye Assembly: {final_filtered_path} already exists", log_file)    
        return final_filtered_path, output_directory
    else:       
        # Construct the Flye command
        flye_command = ["flye",
                        "--nano-hq",
                        input_fastq,
                        "--genome-size",
                        str(GENOME_SIZE),
                        "--out-dir",
                        output_directory,
                        "--threads",
                        str(cpu_threads),
                        "--keep-haplotypes"]
        
        # Run the command
        log_print(f"CMD:\t{' '.join(flye_command)}", log_file)
        flye_result = subprocess.run(flye_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)    
            
        if flye_result.returncode != 0:
            log_print(f"ERROR: {flye_result.stderr}", log_file)
            return None
        else:
            os.rename(os.path.join(output_directory, "assembly.fasta"), assembly_file_path)
            final_assembly_path = move_file_up(assembly_file_path, log_file, move_bool = False)
            log_print("PASS: Successfully generated ONT Flye de novo Assembly", log_file)
        
        # Return the path to the assembled file    
        return final_assembly_path, output_directory

# Main ONT Folder Processing Function
def process_ONT(ONT_FOLDER, CURRENT_ORGANISM_KINGDOM, GENOME_SIZE, PERCENT_RESOURCES, busco_db_dict, log_file):
    """
    Processes basecalled ONT reads into a cleaned and indexed de novo assembly.
    
    Args:
        ONT_FOLDER (str): Path to the folder containing a sub-folder with the ONT reads.
        CURRENT_ORGANISM_KINGDOM (str): Single word Kingdom Description: Archaea, Bacteria, Fauna, Flora, Funga, or Protista.
        GENOME_SIZE (int): Expected size in Mega-Bytes/Bases for the genome (Fungi are about 60-80).
        PERCENT_RESOURCES (int): Percentage of resources that can be utilized for processing (default = 80).
        busco_db_dict (dict): Dictionary of all Kingdoms and their respective BUSCO databases.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.
    
    Returns:
        cleaned_ont_assembly (str): Path to the final cleaned ONT assembly FASTA file.
    """    
    # Get the number of CPUs available on the system
    num_cpus = multiprocessing.cpu_count()

    # Get the amount of RAM (GB) currently available
    mem_info = psutil.virtual_memory()

    # Calculate the number of threads based on PERCENT_RESOURCES for available CPUs & RAM
    cpu_threads = int(math.floor(num_cpus * PERCENT_RESOURCES))
    ram_gb = int(mem_info.total / (1024.0 ** 3) * PERCENT_RESOURCES)
    
    # Take ONT Folder with FASTQ Files, extract, combine
    ont_raw_data_dir, combined_ont_fastq = ont_combine_fastq_gz(ONT_FOLDER, log_file)

## QUALITY CONTROL CHECK AREA
    # Quality Control Check of Combined ONT Reads with NanoStat
    nanostat_dir = assess_with_nanostat(combined_ont_fastq, cpu_threads, log_file)

    # Flye Assembly of combined ONT FASTQs
    ont_flye_assembly, flye_dir = assemble_ont_flye(combined_ont_fastq, cpu_threads, log_file, GENOME_SIZE)
    
    # Decontamination of ONT Flye Assembly
    cleaned_ont_assembly, removed_csv = clean_dirty_fasta(ont_flye_assembly, ONT_FOLDER, CURRENT_ORGANISM_KINGDOM, log_file)                
    base_name = cleaned_ont_assembly.split('/')[-1].split('_ont')[0]
    
    # Cleanup ONT assembly bulk files and keep only the items_to_keep found in the ONT_FOLDER
    items_to_keep = [nanostat_dir.replace('_ont_combined_NanoStat',''), nanostat_dir,
                     nanostat_dir.replace('_ont_combined_NanoStat',f'_ont_flye_filtered_{busco_db_dict[CURRENT_ORGANISM_KINGDOM][0].split("/")[-1]}_busco'),
                     nanostat_dir.replace('_ont_combined_NanoStat',f'_ont_flye_filtered_{busco_db_dict[CURRENT_ORGANISM_KINGDOM][1].split("/")[-1]}_busco'),
                     nanostat_dir.replace('_ont_combined_NanoStat','_ont_flye_filtered_quast'),
                     combined_ont_fastq, ont_flye_assembly,
                     cleaned_ont_assembly, removed_csv,
                     f'{ONT_FOLDER}/{base_name}_ont_flye_bwa_aligned.bam']
    
    chopping_block = os.listdir(ONT_FOLDER)
    chopping_block = [os.path.join(ONT_FOLDER, item) for item in chopping_block]
        
    for item in chopping_block:
        item_path = os.path.join(ONT_FOLDER, item)
        if item_path not in items_to_keep:
            if os.path.isfile(item_path):
                os.remove(item_path)
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)

## QUALITY CONTROL CHECK AREA
    if __name__ != "__main__":
        pass
    else:
        # Quality Control Check Cleaned ONT Assembly with QUAST
        quast_thread = Thread(target = assess_with_quast, args = (cleaned_ont_assembly, log_file, cpu_threads))
        quast_thread.start()
          
        # Quality Control Check Pilon Polished Assembly with BUSCO agasint first database
        first_busco_thread = Thread(target = assess_with_busco, args = (cleaned_ont_assembly, log_file, busco_db_dict[CURRENT_ORGANISM_KINGDOM][0]))
        first_busco_thread.start()
        
        # Quality Control Check Pilon Polished Assembly with BUSCO agasint second database
        second_busco_thread = Thread(target = assess_with_busco, args = (cleaned_ont_assembly, log_file, busco_db_dict[CURRENT_ORGANISM_KINGDOM][1]))
        second_busco_thread.start()
            
        # Wait for all QC threads to finish
        quast_thread.join()
        first_busco_thread.join()
        second_busco_thread.join()
    
    return cleaned_ont_assembly

## Debuging Main Space & Example
if __name__ == "__main__":
    print('EGAP Oxford Nanopore Technologies (ONT) Pipeline')    
    # Get working environment information
    environment_dir = get_env_dir()
    
    # Argument Parsing
    parser = argparse.ArgumentParser(description='Process ONT Folder and Genome Size')
    
    # Default values
    default_folder = f'{environment_dir}/Entheome/Ps_aff_hopii/MODULAR_TEST/ONT_MinION/'
    default_organism_kingdom = 'Funga'
    default_genome_size = 60
    default_percent_resources = 80
    
    # Add arguments with default values
    parser.add_argument('--ont_folder', default = default_folder,
                        help = f'Path to the ONT Folder. (default: {default_folder})')
    parser.add_argument('--organism_kingdom',default = default_organism_kingdom,
                        help = f'Kingdom the current organism data belongs to. (default: {default_organism_kingdom})')
    parser.add_argument('--genome_size', type = int, default = default_genome_size,
                        help = f'Genome Size. (default: {default_genome_size})')
    parser.add_argument('--resource_use', type = int, default = default_percent_resources,
                        help = f'Percent of Resources to use. (default: {default_percent_resources})')
    
    # Parse the arguments
    args = parser.parse_args()
    
    ONT_FOLDER = args.ont_folder
    CURRENT_ORGANISM_KINGDOM = args.organism_kingdom
    GENOME_SIZE = args.genome_size
    PERCENT_RESOURCES = (args.resource_use/100)
    
    # Generate log file with the desired behavior
    debug_log = f'{ONT_FOLDER}ONT_log.tsv'
    log_file = generate_log_file(debug_log, use_numerical_suffix=False)
    
    # Generate BUSCO Database Dictionary
    busco_db_dict = {'Archaea':  [f'{environment_dir}/EGAP/BUSCO_Databases/archaea_odb10',
                                  f'{environment_dir}/EGAP/BUSCO_Databases/euryarchaeota_odb10',],
                     'Bacteria': [f'{environment_dir}/EGAP/BUSCO_Databases/actinobacteria_phylum_odb10',
                                  f'{environment_dir}/EGAP/BUSCO_Databases/proteobacteria_odb10',],
                     'Fauna':    [f'{environment_dir}/EGAP/BUSCO_Databases/vertebrata_odb10',
                                  f'{environment_dir}/EGAP/BUSCO_Databases/arthropoda_odb10',],
                     'Flora':    [f'{environment_dir}/EGAP/BUSCO_Databases/eudicots_odb10',
                                  f'{environment_dir}/EGAP/BUSCO_Databases/liliopsida_odb10'],
                     'Funga':    [f'{environment_dir}/EGAP/BUSCO_Databases/basidiomycota_odb10',
                                  f'{environment_dir}/EGAP/BUSCO_Databases/agaricales_odb10'],
                     'Protista': [f'{environment_dir}/EGAP/BUSCO_Databases/alveolata_odb10',
                                  f'{environment_dir}/EGAP/BUSCO_Databases/euglenozoa_odb10']}
    
    # Run main ONT Cleaning function
    cleaned_ont_assembly = process_ONT(ONT_FOLDER, CURRENT_ORGANISM_KINGDOM, GENOME_SIZE, PERCENT_RESOURCES, busco_db_dict, log_file)
