# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 12:20:03 2023

@author: ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting

Command Line Example:
    python EGAP_Illumina_Only.py -i /path/to/folder -r FLOAT

The -i, --input_folder must have input FASTQ.GZ file, matching primers.txt and Index.txt
The -r, --resources Percentage of resources to use (0.01-1.00; default: 0.2) 

"""
import os, argparse, gzip, shutil

from EGAP_QC import illu_only_qc_checks
from EGAP_Tools import log_print, initialize_logging_environment, run_subprocess_cmd, find_file, get_resource_values

# Global output_area variable
CPU_THREADS = 1
DEFAULT_LOG_FILE = None
ENVIRONMENT_TYPE = None

def parse_fq_files(input_folder):
    """
    Searches for and returns paths of FASTQ files in a given folder.

    Args:
        input_folder (str): Path to the folder containing FASTQ files.

    Returns:
        list of str: Paths of found "_1.fq" and "_2.fq" files, or None if not found.
    """
    log_print("Parsing input folder for .fq files...")
    fq_files = {'_1.fq': None, '_2.fq': None}

    # First, search for .fq files
    for end in ['_1.fq', '_2.fq']:
        for file in os.listdir(input_folder):
            if file.endswith(end):
                fq_files[end] = os.path.join(input_folder, file)

    # Check if both .fq files are found
    if all(fq_files.values()):
        return [fq_files['_1.fq'], fq_files['_2.fq']]

    # If not, search for .fq.gz files and unzip them
    for end in ['_1.fq', '_2.fq']:
        for file in os.listdir(input_folder):
            if file.endswith(end + '.gz'):
                gz_file_path = os.path.join(input_folder, file)
                unzipped_file_path = gz_file_path.rstrip('.gz')

                # Unzip the file
                with gzip.open(gz_file_path, 'rb') as f_in:
                    with open(unzipped_file_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                fq_files[end] = unzipped_file_path
                break

    return [fq_files['_1.fq'], fq_files['_2.fq']]

def trimmomatic_prep(input_fq_list, CPU_THREADS):
    """
    Prepares and executes Trimmomatic for sequence trimming.

    Args:
        input_fq_list (list of str): Paths to input FASTQ files.
        CPU_THREADS (int): Number of CPU threads to use.

    Returns:
        tuple of str: Paths to trimmed forward and reverse FASTQ files.
    """  
    log_print("Trimming Illumina reads with Trimmomatic...")
    # Build output file path names based on input_fq_list
    trimmo_f_pair_path = input_fq_list[0].replace('_1.fq','_forward_paired.fq')
    trimmo_f_unpair_path = input_fq_list[0].replace('_1.fq','_forward_unpaired.fq')
    trimmo_r_pair_path = input_fq_list[1].replace('_2.fq','_reverse_paired.fq')
    trimmo_r_unpair_path = input_fq_list[1].replace('_2.fq','_reverse_unpaired.fq')
    
    # Default path to Trimmomatic paths
    default_trimmo_path = "trimmomatic-0.39.jar"
    default_trimmo_adapters_path = "TruSeq3-PE.fa"
    trimmo_path = find_file(default_trimmo_path)
    trimmo_adapters_path = find_file(default_trimmo_adapters_path)   
    
    # Make variables for other commands as needed: {trimmo_adapters_path}:{}:{}:{}:{}, SLIDINGWINDOW, MINLEN
    # Joining command list into a single string
    trimmo_cmd =["java", "-jar", trimmo_path,
                 "PE","-phred33","-threads", str(CPU_THREADS),
                 input_fq_list[0], input_fq_list[1],
                 trimmo_f_pair_path, trimmo_f_unpair_path,
                 trimmo_r_pair_path, trimmo_r_unpair_path,
                 f"ILLUMINACLIP:{trimmo_adapters_path}:2:30:10:11",
                 # 'HEADCROP:10', 'CROP:145',
                 "SLIDINGWINDOW:50:25", "MINLEN:125"]

    if os.path.exists(trimmo_f_pair_path) and os.path.exists(trimmo_r_pair_path):
        log_print("PASS:\tSkipping Trimmomatic, output files already exist")
    else:
        # Executing the command
        _ = run_subprocess_cmd(trimmo_cmd, False)
    
    return trimmo_f_pair_path, trimmo_r_pair_path

def bbduk_map(trimmo_f_pair_path, trimmo_r_pair_path):
    """
    Runs BBDuk for quality trimming and adapter removal.

    Args:
        trimmo_f_pair_path (str): Path to forward paired FASTQ file.
        trimmo_r_pair_path (str): Path to reverse paired FASTQ file.

    Returns:
        tuple of str: Paths to forward and reverse mapped FASTQ files.
    """
    bbduk_f_map_path = trimmo_f_pair_path.replace('_forward_paired.fq','_forward_mapped.fq')
    bbduk_r_map_path = trimmo_r_pair_path.replace('_reverse_paired.fq','_reverse_mapped.fq')
    
    # Default path to BBDuk adapters
    default_adapters_path = "adapters.fa"
    default_bbduk_path = "bbduk.sh"
    adapters_path = find_file(default_adapters_path)
    bbduk_path = find_file(default_bbduk_path)
        
    # bbduk Map Illumina Reads with run_subprocess_cmd shell_check = True
    bbduk_cmd = [bbduk_path,
                f"in1={trimmo_f_pair_path}", f"in2={trimmo_r_pair_path}",
                f"out1={bbduk_f_map_path}", f"out2={bbduk_r_map_path}",
                f"ref={adapters_path}", "ktrim=r", "k=23", "mink=11", "hdist=1",
                "tpe","tbo", "qtrim=rl", "trimq=20"]
    if os.path.exists(bbduk_f_map_path) and os.path.exists(bbduk_r_map_path):
        log_print("PASS:\tSkipping bbduk, Mapped outputs alredy exist")
    else:
        _ = run_subprocess_cmd(bbduk_cmd, False)
    
    return bbduk_f_map_path, bbduk_r_map_path

def clumpify_dedup(bbduk_f_map_path, bbduk_r_map_path):
    """
    De-duplicates FASTQ files using Clumpify.

    Args:
        bbduk_f_map_path (str): Path to forward mapped FASTQ file.
        bbduk_r_map_path (str): Path to reverse mapped FASTQ file.

    Returns:
        tuple of str: Paths to forward and reverse de-duplicated FASTQ files.
    """
    clump_f_dedup_path = bbduk_f_map_path.replace('_forward_mapped.fq','_forward_dedup.fq')
    clump_r_dedup_path = bbduk_r_map_path.replace('_reverse_mapped.fq','_reverse_dedup.fq')
    
    default_clumpify_path = "clumpify.sh"
    
    clumpify_path = find_file(default_clumpify_path)
    
    clumpify_cmd = [clumpify_path,
                    f"in={bbduk_f_map_path}",
                    f"in2={bbduk_r_map_path}",
                    f"out={clump_f_dedup_path}",
                    f"out2={clump_r_dedup_path}",
                    "dedupe"]
    
    if os.path.exists(clump_f_dedup_path) and os.path.exists(clump_r_dedup_path):
        log_print("PASS:\tSkipping clumpify, deduplicated outputs already exist")
    else:
        _ = run_subprocess_cmd(clumpify_cmd, False)
    
    return clump_f_dedup_path, clump_r_dedup_path

def parse_bbmerge_output(output):
    """
    Extracts average insert size and standard deviation from BBMerge output.

    Args:
        output (str): Output string from BBMerge command.

    Returns:
        tuple: Average insert size and standard deviation.
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

def masurca_config_gen(input_folder, input_fq_list, clump_f_dedup_path, clump_r_dedup_path, CPU_THREADS):
    """
    Generates a MaSuRCA config and runs genome assembly.

    Args:
        input_folder (str): Path for config file generation.
        input_fq_list (list of str): Paths to input FASTQ files.
        clump_f_dedup_path (str): Path to forward deduplicated file.
        clump_r_dedup_path (str): Path to reverse deduplicated file.
        CPU_THREADS (int): Number of CPU threads to use.

    Returns:
        str: Path to the scaffolded assembly file.
    """
    bbmap_out_path = clump_f_dedup_path.replace('_forward_dedup.fq','_bbmap')
    
    default_bbmerge_path = "bbmerge.sh"
    bbmerge_path = find_file(default_bbmerge_path)
    
    # bbmerge to Determine average insert size and standard deviation
    bbmerge_cmd = [bbmerge_path,
                   f"in1={clump_f_dedup_path}",
                   f"in2={clump_r_dedup_path}",
                   f"out={bbmap_out_path}",
                   "ihist=insert_size_histogram.txt"]
    
    if os.path.exists(bbmap_out_path):
        log_print("PASS:\tSkipping bbmap, insert size histogram output already exists")
        
        # Open the file for writing
        with open(f"{input_folder}/insert_size_histogram.txt", 'r') as file:
            for line in file:
                if "#Mean	" in line:
                    avg_insert = round(float(line.replace("#Mean	","").replace("\n","")),0)
                if "#STDev	" in line:
                    std_dev = round(float(line.replace("#STDev	","").replace("\n","")),0)
    else:
        bbmerge_output = run_subprocess_cmd(bbmerge_cmd, False)
   
        # Parse standard out text for avg_insert and std_dev
        try:
            avg_insert, std_dev = parse_bbmerge_output(bbmerge_output)
        except:
            avg_insert = 251
            std_dev = 30
        
    # Define the content to be written to the file
    # TODO: Make variables for other commands as needed: USE_LINKING_MATES, LIMIT_JUMP_COVERAGE, CA_PARAMETERS, MEGA_READS_ONE_PASS, LHE_COVERAGE, KMER_COUNT_THRESHOLD, JF_SIZE, DO_HOMOPOLYMER_TRIM
    config_content = ["DATA\n",
                      f"PE= pe {avg_insert} {std_dev} {clump_f_dedup_path} {clump_r_dedup_path}\n",
                      "END\n",
                      "PARAMETERS\n",
                      "GRAPH_KMER_SIZE = auto\n",
                      "USE_LINKING_MATES = 1\n",
                      "LIMIT_JUMP_COVERAGE = 300\n",
                      "CA_PARAMETERS = cgwErrorRate=0.15\n",
                      "MEGA_READS_ONE_PASS=0\n",
                      "LHE_COVERAGE=35\n",
                      "KMER_COUNT_THRESHOLD = 1\n",
                      f"NUM_THREADS = {CPU_THREADS}\n",
                      "JF_SIZE = 2500000000\n",
                      "DO_HOMOPOLYMER_TRIM = 0\n",
                      "END"]
    
    # Specify the file name
    config_path = f"{input_folder}/masurca_config_file.txt"
    
    # Open the file for writing
    with open(config_path, 'w') as file:
        for entry in config_content:
            file.write(entry)

    # Masurca Assembly to generate config assembly.sh
    masurca_config_cmd = ["masurca", "masurca_config_file.txt"]
    _ = run_subprocess_cmd(masurca_config_cmd, False)

    # Paths for MaSuRCA output
    default_scaffolded_assembly_path = os.path.join(input_folder, "CA", "primary.genome.scf.fasta")
    scaffolded_assembly_path = os.path.join(input_folder, "primary.genome.scf.fasta")

    # Check for MaSuRCA output
    if os.path.exists(default_scaffolded_assembly_path):
        log_print("PASS:\tSkipping MaSuRCA, moving output files")
        shutil.move(default_scaffolded_assembly_path, scaffolded_assembly_path)
    elif os.path.exists(scaffolded_assembly_path):
        log_print("PASS:\tSkipping MaSuRCA Assembly; scaffolded assembly already exists")    
    else:
        masurca_assemble_cmd = ["bash",f"{input_folder}/assemble.sh"]
        _ = run_subprocess_cmd(masurca_assemble_cmd, False)

    return scaffolded_assembly_path

def illumina_only_main(INPUT_FOLDER, PERCENT_RESOURCES):
    """
    Main function to run the Illumina-only bioinformatics pipeline.
    
    Args:
        INPUT_FOLDER (str): Path to the input folder with FASTQ files.
        PERCENT_RESOURCES (float): Percentage of resources to use.
    
    Returns:
        str: Path to the completed scaffolded assembly.
    """
    # CPU Threads count setup
    CPU_THREADS, _ = get_resource_values(PERCENT_RESOURCES)
    
    # Generate list of raw input fq files
    input_fq_list = parse_fq_files(INPUT_FOLDER)
    initialize_logging_environment(INPUT_FOLDER)
    
    # Run Trimmomatic on the raw input files
    trimmo_f_pair_path, trimmo_r_pair_path = trimmomatic_prep(input_fq_list, CPU_THREADS)
    
    # Run bbduk on the trimmed files
    bbduk_f_map_path, bbduk_r_map_path = bbduk_map(trimmo_f_pair_path, trimmo_r_pair_path)
    
    # Run clumpify on the mapped files
    clump_f_dedup_path, clump_r_dedup_path = clumpify_dedup(bbduk_f_map_path, bbduk_r_map_path)
    
    # Run Masurca to generate a Scaffolded Assembly
    current_working_dir = os.getcwd()
    os.chdir(INPUT_FOLDER)
    scaffolded_assmebly_path = masurca_config_gen(INPUT_FOLDER, input_fq_list, clump_f_dedup_path, clump_r_dedup_path, CPU_THREADS)
    
    # Run QC Checks on Scaffolded Assembly
    illu_only_qc_checks(scaffolded_assmebly_path)
    
    os.chdir(current_working_dir)
    
    return scaffolded_assmebly_path

if __name__ == "__main__":
    # Set Default Values
    default_input_folder = "/mnt/d/ENTHEOME/Ps_zapatecorum/RECREATION"
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
