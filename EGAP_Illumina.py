# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 12:20:03 2023

@author: ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting with the help of ChatGPT 4.0

EGAP Illumina Pipeline

Command Line Example:
    python EGAP_Illumina.py -i /path/to/illumina/folder -d READS_DATA_STRING
                            -k ORGANISM_KINGDOM_STRING  -es ESTIMATED_GENOME_SIZE_MBP
                            -r PERCENT_RESOURCES_FLOAT  -rf /path/to/reference_sequence
Arguments:
    -i, --input_folder: Path to the input folder containing FASTQ files.
                        (default: /mnt/d/ENTHEOME/Ps_semilanceata/Illumina_PE150)
    -d, --reads_data: Specify the type of sequencing data. Must be one of:
                      'illu' for Illumina, 'ont' for Nanopore, or 'hybrid' for combined data.
                      (default: illu)
    -k, --organism_kingdom: The kingdom to which the current organism belongs.
                            Options: Archaea, Bacteria, Fauna, Flora, Funga, Protista.
                            (default: Funga)
    -es, --est_size: Estimated genome size in Megabase pairs (Mbp). For example, '60' will
                     be converted to 60,000,000 base pairs. (default: 60)
    -r, --resources: Percentage of available system resources (CPU and memory) to use.
                     Must be between 0.01 and 1.00. (default: 0.45)
    -rf, --ref_seq: Optional path to the reference genome for guiding the assembly.
                    (default: None)
"""
# Base Python Imports
import os, argparse, shutil, glob, gzip
from threading import Thread


# Custom Python Imports
from EGAP_QC import assembly_qc_checks, assess_with_fastqc
from EGAP_Tools import log_print, initialize_logging_environment, get_resource_values, cleanup
from EGAP_Masurca import masurca_de_novo, masurca_ref_seq


# Global output_area variable
CPU_THREADS = 1
RAM_GB = 1
DEFAULT_LOG_FILE = None
ENVIRONMENT_TYPE = None


def combine_files(file_list, output_file):
    """
    Combines multiple FASTQ or FASTQ.GZ files into a single file.

    Args:
        file_list (list): List of file paths to be combined.
        output_file (str): Path of the output file.
    """
    log_print(f"Combining Illumina FASTQ files: {file_list} into {output_file}")
    with open(output_file, 'wb') as wfd:
        for f in file_list:
            if f.endswith(".gz"):
                with gzip.open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
            else:
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
                    

def parse_fq_files(input_folder):
    """
    Searches for and returns paths of paired-end FASTQ or FASTQ.GZ files in a given folder and its sub-folders.
    
    Args:
        input_folder (str): Path to the folder containing FASTQ files.

    Returns:
        list of str: Paths of combined "_1.fastq" and "_2.fastq" files.
    """
    log_print('Parsing current folder for separated Illumina reads files...')
    # If multiple pairs exist, combine them into a single file per direction (_1, _2)
    combined_1 = os.path.join(input_folder, f"{input_folder.split('/')[-2]}_1.fastq")
    combined_2 = os.path.join(input_folder, f"{input_folder.split('/')[-2]}_2.fastq")

    log_print(f"Parsing input folder for FASTQ files: {combined_1}, {combined_2}...")

    if os.path.isfile(combined_1) and os.path.isfile(combined_2):
        log_print("NOTE:\tCombined files already exist.")
        return [combined_1, combined_2]
    else:

        # Look for paired-end files with various extensions in the input folder and sub-folders
        pattern_1 = ['*_1.fq', '*_1.fastq', '*_1.fq.gz', '*_1.fastq.gz']
        pattern_2 = ['*_2.fq', '*_2.fastq', '*_2.fq.gz', '*_2.fastq.gz']
    
        # Find files in the given folder and sub-folders
        found_1 = []
        found_2 = []
        for pattern in pattern_1:
            found_1.extend(glob.glob(os.path.join(input_folder, '**', pattern), recursive=True))
        for pattern in pattern_2:
            found_2.extend(glob.glob(os.path.join(input_folder, '**', pattern), recursive=True))
    
        # Ensure we have matching pairs based on naming convention
        paired_1 = {}
        paired_2 = {}
    
        # Populate dictionaries based on base filenames (without _1 or _2)
        for f in found_1:
            base_name = os.path.basename(f).replace("_1.fq", "").replace("_1.fastq", "").replace("_1.fq.gz", "").replace("_1.fastq.gz", "")
            paired_1[base_name] = paired_1.get(base_name, []) + [f]
    
        for f in found_2:
            base_name = os.path.basename(f).replace("_2.fq", "").replace("_2.fastq", "").replace("_2.fq.gz", "").replace("_2.fastq.gz", "")
            paired_2[base_name] = paired_2.get(base_name, []) + [f]
    
        # Find complete pairs
        complete_pairs = []
        for base_name in paired_1:
            if base_name in paired_2:
                complete_pairs.append((paired_1[base_name], paired_2[base_name]))
    
        if not complete_pairs:
            log_print("ERROR: No complete paired-end FASTQ files found.")
            raise FileNotFoundError("Required FASTQ files are missing in the folder or sub-folders.")
            
        # Flatten and combine files
        files_1 = [item for sublist in [pair[0] for pair in complete_pairs] for item in sublist]
        files_2 = [item for sublist in [pair[1] for pair in complete_pairs] for item in sublist]
    
        combine_files(files_1, combined_1)
        combine_files(files_2, combined_2)
    
        return [combined_1, combined_2]


def illumina_prep(ILLU_FOLDER, READS_DATA, CURRENT_ORGANISM_KINGDOM, EST_SIZE, PERCENT_RESOURCES, CPU_THREADS, RAM_GB, REF_SEQ):
    """
    Main pipeline for preparing Illumina data for MaSuRCA hybrid genome assembly.

    Args:
        ILLU_FOLDER (str): Path to the folder containing input data.
        READS_DATA (str): Type of sequencing reads ("illu", "ont", or "hybrid").
        CURRENT_ORGANISM_KINGDOM (str): Kingdom to which the organism belongs (e.g., "Funga").
        EST_SIZE (str): Estimated genome size in Megabase pairs (Mbp) (e.g. "60").
        PERCENT_RESOURCES (float): Percentage of system resources (CPU and memory) to use.
        REF_SEQ (str): Path to the reference genome (optional).

    Returns:
        input_fq_list (list): List of paths to the combined raw Illumina reads in FASTQ format.
        fastqc_output_dirs (list): List of paths to the FastQC folders containing data on the combined raw Illumina reads.
    """
    initialize_logging_environment(ILLU_FOLDER)
    
    # Generate list of raw input fq files
    input_fq_list = parse_fq_files(ILLU_FOLDER)

    # Run QC Checks on Illumina Reads
    r_fastqc_cpus = int(round(CPU_THREADS/2, 0))
    f_fastqc_cpus = CPU_THREADS - r_fastqc_cpus
    
    # Establish FastQC output dirs
    fastqc_output_dirs = [fq_file.replace('.fastq','_fastQC') for fq_file in input_fq_list]
    
    # Quality Control Check A FastQC on Trimmed Illumina Reads 
    f_fastqc_thread = Thread(target = assess_with_fastqc, args = (input_fq_list[0], fastqc_output_dirs[0], f_fastqc_cpus))
    f_fastqc_thread.start()            
    r_fastqc_thread = Thread(target = assess_with_fastqc, args = (input_fq_list[1], fastqc_output_dirs[1], r_fastqc_cpus))
    r_fastqc_thread.start()
    
    # Wait for all QC threads to finish
    f_fastqc_thread.join()
    r_fastqc_thread.join()

    os.chdir(ILLU_FOLDER)

    # # Create a set of absolute paths to keep
    # keep_paths = set(input_fq_list + fastqc_output_dirs + [DEFAULT_LOG_FILE])

    # log_print(f"NOTE:\tDESIRED FILES TO KEEP: {', '.join(keep_paths)}")
    
    # # Cleanup files    
    # cleanup(keep_paths, ILLU_FOLDER, DEFAULT_LOG_FILE)

    log_print("PASS:\tProcessed Illumina Raw Data; ready for next steps.\n")

    return input_fq_list, fastqc_output_dirs


def illumina_only_main(ILLU_FOLDER, READS_DATA, CURRENT_ORGANISM_KINGDOM, EST_SIZE, PERCENT_RESOURCES, REF_SEQ):
    """
    Main pipeline for performing Illumina-only genome assembly, produces de novo and referenced assemblys (if REF_SEQ is not None).

    Args:
        ILLU_FOLDER (str): Path to the folder containing input data.
        READS_DATA (str): Type of sequencing reads ("illu", "ont", or "hybrid").
        CURRENT_ORGANISM_KINGDOM (str): Kingdom to which the organism belongs (e.g., "Funga").
        EST_SIZE (str): Estimated genome size in Megabase pairs (Mbp) (e.g. "60").
        PERCENT_RESOURCES (float): Percentage of system resources (CPU and memory) to use.
        REF_SEQ (str): Path to the reference genome (optional).

    Returns:
        ref_seq_assembly_path (str): Path to the MaSuRCA referenced assembly FASTA file.
        de_novo_assembly_path (str): Path to the MaSuRCA de novo assembly assembly FASTA file.
    """
    if READS_DATA != "illu":
        log_print(f"ERROR:\tThe input data should be 'illu' and is '{READS_DATA}'. Please check input.")
        return None
    
    # CPU Threads count setup
    CPU_THREADS, RAM_GB = get_resource_values(PERCENT_RESOURCES)

    input_fq_list, fastqc_output_dirs = illumina_prep(ILLU_FOLDER, READS_DATA, CURRENT_ORGANISM_KINGDOM, PERCENT_RESOURCES, EST_SIZE, REF_SEQ)

    # Create unique folders for De-Novo and Ref-Seq Assemblies
    ref_seq_folder = os.path.join(ILLU_FOLDER, "Ref-Seq-Assembly")
    de_novo_folder = os.path.join(ILLU_FOLDER, "De-Novo-Assembly")

    # If reference sequence exists, generate a referenced assembly
    ref_assembly_path = None
    if REF_SEQ:
        log_print('Generating MaSuRCA referenced assembly...')
        os.makedirs(ref_seq_folder, exist_ok=True)
        ref_assembly_path = masurca_ref_seq(ILLU_FOLDER, ref_seq_folder, input_fq_list, EST_SIZE, CPU_THREADS, RAM_GB, REF_SEQ)
        
        # Run QC Checks on Reference-Guided Assembly
        ref_seq_quast_output_dir, ref_seq_compleasm_output_dir_1, ref_seq_compleasm_output_dir_2 = assembly_qc_checks(ref_assembly_path, CURRENT_ORGANISM_KINGDOM, READS_DATA, CPU_THREADS)

    # Generate a de novo assembly
    log_print('Generating MaSuRCA de novo assembly...')
    os.makedirs(de_novo_folder, exist_ok=True)
    de_novo_assembly_path = masurca_de_novo(ILLU_FOLDER, de_novo_folder, input_fq_list, EST_SIZE, CPU_THREADS, RAM_GB, REF_SEQ)

    # Run QC Checks on Scaffolded Assembly # UPDATE TO USE _clean.fasta when fully implemented
    de_novo_quast_output_dir, de_novo_compleasm_output_dir_1, de_novo_compleasm_output_dir_2 = assembly_qc_checks(de_novo_assembly_path, READS_DATA, CURRENT_ORGANISM_KINGDOM, REF_SEQ, CPU_THREADS)
 
    os.chdir(ILLU_FOLDER)

    # Create a set of absolute paths to keep
    keep_paths = [ILLU_FOLDER, DEFAULT_LOG_FILE, de_novo_assembly_path,
                  de_novo_quast_output_dir,
                  de_novo_compleasm_output_dir_1, de_novo_compleasm_output_dir_2]
    if REF_SEQ:
        keep_paths.append(ref_assembly_path, ref_seq_quast_output_dir,
                          ref_seq_compleasm_output_dir_1, ref_seq_compleasm_output_dir_2)
    
    # # add to keep_paths if REF_SEQ != None
    # if REF_SEQ != None:
    #     keep_paths = set(keep_paths, [ref_assembly_path, ref_seq_quast_output_dir, ref_seq_compleasm_output_dir_1, ref_seq_compleasm_output_dir_2])

    # Cleanup files    
    # cleanup(keep_paths, ILLU_FOLDER, DEFAULT_LOG_FILE) # Gets hung up removing CA folder specifically -7

    log_print("PASS:\Illumina reads Sucessfully Assembled; ready for next steps.\n")

    return keep_paths


## Debuging Main Space & Example
if __name__ == "__main__":
    print('EGAP Illumina Pipeline')       

    # Create the parser
    parser = argparse.ArgumentParser(description="Run EGAP Illumina Pipelie")

    # Set Default Values
    default_input_folder = "/mnt/d/ENTHEOME/Ps_semilanceata/Illumina_PE150" # "/mnt/d/Tryptomics/Genomics/AM_Ps_cubensis/Tallebudgera-Valley-20220505" #"/mnt/d/Tryptomics/Genomics/MG_Ps_cubensis_GT/Illumina_MiSeq_PE" # "/mnt/d/Tryptomics/Example/"
    default_reads_data = "illu"
    default_organism_kingdom = "Funga"
    default_estimated_genome_size = 60
    default_percent_resources = 0.45
    default_reference_sequence = None # Psilocybe cubensis reference sequence: "/mnt/d/Tryptomics/Genomics/MG_Ps_cubensis_PE_REF_SEQ/GCA_017499595.gbff" 

    # Add arguments
    parser.add_argument("-i", "--input_folder", type=str, default=default_input_folder,
                        help="Path to the input folder containing FASTQ files")
    parser.add_argument('--reads_data', '-d',
                        type = str, default = default_reads_data,
                        choices = ["illu","ont","hybrid"],
                        help = f'Indicate if the provided data are generated from the same organism or different organisms (default: {default_reads_data})')
    parser.add_argument('--organism_kingdom', '-k',
                        type=str, default=default_organism_kingdom,
                        help = f'Kingdom the current organism data belongs to. (default: {default_organism_kingdom})')
    parser.add_argument("--est_size", "-es",
                        type=int, default=default_estimated_genome_size,
                        help="Estimaged size of the genome in Mbp (aka million-base-pairs). (default: {default_estimated_genome_size}Mbp => {int(default_estimated_genome_size)*1000000})")
    parser.add_argument("-r", "--resources",
                        type=float, default=default_percent_resources,
                        help=f"Percentage of resources to use. (0.01-1.00; default: {default_percent_resources})")
    parser.add_argument("--ref_seq", "-rf",
                        type=str, default=default_reference_sequence,
                        help="Path to the reference genome for assembly. (default: {default_reference_sequence})")

    # Parse arguments
    args = parser.parse_args()
    INPUT_FOLDER = args.input_folder
    READS_DATA = args.reads_data
    CURRENT_ORGANISM_KINGDOM = args.organism_kingdom
    EST_SIZE = args.est_size * 1000000
    PERCENT_RESOURCES = args.resources
    REF_SEQ = args.ref_seq
        
    # Hybrid Prep Main function
    input_fq_list, fastqc_output_dirs = illumina_prep(INPUT_FOLDER, READS_DATA,
                                                      CURRENT_ORGANISM_KINGDOM,
                                                      EST_SIZE, PERCENT_RESOURCES,
                                                      REF_SEQ)

    # Illumina Only Main function
    illu_final_paths_list = illumina_only_main(INPUT_FOLDER,
                                               READS_DATA,
                                               CURRENT_ORGANISM_KINGDOM,
                                               EST_SIZE,
                                               PERCENT_RESOURCES, 
                                               REF_SEQ)