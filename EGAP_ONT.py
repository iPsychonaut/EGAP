# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 11:49:17 2023

@author: ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting with the help of ChatGPT 4.0

EGAP Oxford Nanopore Technologies (ONT) Pipeline

Command Line Example:
    python EGAP_ONT.py -i /path/to/ont/folder -d READS_DATA_STRING
                       -k ORGANISM_KINGDOM_STRING  -es ESTIMATED_GENOME_SIZE_MBP
                       -r PERCENT_RESOURCES_FLOAT  -rf /path/to/reference_sequence
Arguments:
    -i, --ont_folder: Path to the ONT Folder.
                      (default: /mnt/d/Entheome/Ps_semilanceata/ONT_MinION/)
    -d, --reads_data: Specify the type of sequencing data. Must be one of:
                      'illu' for Illumina, 'ont' for Nanopore, or 'hybrid' for combined data.
                      (default: ont)
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
import os, glob, gzip, argparse


# Required Python Imports
from tqdm import tqdm
from Bio import SeqIO


# Custom Python Imports
from EGAP_Tools import log_print, initialize_logging_environment, get_resource_values, cleanup, run_subprocess_cmd
from EGAP_QC import assess_with_nanoplot, assembly_qc_checks
from EGAP_Flye import assemble_ont_flye


# Global output_area variable
CPU_THREADS = 1
RAM_GB = 1
DEFAULT_LOG_FILE = None
ENVIRONMENT_TYPE = None


# Function to extract and combine multiple FASTQ.GZ files into a single FASTQ file
def ont_combine_fastq_gz(ONT_FOLDER):
    """
    Combine multiple ONT FASTQ.GZ files into a single FASTQ file.

    Args:
        ONT_FOLDER (str): Path to the folder containing ONT FASTQ.GZ files.

    Returns:
        combined_ont_fastq_path (str): Path to the combined FASTQ file.
    """
    log_print('Combining ONT FASTQ.GZ files...')

    # Search for the folder containing multiple .gz files
    ont_raw_data_dir = next((subdir for subdir in glob.glob(os.path.join(ONT_FOLDER, "*"))
                            if os.path.isdir(subdir) and any(file.endswith(".gz") for file in os.listdir(subdir))), None)
    if ont_raw_data_dir is None:
        log_print(f"NOTE: No directory containing '.gz' files found within '{ONT_FOLDER}'")
        log_print(f"Attempting with ONT_FOLDER '{ONT_FOLDER}'")
        ont_raw_data_dir = ONT_FOLDER

    # Get the base name from the path and append the '.fastq' extension
    base_name = ONT_FOLDER.split("/")[-3]
    if base_name == "ENTHEOME":
        base_name  = ONT_FOLDER.split("/")[-2]
    combined_ont_fastq_path = os.path.join(ONT_FOLDER, f'{base_name}_ont_combined.fastq')
    return_raw_data_dir = os.path.join(ONT_FOLDER, base_name)
    raw_file_list = glob.glob(os.path.join(ont_raw_data_dir, '*.fastq.gz'))
    
    # Check if the combined file already exists
    print(combined_ont_fastq_path)
    if os.path.isfile(combined_ont_fastq_path):
        log_print(f"NOTE: Skipping extraction & combination: Combined fastq file: {combined_ont_fastq_path} already exists")
        return return_raw_data_dir, raw_file_list, combined_ont_fastq_path

    # Combine the FASTQ files
    with open(combined_ont_fastq_path, 'w') as combined_file:
        for filename in tqdm(raw_file_list, desc='Combining files'):
            with gzip.open(filename, 'rt') as gz_file:
                for record in SeqIO.parse(gz_file, "fastq"):
                    try:
                        SeqIO.write(record, combined_file, "fastq")
                    except Exception as e:
                        log_print(f"ERROR: in FASTQ record: {e}")
                        raise e

    log_print(f"PASS: Successfully created combined fastq file: {combined_ont_fastq_path}")

    return return_raw_data_dir, raw_file_list, combined_ont_fastq_path


def run_porechop(combined_ont_fastq, length_threshold, quality_threshold, CPU_THREADS):
    log_print('Removing Adapter sequences from ONT reads with Porechop...')
    
    chop_filtered_out = combined_ont_fastq.replace(".fq",".fastq").replace(".fastq","_chopped_filtered.fastq")
    
    if os.path.isfile(chop_filtered_out):
        log_print(f"NOTE: Skipping Porechop: Chopped & Filtered fastq file {chop_filtered_out} already exists")
    else:
        porechop_cmd = ["porechop", "-i", combined_ont_fastq,
                        "-o", chop_filtered_out,
                        "-v", "2",
                        "-t", str(CPU_THREADS)]
        
        _ = run_subprocess_cmd(porechop_cmd, False)
        
        log_print('PASS:\tAdapter sequences removed from ONT reads.')
    
        run_nanofilt(chop_filtered_out, length_threshold, quality_threshold, CPU_THREADS)
    
    return chop_filtered_out


def run_nanofilt(chopped_ont_fastq, length_threshold, quality_threshold, CPU_THREADS):
    log_print('Filtering out Low-Quality ONT reads with NanoFilt...')
    
    nanofilt_cmd = ["NanoFilt",
                    "-l", str(length_threshold),
                    "-q", str(quality_threshold),
                    chopped_ont_fastq]
    
    _ = run_subprocess_cmd(nanofilt_cmd, False)
    
    log_print('PASS:\tLow-Quality reads filtered from ONT reads.')
    

def ont_prep(ONT_FOLDER, READS_DATA, CURRENT_ORGANISM_KINGDOM, EST_SIZE, PERCENT_RESOURCES, REF_SEQ):
    """
    Main pipeline for preparing ONT data for MaSuRCA hybrid genome assembly.
    
    Args:
        ONT_FOLDER (str): Path to the folder containing input ONT data.
        READS_DATA (str): Type of sequencing reads ("illu", "ont", or "hybrid").
        CURRENT_ORGANISM_KINGDOM (str): Kingdom to which the organism belongs (e.g., "Funga").
        EST_SIZE (str): Estimated genome size in Megabase pairs (Mbp) (e.g. "60").
        PERCENT_RESOURCES (float): Percentage of system resources (CPU and memory) to use.
        REF_SEQ (str): Path to the reference genome (optional).
        
    Returns:
        trimmed_ont_fastq (str): Path to the trimmed ONT FASTQ file.
    """
    initialize_logging_environment(ONT_FOLDER)
    
    CPU_THREADS, ram_gb = get_resource_values(PERCENT_RESOURCES)
    
    # Take ONT Folder with FASTQ Files, extract, combine
    ont_raw_data_dir, raw_file_list, combined_ont_fastq = ont_combine_fastq_gz(ONT_FOLDER)

    # Quality Control Check of Combined ONT Reads with NanoPlot
    raw_nanoplot_dir, length_threshold, quality_threshold = assess_with_nanoplot(combined_ont_fastq, CPU_THREADS)
    
    # RUN Porechop & NanoFilt ON combined_ont_fastq
    filtered_ont_fastq = run_porechop(combined_ont_fastq, length_threshold, quality_threshold, CPU_THREADS)
    
    # Quality Control Check of Combined ONT Reads with NanoPlot
    filtered_nanoplot_dir, _, _ = assess_with_nanoplot(filtered_ont_fastq, CPU_THREADS)

    log_print("PASS:\tProcessed ONT Raw Data; ready for next steps.\n")

    return filtered_ont_fastq


def ont_only_main(ONT_FOLDER, READS_DATA, CURRENT_ORGANISM_KINGDOM, EST_SIZE, PERCENT_RESOURCES, REF_SEQ):
    """
    Main pipeline for performing ONT-only genome assembly, produces de novo and referenced assemblys (if REF_SEQ is not None).
    
    Args:
        ONT_FOLDER (str): Path to the folder containing input ONT data.
        READS_DATA (str): Type of sequencing reads ("illu", "ont", or "hybrid").
        CURRENT_ORGANISM_KINGDOM (str): Kingdom to which the organism belongs (e.g., "Funga").
        EST_SIZE (str): Estimated genome size in Megabase pairs (Mbp) (e.g. "60").
        PERCENT_RESOURCES (float): Percentage of system resources (CPU and memory) to use.
        REF_SEQ (str): Path to the reference genome (optional).
        
    Returns:
        ref_seq_assembly_path (str): Path to the MaSuRCA referenced assembly FASTA file.
        de_novo_assembly_path (str): Path to the MaSuRCA de novo assembly assembly FASTA file.
    """
    initialize_logging_environment(ONT_FOLDER)

    if READS_DATA != "ont":
        log_print(f"ERROR:\tThe input data should be 'ont' and is '{READS_DATA}'. Please check input.")
        return None
    
    trimmed_ont_fastq = ont_prep(ONT_FOLDER, READS_DATA, CURRENT_ORGANISM_KINGDOM, EST_SIZE, PERCENT_RESOURCES, REF_SEQ)
    
    # Create unique folders for De-Novo and Ref-Seq Assemblies
    de_novo_folder = os.path.join(ONT_FOLDER, "De-Novo-Assembly")
    
    os.makedirs(de_novo_folder, exist_ok=True)
    ont_flye_assembly, ont_flye_indexed_bam, flye_dir = assemble_ont_flye(trimmed_ont_fastq, READS_DATA, CURRENT_ORGANISM_KINGDOM, EST_SIZE, CPU_THREADS, REF_SEQ)
    
    # Quality Control Check of the ONT de novo Assembly with QUAST & 2x Compleasm BUSCO # UPDATE TO USE _clean.fasta when fully implemented
    quast_output_dir, compleasm_output_dir_1, compleasm_output_dir_2 = assembly_qc_checks(ont_flye_assembly, READS_DATA, CURRENT_ORGANISM_KINGDOM, REF_SEQ, CPU_THREADS)

    
    # Create a set of absolute paths to keep
    keep_paths = [ONT_FOLDER, DEFAULT_LOG_FILE, ont_flye_assembly,
                  quast_output_dir,
                  compleasm_output_dir_1, compleasm_output_dir_2]
    
    # # add to keep_paths if REF_SEQ != None
    # if REF_SEQ != None:
    #     keep_paths = set(keep_paths, [ref_assembly_path, ref_seq_quast_output_dir, ref_seq_compleasm_output_dir_1, ref_seq_compleasm_output_dir_2])

    # Cleanup files    
    # cleanup(keep_paths, ILLU_FOLDER, DEFAULT_LOG_FILE) # Gets hung up removing CA folder specifically -7

    log_print("PASS:\ONT reads Sucessfully Assembled; ready for next steps.\n")

    return keep_paths


## Debuging Main Space & Example
if __name__ == "__main__":
    print('EGAP Oxford Nanopore Technologies (ONT) Pipeline')       
    
    # Argument Parsing
    parser = argparse.ArgumentParser(description='RUN EGAP ONT Pipeline')
    
    # Default values
    default_folder = '/mnt/d/Entheome/Ps_semilanceata/ONT_MinION/'
    default_reads_data = "ont"
    default_organism_kingdom = 'Funga'
    default_estimated_genome_size = 60
    default_percent_resources = 0.45
    default_reference_sequence = None
    
    # Add arguments with default values
    parser.add_argument('--ont_folder', '-i',
                        type = str, default = default_folder,
                        help = f'Path to the ONT Folder. (default: {default_folder})')
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
    
    # Parse the arguments
    args = parser.parse_args()
    ONT_FOLDER = args.ont_folder
    READS_DATA = args.reads_data
    CURRENT_ORGANISM_KINGDOM = args.organism_kingdom
    EST_SIZE = args.est_size
    PERCENT_RESOURCES = args.resource_use
    REF_SEQ = args.ref_seq

    # Run main ONT Cleaning function
    combined_ont_fastq = ont_prep(ONT_FOLDER, READS_DATA,
                                  CURRENT_ORGANISM_KINGDOM, EST_SIZE,
                                  PERCENT_RESOURCES, REF_SEQ)

    
    ont_final_paths_list = ont_only_main(ONT_FOLDER, READS_DATA,
                                         CURRENT_ORGANISM_KINGDOM,
                                         EST_SIZE, PERCENT_RESOURCES,
                                         REF_SEQ)
