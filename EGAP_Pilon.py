# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 13:49:51 2023

@author: ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting with the help of ChatGPT 4.0



"""
# Base Python Imports
import os, argparse, subprocess


# Custom Python Imports
from EGAP_Tools import log_print, initialize_logging_environment, run_subprocess_cmd, get_resource_values, find_file, generate_sam, convert_sam_to_bam


# Global output_area variable
CPU_THREADS = 1
RAM_GB = 1
DEFAULT_LOG_FILE = None
ENVIRONMENT_TYPE = None


# Function to trim Illumina FASTQ files and generate Forward Paired and Reverse Paired FASTQ Reads as a list
def trim_with_trimmomatic(folder_name, combined_files, option1, option2):
    """
    Trim Illumina FASTQ files using Trimmomatic.

    Args:
        folder_name (str): Path to the folder containing Illumina data.
        combined_files (list): List of paths to the files for raw Illumina Reads for trimming.
        data_type (str): Either 'PE' for Paired-End reads or 'SE' Single-End reads.
        option1 (str): Additional Trimmomatic option.
        option2 (str): Additional Trimmomatic option.

    Returns:
        fq_paired_list (list): List containing paths to Forward Paired and Reverse Paired FASTQ files.
        fastqc_output_dirs (list): List containing paths to FastQC output directories for the Forward Paired and Reverse Paired FASTQ files.
    """    
    # Check for trimmomatic jar file
    log_print(f"Running Trimmomatic on {combined_files}...")

    # Default path to Trimmomatic paths
    default_trimmo_path = "trimmomatic-0.39.jar"
    default_trimmo_adapters_path = "TruSeq3-PE.fa"
    trimmo_path = find_file(default_trimmo_path)
    trimmo_adapters_path = find_file(default_trimmo_adapters_path)   

    # Generate directories and lists
    fq_paired_list = []
    fastqc_output_dirs = []    
    illumina_dir = '/'.join(combined_files[0].split('/')[:-1])
        
    # Check for input files
    fwd_file = combined_files[0]
    rev_file = combined_files[1]
    
    # Set up output files
    fwd_paired_out = os.path.join(illumina_dir, os.path.basename(fwd_file).replace('_combined_1.fq', '_trimmomatic_forward_paired.fq'))
    fwd_unpaired_out = os.path.join(illumina_dir, os.path.basename(fwd_file).replace('_combined_1.fq', '_trimmomatic_forward_unpaired.fq'))
    rev_paired_out = os.path.join(illumina_dir, os.path.basename(rev_file).replace('_combined_2.fq', '_trimmomatic_reverse_paired.fq'))
    rev_unpaired_out = os.path.join(illumina_dir, os.path.basename(rev_file).replace('_combined_2.fq', '_trimmomatic_reverse_unpaired.fq'))
    
    # Check if the final combined file exists
    if os.path.isfile(fwd_paired_out) and os.path.isfile(rev_paired_out):
        log_print(f"PASS:\tSkipping Trimmomtaic trimming: Trimmed files already exist {fwd_paired_out}, {rev_paired_out}")
        fq_paired_list.append(fwd_paired_out)
        fq_paired_list.append(rev_paired_out)     
    else:
        
        # Build command for Trimmomatic
        trimmomatic_cmd = ["java", "-jar", trimmo_path, "PE", "-phred33", 
                           "-threads", str(CPU_THREADS),
                           fwd_file, rev_file, 
                           fwd_paired_out, fwd_unpaired_out, 
                           rev_paired_out, rev_unpaired_out, 
                           f"ILLUMINACLIP:{trimmo_adapters_path}:2:30:10:11", 
                           option1, option2, "SLIDINGWINDOW:50:25", "MINLEN:125"]
        
        # Execute command
        log_print(f"CMD:\t{' '.join(trimmomatic_cmd)}")
        trimmomatic_result = subprocess.Popen(trimmomatic_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = trimmomatic_result.communicate()
        return_code = trimmomatic_result.returncode
        if return_code != 0:
            log_print(f"ERROR:\tTrimmomatic failed for:\t{fwd_file}; {rev_file}; {stderr.decode('utf-8')}")
        else:
            log_print("PASS:\tGenerated all Trimmomatic output files successfully")

    # Return the paths to the trimmed files
    fq_paired_list = [fwd_paired_out, rev_paired_out]
    fastqc_output_dirs = [os.path.join(folder_name, 'fastqc_forward_paired'), os.path.join(folder_name, 'fastqc_reverse_paired')]
    
    return fq_paired_list, fastqc_output_dirs


def pilon_polish_assembly(assembly_file, bam_file, CPU_THREADS, RAM_GB):    
    """
    Uses Pilon to polish the Cleaned ONT Flye Assembly with Illiumna Binary Alignment Map
    
    Args:
        assembly_file (str): Path to the Flye Assembly FASTA file.
        bam_file (str): A Binary Alignment Map (BAM file) from the Illumina Sequence Alignment Map with BWA. 
        CPU_THREADS (int): Number of threads available.
        RAM_GB (int): Number of Gigabytes of available RAM.

    Returns:
        pilon_fasta (str): Path to the final polished Pilon FASTA file output.
    """
    log_print(f"Pilon Polishing {assembly_file} with {bam_file}...")

    pilon_output = assembly_file.replace(".fasta","_polished")
    pilon_fasta = pilon_output + ".fasta"
    
    # Check if the output file already exists
    if os.path.isfile(pilon_fasta):
        log_print(f"PASS:\tSkipping Pilon Polish, {pilon_fasta} already exists")
        return pilon_fasta
    else:
        os.makedirs(pilon_output, exist_ok=True)  # add exist_ok=True


    # Default path to Trimmomatic paths
    default_pilon_path = "pilon-1.24.jar"
    pilon_path = find_file(default_pilon_path)

    # If the file doesn't exist, run Pilon
    pilon_cmd = ["java", f"-Xmx{RAM_GB}G",
                  "-jar", pilon_path,
                  "--genome", assembly_file,
                  "--frags", bam_file,
                  "--output", pilon_output,
                  "--changes", "--vcf",
                  "--threads", str(CPU_THREADS)]
    
    _ = run_subprocess_cmd(pilon_cmd, False)   

    return pilon_fasta


def final_hybrid_pilon(ASSEMBLY_FASTA, illumina_reads_list, READS_DATA, CURRENT_ORGANISM_KINGDOM, PERCENT_RESOURCES, REF_SEQ):
    """
    Main function to call for the Pilon polishing of Cleaned ONT Reads with Trimmed Illumina Reads.
    
    Args:
        ASSEMBLY_FASTA (str): Path to the folder containing a sub-folder with the ONT reads.
        illumina_reads_list (list): 
        CURRENT_ORGANISM_KINGDOM (str): Single word Kingdom Description: Archaea, Bacteria, Fauna, Flora, Funga, or Protista.
        PERCENT_RESOURCES (int): Amount of resources availalbe for use during processing.
    
    Returns:
        pilon_output (str): Path to the final polished Pilon FASTA file output.
    """
    log_print(f"Running Pilon Polish of Hybrid Assembly {ASSEMBLY_FASTA}; with trimmed Illumina Reads: {illumina_reads_list}")
    
    # Generate log file with the desired behavior
    base_folder = '/'.join(ASSEMBLY_FASTA.split('/')[:-2])
    initialize_logging_environment(base_folder)
    
    CPU_THREADS, RAM_GB = get_resource_values(PERCENT_RESOURCES)
    
    # Process Illumina Reads with Trimmomatic
    fq_paired_list, fastqc_output_dirs = trim_with_trimmomatic(base_folder, illumina_reads_list, 'HEADCROP:10', 'CROP:145')
    
    illu_bam_file = ASSEMBLY_FASTA.replace(".fasta", "_bwa_aligned_sorted.bam")
    
    if os.path.isfile(illu_bam_file):
        bam_file = illu_bam_file
    else:
        # Generate the SAM file based on the indexed ont assembly and trimmed paired illumina reads
        sam_file = generate_sam(ASSEMBLY_FASTA, fq_paired_list, CPU_THREADS)
        
        # Generate an indexed BAM file from the SAM file
        bam_file = convert_sam_to_bam(sam_file, CPU_THREADS)
    
    # Perfrom final Pilon Polish on the cleaned ont assembly with the generated BAM file
    pilon_polished_assembly = pilon_polish_assembly(ASSEMBLY_FASTA, bam_file, CPU_THREADS, RAM_GB)

    return pilon_polished_assembly


## Debuging Main Space & Example
if __name__ == "__main__":
    print('EGAP Pilon Pilosher')
    
    # Argument Parsing
    parser = argparse.ArgumentParser(description='Run EGAP Pilon Pipeline')
   
    # Default values
    default_assembly = '/mnt/d/Entheome/Ps_semilanceata/ENTHEOME_Ps_semilanceata_de_novo.fasta'
    default_illumina_f = '/mnt/d/Entheome/Ps_semilanceata/Illumina_PE150/Ps_semilanceata__1.fq'
    default_illumina_r = '/mnt/d/Entheome/Ps_semilanceata/Illumina_PE150/Ps_semilanceata__2.fq'
    default_reads_data = 'hybrid'
    default_organism_kingdom = 'Funga'
    default_estimated_genome_size = 60
    default_percent_resources = 0.45
    default_reference_sequence= None
    
    # Add arguments with default values
    parser.add_argument('--input_assembly', '-ia',
                        default = default_assembly,
                        help = f'Path to the Assembly fasta file. (default: {default_assembly})')
    parser.add_argument('--illumina_f_input', '-if',
                        default = default_illumina_f,
                        help = f'Path to the Trimmed Forward Illumina Reads. (default: {default_illumina_f})')
    parser.add_argument('--illumina_r_input', '-ir',
                        default = default_illumina_r,
                        help = f'Path to the Trimmed Reverse Illumina Reads. (default: {default_illumina_r})')
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
    ASSEMBLY_FASTA = args.input_assembly
    ILLUM_F_READS = args.illumina_f_input
    ILLUM_R_READS = args.illumina_r_input
    READS_DATA = args.reads_data
    CURRENT_ORGANISM_KINGDOM = args.organism_kingdom
    EST_SIZE = args.est_size * 1000000
    PERCENT_RESOURCES = args.resources
    REF_SEQ = args.ref_seq
        
    # Hybrid Pilon Polish
    pilon_polished_assembly = final_hybrid_pilon(ASSEMBLY_FASTA, [ILLUM_F_READS, ILLUM_R_READS],
                                                 READS_DATA, CURRENT_ORGANISM_KINGDOM,
                                                 PERCENT_RESOURCES, REF_SEQ)