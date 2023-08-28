# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 13:49:51 2023

@author: ian.michael.bollinger@gmail.com with the help of ChatGPT 4.0

Command Line Example:
    python EGAP_pilon_polish.py --ont_input /path/to/cleaned_ont_reads.fastq --illumina_f_input /path/to/forward_reads.fq --illumina_r_input /path/to/reverse_reads.fq --organism_kingdom STRING --resource_use INTEGER

The --organism_kingdom must be from the following: Archaea, Bacteria, Fauna, Flora, Funga, or Protista
"""
import os, subprocess, argparse, multiprocessing, psutil, math, shutil
from check_tools import get_env_dir, check_for_jars
from EGAP_qc import assess_with_quast, assess_with_busco
from log_print import log_print, generate_log_file
from threading import Thread

# Function to Index Cleaned ONT Flye Assembly with BWA
def index_cleaned_ont_assembly(input_assembly, log_file):
    """
    Index Cleaned ONT Flye Assembly with BWA.

    Args:
        input_assembly (str): Path to an cleaned FASTA assembly.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.
    """
    log_print(f'Indexing {input_assembly} with BWA-MEM Algorithm...', log_file)
    # Check if the index files already exist
    index_files = [f"{input_assembly}.{ext}" for ext in ["amb", "ann", "bwt", "pac", "sa"]]
    if all(os.path.exists(path) for path in index_files):
        log_print(f"PASS:\tSkipping Indexing, Index files for {input_assembly} already exist", log_file)
        return
    
    # Run BWA-MEM Index
    bwa_command = ["bwa",
                   "index",
                   input_assembly]
    # Run the command
    log_print(f"CMD:\t{' '.join(bwa_command)}", log_file)
    bwa_result = subprocess.run(bwa_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)    
    if bwa_result.returncode != 0:
        log_print(f"ERROR:\t{bwa_result.stderr}", log_file)
        return
    else:
        log_print("PASS:\tSuccessfully indexed Cleaned ONT Flye Assembly", log_file)

# Create a SAM Map of the Decontaminated Flye with the Trimmed Illumina Paired Forward & Reverse Reads
def generate_illumina_sam(cleaned_ont_assembly, forward_reads, reverse_reads, cpu_threads, log_file):
    """
    Generates a SAM Map of the ONT Cleaned Flye de novo Assembly with the Trimmed Illumina Paired Forward & Reverse Reads.

    Args:
        cleaned_ont_assembly (str): Path to a cleaned ONT Assembly FASTQ file.
        forward_reads (str): Path to the Trimmed Paired Forward ILlumina FASTQ file.
        reverse_reads (str): Path to the Trimmed Paired Reverse ILlumina FASTQ file.
        cpu_threads (int): Number of threads available for processing.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.

    Returns:
        output_sam (str): A SAM Map File of the ONT Cleaned Flye de novo Assembly with the Trimmed Illumina Paired Forward & Reverse Reads.
    """
    log_print(f'Generating SAM from {cleaned_ont_assembly}, {forward_reads}, and {reverse_reads} with BWA-MEM Algorithm...', log_file)
        
    # Check if the SAM file already exists
    output_sam = cleaned_ont_assembly.replace("_filtered.fasta", "_bwa_aligned.sam")
    if os.path.isfile(output_sam):
        log_print(f"PASS:\tSkipping SAM generation, SAM file {output_sam} already exists", log_file)
        return output_sam
    
    # Open the output file
    with open(output_sam, 'w') as sam_output_file:
        # Run the BWA-MEM Mapping
        sam_cmd = ["bwa", "mem",
                   "-t", str(cpu_threads),
                   cleaned_ont_assembly,
                   forward_reads, reverse_reads]
        log_print(f"CMD:\t{' '.join(sam_cmd)}", log_file)
        
        # Run the BWA-MEM Mapping command
        sam_result = subprocess.run(sam_cmd, stdout=sam_output_file, stderr=subprocess.PIPE)
        
        # Check if BWA-MEM ran successfully
        if sam_result.returncode != 0:
            log_print(f"ERROR:\t{sam_result.stderr.decode()}", log_file)
            return None
        else:
            log_print("PASS:\tSuccessfully generated SAM file", log_file)

    return output_sam

# Generate Binary Alignment Map from the Illumina Sequence Alignment Map with BWA
def convert_sam_to_bam(sam_file, cpu_threads, log_file):
    """
    Generates Binary Alignment Map (BAM file) from the Illumina Sequence Alignment Map with BWA.

    Args:
        sam_file (str): A SAM Map File of the ONT Cleaned Flye de novo Assembly with the Trimmed Illumina Paired Forward & Reverse Reads.
        cpu_threads (int): Number of threads available for processing.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.

    Returns:
        output_bam (str): A Binary Alignment Map (BAM file) from the Illumina Sequence Alignment Map with BWA. 
    """
    log_print(f'Converting {sam_file} into BAM file with Samtools...', log_file)

    # Check if the BAM file already exists
    output_bam = sam_file.replace(".sam", ".bam")
    output_bai = output_bam + ".bai"  # Index file
    if os.path.isfile(output_bam) and os.path.isfile(output_bai):
        log_print(f"PASS:\tSkipping BAM Generation, BAM file {output_bam} and index {output_bai} already exist", log_file)
        return output_bam

    # First command: samtools view
    view_cmd = ["samtools",
                "view",
                "-S",
                "-b", sam_file]
    log_print(f"CMD:\t{' '.join(view_cmd)}", log_file)
    view_process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE)

    # Second command: samtools sort
    sort_cmd = ["samtools",
                "sort",
                "-@", str(cpu_threads),
                "-o", output_bam]
    log_print(f"CMD:\t{' '.join(sort_cmd)}", log_file)
    sort_process = subprocess.Popen(sort_cmd, stdin=view_process.stdout, stdout=subprocess.PIPE)
    
    # Close the stdout of view_process to allow view_process to receive a SIGPIPE if sort_process exits.
    view_process.stdout.close()

    # Wait for sort_process to finish
    sort_output, sort_error = sort_process.communicate()

    # Check if Sam to Bam ran successfully
    if sort_process.returncode != 0:
        log_print(f"ERROR:\t{sort_error.decode('utf-8')}", log_file)
        return None

    # Third command: samtools index
    index_cmd = ["samtools",
                 "index",
                 "-@", str(cpu_threads),
                 output_bam]
    log_print(f"CMD:\t{' '.join(index_cmd)}", log_file)
    index_process = subprocess.Popen(index_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    index_output, index_error = index_process.communicate()

    # Check if indexing was successful
    if index_process.returncode != 0:
        log_print(f"ERROR:\t{index_error.decode('utf-8')}", log_file)
        return None

    log_print("PASS:\tSuccessfully generated BAM file", log_file)
    return output_bam

# Polish the Cleaned ONT Flye Assembly with Illiumna Binary Alignment Map
def pilon_polish_assembly(assembly_file, bam_file, ram_gb, cpu_threads, log_file):    
    """
    Uses Pilon to polish the Cleaned ONT Flye Assembly with Illiumna Binary Alignment Map
    
    Args:
        assembly_file (str): Path to the Flye Assembly FASTA file.
        bam_file (str): A Binary Alignment Map (BAM file) from the Illumina Sequence Alignment Map with BWA. 
        ram_gb (int): Number of Gigabytes of available RAM.
        cpu_threads (int): Number of threads available for processing.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.

    Returns:
        pilon_fasta (str): Path to the final polished Pilon FASTA file output.
    """
    log_print(f'Pilon Polishing {assembly_file} with {bam_file}...', log_file)
    # Create the output directory if it does not exist
    ont_remove_list = assembly_file.split('/')

    # Filter the list for items containing 'ont'
    filtered_items = [item for item in ont_remove_list if 'ONT_MinION' in item]

    # Return the only item, if it exists
    ont_remove = filtered_items[0] if len(filtered_items) == 1 else None

    # Build out pilon output folder and file
    pilon_output = assembly_file.replace('_ont_flye_filtered.fasta','_polished_pilon')
    pilon_output = pilon_output.replace(f'{ont_remove}/','')
    pilon_fasta = f'{pilon_output}.fasta'
    
    # Check if the output file already exists
    if os.path.isfile(pilon_fasta):
        log_print(f"PASS:\tSkipping Pilon Polish, {pilon_fasta} already exists", log_file)
        return pilon_fasta
    else:
        os.makedirs(pilon_output, exist_ok=True)  # add exist_ok=True

    # Check for specific Third-Party JAR Files with adjusted program_dict
    program_dict = {"pilon": ("https://github.com/broadinstitute/pilon", "pilon*.jar")}
    pilon_jar_path_dict, base_dir = check_for_jars(program_dict, log_file)
    pilon_jar_path = pilon_jar_path_dict['pilon']

    # If the file doesn't exist, run Pilon
    pilon_cmd = ["java", f"-Xmx{ram_gb}G",
                 "-jar", pilon_jar_path,
                 "--genome", assembly_file,
                 "--frags", bam_file,
                 "--output", pilon_output,
                 "--changes", "--vcf", "--diploid",
                 "--threads", str(cpu_threads)]
    log_print(f"CMD:\t{' '.join(pilon_cmd)}", log_file)
    pilon_result = subprocess.run(pilon_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Check if Pilon ran successfully
    if pilon_result.returncode != 0:
        log_print(f"ERROR:\t{pilon_result.stderr}", log_file)
        return None
    else:
        log_print("PASS:\tSuccessfully polished assembly with Pilon", log_file)
        return pilon_fasta

# Function to polish Cleaned ONT Reads with Trimmed Illumina Reads
def final_pilon_polish(cleaned_ont_assembly, illumina_reads_list, CURRENT_ORGANISM_KINGDOM, PERCENT_RESOURCES, busco_db_dict, log_file):
    """
    Main function to call for the Pilon polishing of Cleaned ONT Reads with Trimmed Illumina Reads.
    
    Args:
        cleaned_ont_assembly (str): Path to the folder containing a sub-folder with the ONT reads.
        forward_reads (int): Expected size in Mega-Bytes/Bases for the genome (Fungi are about 60-80).
        reverse_reads (int) : Percentage of resources that can be utilized for processing (default = 80).
        PERCENT_RESOURCES (int): Amount of resources availalbe for use during processing.
        busco_db_dict (dict): Dictionary of all Kingdoms and their respective BUSCO databases.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.
    
    Returns:
        pilon_output (str): Path to the final polished Pilon FASTA file output.
    """    
    # Get the number of CPUs available on the system
    num_cpus = multiprocessing.cpu_count()
    
    # Get the amount of RAM (GB) currently available
    mem_info = psutil.virtual_memory()
    
    # Calculate the number of threads based on PERCENT_RESOURCES for available CPUs & RAM
    cpu_threads = int(math.floor(num_cpus * PERCENT_RESOURCES))
    ram_gb = int(mem_info.total / (1024.0 ** 3) * PERCENT_RESOURCES)
    
    output_sam = cleaned_ont_assembly.replace("_filtered.fasta", "_bwa_aligned.sam")
    output_bam = output_sam.replace(".sam", ".bam")
    
    if os.path.isfile(output_bam):
        log_print(f"PASS:\tSkipping Indexing, SAM, and BAM file generation, {output_bam} already exists", log_file)
        bam_file = output_bam
    else:
        # Index Cleaned ONT Flye Assembly with BWA
        index_cleaned_ont_assembly(cleaned_ont_assembly, log_file)
        
        # Establish data_type
        if len(illumina_reads_list) == 2:
            data_type = 'PE'
        elif len(illumina_reads_list) == 1:
            data_type = 'SE'
        
        print(data_type)
        
        if data_type == 'PE':
            # Generate the SAM file based on the indexed cleaned ont assembly and trimmed paired illumina reads
            sam_file = generate_illumina_sam(cleaned_ont_assembly, illumina_reads_list[0], illumina_reads_list[1], cpu_threads, log_file)
        elif data_type == 'SE':
            # Generate the SAM file based on the indexed cleaned ont assembly and trimmed paired illumina reads
            print('UNLOGGED:\t DEVELOP SAM FILE BASED ON SINGLE-END READS ONLY')
            # sam_file = generate_illumina_sam(cleaned_ont_assembly, illumina_reads_list[0], cpu_threads, log_file)
        
        # Generate the BAM file from the SAM file
        bam_file = convert_sam_to_bam(sam_file, cpu_threads, log_file)    
    
    # Perfrom final Pilon Polish on the cleaned ont assembly with the generated BAM file
    pilon_output = pilon_polish_assembly(cleaned_ont_assembly, bam_file, ram_gb, cpu_threads, log_file)    
    
    # Clean up after final pilon polish
    ONT_FOLDER = os.path.dirname(cleaned_ont_assembly)
    base_name = cleaned_ont_assembly.split('/')[-1].split('_ont')[0]
    nanostat_dir = f'{ONT_FOLDER}/{base_name}_ont_combined_NanoStat'
    # Cleanup ONT assembly bulk files and keep only the items_to_keep found in the ONT_FOLDER
    items_to_keep = [nanostat_dir.replace('_ont_combined_NanoStat',''), nanostat_dir,
                     nanostat_dir.replace('_ont_combined_NanoStat',f'_ont_flye_filtered_{busco_db_dict[CURRENT_ORGANISM_KINGDOM][0].split("/")[-1]}_busco'),
                     nanostat_dir.replace('_ont_combined_NanoStat',f'_ont_flye_filtered_{busco_db_dict[CURRENT_ORGANISM_KINGDOM][1].split("/")[-1]}_busco'),
                     nanostat_dir.replace('_ont_combined_NanoStat','_ont_flye_filtered_quast'),
                     cleaned_ont_assembly.replace('flye_filtered.fasta','combined.fastq'),
                     cleaned_ont_assembly.replace('_filtered',''), cleaned_ont_assembly, 
                     f'{ONT_FOLDER}/cumulative_removed_sequences.csv',
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
        quast_thread = Thread(target = assess_with_quast, args = (pilon_output, log_file, cpu_threads))
        quast_thread.start()
          
        # Quality Control Check Pilon Polished Assembly with BUSCO agasint first database
        first_pilon_busco_thread = Thread(target = assess_with_busco, args = (pilon_output, log_file, busco_db_dict[CURRENT_ORGANISM_KINGDOM][0]))
        first_pilon_busco_thread.start()

        # Quality Control Check Pilon Polished Assembly with BUSCO agasint the second database
        second_pilon_busco_thread = Thread(target = assess_with_busco, args = (pilon_output, log_file, busco_db_dict[CURRENT_ORGANISM_KINGDOM][1]))
        second_pilon_busco_thread.start()
            
        # Wait for all QC threads to finish
        quast_thread.join()
        first_pilon_busco_thread.join()
        second_pilon_busco_thread.join()
    
    return pilon_output

## Debuging Main Space & Example
if __name__ == "__main__":
    print('EGAP Pilon Pilosher')   
    # Argument Parsing
    parser = argparse.ArgumentParser(description='Pilon Polish Cleaned ONT Read with Trimmed Illumina Reads')

    # Get working environment information
    environment_dir = get_env_dir()
    
    # Default values
    default_ont = f'{environment_dir}/Entheome/Ps_aff_hopii/MODULAR_TEST/ONT_MinION/B1_3_ont_flye_filtered.fasta'
    default_illumina_f = f'{environment_dir}/Entheome/Ps_aff_hopii/MODULAR_TEST/Illumina_PE150/B1_3_trimmomatic_forward_paired.fq'
    default_illumina_r = f'{environment_dir}/Entheome/Ps_aff_hopii/MODULAR_TEST/Illumina_PE150/B1_3_trimmomatic_reverse_paired.fq'
    default_organism_kingdom = 'Funga'
    default_percent_resources = 80
    
    # Add arguments with default values
    parser.add_argument('--ont_input', default = default_ont,
                        help = f'Path to the Cleaned ONT Reads. (default: {default_ont})')
    parser.add_argument('--illumina_f_input', default = default_illumina_f,
                        help = f'Path to the Trimmed Forward Illumina Reads. (default: {default_illumina_f})')
    parser.add_argument('--illumina_r_input', default = default_illumina_r,
                        help = f'Path to the Trimmed Reverse Illumina Reads. (default: {default_illumina_r})')
    parser.add_argument('--organism_kingdom',default = default_organism_kingdom,
                        help = f'Kingdom the current organism data belongs to. (default: {default_organism_kingdom})')
    parser.add_argument('--resource_use', type = int, default = default_percent_resources,
                        help = f'Percent of Resources to use. (default: {default_percent_resources})')
    
    # Parse the arguments
    args = parser.parse_args()
    ONT_READS = args.ont_input
    ILLUM_F_READS = args.illumina_f_input
    ILLUM_R_READS = args.illumina_r_input
    CURRENT_ORGANISM_KINGDOM = args.organism_kingdom
    PERCENT_RESOURCES = (args.resource_use/100)
    
    # Generate Main Logfile
    illumina_reads_list = [ILLUM_F_READS, ILLUM_R_READS]
    base_folder = '/'.join(ONT_READS.split('/')[:-2])
    debug_log = f'{base_folder}Pilon_log.tsv'
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

    # Run main Pilon Polish
    pilon_output = final_pilon_polish(ONT_READS, illumina_reads_list, CURRENT_ORGANISM_KINGDOM, PERCENT_RESOURCES, busco_db_dict, log_file)
