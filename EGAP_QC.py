# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 11:46:13 2023

@author: ian.michael.bollinger@gmail.com with the help of ChatGPT 4.0
"""
import os, subprocess, csv, re, zipfile, multiprocessing, math, psutil, argparse
import pandas as pd
from EGAP_Tools import log_print, initialize_logging_environment, get_env_dir, find_file, run_subprocess_cmd

CPU_THREADS = 1
DEFAULT_LOG_FILE = None
ENVIRONMENT_TYPE = None

# Function to run QUAST on an assembly
def assess_with_quast(assembly_file_path, assembly_type):
    """
    Quality Control Check an Assembly FASTA File with QUAST.
    
    Args:
        assembly_file_path (str): Path to the Assembly FASTA file.
        assembly_type (str): Description of assmelby type can only be: 'illu_only', 'ont_only', 'hybrid_same', 'hybrid_different'
        
    Returns:
        quast_dir (str): Path to the folder containing the QUAST output data.
    """
    global CPU_THREADS
    log_print(f'QUAST Analysis of {os.path.basename(assembly_file_path)} without reference...')
    quast_dir = assembly_file_path.replace('.fasta','_quast')
    summary_file_name = os.path.join(quast_dir, "report.tsv")

    # Check if the output directory already exists
    if os.path.isfile(summary_file_name):
        log_print(f"PASS:\tOutput summary file already exists: {summary_file_name}")
    else:
        # Create output directory
        os.makedirs(quast_dir, exist_ok=True)
    
        # Check if report file exists
        if not os.path.isfile(summary_file_name):
            # Generate quast command
            quast_cmd = ["quast",
                         "-o", quast_dir,
                         "-t", str(CPU_THREADS),
                         assembly_file_path]
            _ = run_subprocess_cmd(quast_cmd, False)
        
    # Check summary file
    quast_stats = {}

    # Open the summary file
    with open(summary_file_name, 'r') as f:
        # Create a CSV reader that splits on tabs
        reader = csv.reader(f, delimiter='\t')
    
        # Skip the header row
        next(reader)
    
        # Process each line
        for row in reader:
            # Remove leading/trailing whitespace and remove any "(>= xxx bp)" suffix
            stat_name = re.sub(r'\s*\(\>= \d+ bp\)\s*', '', row[0]).strip()
            # Convert stat_value to float if possible, otherwise leave as string
            try:
                stat_value = float(row[1])
            except ValueError:
                stat_value = row[1]
            # Add to dictionary
            quast_stats[stat_name] = stat_value
    
    # Check quality of assembly
    total_contigs = quast_stats['# contigs']
    contigs_10k = quast_stats['# contigs']
    N50 = quast_stats['N50']
    
    # Quast Logic Check for Illumin only VS ONT Only VS Hybrid Assembly
    if assembly_type == "illu_only":
        # Logic for QUAST Illumina Only Stats PASS-WARN-ERROR
        if contigs_10k > 0.8 * total_contigs and N50 > 1000:
            log_print(f"PASS:\tHigh-quality assembly based on QUAST stats: {contigs_10k} contigs > 80% of total contigs; N50: {N50} > 1000")
        elif contigs_10k > 0.5 * total_contigs and N50 > 500:
            log_print(f"NOTE:\tModerate-quality assembly based on QUAST stats: {contigs_10k} contigs > 50% of total contigs; N50: {N50} > 500")
        else:
            log_print(f"ERROR:\tLow-quality assembly based on QUAST stats: {contigs_10k} contigs; N50: {N50}")    
    elif assembly_type == "ont_only":
        # TODO: Logic for QUAST ONT Only Stats PASS-WARN-ERROR
        log_print("ERROR:\tNO LOGIC IN PLACE FOR ONT ONLY ASSEMBLY")
        
    elif assembly_type == "hybrid_same" or "hybrid_different":
        # Logic for QUAST Hybrid Stats PASS-WARN-ERROR
        if contigs_10k > 0.8 * total_contigs and N50 > 100000:
            log_print(f"PASS:\tHigh-quality assembly based on QUAST stats: {contigs_10k} contigs > 80% of total contigs; N50: {N50} > 100000")
        elif contigs_10k > 0.5 * total_contigs and N50 > 50000:
            log_print(f"NOTE:\tModerate-quality assembly based on QUAST stats: {contigs_10k} contigs > 50% of total contigs; N50: {N50} > 50000")
        else:
            log_print(f"ERROR:\tLow-quality assembly based on QUAST stats: {contigs_10k} contigs; N50: {N50}")    
    else:
        log_print(f"ERROR:\tInvalid assembly type: {assembly_type}")
    return quast_dir

# Function to run FastQC on a provided sequence file an produce a fastqc_df from the summary.txt
def assess_with_fastqc(sequence_file, fastqc_dir, log_file, cpu_threads):
    """
    Processes basecalled ONT reads into a cleaned and indexed de novo assembly.
    
    Args:
        sequence_file (str): Path to a FASTQ file to analyze.
        fastqc_dir (str): Path to the directory to store output data.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.
        cpu_threads (int): Number of CPU threads to utilize.
    
    Returns:
        fastqc_dir (str): Path to the folder contaiing the FASTQC ouput data.
    """
    log_print(f"Running FastQC Analysis on {os.path.basename(sequence_file)}...")
    if not os.path.exists(fastqc_dir):
        os.makedirs(fastqc_dir)
    
    # Construct the expected summary file path
    summary_file_path = os.path.join(fastqc_dir, os.path.basename(sequence_file).replace('.fq',"_fastqc_summary.txt"))
    
    # Check if the summary file already exists
    if not os.path.exists(summary_file_path):
        # Construct the expected zip file path
        zip_file = os.path.join(fastqc_dir, os.path.basename(sequence_file).replace('.fq',"_fastqc.zip"))

        # Check if the zip file already exists
        if not os.path.exists(zip_file):
            # Construct the FastQC command
            fastqc_cmd = ["fastqc",
                          sequence_file,
                          "--threads",
                          str(cpu_threads),
                          "-o",
                          fastqc_dir]
            _ = run_subprocess_cmd(fastqc_cmd, False)

        else:
            log_print(f"PASS:\tSkipping FastQC Analysis: FastQC folder for {sequence_file} already exists")
        
        # Get the expected directory name within the zip file
        zip_dir_name = os.path.basename(sequence_file).replace('.fq',"_fastqc")
        
        # Extract the summary.txt file from the zip file
        with zipfile.ZipFile(zip_file, "r") as zip_ref:
            summary_file_in_zip = os.path.join(zip_dir_name, "summary.txt")
            if summary_file_in_zip in zip_ref.namelist():
                zip_ref.extract(summary_file_in_zip, fastqc_dir)
            else:
                log_print(f"ERROR:\tUnable to find {summary_file_in_zip} in {zip_file}")
                return None
        
        # Rename the summary.txt file to include the sequence file base name
        os.rename(os.path.join(fastqc_dir, summary_file_in_zip), summary_file_path)
    
    # Read the summary.txt file into a pandas DataFrame
    fastqc_df = pd.read_csv(summary_file_path, sep="\t", header=None)
    fastqc_df.columns = ["status", "module", "filename"]
    
    # Logic for FastQC PASS-WARN-ERROR
    if not fastqc_df.empty:
        log_print(f"PASS:\tFastQC completed for {sequence_file}:")
        status_counts = fastqc_df['status'].value_counts()
        warn_count = status_counts.get('WARN', 0)
        
        # TODO: If warning is in Overrepresented sequences AND?/OR? Adapter Content then there is likely a Failed Primers Trimming Error
        # Rules for Three (3) Buckets: PASS, WARN, ERROR
        if warn_count == 1:
             log_print(f'PASS:\tFastQC Count: {warn_count}. {sequence_file} is acceptable')
        elif warn_count < 3:
             log_print(f'WARN:\tFastQC Count Approaching limit (3): {warn_count}. {sequence_file} is acceptable; should be reviewed')
        elif warn_count > 3:
            log_print(f'ERROR:\tFastQC Count Exceeded limit (3): {warn_count} data requires review')
        else:
            log_print(f'PASS:\tFastQC: {sequence_file} is acceptable')
    
    return fastqc_dir

# Function to analyze a ONT Reads as a FASTQ file using NanoStat
def assess_with_nanostat(input_fastq, cpu_threads):
    """
    Analyze a FASTQ file using NanoStat to retrieve statistics about the sequencing run.
    
    Args:
        input_fastq (str): Path to the FASTQ file to be analyzed.
        cpu_threads (int): Number of CPU threads to be used by NanoStat.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.
    """
    log_print(f"Running NanoStat Analysis on {os.path.basename(input_fastq)}...")
    nanostat_dir = input_fastq.replace('.fastq','_NanoStat')
    csv_path = os.path.join(nanostat_dir, "nanostat_results.csv")
    
    if not os.path.exists(csv_path):
        # Generate nanostat_dir if it doesn't exist
        if not os.path.exists(nanostat_dir):
            os.makedirs(nanostat_dir)
        
        # Construct the NanoStat command
        nanostat_command = ["NanoStat",
                            "--fastq",
                            input_fastq,
                            "--outdir",
                            nanostat_dir,
                            "--threads",
                            str(cpu_threads)]
        
        # Run the command
        log_print(f"CMD:\t{' '.join(nanostat_command)}")
        nanostat_result = subprocess.run(nanostat_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)    
        
        # Check if NanoStat ran successfully
        if nanostat_result.returncode != 0:
            log_print(f"ERROR:\tNanoStat failed: {nanostat_result.stdout.decode()} {nanostat_result.stderr.decode()}")
            return None
        
        # Extract and determine qualities from NanoStat output 
        nanostat_output = nanostat_result.stdout.decode().splitlines()
        output_dict = {}
        for line in nanostat_output:
            parts = line.split(":")
            if len(parts) == 2:
                key, value = parts
                output_dict[key.strip()] = value.strip()

        # Convert dictionary to DataFrame
        nanostat_df = pd.DataFrame([output_dict])
        
        # Save the DataFrame to a CSV file
        nanostat_df.to_csv(csv_path, index=False)
    
    # Load the saved CSV for data analysis
    loaded_df = pd.read_csv(csv_path)
   
    # Check for Mean read quality
    mean_quality = loaded_df["Mean read quality"].iloc[0]
    if mean_quality >= 10:
        log_print(f'PASS:\tMean read quality is greater than ten (>10): {mean_quality}')
    elif mean_quality >= 7:
        log_print(f'WARN:\tMean read quality is in the warning range, less than ten, greater than seven (>7): {mean_quality}')
    else:
        log_print(f'ERROR:\tMean read quality does not meet the minimum criteria of at least seven (>7): {mean_quality}')

    # Check for Read length N50
    n50_length = float(loaded_df["Read length N50"].iloc[0].replace(',', ''))
    if n50_length >= 7000:
        log_print(f'PASS:\tRead length N50 meets the criteria of at least 7000bp: {n50_length}')
    elif n50_length >= 3500:
        log_print(f'WARN:\tRead length N50 is in the warning range of at least 3500bp {n50_length}')
    else:
        log_print(f'ERROR:\tRead length N50 does not meet the minimum criteria of at least 3500bp: {n50_length}')

    # Check for Percentage of reads above Q10
    percentage_q10 = float(loaded_df[">Q10"].iloc[0].split(" ")[1].strip("()%"))
    if percentage_q10 >= 75:
        log_print(f'PASS:\tPercentage of reads above Q10 meets the criteria of at least seventy-five (>75): {percentage_q10}')
    elif percentage_q10 >= 60:
        log_print(f'WARN:\tPercentage of reads above Q10 is in the warning range of at least sixty (>60): {percentage_q10}')
    else:
        log_print(f'ERROR:\tPercentage of reads above Q10 does not meet the minimum criteria of at least sixty (>60): {percentage_q10}')
        
    return nanostat_dir

def assess_with_compleasm(input_fasta, lineage_name, compleasm_threads):
    """
    Quality Control Check of an organism sequence with a provided compleasm Database.
    
    Args:
        input_fasta (str): Path to the FASTA file to be analyzed.
        database_dir (str): Path to Database (lineage) for comparison against the input FASTA file.
    
    Returns:
        compleasm_output (str): Path to the folder containing the compleasm output data.
    """
    log_print(f'Compleasm Analysis of {os.path.basename(input_fasta)} against {lineage_name}...')
    
    input_dir = os.path.dirname(input_fasta)
    base_name = os.path.basename(input_fasta)
    compleasm_output = base_name.replace('.fasta',f"_{lineage_name}_compleasm")  
    
    # Default path to Compleasm
    default_compleasm_path = "compleasm.py" 
    compleasm_path = find_file(default_compleasm_path)
        
    # Generate compleasm command
    compleasm_cmd = ["python", compleasm_path, "run",
                     "-a", input_fasta,
                     "-o", compleasm_output,
                     "-l", lineage_name,
                     "-t", str(compleasm_threads)]
    _ = run_subprocess_cmd(compleasm_cmd, False)
 
    # Check summary file
    summary_file_name = os.path.join(input_dir, compleasm_output, "summary.txt")
    if os.path.isfile(summary_file_name):
        log_print(f"PASS:\tOutput summary file already exists: {summary_file_name}")
    else:
        # Run compleasm
        compleasm_result = subprocess.run(compleasm_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Check if compleasm ran successfully
        if compleasm_result.returncode != 0:
            log_print(f"ERROR:\t{compleasm_result.stderr.decode()}")
            return compleasm_output
        
    # Establish Summary file and key varibles
    compleasm_stats = {}
    
    # Define the regex patterns for the values we're interested in
    single_copy_pattern = re.compile(r"S:(\d+.\d+)%, \d+")
    duplicated_pattern = re.compile(r"D:(\d+.\d+)%, \d+")
    fragmented_pattern = re.compile(r"F:(\d+.\d+)%, \d+")
    incomplete_pattern = re.compile("I:(\d+.\d+)%, \d+")
    missing_pattern = re.compile("M:(\d+.\d+)%, \d+")
    n_count_pattern = re.compile("N:(\d+)")
    
    with open(summary_file_name, "r") as f:
        for line in f:
            # Parse the line to extract the values
            single_copy = re.search(single_copy_pattern, line)
            duplicated = re.search(duplicated_pattern, line)
            fragmented = re.search(fragmented_pattern, line)
            incomplete = re.search(incomplete_pattern, line)
            missing = re.search(missing_pattern, line)
            n_count = re.search(n_count_pattern, line)
            
            # If we found the values, add them to our stats dictionary
            if single_copy:
                compleasm_stats['Single copy'] = float(single_copy.group(1))
            if duplicated:
                compleasm_stats['Duplicated'] = float(duplicated.group(1))
            if fragmented:
                compleasm_stats['Fragmented'] = float(fragmented.group(1))
            if incomplete:
                compleasm_stats['Incomplete'] = float(incomplete.group(1))
            if missing:
                compleasm_stats['Missing'] = float(missing.group(1))
            if n_count:
                compleasm_stats['N_count'] = int(n_count.group(1))
        compleasm_stats['Completeness'] = round(compleasm_stats['Single copy'] + compleasm_stats['Duplicated'],2)
    
        # Logic for COMPLEASM PASS-WARN-ERROR
        if compleasm_stats['Completeness'] < 75:
            log_print(f"ERROR:\t{lineage_name} compleasm BUSCO completeness (n={compleasm_stats['N_count']}) {compleasm_stats['Completeness']}% < 75%; S: {compleasm_stats['Single copy']}%; D: {compleasm_stats['Duplicated']}%, F: {compleasm_stats['Fragmented']}; I: {compleasm_stats['Incomplete']}; M: {compleasm_stats['Missing']}")
        elif compleasm_stats['Completeness'] >= 75 and compleasm_stats['Completeness'] < 85:
            log_print(f"NOTE:\t{lineage_name} compleasm BUSCO completeness (n={compleasm_stats['N_count']}) {compleasm_stats['Completeness']}% < 75%; S: {compleasm_stats['Single copy']}%; D: {compleasm_stats['Duplicated']}%, F: {compleasm_stats['Fragmented']}; I: {compleasm_stats['Incomplete']}; M: {compleasm_stats['Missing']}")
        else:
            log_print(f"PASS:\t{lineage_name} compleasm BUSCO completeness (n={compleasm_stats['N_count']}) {compleasm_stats['Completeness']}% < 75%; S: {compleasm_stats['Single copy']}%; D: {compleasm_stats['Duplicated']}%, F: {compleasm_stats['Fragmented']}; I: {compleasm_stats['Incomplete']}; M: {compleasm_stats['Missing']}")
        return compleasm_output

def illu_only_qc_checks(scaffolded_assembly_path, CPU_THREADS):
    """
    Performs quality control checks on the Illumina Only scaffolded assembly using QUAST and CompleAsm.
    
    This function runs two quality control tools: QUAST and CompleAsm. QUAST is used for assessing the quality of 
    the assembled sequences, while CompleAsm is run twice with different lineage specifications (Basidiomycota and 
    Agaricales) for completeness analysis of the assembly.
    
    Parameters:
        scaffolded_assembly_path (str): The path to the scaffolded genome assembly file.
    
    Notes:
        - The output of QUAST is stored in a directory named after the scaffolded assembly file with '_quast' appended.
        - CompleAsm is used to analyze the completeness of the assembly relative to specific lineages, namely 
          Basidiomycota and Agaricales. The output for each lineage is stored in separate directories.
        - The function uses the global `CPU_THREADS` variable to set the number of threads for both QUAST and CompleAsm.
    """
    quast_out = assess_with_quast(scaffolded_assembly_path, "illu_only")

    compleasm_threads = int(round(CPU_THREADS/2,0))
    
    basidio_out = assess_with_compleasm(scaffolded_assembly_path, "basidiomycota", compleasm_threads)
    agaricales_out = assess_with_compleasm(scaffolded_assembly_path, "agaricales", compleasm_threads)
    
# Debuging Main Space & Example
if __name__ == "__main__":
    print('EGAP_qc.py DEBUG')    
    # Argument Parsing
    parser = argparse.ArgumentParser(description='Quality Control Checks')
    
    # Default values
    default_file = '/mnt/e/Entheome/Ps_aff_hopii/MODULAR_TEST/ONT_MinION/B1_3_ont_combined.fastq'
    default_organism_kingdom = 'Funga'
    default_percent_resources = 40
    
    # Add arguments with default values
    parser.add_argument('--input_file', default = default_file,
                        help = f'Path to the File to provide QC on. (default: {default_file})')
    parser.add_argument('--organism_kingdom',default = default_organism_kingdom,
                        help = f'Kingdom the current organism data belongs to. (default: {default_organism_kingdom})')
    parser.add_argument('--resource_use', type = int, default = default_percent_resources,
                        help = f'Percent of Resources to use. (default: {default_percent_resources})')
    
    # Parse the arguments
    args = parser.parse_args()
    INPUT_FILE = args.input_file
    CURRENT_ORGANISM_KINGDOM = args.organism_kingdom
    PERCENT_RESOURCES = args.resource_use
    
    # Get working environment information
    environment_dir = get_env_dir(INPUT_FILE)
    
    # Get the number of CPUs available on the system
    num_cpus = multiprocessing.cpu_count()

    # Get the amount of RAM (GB) currently available
    mem_info = psutil.virtual_memory()

    # Calculate the number of threads based on PERCENT_RESOURCES for available CPUs & RAM
    cpu_threads = int(math.floor(num_cpus * PERCENT_RESOURCES))
    ram_gb = int(mem_info.total / (1024.0 ** 3) * PERCENT_RESOURCES)

    # Generate log file with the desired behavior
    input_folder = os.path.dirname(INPUT_FILE)
    initialize_logging_environment(input_folder)
    
    # # Generate BUSCO Database Dictionary
    # busco_db_dict = {'Archaea': [f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/archaea_odb10',
    #                              f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/euryarchaeota_odb10',],
    #                  'Bacteria': [f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/actinobacteria_phylum_odb10',
    #                               f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/proteobacteria_odb10',],
    #                  'Fauna': [f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/vertebrata_odb10',
    #                            f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/arthropoda_odb10',],
    #                  'Flora': [f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/eudicots_odb10',
    #                            f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/liliopsida_odb10'],
    #                  'Funga': [f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/basidiomycota_odb10',
    #                            f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/agaricales_odb10'],
    #                  'Protista': [f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/alveolata_odb10',
    #                               f'{environment_dir}/EGAP/EGAP_Databases/BUSCO_Databases/euglenozoa_odb10']}
    
    # cleaned_ont_assembly = "/mnt/e/Entheome/Ps_aff_hopii/MODULAR_TEST/Ps_aff_hopii_B1_3_polished_pilon.fasta"
    # compleasm_output = assess_with_compleasm(cleaned_ont_assembly, log_file, busco_db_dict[CURRENT_ORGANISM_KINGDOM][0])
    # compleasm_output = assess_with_compleasm(cleaned_ont_assembly, log_file, busco_db_dict[CURRENT_ORGANISM_KINGDOM][1])
