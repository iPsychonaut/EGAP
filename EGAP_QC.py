# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 11:46:13 2023

@author: ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting
"""
# Base Python Imports
import os, re, zipfile, argparse


# Required Python Imports
import pandas as pd


# Custom Python Imports
from EGAP_Tools import log_print, initialize_logging_environment, run_subprocess_cmd, get_resource_values


# Global variables
CPU_THREADS = 1
DEFAULT_LOG_FILE = None
ENVIRONMENT_TYPE = None


def assess_with_quast(assembly_file_path, assembly_type, quast_threads, CURRENT_ORGANISM_KINGDOM, REF_SEQ=None):
    """
    Quality Control Check an Assembly FASTA File with QUAST.
    
    Args:
        assembly_file_path (str): Path to the Assembly FASTA file.
        assembly_type (str): Description of assmelby type can only be: 'illu_only', 'ont_only', 'hybrid_same', 'hybrid_different'
        CURRENT_ORGANISM_KINGDOM
        REF_SEQ
        quast_threads (int): Number of available threads.
        
    Returns:
        quast_dir (str): Path to the folder containing the QUAST output data.

    Notes:
        - This function uses QUAST (Quality Assessment Tool for Genome Assemblies) to evaluate the quality of an assembly file.
        - It processes the assembly file without a reference genome.

    Considerations:
        - Ensure QUAST is properly installed and accessible in the environment where the script is run.
        - The function assumes the assembly file is in FASTA format.

    Examples:
        assess_with_quast("/path/to/assembly.fasta", "illu_only", 4)
    """
    log_print(f'QUAST Analysis of {os.path.basename(assembly_file_path)} with reference...')
    quast_dir = assembly_file_path.replace('.fasta','_quast')
    summary_file_name = os.path.join(quast_dir, "report.tsv")

    if CURRENT_ORGANISM_KINGDOM == "Funga":
        gene_finding_param = ["--glimmer","--fungus"]
    elif CURRENT_ORGANISM_KINGDOM == "Bacteria":
        gene_finding_param = ["--glimmer"]
    elif CURRENT_ORGANISM_KINGDOM == "Archea":
        gene_finding_param = ["--glimmer"]
    elif CURRENT_ORGANISM_KINGDOM == "Flora" or CURRENT_ORGANISM_KINGDOM == "Fauna":
        gene_finding_param = ["--glimmer", "--eukaryote"]

    # Check if the output directory already exists
    if os.path.isfile(summary_file_name):
        log_print(f"PASS:\tOutput summary file already exists: {summary_file_name}")
    else:
        # Create output directory
        os.makedirs(quast_dir, exist_ok=True)
    
        # Check if report file exists
        if not os.path.isfile(summary_file_name):
            if REF_SEQ != None:                   
                REFERENCE_GFF = REF_SEQ.replace("fna.gz","gff")
                # Generate quast command with updated parameters
                quast_cmd = [
                    "quast",
                    "-o", quast_dir,
                    "-t", str(quast_threads),
                    "-r", REF_SEQ,
                    "-g", REFERENCE_GFF,
                    "--operons", REFERENCE_GFF,
                    "--circos",
                    "--plots-format", "svg",
                    "--min-alignment", "100",         # Increase minimum alignment length
                    "--min-identity", "97.0",         # Increase alignment identity threshold
                    "--ambiguity-usage", "all",       # Adjust ambiguity resolution
                    "--min-contig", "1000",           # Filter short contigs
                    "--extensive-mis-size", "1000",    # Adjust extensive misassembly size
                    "--fragmented"
                ] + gene_finding_param + [assembly_file_path]
            else:
                # Generate quast command with updated parameters
                quast_cmd = [
                    "quast",
                    "-o", quast_dir,
                    "-t", str(quast_threads),
                    "--circos",
                    "--plots-format", "svg",
                    "--min-alignment", "100",         # Increase minimum alignment length
                    "--min-identity", "97.0",         # Increase alignment identity threshold
                    "--ambiguity-usage", "all",       # Adjust ambiguity resolution
                    "--min-contig", "1000",           # Filter short contigs
                    "--extensive-mis-size", "1000",    # Adjust extensive misassembly size
                    "--fragmented"
                    ] + gene_finding_param + [assembly_file_path]
                
            _ = run_subprocess_cmd(quast_cmd, False)
    
    # # TODO: REWORK QUAST PASS/WARN/ERROR STATS ANALYSIS OF DATA SUMMARY FILE
    # # Check summary file
    # quast_stats = {}

    # # Open the summary file
    # with open(summary_file_name, 'r') as f:
    #     # Create a CSV reader that splits on tabs
    #     reader = csv.reader(f, delimiter='\t')
    
    #     # Skip the header row
    #     next(reader)
    
    #     # Process each line
    #     for row in reader:
    #         # Remove leading/trailing whitespace and remove any "(>= xxx bp)" suffix
    #         stat_name = re.sub(r'\s*\(\>= \d+ bp\)\s*', '', row[0]).strip()
    #         # Convert stat_value to float if possible, otherwise leave as string
    #         try:
    #             stat_value = float(row[1])
    #         except ValueError:
    #             stat_value = row[1]
    #         # Add to dictionary
    #         quast_stats[stat_name] = stat_value
    
    # # Check quality of assembly
    # total_contigs = quast_stats['# contigs']
    # contigs_10k = quast_stats['# contigs']
    # N50 = quast_stats['N50']
    
    # # Quast Logic Check for Illumina only VS ONT Only VS Hybrid Assembly
    # if assembly_type == "illu_only":
    #     # Logic for QUAST Illumina Only Stats PASS-WARN-ERROR
    #     if contigs_10k > 0.50 * total_contigs and N50 > 50000:
    #         log_print(f"PASS:\tHigh-quality assembly based on QUAST stats: {contigs_10k} contigs > 80% of total contigs; N50: {N50} > 1000")
    #     elif contigs_10k > 0.25 * total_contigs and N50 > 10000:
    #         log_print(f"NOTE:\tModerate-quality assembly based on QUAST stats: {contigs_10k} contigs > 50% of total contigs; N50: {N50} > 500")
    #     else:
    #         log_print(f"ERROR:\tLow-quality assembly based on QUAST stats: {contigs_10k} contigs; N50: {N50}")    
    # elif assembly_type == "ont_only":
    #     # TODO: Logic for QUAST ONT Only Stats PASS-WARN-ERROR
    #     log_print("ERROR:\tNO LOGIC IN PLACE FOR ONT ONLY ASSEMBLY")
        
    # elif assembly_type == "hybrid_same" or "hybrid_different":
    #     # Logic for QUAST Hybrid Stats PASS-WARN-ERROR
    #     if contigs_10k > 0.75 * total_contigs and N50 > 100000:
    #         log_print(f"PASS:\tHigh-quality assembly based on QUAST stats: {contigs_10k} contigs > 80% of total contigs; N50: {N50} > 100000")
    #     elif contigs_10k > 0.50 * total_contigs and N50 > 50000:
    #         log_print(f"NOTE:\tModerate-quality assembly based on QUAST stats: {contigs_10k} contigs > 50% of total contigs; N50: {N50} > 50000")
    #     else:
    #         log_print(f"ERROR:\tLow-quality assembly based on QUAST stats: {contigs_10k} contigs; N50: {N50}")    
    # else:
    #     log_print(f"ERROR:\tInvalid assembly type: {assembly_type}")
    
    return quast_dir


def assess_with_fastqc(sequence_file, fastqc_dir, fastqc_threads):
    """
    Processes basecalled ONT reads into a cleaned and indexed de novo assembly.
    
    Args:
        sequence_file (str): Path to a FASTQ file to analyze.
        fastqc_dir (str): Path to the directory to store output data.
        fastqc_threads (int): Number of available threads.
    
    Returns:
        fastqc_dir (str): Path to the folder contaiing the FASTQC ouput data.

    Notes:
        - FastQC is used for quality control checks on sequence data stored in FASTQ format.
        - The function generates a directory with FastQC reports.

    Considerations:
        - FastQC should be installed and accessible in the script's running environment.
        - Handles the creation of the output directory if it does not exist.

    Examples:
        assess_with_fastqc("/path/to/sequence.fq", "/path/to/fastqc_output", "/path/to/log.txt", 4)
    """
    log_print(f"Running FastQC Analysis on {os.path.basename(sequence_file)}...")
    if not os.path.isfile(fastqc_dir):
        os.makedirs(fastqc_dir)
    
    # Construct the expected summary file path
    summary_file_path = os.path.join(fastqc_dir, os.path.basename(sequence_file).replace('.fq','.fastq').replace('.fastq',"_fastqc_summary.txt"))
    
    # Check if the summary file already exists
    if not os.path.isfile(summary_file_path):
        # Construct the expected zip file path
        zip_file = os.path.join(fastqc_dir, os.path.basename(sequence_file).replace('.fq','.fastq').replace('.fastq',"_fastqc.zip"))

        # Check if the zip file already exists
        if not os.path.isfile(zip_file):
            # Construct the FastQC command
            fastqc_cmd = ["fastqc",
                          sequence_file,
                          "--threads",
                          str(fastqc_threads),
                          "-o",
                          fastqc_dir]
            _ = run_subprocess_cmd(fastqc_cmd, False)

        else:
            log_print(f"NOTE:\tSkipping FastQC Analysis: FastQC folder for {sequence_file} already exists")

    return fastqc_dir

    # # TODO: REWORK FASTQC PASS/WARN/ERROR STATS ANALYSIS OF DATA IN ZIP FILE        
    #     # Get the expected directory name within the zip file
    #     zip_dir_name = os.path.basename(sequence_file).replace('.fq',"_fastqc")
        
    #     # Extract the summary.txt file from the zip file
    #     with zipfile.ZipFile(zip_file, "r") as zip_ref:
    #         summary_file_in_zip = os.path.join(zip_dir_name, "summary.txt")
    #         if summary_file_in_zip in zip_ref.namelist():
    #             zip_ref.extract(summary_file_in_zip, fastqc_dir)
    #         else:
    #             log_print(f"ERROR:\tUnable to find {summary_file_in_zip} in {zip_file}")
    #             return None
        
    #     # Rename the summary.txt file to include the sequence file base name
    #     os.rename(os.path.join(fastqc_dir, summary_file_in_zip), summary_file_path)
    
    # # Read the summary.txt file into a pandas DataFrame
    # fastqc_df = pd.read_csv(summary_file_path, sep="\t", header=None)
    # fastqc_df.columns = ["status", "module", "filename"]
    
    # # Logic for FastQC PASS-WARN-ERROR
    # if not fastqc_df.empty:
    #     log_print(f"PASS:\tFastQC completed for {sequence_file}:")
    
    # return fastqc_dir
    
    
def parse_nanostats(file_path):
    """
    Parses a NanoStat TSV file using pandas to extract the median read length and median read quality.

    Parameters:
    -----------
    file_path : str
        Path to the NanoStats TSV file.

    Returns:
    --------
    tuple[int, float]
        A tuple containing:
        - length_threshold: int -> Median read length.
        - quality_threshold: int -> Median read quality.

    Raises:
    -------
    FileNotFoundError:
        If the file cannot be found at the specified path.
    KeyError:
        If the required fields are not present in the TSV.

    Example:
    --------
    >>> length, quality = parse_nanostats('/mnt/data/NanoStats.tsv')
    >>> print(f"Length Threshold: {length}, Quality Threshold: {quality}")
    Length Threshold: 2355, Quality Threshold: 13
    """
    try:
        df = pd.read_csv(file_path, sep='\t', index_col=0)
        length_threshold = int(round(float(df.loc["mean_read_length"].values[0]),0))
        quality_threshold = int(round(float(df.loc["median_qual"].values[0]), 0))

    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {file_path}")
    except KeyError as e:
        raise KeyError(f"Missing required data: {e}")

    # PASS/NOTE/WARN/ERROR STATS ANALYSIS OF DATA SUMMARY FILE
    if length_threshold < 2000:
        log_print(f"NOTE:\tMedian Length Threshold is < 2000: {length_threshold}")
    if length_threshold < 1000:
        log_print(f"WARN:\tMedian Length Threshold is < 1000: {length_threshold}")
    elif length_threshold < 500:
        log_print(f"ERROR:\tPoor Quality Reads; Median Length Threshold is < 500: {length_threshold}")
    elif length_threshold >= 2000:
        log_print(f"PASS:\tGood Quality Reads; Median Length Threshold is > 2000: {length_threshold}")
    if quality_threshold < 17:
        log_print(f"NOTE:\tMedian Quality Threshold is < 15: {quality_threshold}")
    if quality_threshold < 12:
        log_print(f"WARN:\tMedian Quality Threshold is < 15: {quality_threshold}")
    elif quality_threshold < 10:
        log_print(f"ERROR:\tPoor Quality Reads; Median Length Threshold is < 10: {quality_threshold}")
    elif quality_threshold >= 17:
        log_print(f"PASS:\tGood Quality Reads; Median Length Threshold is => 17: {quality_threshold}")

    return length_threshold, quality_threshold


def assess_with_nanoplot(input_fastq, nanoplot_threads):
    """
    Analyze a FASTQ file using NanoPlot to retrieve statistics about the sequencing run.
    
    Args:
        input_fastq (str): Path to the FASTQ file to be analyzed.
        nanoplot_threads (int): Number of available threads.

    Notes:
        - NanoPlot is utilized for statistical analysis of Nanopore sequencing data.
        - The function generates a directory with NanoPlot results in CSV format.

    Considerations:
        - Ensure NanoPlot is installed and the input FASTQ file is from Nanopore sequencing.
        - Output directory is created if it does not exist.

    Examples:
        assess_with_nanoplot("/path/to/reads.fastq", 4)
    """
    log_print(f"Running NanoPlot Analysis on {os.path.basename(input_fastq)}...")
    nanoplot_dir = input_fastq.replace('.fastq','_NanoPlot')
    tsv_path = os.path.join(nanoplot_dir, "NanoStats.txt")
    
    if not os.path.isfile(tsv_path):
        # Generate nanostat_dir if it doesn't exist
        if not os.path.isfile(nanoplot_dir):
            os.makedirs(nanoplot_dir)
        
        # Construct the NanoStat command
        nanoplot_command = ["NanoPlot",
                            "--verbose",
                            "-t", str(nanoplot_threads),
                            "--fastq",
                            input_fastq,
                            "--outdir",
                            nanoplot_dir,
                            "--tsv_stats"]
        
        # Run the command
        _ = run_subprocess_cmd(nanoplot_command, False)
    
    length_threshold, quality_threshold = parse_nanostats(tsv_path)
      
    return nanoplot_dir, length_threshold, quality_threshold


def assess_with_compleasm(input_fasta, lineage_name, compleasm_threads):
    """
    Quality Control Check of an organism sequence with a provided compleasm Database.
    
    Args:
        input_fasta (str): Path to the FASTA file to be analyzed.
        lineage_name (str): name of the organism's lineage database for BUSCO.
        compleasm_threads (int): Number of available threads.
    
    Returns:
        compleasm_output (str): Path to the folder containing the compleasm output data.

    Notes:
        - CompleAsm is used for completeness assessment of the assembly against a specific lineage database.
        - Generates output directory with CompleAsm results.

    Considerations:
        - Ensure CompleAsm is installed and accessible, along with the required lineage database.
        - The function assumes the assembly file is in FASTA format.

    Examples:
        assess_with_compleasm("/path/to/assembly.fasta", "lineage_name", 4)
    """
    log_print(f'Compleasm Analysis of {os.path.basename(input_fasta)} against {lineage_name}...')
    
    input_dir = os.path.dirname(input_fasta)
    base_name = os.path.basename(input_fasta)
    compleasm_output = os.path.join(input_dir, base_name.replace('.fasta',f"_{lineage_name}_compleasm"))
    
    # Check summary file
    summary_file_name = os.path.join(input_dir, compleasm_output, "summary.txt")
    if os.path.isfile(summary_file_name):
        log_print(f"PASS:\tOutput summary file already exists: {summary_file_name}")
    else:
        log_print("Running Compleasm command...")

        # Generate compleasm command
        compleasm_cmd = ["compleasm", "run",
                         "-a", input_fasta,
                         "-o", compleasm_output,
                         "-l", lineage_name,
                         "-t", str(compleasm_threads)]
         
        # Run compleasm
        _ = run_subprocess_cmd(compleasm_cmd, shell_check=False)
        
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


def assembly_qc_checks(input_assembly_path, READS_DATA, CURRENT_ORGANISM_KINGDOM, REF_SEQ, CPU_THREADS):
    """
    Performs quality control checks on the Illumina Only scaffolded assembly using QUAST and CompleAsm.
    
    This function runs two quality control tools: QUAST and CompleAsm. QUAST is used for assessing the quality of 
    the assembled sequences, while CompleAsm is run twice with different lineage specifications (Basidiomycota and 
    Agaricales) for completeness analysis of the assembly.
    
    Parameters:
        input_assembly_path (str): The path to the scaffolded genome assembly file.
        READS_DATA (str): Indicate if the provided data are generated from the single source: 'illu' or 'ont'; or 'hybrid'.
        CURRENT_ORGANISM_KINGDOM (str): Kingdom the current organism data belongs to.
        REF_SEQ (str): Path to the reference gzipped fasta to use.
        CPU_THREADS (int): Number of available threads.
    
    Notes:
        - The output of QUAST is stored in a directory named after the scaffolded assembly file with '_quast' appended.
        - CompleAsm is used to analyze the completeness of the assembly relative to specific lineages, namely 
          Basidiomycota and Agaricales. The output for each lineage is stored in separate directories.
        - This function orchestrates running multiple quality control checks (QUAST, CompleAsm) on genome assemblies.
        - Utilizes threading to run the checks concurrently for efficiency.

    Considerations:
        - Ensure all dependent tools (QUAST, CompleAsm) are installed and properly configured.
        - Adjust the number of CPU threads based on the available system resources.

    Examples:
        assembly_qc_checks("/path/to/assembly.fasta", "illu", "Funga", 4)
    """    
    # Generate BUSCO Database Dictionary
    busco_db_dict = {'Archaea': ['archaea',
                                  'euryarchaeota',],
                      'Bacteria': ['actinobacteria',
                                  'proteobacteria',],
                      'Fauna': ['vertebrata',
                                'arthropoda',],
                      'Flora': ['eudicots',
                                'liliopsida'],
                      'Funga': ['basidiomycota',
                                'agaricales'],
                      'Protista': ['alveolata',
                                  'euglenozoa']}
            
    quast_dir = assess_with_quast(input_assembly_path, READS_DATA, CPU_THREADS, CURRENT_ORGANISM_KINGDOM, REF_SEQ)
    compleasm_ouput_1 = assess_with_compleasm(input_assembly_path, busco_db_dict[CURRENT_ORGANISM_KINGDOM][0], CPU_THREADS)
    compleasm_outupt_2 = assess_with_compleasm(input_assembly_path, busco_db_dict[CURRENT_ORGANISM_KINGDOM][1], CPU_THREADS)
    
    return quast_dir, compleasm_ouput_1, compleasm_outupt_2


## Debuging Main Space & Example
if __name__ == "__main__":
    print('EGAP Quality Control Checks')    
    # Argument Parsing
    parser = argparse.ArgumentParser(description='EGAP Quality Control Checks')
    
    # Default values
    default_file = '/mnt/d/Entheome/Ps_semilanceata/Ps_semilanceata_A1_3_polished_pilon.fasta'
    default_reads_data = "hybrid"
    default_organism_kingdom = 'Funga'
    default_percent_resources = 0.9
    default_reference_fasta = "/mnt/d/Tryptomics/Genomics/MG_Ps_cubensis_PE_REF_SEQ/GCA_017499595.fna.gz"
    
    # Add arguments with default values
    parser.add_argument('--input_file', default = default_file,
                        help = f'Path to the File to provide QC on. (default: {default_file})')
    parser.add_argument('--reads_data', '-d',
                        type = str, default = default_reads_data,
                        choices = ["illu","ont","hybrid"],
                        help = f"Indicate if the provided data are generated from the single source: 'illu' or 'ont'; or 'hybrid'. (default: {default_reads_data})")
    parser.add_argument('--organism_kingdom',default = default_organism_kingdom,
                        help = f'Kingdom the current organism data belongs to. (default: {default_organism_kingdom})')
    parser.add_argument('--resource_use', type = float, default = default_percent_resources,
                        help = f'Percent of Resources to use. (default: {default_percent_resources})')
    parser.add_argument('--ref', "-r",
                        type = str, default = default_reference_fasta,
                        help = f'Path to the reference fasta to use. (default: {default_reference_fasta})')
    # Parse the arguments
    args = parser.parse_args()
    INPUT_FILE = args.input_file
    READS_DATA = args.reads_data
    CURRENT_ORGANISM_KINGDOM = args.organism_kingdom
    PERCENT_RESOURCES = args.resource_use
    REFERENCE_FASTA = args.ref
    
    CPU_THREADS, ram_gp = get_resource_values(PERCENT_RESOURCES)

    # Generate log file with the desired behavior
    input_folder = os.path.dirname(INPUT_FILE)
    initialize_logging_environment(input_folder)

    # Assembly (QUAST, 2x Compleasm BUSCOs) QC Debug    
    quast_output, compleasm_ouput_1, compleasm_outupt_2 = assembly_qc_checks(INPUT_FILE, READS_DATA, CURRENT_ORGANISM_KINGDOM, REFERENCE_FASTA, CPU_THREADS)

    ## ONT NanoStat QC Debug
    # raw_ont_fastq = "/mnt/d/Entheome/Ps_semilanceata/ONT_MinION/*.fastq"
    # nano_stat_dir = assess_with_nanostat(raw_ont_fastq, CPU_THREADS)

    
    ## Illumina FastQC Debug
    # fq_paired_list = ["/mnt/d/Entheome/Ps_semilanceata/Illumina_PE150/*_1.fastq", "/mnt/d/Entheome/Ps_semilanceata/Illumina_PE150/*_2.fastq"]
    # fastqc_output_dirs = [fq_file.replace('trimmomatic','fastqc') for fq_file in fq_paired_list]
    # fastqc_output_dirs = [fq_file.replace('.fq','') for fq_file in fq_paired_list]
    # f_fastqc_cpus = int(round(CPU_THREADS/2, 0))
    # r_fastqc_cpus = CPU_THREADS-f_fastqc_cpus
    
    # # Quality Control Check A FastQC on Trimmed Illumina Reads 
    # f_fastqc_thread = Thread(target = assess_with_fastqc, args = (fq_paired_list[0], fastqc_output_dirs[0], f_fastqc_cpus))
    # f_fastqc_thread.start()            
    # r_fastqc_thread = Thread(target = assess_with_fastqc, args = (fq_paired_list[1], fastqc_output_dirs[1], r_fastqc_cpus))
    # r_fastqc_thread.start()
    
    # # Wait for all QC threads to finish
    # f_fastqc_thread.join()
    # r_fastqc_thread.join()
