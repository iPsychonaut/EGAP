# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 11:48:13 2023

@author: ian.michael.bollinger@gmail.com with the help of ChatGPT 4.0

Command Line Example:
    python EGAP_illumina.py --illumina_folder /path/to/folder --primer_type STRING --resource_use INTEGER

The --primer_type must be a string similar to 'TruSeq3-PE' to represent the Illumina primer type to use with trimmomatic.
The --illumina_folder must have a sub-folder with raw Illumina .fq.gz files and their matching MD5.txt file.

Note:
1. The script supports both Single-End (SE) and Paired-End (PE) datasets.
2. For SE datasets, the script expects files with a '_combined.fq' suffix.
3. For PE datasets, the script expects files with '_combined_1.fq' and '_combined_2.fq' suffixes.
"""  
import os, subprocess, glob, fnmatch, hashlib, shutil, gzip, argparse, multiprocessing, queue, psutil, math, sys
import pandas as pd
from threading import Thread
from log_print import log_print, generate_log_file
from check_tools import get_md5, search_directory_for_file, check_for_jars
from EGAP_qc import assess_with_fastqc

# Function to run MD5 checksums on all .FQ.GZ files in the provided directory given there is also an MD5.txt file
def md5_check(illumina_raw_data_dir, trimmomatic_output_dir, illumina_df, log_file):
    """
    Run MD5 checksums on all .FQ.GZ files in the provided directory given there is also an MD5.txt file.
    
    Args:
        illumina_raw_data_dir (str): Path to the main Illumina Raw Reads Data Directory.
        illumina_df (DataFrame): A dataframe that will contain all the data.
    """
    md5_file_path = os.path.join(illumina_raw_data_dir, 'MD5.txt')
    if os.path.exists(md5_file_path):
        with open(md5_file_path, 'r') as f:
            md5_data = f.read().splitlines()
        md5_dict = dict(line.split() for line in md5_data)
        md5_df = pd.DataFrame(md5_dict.items(), columns=['MD5', 'Filename'])
        illumina_df = pd.concat([illumina_df, md5_df], ignore_index=True)
        
        # For each gz file, extract and verify MD5 checksum
        for file_name in os.listdir(illumina_raw_data_dir):
            if file_name.endswith('.gz'):
                gz_file_path = os.path.join(illumina_raw_data_dir, file_name)
                out_file_path = os.path.join(trimmomatic_output_dir, file_name.rsplit('.gz', 1)[0])
        
                # Check if extracted file already exists, if it does, skip the unzipping and MD5 check process
                if os.path.isfile(out_file_path):
                    log_print(f"PASS:\tSkipping extraction and MD5 check: {out_file_path} already exists", log_file)
                    continue
                
                # Get original MD5 from the DataFrame
                matching_row = illumina_df[illumina_df['Filename'].str.contains(file_name)].head(1)
                original_md5 = matching_row['MD5'].iloc[0] if not matching_row.empty else 'None'
                
                # Verify MD5 checksums
                if original_md5 != 'None':
                    new_md5 = get_md5(gz_file_path)
                    if original_md5 != new_md5:
                        log_print(f"ERROR:\tMD5 checksum mismatch for {out_file_path}: original {original_md5}, new {new_md5}", log_file)
                        break
                    else:
                        log_print(f'PASS:\tOriginal MD5 MATCHES for {file_name.split("/")[-1]}', log_file)
                        # If file not already extracted, extract
                        if not os.path.isfile(out_file_path):
                            with gzip.open(gz_file_path, 'rb') as f_in:
                                with open(out_file_path, 'wb') as f_out:
                                    shutil.copyfileobj(f_in, f_out)
                else:
                    log_print(f"ERROR:\tOriginal MD5 checksum not found for {file_name.split('/')[-1]}", log_file)
                    break

# Function to extract all Illumina reads and check them against their provided md5 checksum
def illumina_extract_and_check(folder_name, log_file):
    """
    Extract and verify Illumina MD5 checksums.

    Args:
        folder_name (str): Path to the folder containing Illumina data.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.

    Returns:
        trimmomatic_output_dir (str): Path to the Trimmomatic output directory.
    """
    # Generate a dataframe to contain md5 checksum information
    illumina_df = pd.DataFrame(columns=['MD5', 'Filename'])
    
    # Search the Raw Reads directory for the MD5.txt file and process it
    illumina_raw_data_dir = search_directory_for_file(folder_name, "MD5.txt")
        
    log_print(f"Running MD5 Checksum Analysis on Raw Illumina FASTQ.GZ files in {illumina_raw_data_dir}...", log_file)
    if illumina_raw_data_dir is None:
        log_print(f"ERROR:\tNo directory containing 'MD5.txt' found within '{folder_name}'", log_file)
        sys.exit(f"ERROR:\tNo directory containing 'MD5.txt' found within '{folder_name}'")
    
    # Generate Trimmotaic Output Directory if it doesn't exist
    base_data_dir = '/'.join(illumina_raw_data_dir.split('/')[:-1])
    trimmomatic_output_dir = f'{illumina_raw_data_dir}_trimmomatic_output/'
    os.makedirs(trimmomatic_output_dir, exist_ok=True)
    
    # Extract base_fname from the first .fq file
    base_fname = f'{illumina_raw_data_dir.split("/")[-1]}'
    
    # Generate List of FQ files in the illumina raw data
    fq_files = [f for f in os.listdir(illumina_raw_data_dir) if f.endswith('.fq.gz')]
    combined_files = []
    
    # Determine if data is Paired-End (PE) or Single-End (SE)
    found_1 = False
    found_2 = False
    data_type = ''
    for file in fq_files:
        if '_1' in file:
            found_1 = True
        elif '_2' in file:
            found_2 = True
    if found_1 and found_2:
        data_type = 'PE'
    else:
        data_type = 'SE'
    
    # Run Paired-End MD5 check and combine
    if data_type == 'PE':
        log_print('Performing MD5 check on Paired-End Illumina reads...', log_file)
        
        # Generate combined file names    
        combined_1_file = os.path.join(base_data_dir, f'{base_fname}_combined_1.fq')
        combined_2_file = os.path.join(base_data_dir, f'{base_fname}_combined_2.fq')

        # Check if the final combined file exists
        if os.path.isfile(combined_1_file) and os.path.isfile(combined_2_file):
            log_print(f"PASS:\tSkipping MD5 check & concatenation: Combined files already exist {combined_1_file}, {combined_2_file}", log_file)
            combined_files.append(combined_1_file)
            combined_files.append(combined_2_file)
            return trimmomatic_output_dir, combined_files, data_type
        else:
            # If MD5.txt file exists, extract MD5 checksums and match with calculated MD5 checksums
            md5_check(illumina_raw_data_dir, trimmomatic_output_dir, illumina_df, log_file)
            log_print(f"Combining '_1.fq' and '_2.fq' files into: {combined_1_file}, {combined_2_file}...", log_file)
            # For each file pair, combine them if present
            for file_name in os.listdir(illumina_raw_data_dir): 
                # Combine files into one master forward and reverse file
                with open(combined_1_file, 'wb') as wfd:
                    for f in [os.path.join(trimmomatic_output_dir, f) for f in os.listdir(trimmomatic_output_dir) if base_fname in f and '_1.fq' in f]:
                        with open(f, 'rb') as fd:
                            shutil.copyfileobj(fd, wfd)
                with open(combined_2_file, 'wb') as wfd:
                    for f in [os.path.join(trimmomatic_output_dir, f) for f in os.listdir(trimmomatic_output_dir) if base_fname in f and '_2.fq' in f]:
                        with open(f, 'rb') as fd:
                            shutil.copyfileobj(fd, wfd)
                
        # Move combined files up one directory
        combined_files.append(combined_1_file)
        combined_files.append(combined_2_file)             
    
    # Run Single-End MD5 check and combine
    elif data_type == "SE":
        log_print('Performing MD5 check on Single-End Illumina reads...', log_file)   
        
        # Generate combined file names
        combined_file = os.path.join(trimmomatic_output_dir, f'{base_fname}_combined.fq')
        combined_paired_file = combined_file.replace('_combined.fq','_combined_trimmomatic_forward_paired.fq')
        if os.path.isfile(combined_file):
            log_print(f"PASS:\tSkipping MD5 check & concatenation: Combined file already exists {combined_file}", log_file)
            combined_files.append(combined_file)
            return trimmomatic_output_dir, combined_files, data_type
        else:
            log_print(f"Combining files into: {combined_file}...", log_file)
            
            # If MD5.txt file exists, extract MD5 checksums and match with calculated MD5 checksums
            md5_check(illumina_raw_data_dir, trimmomatic_output_dir, illumina_df, log_file)

            # For each file pair, combine them if present
            for file_name in os.listdir(illumina_raw_data_dir):
                # Combine files into one master forward and reverse file
                with open(combined_file, 'wb') as wfd:
                    for f in [os.path.join(trimmomatic_output_dir, f) for f in os.listdir(trimmomatic_output_dir) if base_fname in f and '.fq' in f]:
                        with open(f, 'rb') as fd:
                            shutil.copyfileobj(fd, wfd)
                
        # Move combined files up one directory
        combined_files.append(combined_file)
            
    else:
        log_print(f"ERROR:\tNo .fq files found in {trimmomatic_output_dir}", log_file)
        sys.exit(f"ERROR:\tNo .fq files found in {trimmomatic_output_dir}")

    return trimmomatic_output_dir, combined_files, data_type

# Function to trim Illumina FASTQ files and generate Forward Paired and Reverse Paired FASTQ Reads as a list
def trim_with_trimmomatic(folder_name, combined_files, data_type, ILLU_PRIMER_TYPE, option1, option2, log_file, cpu_threads):
    """
    Trim Illumina FASTQ files using Trimmomatic.

    Args:
        folder_name (str): Path to the folder containing Illumina data.
        data_type (str): Either 'PE' for Paired-End reads or 'SE' Single-End reads.
        ILLU_PRIMER_TYPE (str): Illumina primer type to use with Trimmomatic.
        option1 (str): Additional Trimmomatic option.
        option2 (str): Additional Trimmomatic option.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.
        cpu_threads (int): Number of CPU threads to use.

    Returns:
        fq_paired_list (list): List containing paths to Forward Paired and Reverse Paired FASTQ files.
        fastqc_output_dirs (list): List containing paths to FastQC output directories for the Forward Paired and Reverse Paired FASTQ files.
    """    
    # Check for trimmomatic jar file
    log_print(f"Running Trimmomatic on {data_type} files {combined_files} using {ILLU_PRIMER_TYPE}...", log_file)
    jar_paths_dict, _ = check_for_jars({"trimmomatic": ("https://github.com/usadellab/Trimmomatic",
                                                        "trimmomatic-*.jar",
                                                        "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip")},
                                       'EGAP')
    trimmomatic_jar_path = jar_paths_dict['trimmomatic']

    # Generate directories and lists
    fq_paired_list = []
    fastqc_output_dirs = []    
    base_name = '_'.join(os.path.basename(combined_files[0]).split('_')[:-2])
    primers_path = f"{'/'.join(trimmomatic_jar_path.split('/')[:-1])}/adapters/{ILLU_PRIMER_TYPE}.fa"        
    illumina_dir = '/'.join(combined_files[0].split('/')[:-1])
    trimmomatic_dir = f'{illumina_dir}/{base_name}_trimmomatic_output'
        
    if data_type == 'PE':
        # Check for input files
        fwd_file = combined_files[0]
        rev_file = combined_files[1]
        
        # Set up output files
        fwd_paired_out = os.path.join(illumina_dir, os.path.basename(fwd_file).replace('_combined_1.fq', '_trimmomatic_forward_paired.fq'))
        fwd_unpaired_out = os.path.join(trimmomatic_dir, os.path.basename(fwd_file).replace('_combined_1.fq', '_trimmomatic_forward_unpaired.fq'))
        rev_paired_out = os.path.join(illumina_dir, os.path.basename(rev_file).replace('_combined_2.fq', '_trimmomatic_reverse_paired.fq'))
        rev_unpaired_out = os.path.join(trimmomatic_dir, os.path.basename(rev_file).replace('_combined_2.fq', '_trimmomatic_reverse_unpaired.fq'))
        
        # Check if the final combined file exists
        if os.path.isfile(fwd_paired_out) and os.path.isfile(rev_paired_out):
            log_print(f"PASS:\tSkipping Trimmomtaic trimming: Trimmed files already exist {fwd_paired_out}, {rev_paired_out}", log_file)
            fq_paired_list.append(fwd_paired_out)
            fq_paired_list.append(rev_paired_out)     
        else:
            
            # Build command for Trimmomatic
            trimmomatic_cmd = ["java", "-jar", trimmomatic_jar_path, data_type, "-phred33", 
                               "-threads", str(cpu_threads),
                               fwd_file, rev_file, 
                               fwd_paired_out, fwd_unpaired_out, 
                               rev_paired_out, rev_unpaired_out, 
                               f"ILLUMINACLIP:{primers_path}:2:30:10:11", 
                               option1, option2, "SLIDINGWINDOW:50:25", "MINLEN:125"]
            
            # Execute command
            log_print(f"CMD:\t{' '.join(trimmomatic_cmd)}", log_file)
            trimmomatic_result = subprocess.Popen(trimmomatic_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = trimmomatic_result.communicate()
            return_code = trimmomatic_result.returncode
            if return_code != 0:
                log_print(f"ERROR:\tTrimmomatic failed for:\t{fwd_file}; {rev_file}; {stderr.decode('utf-8')}", log_file)
            else:
                log_print(f"PASS:\tGenerated all Trimmomatic output files successfully", log_file)

        # Return the paths to the trimmed files
        fq_paired_list = [fwd_paired_out, rev_paired_out]
        fastqc_output_dirs = [os.path.join(folder_name, 'fastqc_forward_paired'), os.path.join(folder_name, 'fastqc_reverse_paired')]
    
    elif data_type == 'SE':
        # Check for input files
        input_file = combined_files[0]

        # Set up output files
        paired_out = os.path.join(trimmomatic_dir, os.path.basename(fwd_file).replace('_combined_1.fq', '_trimmomatic_forward_paired.fq'))
        unpaired_out = os.path.join(trimmomatic_dir, os.path.basename(fwd_file).replace('_combined_1.fq', '_trimmomatic_forward_unpaired.fq'))

        # Check if the final combined file exists
        if os.path.isfile(paired_out):
            log_print(f"PASS:\tSkipping Trimmomtaic trimming: Trimmed files already exist {paired_out}", log_file)
            fq_paired_list.append(paired_out)  
        else:
            print(f'output: {paired_out}\n{unpaired_out}')        
            
            # Build command for Trimmomatic
            trimmomatic_cmd = ["java", "-jar", trimmomatic_jar_path, data_type, "-phred33", 
                               "-threads", str(cpu_threads),
                               input_file, paired_out, unpaired_out, 
                               f"ILLUMINACLIP:{primers_path}:2:30:10:11", 
                               option1, option2, "SLIDINGWINDOW:50:25", "MINLEN:125"]
            
            # Execute command
            log_print(f"CMD:\t{' '.join(trimmomatic_cmd)}", log_file)
            result = subprocess.run(trimmomatic_cmd, capture_output=True, text=True)
            log_print(result.stdout, log_file)
            if return_code != 0:
                log_print(f"ERROR:\tTrimmomatic failed for:\t{input_file}; {stderr.decode('utf-8')}", log_file)
            else:
                log_print(f"PASS:\tGenerated all Trimmomatic output files successfully", log_file)
                
        # Return the paths to the trimmed files
        fq_paired_list = [paired_out]
        fastqc_output_dirs = [os.path.join(folder_name, 'fastqc_forward_paired'), os.path.join(folder_name, 'fastqc_reverse_paired')]
    
    return fq_paired_list, fastqc_output_dirs

# Main Illumina Folder Processing Function
def process_illumina(ILLUMINA_FOLDER, ILLU_PRIMER_TYPE, PERCENT_RESOURCES, log_file):
    """
    Process the Illumina folder by extracting and verifying the data, trimming with Trimmomatic, and conducting quality control.

    Args:
        ILLUMINA_FOLDER (str): Path to the main Illumina folder.
        ILLU_PRIMER_TYPE (str): Illumina primer type to use with Trimmomatic.
        PERCENT_RESOURCES (int): Percentage of system resources to use.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.
    """
    # Extract and verify Illumina data
    trimmomatic_output_dir, combined_files, data_type = illumina_extract_and_check(ILLUMINA_FOLDER, log_file)
    
    # Check the number of available CPU cores and allocate based on the user's resource percentage
    num_cpus = multiprocessing.cpu_count()
    cpu_threads = int(round(num_cpus * (PERCENT_RESOURCES / 100.0),0))
    
    # Trim the Illumina data
    fq_paired_list, fastqc_output_dirs = trim_with_trimmomatic(trimmomatic_output_dir,
                                                               combined_files, data_type, ILLU_PRIMER_TYPE,
                                                               'SLIDINGWINDOW:4:15', 'MINLEN:36',
                                                               log_file, cpu_threads)

## QUALITY CONTROL CHECK AREA
    if __name__ != "__main__":
        pass
    else:
        r_fastqc_cpus = int(round(num_cpus/2, 0))
        f_fastqc_cpus = cpu_threads - r_fastqc_cpus
        
        # Establish FastQC output dirs
        fastqc_output_dirs = [fq_file.replace('trimmomatic','fastqc') for fq_file in fq_paired_list]
        fastqc_output_dirs = [fq_file.replace('.fq','') for fq_file in fq_paired_list]
        
        if data_type == 'PE':
            # Quality Control Check A FastQC on Trimmed Illumina Reads 
            f_fastqc_thread = Thread(target = assess_with_fastqc, args = (fq_paired_list[0], fastqc_output_dirs[0], log_file, f_fastqc_cpus))
            f_fastqc_thread.start()            
            r_fastqc_thread = Thread(target = assess_with_fastqc, args = (fq_paired_list[1], fastqc_output_dirs[1], log_file, r_fastqc_cpus))
            r_fastqc_thread.start()
            
            # Wait for all QC threads to finish
            f_fastqc_thread.join()
            r_fastqc_thread.join()
        elif data_type == 'SE':
            # Quality Control Check A FastQC on Trimmed Illumina Reads 
            fastqc_thread = Thread(target = assess_with_fastqc, args = (fq_paired_list[0], fastqc_output_dirs[0], log_file, f_fastqc_cpus))
            fastqc_thread.start()            
            
            # Wait for all QC threads to finish
            fastqc_thread.join()
        
    return fq_paired_list, fastqc_output_dirs, data_type

## Debuging Main Space & Example
if __name__ == "__main__":
    print('EGAP Illumina Pipeline')   
    # Argument Parsing
    parser = argparse.ArgumentParser(description='Process Illumina Folder')
    
    # Default values
    default_folder = f'/mnt/e/Entheome/Ps_aff_hopii/MODULAR_TEST/Illumina_PE150/B1_3'
    default_primer = 'TruSeq3-PE'
    default_percent_resources = 80
    
    # Add arguments with default values
    parser.add_argument('--illumina_folder', default = default_folder,
                        help = f'Path to the Illumina (PE or SE) Folder. (default: {default_folder})')
    parser.add_argument('--primer_type', default = default_primer,
                        help = f'Type of Illumina Primers used. (default: {default_primer})')
    parser.add_argument('--resource_use', type = int, default = default_percent_resources,
                        help = f'Percent of Resources to use. (default: {default_percent_resources})')
    
    # Parse the arguments
    args = parser.parse_args()
    ILLUMINA_FOLDER = args.illumina_folder
    ILLU_PRIMER_TYPE = args.primer_type
    PERCENT_RESOURCES = args.resource_use
    
    # Generate log file with the desired behavior
    debug_log = f'{ILLUMINA_FOLDER}/Illumina_log.tsv'
    log_file = generate_log_file(debug_log, use_numerical_suffix=False)
    
    # Run main Illumina Trimming function
    fq_paired_list, fastqc_output_dirs, data_type = process_illumina(ILLUMINA_FOLDER, ILLU_PRIMER_TYPE, PERCENT_RESOURCES, log_file)
