#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
process_inputs.py

Updated on Tue Apr 01 2025

This script processes input parameters for the Entheome Genome Assembly Pipeline,
either from a CSV file or command-line parameters, and calculates resource values.
It also determines the Illumina input type (concatenated, directory, or SRA),
verifies MD5 checksums for Illumina raw reads and concatenates paired-end files,
downloads test data from SRA/NCBI, and combines ONT FASTQ.GZ files from a directory.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""

import os, sys, math, json, hashlib, subprocess, glob, shutil, gzip, zipfile, traceback
import pandas as pd
from Bio import SeqIO


# --------------------------------------------------------------
# Calculate resource allocation
# --------------------------------------------------------------
def get_resource_values(percent_resources, total_cpu, total_ram):
    """Calculate CPU threads and RAM based on a percentage of total resources.

    Args:
        percent_resources (float): Percentage of resources to allocate (0.0 to 1.0).
        total_cpu (int): Total available CPU threads.
        total_ram (int): Total available RAM in GB.

    Returns:
        tuple: (number of CPU threads, RAM in GB).
    """
    cpu_threads = int(math.floor(total_cpu * percent_resources))
    ram_gb = int(total_ram * percent_resources)
    return cpu_threads, ram_gb


# --------------------------------------------------------------
# Check existence of raw data files
# --------------------------------------------------------------
def all_files_exist(df, output_dir):
    """Verify if all expected raw data files exist in the base folder.

    Checks for Illumina, ONT, PacBio, and reference sequence files based on DataFrame entries.

    Args:
        df (pandas.DataFrame): DataFrame containing sample metadata.
        output_dir (str): Base directory for output files.

    Returns:
        bool: True if all expected files exist, False otherwise.
    """
    all_exist = True
    base_folders = set(row["SAMPLE_ID"].split("-")[0] for _, row in df.iterrows())  # Unique base folders
    for base_folder in base_folders:
        base_dir = os.path.join(output_dir, "EGAP_Test_Data", base_folder)
        for index, row in df.iterrows():
            sample_id = row["SAMPLE_ID"]
            if row["ILLUMINA_SRA"] != "None":
                illu_f = os.path.join(base_dir, "Illumina", f"{row['ILLUMINA_SRA']}_1.fastq.gz")
                illu_r = os.path.join(base_dir, "Illumina", f"{row['ILLUMINA_SRA']}_2.fastq.gz")
                if not (os.path.exists(illu_f) and os.path.exists(illu_r)):
                    print(f"Missing Illumina files for {sample_id} in base dir: {illu_f} or {illu_r}")
                    all_exist = False
                else:
                    print(f"Found Illumina files for {sample_id} in base dir: {illu_f}, {illu_r}")
            
            if row["ONT_SRA"] != "None":
                ont = os.path.join(base_dir, "ONT", f"{row['ONT_SRA']}.fastq.gz")
                if not os.path.exists(ont):
                    print(f"Missing ONT file for {sample_id} in base dir: {ont}")
                    all_exist = False
                else:
                    print(f"Found ONT file for {sample_id} in base dir: {ont}")
            
            if row["PACBIO_SRA"] != "None":
                pacbio = os.path.join(base_dir, "PacBio", f"{row['PACBIO_SRA']}.fastq.gz")
                pacbio_pass = os.path.join(base_dir, "PacBio", f"{row['PACBIO_SRA']}_pass.fastq.gz")
                if not (os.path.exists(pacbio) or os.path.exists(pacbio_pass)):
                    print(f"Missing PacBio file for {sample_id} in base dir: {pacbio} or {pacbio_pass}")
                    all_exist = False
                else:
                    found_file = pacbio if os.path.exists(pacbio) else pacbio_pass
                    print(f"Found PacBio file for {sample_id} in base dir: {found_file}")
            
            if row["REF_SEQ_GCA"] != "None":
                ref = os.path.join(base_dir, "RefSeq", f"{row['REF_SEQ_GCA']}.fasta")
                if not os.path.exists(ref):
                    print(f"Missing GCA file for {sample_id} in base dir: {ref}")
                    all_exist = False
                else:
                    print(f"Found GCA file for {sample_id} in base dir: {ref}")
    
    print(f"All files exist check result: {all_exist}")
    return all_exist


# --------------------------------------------------------------
# Process input parameters and resources
# --------------------------------------------------------------
def process_inputs(input_csv, params_dict, percent, cpu, ram, total_cpu, total_ram):
    """Process input parameters and calculate resources for the pipeline.

    Handles input from CSV or command-line parameters, determines input types,
    downloads data if needed, and saves processed metadata.

    Args:
        input_csv (str): Path to input CSV file or 'None'.
        params_dict (dict): Dictionary of pipeline parameters.
        percent (str or float): Percentage of resources to use.
        cpu (str): Number of CPU threads or 'None'.
        ram (str): Amount of RAM in GB or 'None'.
        total_cpu (int): Total available CPU threads.
        total_ram (int): Total available RAM in GB.
    """
    initial_dir = os.getcwd()
    df = None
    output_dir = params_dict.get("output_dir", initial_dir)
    completed_flag = os.path.join(initial_dir, "completed.flag")
    input_csv = os.path.join(initial_dir, os.path.basename(input_csv))
    
    try:
        # Flush stdout to ensure immediate logging
        sys.stdout.flush()
        print(f"DEBUG: Starting process_inputs with input_csv={input_csv}, params_dict={params_dict}")
        if percent is None or percent == "null" or percent == "None":
            percent = float(params_dict.get("percent_resources", 1.0) or 1.0)
        else:
            percent = float(percent)

        no_file = params_dict.get("no_file", "None")
        has_reads = any(param is not None and param != "None" for param in [
            params_dict.get("raw_ont_reads"), params_dict.get("raw_illu_reads_1"),
            params_dict.get("raw_illu_reads_2"), params_dict.get("raw_pacbio_reads")
        ])
        has_ref = params_dict.get("reference_sequence") is not None and params_dict.get("reference_sequence") != "None"

        print(f"Input CSV: {input_csv}, Exists: {os.path.exists(input_csv)}")
        if (not input_csv or input_csv == "null" or not os.path.exists(input_csv)) and not has_reads and not has_ref:
            input_csv = no_file

        if input_csv and input_csv != "None" and os.path.exists(input_csv):
            df = pd.read_csv(input_csv, na_values=["None"])
            df.fillna("None", inplace=True)
            print(f"Loaded CSV with columns: {df.columns.tolist()}")
            print(f"CSV data: {df.to_dict(orient='records')}")
        else:
            sample_dict = {
                "ONT_SRA": [params_dict.get("ont_sra", None)],
                "ONT_RAW_DIR": [params_dict.get("raw_ont_dir", None)],
                "ONT_RAW_READS": [params_dict.get("raw_ont_reads", None)],
                "ILLUMINA_SRA": [params_dict.get("illu_sra", None)],
                "ILLUMINA_RAW_DIR": [params_dict.get("raw_illu_dir", None)],
                "ILLUMINA_RAW_F_READS": [params_dict.get("raw_illu_reads_1", None)],
                "ILLUMINA_RAW_R_READS": [params_dict.get("raw_illu_reads_2", None)],
                "PACBIO_SRA": [params_dict.get("pacbio_sra", None)],
                "PACBIO_RAW_DIR": [params_dict.get("raw_pacbio_dir", None)],
                "PACBIO_RAW_READS": [params_dict.get("raw_pacbio_reads", None)],
                "SAMPLE_ID": [params_dict.get("sample_id", None)],
                "ORGANISM_KINGDOM": [params_dict.get("organism_kingdom", None)],
                "ORGANISM_KARYOTE": [params_dict.get("organism_karyote", None)],
                "COMPLEASM_1": [params_dict.get("compleasm_1", None)],
                "COMPLEASM_2": [params_dict.get("compleasm_2", None)],
                "EST_SIZE": [params_dict.get("estimated_genome_size", None)],
                "REF_SEQ_GCA": [params_dict.get("reference_sequence_gca", None)],
                "REF_SEQ": [params_dict.get("reference_sequence", None)]
            }
            df = pd.DataFrame(sample_dict)
            print(f"Created DataFrame from params: {df.to_dict(orient='records')}")

        if percent != 1.0:
            if cpu == "None" and ram == "None":
                cpu, ram = get_resource_values(percent, total_cpu, total_ram)
            elif cpu != "None" and ram == "None":
                _, ram = get_resource_values(percent, total_cpu, total_ram)
            elif ram != "None" and cpu == "None":
                cpu, _ = get_resource_values(percent, total_cpu, total_ram)

        cpu = str(cpu)
        ram = str(ram)

        df.replace("", "None", inplace=True)
        

        def determine_illumina_input(row):
            if row["ILLUMINA_RAW_F_READS"] != "None" and row["ILLUMINA_RAW_R_READS"] != "None":
                return "concatenated"
            elif row["ILLUMINA_RAW_DIR"] != "None":
                return "directory"
            elif row["ILLUMINA_SRA"] != "None":
                return "sra"
            else:
                return "None"
        df["ILLUMINA_INPUT_TYPE"] = df.apply(determine_illumina_input, axis=1)
        
        def determine_ont_input(row):
            if row["ONT_RAW_READS"] != "None":
                return "concatenated"
            elif row["ONT_RAW_DIR"] != "None":
                return "directory"
            elif row["ONT_SRA"] != "None":
                return "sra"
            else:
                return "None"
        df["ONT_INPUT_TYPE"] = df.apply(determine_ont_input, axis=1)
        
        def determine_pacbio_input(row):
            if row["PACBIO_RAW_READS"] != "None":
                return "concatenated"
            elif row["PACBIO_RAW_DIR"] != "None":
                return "directory"
            elif row["PACBIO_SRA"] != "None":
                return "sra"
            else:
                return "None"
        df["PACBIO_INPUT_TYPE"] = df.apply(determine_pacbio_input, axis=1)
        
        def determine_ref_seq_input(row):
            if row["REF_SEQ"] != "None":
                return "concatenated"
            elif row["REF_SEQ_GCA"] != "None":
                return "gca"
            else:
                return "None"
        df["REF_SEQ_INPUT_TYPE"] = df.apply(determine_ref_seq_input, axis=1)

        files_exist = all_files_exist(df, initial_dir)  # Use initial_dir since outputs stay in work dir
        print(f"DEBUG: Files exist check: {files_exist}, completed_flag exists: {os.path.exists(completed_flag)}")
        if os.path.exists(completed_flag) and files_exist:
            print(f"SKIP: All expected output files exist in {initial_dir} and completed.flag found. Skipping processing.")
        else:
            print(f"Processing required: completed.flag exists={os.path.exists(completed_flag)}, all files exist={files_exist}")
            for index, row in df.iterrows():
                sample_id = row["SAMPLE_ID"]
                base_folder = sample_id.split("-")[0]
                print(f"START: Processing row {index} with SAMPLE_ID={sample_id}, Base Folder={base_folder}")
                print(f"SRA Inputs: ILLUMINA_SRA={row['ILLUMINA_SRA']}, ONT_SRA={row['ONT_SRA']}, PACBIO_SRA={row['PACBIO_SRA']}, REF_SEQ_GCA={row['REF_SEQ_GCA']}")
                if (row["ILLUMINA_SRA"] != "None" or row["ONT_SRA"] != "None" or 
                    row["PACBIO_SRA"] != "None" or row["REF_SEQ_GCA"] != "None"):
                    illu_sra_f, illu_sra_r, ont_sra, pacbio_sra, ref_seq_gca, _ = download_test_data(
                        sample_id, row["ILLUMINA_SRA"], row["ONT_SRA"], row["PACBIO_SRA"], row["REF_SEQ_GCA"], output_dir, cpu
                    )
                    print(f"Download Results: ILLUMINA_F={illu_sra_f}, ILLUMINA_R={illu_sra_r}, ONT={ont_sra}, PACBIO={pacbio_sra}, REF_SEQ={ref_seq_gca}")
                    # Update existing columns with downloaded file paths
                    if illu_sra_f and illu_sra_r:
                        df.at[index, "ILLUMINA_SRA"] = f"{illu_sra_f},{illu_sra_r}"  # Store as comma-separated pair
                    if ont_sra:
                        df.at[index, "ONT_SRA"] = ont_sra
                    if pacbio_sra:
                        df.at[index, "PACBIO_SRA"] = pacbio_sra
                    if ref_seq_gca:
                        df.at[index, "REF_SEQ_GCA"] = ref_seq_gca
                else:
                    print(f"SKIP: No SRA or GCA data to download for {sample_id}")
                if row["ILLUMINA_INPUT_TYPE"] == "directory" and row["ILLUMINA_RAW_DIR"] != "None" and os.path.exists(row["ILLUMINA_RAW_DIR"]):
                    combined_gz_list = illumina_extract_and_check(row["ILLUMINA_RAW_DIR"], sample_id)
                    if combined_gz_list:
                        df.at[index, "ILLUMINA_SRA"] = f"{combined_gz_list[0]},{combined_gz_list[1]}"  # Update ILLUMINA_SRA
                    else:
                        print(f"ERROR: Failed to process ILLUMINA_RAW_DIR for {sample_id}")
                else:
                    print(f"SKIP: No valid ILLUMINA_RAW_DIR for {sample_id}")
                if row["ONT_INPUT_TYPE"] == "directory" and row["ONT_RAW_DIR"] != "None":
                    ont_combined = ont_combine_fastq_gz(row["ONT_RAW_DIR"], total_cpu, sample_id)
                    if ont_combined:
                        df.at[index, "ONT_SRA"] = ont_combined  # Update ONT_SRA
                if row["PACBIO_INPUT_TYPE"] == "directory" and row["PACBIO_RAW_DIR"] != "None":
                    fasta_files = glob.glob(os.path.join(row["PACBIO_RAW_DIR"], "*.fasta"))
                    if fasta_files:
                        largest_fasta = max(fasta_files, key=os.path.getsize)
                        print(f"PASS:\tSelected largest PacBio FASTA file: {largest_fasta}")
                        dest_fasta = os.path.join(output_dir, base_folder, "PacBio", f"{sample_id}_pacbio.fasta")
                        os.makedirs(os.path.dirname(dest_fasta), exist_ok=True)
                        shutil.copy(largest_fasta, dest_fasta)
                        df.at[index, "PACBIO_SRA"] = dest_fasta  # Update PACBIO_SRA
                    else:
                        print(f"ERROR:\tNo .fasta files found in PacBio directory: {row['PACBIO_RAW_DIR']}")
            
            # Save updated CSV to output_dir
            output_csv = "processed_inputs.csv"  # Keep it in work dir
            df.to_csv(output_csv, index=False)
            print(f"PASS: Saved updated CSV to {output_csv}")
        
            # Write resource files to work directory
            with open("cpu.txt", "w") as f:
                f.write(cpu)
            with open("ram.txt", "w") as f:
                f.write(ram)
            with open("completed.flag", "w") as f:
                f.write("Processing completed successfully")
            print(f"PASS: Created resource files in {os.getcwd()}")
            sys.stdout.flush()  # Ensure final log is written

    except Exception as e:
        print(f"ERROR: process_inputs failed: {str(e)}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


# --------------------------------------------------------------
# Download test data from SRA/NCBI
# --------------------------------------------------------------
def download_test_data(SAMPLE_ID, ILLUMINA_SRA, ONT_SRA, PACBIO_SRA, REF_SEQ_GCA, output_dir, cpu):
    """Download sequencing and reference data from SRA or NCBI.

    Fetches Illumina, ONT, PacBio, and reference sequence data, organizing them
    in the output directory.

    Args:
        SAMPLE_ID (str): Sample identifier.
        ILLUMINA_SRA (str): Illumina SRA accession or 'None'.
        ONT_SRA (str): ONT SRA accession or 'None'.
        PACBIO_SRA (str): PacBio SRA accession or 'None'.
        REF_SEQ_GCA (str): Reference sequence GCA accession or 'None'.
        output_dir (str): Base directory for output files.
        cpu (str): Number of CPU threads to use.

    Returns:
        tuple: Paths to downloaded files (Illumina forward, reverse, ONT, PacBio, reference, base directory).
    """
    initial_dir = os.getcwd()
    base_folder = SAMPLE_ID.split("-")[0]  # e.g., "Es_coli"
    base_dir = os.path.join(output_dir, base_folder)  # Use output_dir directly
    print(f"Base directory for raw data: {base_dir}")
    try:
        os.makedirs(base_dir, exist_ok=True)
        print(f"PASS: Base directory created or exists: {base_dir}")
    except Exception as e:
        print(f"ERROR: Failed to create directory {base_dir}: {str(e)}")
        return None, None, None, None, None, base_dir

    illu_sra_f = None
    illu_sra_r = None
    ont_sra = None
    pacbio_sra = None
    ref_seq_gca = None

    # Check and download PacBio
    if pd.notna(PACBIO_SRA) and PACBIO_SRA != "None":
        pacbio_base_dir = os.path.join(base_dir, "PacBio")
        os.makedirs(pacbio_base_dir, exist_ok=True)
        pacbio_base = os.path.join(pacbio_base_dir, f"{PACBIO_SRA}.fastq.gz")
        pacbio_base_pass = os.path.join(pacbio_base_dir, f"{PACBIO_SRA}_pass.fastq.gz")
        if os.path.exists(pacbio_base):
            print(f"SKIP:\tPacBio SRA already exists in base dir: {pacbio_base}")
            pacbio_sra = pacbio_base
        elif os.path.exists(pacbio_base_pass):
            print(f"SKIP:\tPacBio SRA already exists in base dir: {pacbio_base_pass}")
            pacbio_sra = pacbio_base_pass
        else:
            print(f"Downloading PacBio reads for {PACBIO_SRA}...")
            os.chdir(pacbio_base_dir)
            prefetch_cmd = ["prefetch", "--max-size", "100G", "--force", "yes", PACBIO_SRA]
            if run_subprocess_cmd(prefetch_cmd, False) != 0:
                print(f"ERROR:\tFailed to prefetch PacBio SRA {PACBIO_SRA}")
                pacbio_sra = None
            else:
                sra_list = glob.glob(f"{PACBIO_SRA}/*.sra")
                if not sra_list:
                    print(f"ERROR:\tNo SRA files found after prefetch in {pacbio_base_dir}/{PACBIO_SRA}")
                    pacbio_sra = None
                else:
                    sra_file = sra_list[0]
                    pacbio_cmd = f"fasterq-dump -e {cpu} -O {pacbio_base_dir} {PACBIO_SRA} && pigz -p {cpu} {PACBIO_SRA}*.fastq"
                    if run_subprocess_cmd(pacbio_cmd, True) != 0:
                        print(f"ERROR:\tFailed to dump PacBio SRA {PACBIO_SRA}")
                        pacbio_sra = None
                    else:
                        if os.path.exists(pacbio_base):
                            print(f"PASS:\tPacBio SRA processed: {pacbio_base}")
                            pacbio_sra = pacbio_base
                        elif os.path.exists(f"{PACBIO_SRA}_pass.fastq.gz"):
                            pacbio_sra = pacbio_base_pass
                            print(f"PASS:\tPacBio SRA processed: {pacbio_sra}")
                        else:
                            print(f"ERROR:\tExpected PacBio SRA file not found after dump")
                            pacbio_sra = None
                    if os.path.exists(sra_file):
                        os.remove(sra_file)
                    sra_dir = os.path.dirname(sra_file)
                    if os.path.exists(sra_dir):
                        shutil.rmtree(sra_dir, ignore_errors=True)
            os.chdir(initial_dir)

    # Check and download Illumina
    if pd.notna(ILLUMINA_SRA) and ILLUMINA_SRA != "None":
        illumina_base_dir = os.path.join(base_dir, "Illumina")
        os.makedirs(illumina_base_dir, exist_ok=True)
        illu_base_f = os.path.join(illumina_base_dir, f"{ILLUMINA_SRA}_1.fastq.gz")
        illu_base_r = os.path.join(illumina_base_dir, f"{ILLUMINA_SRA}_2.fastq.gz")
        if os.path.exists(illu_base_f) and os.path.exists(illu_base_r):
            print(f"SKIP:\tIllumina SRAs already exist in base dir: {illu_base_f}; {illu_base_r}")
            illu_sra_f = illu_base_f
            illu_sra_r = illu_base_r
        else:
            print(f"Downloading Illumina Reads for: {ILLUMINA_SRA}...")
            os.chdir(illumina_base_dir)
            illu_cmd = f"prefetch --force yes {ILLUMINA_SRA} && fasterq-dump -e {cpu} -O {illumina_base_dir} {ILLUMINA_SRA} && pigz -p {cpu} {illumina_base_dir}/{ILLUMINA_SRA}*.fastq"
            if run_subprocess_cmd(illu_cmd, True) != 0:
                print(f"ERROR:\tFailed to download Illumina SRA {ILLUMINA_SRA}")
                illu_sra_f = None
                illu_sra_r = None
            else:
                print(f"PASS:\tIllumina SRA processed: {illu_base_f}, {illu_base_r}")
                illu_sra_f = illu_base_f
                illu_sra_r = illu_base_r
            if os.path.exists(ILLUMINA_SRA):
                shutil.rmtree(ILLUMINA_SRA, ignore_errors=True)
            os.chdir(initial_dir)

    # Check and download ONT
    if pd.notna(ONT_SRA) and ONT_SRA != "None":
        ont_base_dir = os.path.join(base_dir, "ONT")
        os.makedirs(ont_base_dir, exist_ok=True)
        ont_base = os.path.join(ont_base_dir, f"{ONT_SRA}.fastq.gz")
        if os.path.exists(ont_base):
            print(f"SKIP:\tONT SRA already exists in base dir: {ont_base}")
            ont_sra = ont_base
        else:
            print(f"Downloading ONT Reads for: {ONT_SRA}...")
            os.chdir(ont_base_dir)
            ont_cmd = f"prefetch --force yes {ONT_SRA} && fasterq-dump -e {cpu} -O {ont_base_dir} {ONT_SRA} && pigz -p {cpu} {ont_base_dir}/{ONT_SRA}*.fastq"
            if run_subprocess_cmd(ont_cmd, True) != 0:
                print(f"ERROR:\tFailed to download ONT SRA {ONT_SRA}")
                ont_sra = None
            else:
                if os.path.exists(ont_base):
                    print(f"PASS:\tONT SRA processed: {ont_base}")
                    ont_sra = ont_base
                else:
                    print(f"ERROR:\tONT SRA file not found after processing: {ont_base}")
                    ont_sra = None
                if os.path.exists(ONT_SRA):
                    shutil.rmtree(ONT_SRA, ignore_errors=True)
            os.chdir(initial_dir)

    # Check and download Reference Sequence
    if pd.notna(REF_SEQ_GCA) and REF_SEQ_GCA != "None":
        ref_base_dir = os.path.join(base_dir, "RefSeq")
        os.makedirs(ref_base_dir, exist_ok=True)
        ref_base = os.path.join(ref_base_dir, f"{REF_SEQ_GCA}.fasta")
        if os.path.exists(ref_base):
            print(f"SKIP:\tREF_SEQ_GCA already exists in base dir: {ref_base}")
            ref_seq_gca = ref_base
        elif "." not in REF_SEQ_GCA:
            print(f"ERROR:\tReference Sequence GCA requires version number: {REF_SEQ_GCA} has no '.#'")
            ref_seq_gca = None
        else:
            print(f"Downloading Reference Sequence Assembly for: {REF_SEQ_GCA}")
            os.chdir(ref_base_dir)
            datasets_cmd = ["datasets", "download", "genome", "accession", REF_SEQ_GCA, "--include", "genome"]
            if run_subprocess_cmd(datasets_cmd, False) != 0:
                print(f"ERROR:\tFailed to download GCA {REF_SEQ_GCA}")
                ref_seq_gca = None
            else:
                zip_file = "ncbi_dataset.zip"
                if os.path.exists(zip_file):
                    with zipfile.ZipFile(zip_file, 'r') as zip_ref:
                        zip_ref.extractall(ref_base_dir)
                    os.remove(zip_file)
                    genomic_files = glob.glob(os.path.join(ref_base_dir, f"ncbi_dataset/data/{REF_SEQ_GCA}/*_genomic.fna"))
                    if genomic_files:
                        ref_seq_gca = genomic_files[0]
                        print(f"PASS:\tSuccessfully extracted the GCA to: {ref_seq_gca}")
                        shutil.move(ref_seq_gca, ref_base)
                        print(f"PASS:\tSuccessfully moved and renamed the GCA to: {ref_base}")
                        ref_seq_gca = ref_base
                    else:
                        print(f"ERROR:\tNo genomic files found after extraction for {REF_SEQ_GCA}")
                        ref_seq_gca = None
                else:
                    print(f"ERROR:\tZIP file not found after download for {REF_SEQ_GCA}")
                    ref_seq_gca = None
            os.chdir(initial_dir)

    os.chdir(initial_dir)
    return illu_sra_f, illu_sra_r, ont_sra, pacbio_sra, ref_seq_gca, base_dir


# --------------------------------------------------------------
# Combine and verify Illumina reads
# --------------------------------------------------------------
def illumina_extract_and_check(folder_name, SAMPLE_ID):
    """Combine paired-end Illumina reads after MD5 verification.

    Verifies MD5 checksums, concatenates forward and reverse FASTQ files,
    and compresses the results.

    Args:
        folder_name (str): Directory containing Illumina FASTQ.GZ files and MD5.txt.
        SAMPLE_ID (str): Sample identifier for naming output files.

    Returns:
        list or None: Paths to compressed combined forward and reverse files, or None if failed.
    """
    print(f"Running MD5 Checksum Analysis on Raw Illumina FASTQ.GZ files in {folder_name}...")
    illumina_df = pd.DataFrame(columns=["MD5", "Filename"])
    base_folder = SAMPLE_ID.split("-")[0]  # e.g., "Es_coli"
    illumina_dir = os.path.join(os.getcwd(), base_folder, "Illumina")  # Use current dir, not hardcoded EGAP_Test_Data
    os.makedirs(illumina_dir, exist_ok=True)
    combined_1_file = os.path.join(illumina_dir, f"{SAMPLE_ID}_combined_1.fq")
    combined_2_file = os.path.join(illumina_dir, f"{SAMPLE_ID}_combined_2.fq")
    combined_list = [combined_1_file, combined_2_file]
    combined_gz_list = [f"{combined_fq}.gz" for combined_fq in combined_list]
    
    if not os.path.isfile(combined_gz_list[0]) or not os.path.isfile(combined_gz_list[1]):
        if not os.path.isfile(combined_1_file) or not os.path.isfile(combined_2_file):
            md5_check(folder_name, illumina_df)
            raw_1_list = []
            raw_2_list = []
            for filename in os.listdir(folder_name):
                if "_1.fastq.gz" in filename:
                    raw_1_list.append(os.path.join(folder_name, filename))
                elif "_2.fastq.gz" in filename:
                    raw_2_list.append(os.path.join(folder_name, filename))
            if not raw_1_list or not raw_2_list:
                print(f"ERROR:\tNo paired Illumina files found in {folder_name}")
                return None
            fwd_cat_cmd = f"cat {' '.join(raw_1_list)} > {combined_1_file}"
            _ = run_subprocess_cmd(fwd_cat_cmd, shell_check=True)
            rev_cat_cmd = f"cat {' '.join(raw_2_list)} > {combined_2_file}"
            _ = run_subprocess_cmd(rev_cat_cmd, shell_check=True)
        else:
            print(f"SKIP:\tCombined FASTQ files already exist: {combined_list[0]}; {combined_list[1]}.")
        for combined_fq in combined_list:
            gzip_cmd = ["pigz", "-p", str(cpu), combined_fq]  # Fixed -e to -p for pigz
            _ = run_subprocess_cmd(gzip_cmd, shell_check=False)
    else:
        print(f"SKIP:\tGzipped Combined FASTQ files already exist: {combined_gz_list[0]}; {combined_gz_list[1]}.")
    return combined_gz_list if os.path.exists(combined_gz_list[0]) and os.path.exists(combined_gz_list[1]) else None


# --------------------------------------------------------------
# Combine ONT FASTQ.GZ files
# --------------------------------------------------------------
def ont_combine_fastq_gz(ONT_FOLDER, CPU_THREADS, SAMPLE_ID):
    """Combine multiple ONT FASTQ.GZ files into a single FASTQ file.

    Merges ONT FASTQ.GZ files from a directory into a single uncompressed FASTQ file.

    Args:
        ONT_FOLDER (str): Directory containing ONT FASTQ.GZ files.
        CPU_THREADS (int): Number of CPU threads (not used in current implementation).
        SAMPLE_ID (str): Sample identifier for naming output file.

    Returns:
        str or None: Path to combined FASTQ file, or None if failed.
    """
    print("Combining ONT FASTQ.GZ files...")
    base_folder = SAMPLE_ID.split("-")[0]  # e.g., "Es_coli"
    ont_dir = os.path.join(os.getcwd(), base_folder, "ONT")  # Use current dir
    os.makedirs(ont_dir, exist_ok=True)
    ont_raw_data_dir = next((subdir for subdir in glob.glob(os.path.join(ONT_FOLDER, "*"))
                             if os.path.isdir(subdir) and any(file.endswith(".gz") for file in os.listdir(subdir))), None)
    if ont_raw_data_dir is None:
        print(f"NOTE:\tNo directory containing '.gz' files found within '{ONT_FOLDER}'...Attempting with ONT_FOLDER '{ONT_FOLDER}'")
        ont_raw_data_dir = ONT_FOLDER
    combined_ont_fastq_path = os.path.join(ont_dir, f"{SAMPLE_ID}_ont_combined.fastq")
    raw_file_list = glob.glob(os.path.join(ont_raw_data_dir, "*.fastq.gz"))
    if os.path.isfile(combined_ont_fastq_path):
        print(f"NOTE:\tSkipping extraction & combination: Combined fastq file already exists: {combined_ont_fastq_path}")
        return combined_ont_fastq_path
    with open(combined_ont_fastq_path, "w") as combined_file:
        for filename in raw_file_list:
            with gzip.open(filename, "rt") as gz_file:
                for record in SeqIO.parse(gz_file, "fastq"):
                    try:
                        SeqIO.write(record, combined_file, "fastq")
                    except Exception as e:
                        print(f"ERROR:\tFound in FASTQ record: {e}")
                        raise e
    print(f"PASS: Successfully created combined fastq file: {combined_ont_fastq_path}")
    return combined_ont_fastq_path if os.path.exists(combined_ont_fastq_path) else None


# --------------------------------------------------------------
# Verify MD5 checksums for Illumina files
# --------------------------------------------------------------
def md5_check(folder_name, illumina_df):
    """Verify MD5 checksums for Illumina files in the specified folder.

    Compares computed MD5 checksums against those listed in MD5.txt.

    Args:
        folder_name (str): Directory containing Illumina files and MD5.txt.
        illumina_df (pandas.DataFrame): DataFrame to store MD5 and filename data.
    """
    md5_file = os.path.join(folder_name, "MD5.txt")
    if not os.path.exists(md5_file):
        print(f"WARNING: MD5.txt not found in {folder_name}. Skipping MD5 check.")
        return
    with open(md5_file, "r") as f:
        for line in f:
            md5, filename = line.strip().split()
            illumina_df = illumina_df.append({"MD5": md5, "Filename": filename}, ignore_index=True)
    
    for index, row in illumina_df.iterrows():
        file_path = os.path.join(folder_name, row["Filename"])
        if os.path.exists(file_path):
            with open(file_path, "rb") as f:
                file_hash = hashlib.md5(f.read()).hexdigest()
            if file_hash == row["MD5"]:
                print(f"PASS: MD5 check passed for {row['Filename']}")
            else:
                print(f"ERROR: MD5 check failed for {row['Filename']}. Expected {row['MD5']}, got {file_hash}")
        else:
            print(f"ERROR: File not found for MD5 check: {file_path}")


# --------------------------------------------------------------
# Search filesystem for a file
# --------------------------------------------------------------
def find_file(filename, folder=None):
    """Search the filesystem for a file by name.

    Walks the filesystem from a specified or default directory, skipping certain paths.

    Args:
        filename (str): Name of the file to search for.
        folder (str, optional): Starting directory for search. Defaults to 'EGAP_Test_Data'.

    Returns:
        str or None: Absolute path to the file if found, else None.
    """
    print(f"Looking for {filename}")
    if pd.notna(folder):
        root_directory = folder
    else:
        root_directory = os.path.join("EGAP_Test_Data")
    for root, dirs, files in os.walk(root_directory):
        if "$RECYCLE.BIN" in root or "ncbi" in dirs:
            continue
        if filename in files:
            return os.path.join(root, filename)
    return None


# --------------------------------------------------------------
# Execute a subprocess command
# --------------------------------------------------------------
def run_subprocess_cmd(cmd_list, shell_check):
    """Run a subprocess command and stream its output.

    Executes the command and logs its success or failure.

    Args:
        cmd_list (str or list): Command to execute, as a string or list of arguments.
        shell_check (bool): If True, execute the command through the shell.

    Returns:
        int: The subprocess return code.
    """
    if isinstance(cmd_list, str):
        cmd_str = cmd_list
        print(f"CMD:\t{cmd_str}")
        process = subprocess.Popen(cmd_list, shell=shell_check, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True)
    else:
        cmd_str = ' '.join(cmd_list)
        print(f"CMD:\t{cmd_str}")
        process = subprocess.Popen(cmd_list, shell=shell_check, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True)
    for line in process.stdout:
        print(line, end="")
    process.wait()
    if process.returncode != 0:
        print(f"NOTE:\tCommand failed with return code {process.returncode}")
    else:
        print(f"PASS:\tSuccessfully processed command: {cmd_str}")
    return process.returncode


if __name__ == "__main__":
    try:
        if len(sys.argv) != 8:
            print("Usage: python process_inputs.py <input_csv> <params_dict> <percent> <cpu> <ram> <total_cpu> <total_ram>", 
                  file=sys.stderr)
            sys.exit(1)
        
        input_csv = sys.argv[1]
        params_dict = json.loads(sys.argv[2])
        percent = sys.argv[3]
        cpu = sys.argv[4]
        ram = sys.argv[5]
        total_cpu = int(sys.argv[6])
        total_ram = int(sys.argv[7])
        
        process_inputs(input_csv, params_dict, percent, cpu, ram, total_cpu, total_ram)
    except Exception as e:
        print(f"CRITICAL ERROR: Script failed to start: {str(e)}", file=sys.stderr)
        sys.exit(1)