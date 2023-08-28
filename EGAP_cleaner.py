# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 18:18:22 2023

@author: ian.michael.bollinger@gmail.com with the help of ChatGPT 4.0

Command Line Example:
    python EGAP_cleaner.py --dirty_fasta /path/to/dirty.fasta --output_dir /path/to/folder --organism_kingdom STRING

The --organism_kingdom must be from the following: Archaea, Bacteria, Fauna, Flora, Funga, or Protista
"""
import os, subprocess, glob, random, tempfile, shutil
import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from tqdm import tqdm
from log_print import log_print, generate_log_file
from check_tools import get_env_dir, move_file_up

# Function to create BLAST databases from .gbff.gz and .fna files
def create_BLAST_db(database_dir, log_file):
    """
    Create BLAST databases from .gbff.gz and .fna files.

    Args:
        database_dir (str): Path to the combined ONT fastq reads.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.

    Returns:
        _ (list): List of organisms in database.
        db_name (str): Path to database directory.
    """  
    # Find all .gbff.gz and .fna files in the directory
    gbff_files = glob.glob(database_dir + "/*.gbff.gz")
    fna_files = glob.glob(database_dir + "/*.fna")
    all_files = gbff_files + fna_files
    if not all_files:
        log_print("ERROR:\tNo .gbff.gz or .fna files found in the directory", log_file)
        return [None], None
    fasta_files = []
    
    # Convert GBFF to FASTA and add .fna files to the list
    for input_filename in tqdm(all_files, desc="Processing files"):
        if input_filename.endswith('.gbff.gz'):
            output_filename = input_filename.replace('.gbff.gz','.fasta')
            with gzip.open(input_filename, "rt") as handle:
                count = SeqIO.write(SeqIO.parse(handle, "genbank"), output_filename, "fasta")
            log_print(f"NOTE:\tConverted {count} records from {input_filename}", log_file)
            fasta_files.append(output_filename)
        else:
            fasta_files.append(input_filename)
    
    base_fasta_filename = os.path.basename(fasta_files[0]).replace('.fasta', '')
    db_name = f"{database_dir}/{base_fasta_filename}_blast_db"        
    
    # Check if there are already BLAST databases in the directory
    files_in_directory = [os.path.join(database_dir, file) for file in os.listdir(database_dir)]
    filtered_files = [file for file in files_in_directory if "blast_db" in file]
    if filtered_files:
        log_print("PASS:\tBLAST databases already exist in the directory", log_file)
        return [file.replace('.fasta', '_blast_db') for file in fasta_files], db_name

    # Create a BLAST database from each FASTA file
    for fasta_file in tqdm(fasta_files, desc="Creating BLAST databases"):
        log_print(f"Attempting to generate BLASTn database from: {fasta_file}", log_file)
        
        # Generate BLAST database command
        blast_db_cmd = ["makeblastdb",
                        "-in",
                        fasta_file,
                        "-dbtype",
                        "nucl",
                        "-out",
                        db_name]

        # Run the command
        log_print(f"CMD:\t{' '.join(blast_db_cmd)}", log_file)
        blast_db_result = subprocess.run(blast_db_cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)    
        if blast_db_result.returncode != 0:
            log_print(f"ERROR:\t{blast_db_result.stderr}", log_file)
            return [file.replace('.fasta', '_blast_db') for file in fasta_files], db_name

    return [file.replace('.fasta', '_blast_db') for file in fasta_files], db_name

# Function to take in a sequence and return three (3) representative chunks
def chunk_seq(record):
    """
    Create three (3) representative chunks of a given sequence

    Args:
        record (str): A sequence that is to be chunked.

    Returns:
        start (str): Chunk in the first third (start) of the sequence.
        middle (str): Chunk in the second third (middle) of the sequence.
        end (str): Chunk in the last third (end) of the sequence.
    """
    # Seed random number generator to get consistent chunks
    rnd = random.Random(25)
    
    # Get the length of the sequence
    length = len(record)
    
    # If the sequence is shorter than 1500 base pairs, return the whole sequence three times
    if length < 1500:
        return record.seq, record.seq, record.seq

    # Define the boundaries of the three equal parts (start, middle, end) of the sequence
    one_end = int(length/3)
    start = record.seq[0:one_end]
    middle = record.seq[one_end:2*one_end]
    end = record.seq[2*one_end:]

    # Define a helper function to handle each chunk
    def handle_chunk(chunk):
        """
        Chunk handler, manages chunk size to at least 1500 base pairs.

        Args:
            chunk (str): A randome size chunk from a sequence.

        Returns:
            str: Final size-managed chunk.
        """
        # If the chunk is longer than 1500 base pairs, select a random 1500 base pair subsequence
        if len(chunk) > 1500:
            start = rnd.randint(0, len(chunk) - 1500)  # using rnd instead of random
            return chunk[start:start+1500]
        else:
            # If the chunk is shorter than 1500 base pairs, return the entire chunk
            return chunk

    # Use the handle_chunk function to possibly shorten each of the three chunks
    start = handle_chunk(start)
    middle = handle_chunk(middle)
    end = handle_chunk(end)
    
    # Return the potentially shortened chunks
    return start, middle, end

# Function to take in a squence and BLASTn it against a given dataframe
def run_blastn(sequence, db_name):
    """
    Takes in a sequence chunk and BLASTn it against a given dataframe.

    Args:
        sequence (str): A randome size chunk from a sequence.
        db_name (str): Path to a given BLAST database.

    Returns:
        alignment_length (int): The length of the alignment.
        percent_id (float): The alignment's percent ID.
        percent_qcoverage (float): Querey coverage of the alignment as percent.
    """
    # Create a temporary file and write the sequence to it
    with tempfile.NamedTemporaryFile(delete=False, mode='w') as temp_file:
        temp_file.write(f">query\n{sequence}")
        temp_file_name = temp_file.name
    
    # Define the blastn command
    blastn_cline = NcbiblastnCommandline(query = temp_file_name,
                                         db = db_name, evalue = 0.001,
                                         outfmt = 6, out = "blastn_out.txt")
    
    # Run the command
    stdout, stderr = blastn_cline()
    
    # Read the output file into a pandas DataFrame
    blast_results = pd.read_csv("blastn_out.txt", sep="\t", names=["qseqid", "sseqid", "pident", "qlen", "qstart", "qend", "length", "gaps"])
    
    # Clean up the temporary file
    os.unlink(temp_file_name)
    
    # Check if the DataFrame is empty
    if blast_results.empty:
        return None, None, None
        
    # Pull out data based on column names
    alignment_length = blast_results['length'].iloc[0]
    percent_id = blast_results['pident'].iloc[0]
    query_len = blast_results['qlen'].iloc[0]
    qstart = blast_results['qstart'].iloc[0]
    qend = blast_results['qend'].iloc[0]
    gaps = blast_results['gaps'].iloc[0]
    
    # Calculate query coverage
    percent_qcoverage = (alignment_length / query_len) * 100

    return alignment_length, percent_id, percent_qcoverage

# Function to clean an assembly with given database
def clean_dirty_fasta(dirty_fasta, OUTPUT_DIR, current_organism_kingdom, log_file):
    """
    Clean an FASTA assembly with built-in databases.

    Args:
        dirty_fasta (str): Path to an unfiltered fasta.
        OUTPUT_DIR (str): Path to the output folder.
        current_organism_kingdom (str): Single word Kingdom Description: Archaea, Bacteria, Fauna, Flora, Funga, or Protista.
        log_file (obj): Path to the log file where the analysis progress and results will be stored.

    Returns:
        cleaned_assembly (str): Path to the cleaned FASTA file.
    """
    log_print(f"Cleaning all non-{current_organism_kingdom} contigs from {dirty_fasta}...", log_file)
    if dirty_fasta is not None:
        log_print(f"NOTE:\tFound FASTA file: {dirty_fasta}", log_file)
    else:
        log_print("ERROR:\tNo FASTA file found in the directory", log_file)
        return None, None

    cleaned_assembly = dirty_fasta
    

    
    # Set-up User Install Path for Database Directory List generation
    output_split = OUTPUT_DIR.split('/')
    install_path = f"/{output_split[1]}/{output_split[2]}"
    
    # Initialize an empty dataframe for cumulative removed sequences
    cumulative_removed_df = pd.DataFrame()
    
    # Generate Database Directories List
    database_directories = [f"{install_path}/EGAP/Assembled_Databases/Archaea",
                            f"{install_path}/EGAP/Assembled_Databases/Bacteria",
                            f"{install_path}/EGAP/Assembled_Databases/Fauna",
                            f"{install_path}/EGAP/Assembled_Databases/Flora",
                            f"{install_path}/EGAP/Assembled_Databases/Funga",
                            f"{install_path}/EGAP/Assembled_Databases/Protista"] 
    
    # Rename final cleaned FASTA & Generate new path to deposit file in
    last_db = database_directories[-1:]
    last_db_name = last_db[0].split('/')[-1]    
    if last_db_name in cleaned_assembly:
        final_ont_assembly = cleaned_assembly.replace(f'_{last_db_name}','')
        destination_path = os.path.join(OUTPUT_DIR, os.path.basename(final_ont_assembly))
    else:
        final_ont_assembly = cleaned_assembly
        destination_path = os.path.join(OUTPUT_DIR, os.path.basename(final_ont_assembly).replace('.fasta','_filtered.fasta'))
    
    print(f'Cleaned Assembly: {cleaned_assembly}')
    print(f'Final ONT Assembly: {final_ont_assembly}')
    print(f'Destination Path Assembly: {destination_path}')

    
    if os.path.isfile(destination_path):
        log_print(f"PASS:\tSkipping Flye Cleaning: {destination_path} already exists ", log_file)
        
        # update the path to the cleaned assembly after each iteration
        cleaned_assembly = destination_path
    
    
    else:
        for index, database_dir in enumerate(database_directories):        
            if current_organism_kingdom in database_dir:
                log_print(f'PASS:\tSkipping {current_organism_kingdom} because the main data for the organism is of that Kingdom',log_file)
            else:
                # Generate paths for output files
                database_name = database_dir.split('/')[-1]
                if '_filtered' in cleaned_assembly:
                    cleaned_fasta = f"{cleaned_assembly.split('_filtered')[0]}_filtered_{database_name}.fasta"
                else:
                    cleaned_fasta = cleaned_assembly.replace('.fasta', f'_filtered_{database_name}.fasta')
                
                # Check if the filtered csv file already exists
                if os.path.isfile(cleaned_fasta):
                    log_print(f"PASS:\tSkipping {database_name} Flye Cleaning: {cleaned_fasta} already exists ", log_file)
                    # update the path to the cleaned assembly after each iteration
                    cleaned_assembly = cleaned_fasta
                else:
                    # Create the BLAST datbase for searching
                    _, db_name = create_BLAST_db(database_dir, log_file)
                
                    # Generate dataframe to store contig information
                    columns = ['contig_name', 'start_chunk', 'mid_chunk', 'end_chunk', 'start_percent_id', 'mid_percent_id', 'end_percent_id', 'start_percent_qcoverage', 'mid_percent_qcoverage', 'end_percent_qcoverage', 'start_alignment_length', 'mid_alignment_length', 'end_alignment_length']
                    df = pd.DataFrame(columns=columns)
                    fasta_sequences = SeqIO.parse(open(cleaned_assembly),'fasta')
                    contig_number = 0
                    for fasta in tqdm(fasta_sequences, desc="Preparing contigs"):
                        name, sequence = fasta.id, str(fasta.seq)
                        try:
                            start_chunk, mid_chunk, end_chunk = chunk_seq(fasta)
                            df.loc[contig_number] = [name, start_chunk, mid_chunk, end_chunk, None, None, None, None, None, None, None, None, None]
                            contig_number += 1
                        except ValueError as e:
                            log_print(f"PASS:\tSkipping contig {name} due to error: {e}", log_file)
                
                    csv_filename = f"{cleaned_assembly.replace('.fasta','.csv')}"
                    df.to_csv(csv_filename, index=False)
    
                    for database_dir, database_cmd in [(database_dir, 'blast_db')]:
                        for index, row in tqdm(df.iterrows(), total=df.shape[0], desc=f"Running local BLAST on {database_dir.split('/')[-1]}"):
                            for chunk in ['start_chunk', 'mid_chunk', 'end_chunk']:
                                chunk_sequence = row[chunk]
                                alignment_length, percent_id, percent_qcoverage = run_blastn(chunk_sequence, db_name)
                                df.at[index, chunk.replace('chunk', 'percent_id')] = percent_id
                                df.at[index, chunk.replace('chunk', 'percent_qcoverage')] = percent_qcoverage
                                df.at[index, chunk.replace('chunk', 'alignment_length')] = alignment_length
                
                        for col in ['start_percent_id', 'mid_percent_id', 'end_percent_id', 'start_percent_qcoverage', 'mid_percent_qcoverage', 'end_percent_qcoverage', 'start_alignment_length', 'mid_alignment_length', 'end_alignment_length']:
                            df[col] = pd.to_numeric(df[col], errors='coerce')
                
                        filtered_df = df[~((df['start_percent_id'] > 95.000) & (df['start_percent_qcoverage'] > 50.000) |
                                            (df['mid_percent_id'] > 95.000) & (df['mid_percent_qcoverage'] > 50.000) |
                                            (df['end_percent_id'] > 95.000) & (df['end_percent_qcoverage'] > 50.000) |
                                            ((df['start_percent_id'] == 100.000) & (df['start_alignment_length'] > 100)) |
                                            ((df['mid_percent_id'] == 100.000) & (df['mid_alignment_length'] > 100)) |
                                            ((df['end_percent_id'] == 100.000) & (df['end_alignment_length'] > 100)))]
                
                        # Generate the Removed Dataframe
                        removed_df = df[~df.index.isin(filtered_df.index)]
                        
                        # Add a column to the removed_df to store the database name
                        removed_df['Database'] = database_dir.split('/')[-1]
                        
                        # Append the removed_df to cumulative_removed_df
                        cumulative_removed_df = pd.concat([cumulative_removed_df, removed_df])                    
                        df = filtered_df
                        
                    contig_list = df['contig_name'].tolist()
                    log_print(f'NOTE:\tKeeping {round((len(contig_list)/contig_number)*100,0)}% of Contigs', log_file)
                
                    with open(cleaned_assembly, 'r') as original, open(cleaned_fasta, 'w') as filtered:
                        fasta_sequences = SeqIO.parse(original, 'fasta')
                        for fasta in fasta_sequences:
                            if fasta.id in contig_list:
                                SeqIO.write(fasta, filtered, 'fasta')
                    
                    log_print(f"PASS:\tFiltered fasta file saved as {cleaned_fasta}", log_file)
                
                # update the path to the cleaned assembly after each iteration
                cleaned_assembly = cleaned_fasta  
   
    # Save the cumulative_removed_df to a single CSV after the loop
    cumulative_removed_csv = os.path.join(OUTPUT_DIR, "cumulative_removed_sequences.csv")
    if not cumulative_removed_df.empty:
        cumulative_removed_df.to_csv(cumulative_removed_csv, index=False)
        cumulative_removed_csv = move_file_up(cumulative_removed_csv, log_file, move_bool = True)
        log_print(f"NOTE:\tAll removed sequences saved in {cumulative_removed_csv}", log_file)
    else:
        log_print(f"PASS:\tNo sequences were removed across all databases", log_file)
         
    # Attempt to Move cleaned assembly up one folder
    try:
        # Copy the cleaned_assembly to destination_path
        shutil.copy2(cleaned_assembly, destination_path)
        log_print(f"PASS:\tCopied {cleaned_assembly} to {destination_path}", log_file)
    except Exception as e:
        log_print(f"NOTE:\tUnable to copy and rename cleaned ONT assembly file. Reason: {e}", log_file)
    cleaned_assembly = destination_path
    
    return cleaned_assembly, cumulative_removed_csv

## Debuging Main Space & Example
if __name__ == "__main__":
    print('EGAP Contig Cleaner')    
    # Get working environment information
    environment_dir = get_env_dir()
    
    # Argument Parsing
    parser = argparse.ArgumentParser(description='Tools to clean a contanminated FASTA file')
    
    # Default values
    default_file = f'{environment_dir}/Entheome/Ps_aff_hopii/MODULAR_TEST/ONT_MinION/B1_3_ont_flye_output/B1_3_ont_flye.fasta'
    default_folder = f'{environment_dir}/Entheome/Ps_aff_hopii/MODULAR_TEST/ONT_MinION/'
    default_organism_kingdom = 'Funga'
    
    # Add arguments with default values
    parser.add_argument('--dirty_fasta', default = default_file,
                        help = f'Path to the Contaminated FASTA file. (default: {default_file})')
    parser.add_argument('--output_dir', default = default_folder,
                        help = f'Path to the desired output folder. (default: {default_folder})')
    parser.add_argument('--organism_kingdom',default = default_organism_kingdom,
                        help = f'Kingdom the current organism data belongs to. (default: {default_organism_kingdom})')
    
    # Parse the arguments
    args = parser.parse_args()
    DIRTY_FASTA = args.dirty_fasta
    OUTPUT_DIR = args.output_dir
    CURRENT_ORGANISM_KINGDOM = args.organism_kingdom
    
    # Generate Main Logfile
    debug_log = f'{OUTPUT_DIR}Contig_Cleaner_log.tsv'
    log_file = generate_log_file(debug_log, use_numerical_suffix=False)
    
    # Run main FASTA Cleaning function
    cleaned_assembly = clean_dirty_fasta(DIRTY_FASTA, OUTPUT_DIR, CURRENT_ORGANISM_KINGDOM, log_file)
