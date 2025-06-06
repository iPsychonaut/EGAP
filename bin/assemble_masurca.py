#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
assemble_masurca.py

Updated on Sat Mar 29 2025

This script runs MaSuRCA assembly with Illumina and optional long reads.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import os, sys, shutil, re
import pandas as pd
from utilities import run_subprocess_cmd, get_current_row_data
from qc_assessment import qc_assessment


# --------------------------------------------------------------
# Locate MaSuRCA CA folder
# --------------------------------------------------------------
def find_ca_folder(current_work_dir):
    """Find the MaSuRCA CA folder in the current working directory.

    Scans subdirectories for a folder starting with 'CA', defaulting to 'CA' if none is found.

    Returns:
        str: Path to the CA folder.
    """
    subfolders = [f.path for f in os.scandir(current_work_dir) if f.is_dir()]
    ca_folder = os.path.join(current_work_dir, "CA")
    for folder in subfolders:
        if os.path.basename(folder).startswith("CA"):
            ca_folder = folder
            break
    return ca_folder


# --------------------------------------------------------------
# Parse BBMerge insert size statistics
# --------------------------------------------------------------
def parse_bbmerge_output(insert_size_histogram_txt):
    """Extract insert size statistics from BBMerge histogram file.

    Parses the histogram to retrieve average insert size and standard deviation.

    Args:
        insert_size_histogram_txt (str): Path to the BBMerge histogram file.

    Returns:
        tuple: (average insert size, standard deviation).

    Raises:
        ValueError: If mean or standard deviation cannot be parsed.
    """
    print(f"Processing insert size histogram: {insert_size_histogram_txt}...")
    avg_insert = None
    std_dev = None
    with open(insert_size_histogram_txt, "r") as file:
        for line in file:
            if "#Mean\t" in line:
                avg_insert = round(float(line.replace("#Mean\t", "").replace("\n", "")), 0)
            if "#STDev\t" in line:
                std_dev = round(float(line.replace("#STDev\t", "").replace("\n", "")), 0)
    if avg_insert is None or std_dev is None:
        raise ValueError("Could not find average insert size and/or standard deviation in the output.")
    return avg_insert, std_dev


# --------------------------------------------------------------
# Compute insert size statistics with BBMerge
# --------------------------------------------------------------
def bbmap_stats(input_folder, reads_list, cpu_threads):
    """Compute insert size statistics using BBMerge.

    Runs BBMerge to generate a histogram or parses an existing one to extract
    insert size statistics.

    Args:
        input_folder (str): Directory for BBMerge output files.
        reads_list (list): List of FASTQ file paths (paired-end or with extra read).

    Returns:
        tuple: (average insert size, standard deviation).
    """
    bbmap_out_path = f"{input_folder}/bbmap_data.fq"
    insert_size_histogram_txt = f"{input_folder}/insert_size_histogram.txt"
    print(f"DEBUG - bbmap_out_path - {bbmap_out_path}")
    avg_insert = 251
    std_dev = 30
    if os.path.isfile(insert_size_histogram_txt):
        print(f"SKIP:\tbbmap insert size histogram output already exists: {insert_size_histogram_txt}")
        avg_insert, std_dev = parse_bbmerge_output(insert_size_histogram_txt)
    else:
        print("Processing fastq files for bbmap stats...")
        # bbmerge_path = find_file("bbmerge.sh",
        #                          folder = os.environ["CONDA_PREFIX"] + "/bin",
        #                          cpu_threads = cpu_threads)
        bbmerge_path = shutil.which("bbmerge.sh") or shutil.which("bbmerge")
        if not bbmerge_path:
            raise FileNotFoundError("bbmerge not found in PATH")
        bbmerge_cmd = [bbmerge_path,
                       f"in1={reads_list[1]}",
                       f"in2={reads_list[2]}",
                       f"out={bbmap_out_path}",
                       f"ihist={input_folder}/insert_size_histogram.txt"]
        _ = run_subprocess_cmd(bbmerge_cmd, False)
        try:
            avg_insert, std_dev = parse_bbmerge_output(insert_size_histogram_txt)
        except:
            print("NOTE:\tUnable to parse avg_insert or std_dev, using default values: 251 and 30 respectively.")
            avg_insert = 251
            std_dev = 30
    return avg_insert, std_dev


# --------------------------------------------------------------
# Modify MaSuRCA script to skip gap closing
# --------------------------------------------------------------
def skip_gap_closing_section(assembly_sh_path):
    """Modify MaSuRCA assemble.sh script to bypass gap closing.

    Edits the script to skip the gap closing step and use the 9-terminator output.

    Args:
        assembly_sh_path (str): Path to the assemble.sh script.

    Returns:
        str: Path to the modified script.
    """
    with open(assembly_sh_path, "r") as f_in:
        original_script = f_in.read()
    pattern = re.compile(r"(if \[ -s \$CA_DIR/9-terminator/genome\.scf\.fasta \];then)"
                         r"(.*?)"
                         r"(else\s+fail 'Assembly stopped or failed, see \$CA_DIR\.log'\nfi)",
                         re.DOTALL)
    replacement_snippet = (r"\1\n"
                           r"  # Force the final terminator to remain '9-terminator'\n"
                           r"  log \"Skipping gap closing step; using 9-terminator as final.\"\n"
                           r"  TERMINATOR='9-terminator'\n"
                           r"\3")
    modified_text = re.sub(pattern, replacement_snippet, original_script)
    modified_script_path = os.path.join(os.path.dirname(assembly_sh_path), "assemble_skip_gap.sh")
    with open(modified_script_path, "w") as f_out:
        f_out.write(modified_text)
    return modified_script_path


# --------------------------------------------------------------
# Generate and run MaSuRCA assembly
# --------------------------------------------------------------
def masurca_config_gen(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """Generate and execute MaSuRCA configuration for genome assembly.

    Configures MaSuRCA for CABOG assembly with Illumina and optional long reads,
    runs the assembly, and performs quality control.

    Args:
        sample_id (str): Sample identifier.
        input_csv (str): Path to metadata CSV file.
        output_dir (str): Directory for output files.
        cpu_threads (int): Number of CPU threads to use.
        ram_gb (int): Available RAM in GB.

    Returns:
        str or None: Path to the gzipped MaSuRCA assembly FASTA, or None if assembly fails.
    """    
    input_df = pd.read_csv(input_csv)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]  # Convert to Series (single row)

    # Identify read paths, reference, and BUSCO lineage info from CSV
    illumina_sra = current_series["ILLUMINA_SRA"]
    illumina_f_raw_reads = current_series["ILLUMINA_RAW_F_READS"]
    illumina_r_raw_reads = current_series["ILLUMINA_RAW_R_READS"]
    ont_sra = current_series["ONT_SRA"]
    ont_raw_reads = current_series["ONT_RAW_READS"]
    pacbio_sra = current_series["PACBIO_SRA"]
    pacbio_raw_reads = current_series["PACBIO_RAW_READS"]
    ref_seq_gca = current_series["REF_SEQ_GCA"]
    ref_seq = current_series["REF_SEQ"]
    species_id = current_series["SPECIES_ID"]
    est_size = current_series["EST_SIZE"]
    species_dir = os.path.join(output_dir, species_id)
    sample_dir = os.path.join(species_dir, sample_id)

    if pd.notna(ont_sra) and pd.isna(ont_raw_reads):
        ont_raw_reads = os.path.join(species_dir, "ONT", f"{ont_sra}.fastq")
    if pd.notna(illumina_sra) and pd.isna(illumina_f_raw_reads) and pd.isna(illumina_r_raw_reads):
        illumina_f_raw_reads = os.path.join(species_dir, "Illumina", f"{illumina_sra}_1.fastq")
        illumina_r_raw_reads = os.path.join(species_dir, "Illumina", f"{illumina_sra}_2.fastq")    
    if pd.notna(pacbio_sra) and pd.isna(pacbio_raw_reads):
        pacbio_raw_reads = os.path.join(species_dir, "PacBio", f"{pacbio_sra}.fastq")
    if pd.notna(ref_seq_gca) and pd.isna(ref_seq):
        ref_seq = os.path.join(species_dir, "RefSeq", f"{species_id}_{ref_seq_gca}_RefSeq.fasta")

    masurca_out_dir = os.path.join(sample_dir, "masurca_assembly")    

    print(f"DEBUG - illumina_sra - {illumina_sra}")
    print(f"DEBUG - illumina_f_raw_reads - {illumina_f_raw_reads}")
    print(f"DEBUG - illumina_r_raw_reads - {illumina_r_raw_reads}")
    print(f"DEBUG - ont_sra - {ont_sra}")
    print(f"DEBUG - ont_raw_reads - {ont_raw_reads}")
    print(f"DEBUG - pacbio_sra - {pacbio_sra}")
    print(f"DEBUG - pacbio_raw_reads - {pacbio_raw_reads}")
    print(f"DEBUG - ref_seq_gca - {ref_seq_gca}")
    print(f"DEBUG - ref_seq - {ref_seq}")   
    print(f"DEBUG - species_id - {species_id}")
    print(f"DEBUG - est_size - {est_size}")
    print(f"DEBUG - masurca_out_dir - {masurca_out_dir}")
    
    # Set Illumina deduplicated read paths only if Illumina reads are present
    illu_dedup_f_reads = None
    illu_dedup_r_reads = None
    if pd.notna(illumina_f_raw_reads) and pd.notna(illumina_r_raw_reads):
        illu_dedup_f_reads = os.path.join(species_dir, "Illumina", f"{species_id}_illu_forward_dedup.fastq")
        illu_dedup_r_reads = os.path.join(species_dir, "Illumina", f"{species_id}_illu_reverse_dedup.fastq")

    print(f"DEBUG - illu_dedup_f_reads - {illu_dedup_f_reads}")
    print(f"DEBUG - illu_dedup_r_reads - {illu_dedup_r_reads}")

    # Set long-read paths (ONT or PacBio), prefer prefiltered, fallback to raw
    highest_mean_qual_long_reads = None
    if pd.notna(ont_raw_reads):
        print("DEBUG - ONT RAW READS EXIST!")
        candidate = os.path.join(species_dir, "ONT", f"{species_id}_ONT_highest_mean_qual_long_reads.fastq")
        highest_mean_qual_long_reads = candidate if os.path.exists(candidate) else ont_raw_reads
    elif pd.notna(pacbio_raw_reads):
        print("DEBUG - PACBIO RAW READS EXIST!")
        candidate = os.path.join(species_dir, "PacBio", f"{species_id}_PacBio_highest_mean_qual_long_reads.fastq")
        highest_mean_qual_long_reads = candidate if os.path.exists(candidate) else pacbio_raw_reads

    print(f"DEBUG - highest_mean_qual_long_reads    - {highest_mean_qual_long_reads}")

    # Check if only reference or no Illumina reads are provided and skip (NEW: Skip PacBio-only runs)
    if pd.isna(illumina_f_raw_reads) and pd.isna(illumina_r_raw_reads):
        print("SKIP:\tNo Illumina paired-end reads provided, required for MaSuRCA assembly")
        return None

    os.makedirs(masurca_out_dir, exist_ok=True)

    # Calculate insert size and standard deviation based on read type
    if pd.notna(illumina_f_raw_reads) and pd.notna(illumina_r_raw_reads):
        avg_insert, std_dev = bbmap_stats(masurca_out_dir,
                                         [ont_raw_reads, illumina_f_raw_reads, illumina_r_raw_reads, pacbio_raw_reads], cpu_threads)
    elif pd.notna(pacbio_raw_reads):
        avg_insert, std_dev = 15000, 4000
    elif pd.notna(ont_raw_reads):
        avg_insert, std_dev = 8000, 3000
    else:
        avg_insert, std_dev = 251, 30  # Fallback for edge cases
        
    # Set the Jellyfish hash size (adjust if available RAM is lower than 62GB)
    est_size_numb = re.match(r"^(\d+(?:\.\d+)?)(\D+)$", est_size).group(1)
    est_size_mult = re.match(r"^(\d+(?:\.\d+)?)(\D+)$", est_size).group(2)
    multipliers = {'m': 10**6, 'g': 10**9}
    print(est_size_mult)
    if est_size_mult in multipliers:
        jf_size = int(float(est_size_numb) * multipliers[est_size_mult])
    else:
        print(f"NOTE:\tUnable to parse input estimated size {est_size}, using default: 25000000")
        jf_size = 25000000
    
    # Ensure work directory output
    starting_work_dir = os.getcwd()
    if "work" not in starting_work_dir:
        current_work_dir = masurca_out_dir
    else:
        current_work_dir = starting_work_dir
    os.chdir(current_work_dir)
    
    # Define the desired assembly file paths
    data_output_folder = find_ca_folder(current_work_dir)
    primary_genome_scf = os.path.join(data_output_folder, "primary.genome.scf.fasta")  # (NEW: Fallback for interrupted assemblies)
    terminator_genome_scf = os.path.join(data_output_folder, "9-terminator", "genome.scf.fasta")  # (NEW: Check for successful assembly output)
    egap_masurca_assembly_path = os.path.join(masurca_out_dir, f"{sample_id}_masurca.fasta")

    # Determine if this is a hybrid assembly
    hybrid_assembly = pd.notna(illumina_f_raw_reads) and pd.notna(illumina_r_raw_reads) and (pd.notna(ont_raw_reads) or pd.notna(pacbio_raw_reads))

    # Set assembler-specific flags
    soap_assembly = 0
    flye_assembly = 0
    if hybrid_assembly:
        use_linking_mates = 0
        close_gaps = 0
    else:
        use_linking_mates = 1
        close_gaps = 1

    # Check for existing assembly and save in work directory (NEW: Handle successful and interrupted assemblies)
    if os.path.exists(terminator_genome_scf):
        print(f"PASS:\tAssembly found at {terminator_genome_scf}, saving to work directory")
        shutil.copy(terminator_genome_scf, egap_masurca_assembly_path)
        egap_masurca_assembly_path, masurca_stats_list, _ = qc_assessment("masurca", input_csv, sample_id, output_dir, cpu_threads, ram_gb)
        return egap_masurca_assembly_path
    elif os.path.exists(primary_genome_scf):
        print(f"PASS:\tInterrupted assembly found at {primary_genome_scf}, saving to work directory")
        shutil.copy(primary_genome_scf, egap_masurca_assembly_path)
        egap_masurca_assembly_path, masurca_stats_list, _ = qc_assessment("masurca", input_csv, sample_id, output_dir, cpu_threads, ram_gb)
        return egap_masurca_assembly_path
    elif os.path.exists(egap_masurca_assembly_path):
        print(f"SKIP:\tMaSuRCA Assembly, scaffolded assembly already exists: {egap_masurca_assembly_path}.")
        egap_masurca_assembly_path, masurca_stats_list, _ = qc_assessment("masurca", input_csv, sample_id, output_dir, cpu_threads, ram_gb)
        return egap_masurca_assembly_path
    else:        
        # Build the configuration file (NEW: Restored original configuration logic)
        config_content = ["DATA\n"]
        
        # Add Illumina paired-end reads if available
        if illu_dedup_f_reads and illu_dedup_r_reads:
            config_content.append(f"PE= pe {int(avg_insert)} {int(std_dev)} "
                                  f"{illu_dedup_f_reads} {illu_dedup_r_reads}\n")

        # Add long-read data (NEW: Use only determined long-read file)
        if highest_mean_qual_long_reads:
            if pd.notna(ont_raw_reads):
                config_content.append(f"NANOPORE={highest_mean_qual_long_reads}\n")
            elif pd.notna(pacbio_raw_reads):
                config_content.append(f"PACBIO={highest_mean_qual_long_reads}\n")

        # Add reference if available (NEW: Use unzipped ref_seq path without transformations)
        if pd.notna(ref_seq):
            config_content.append(f"REFERENCE={ref_seq}\n")
        config_content.append("END\n")

        # Add PARAMETERS section
        config_content.append("PARAMETERS\n")
        config_content.append("GRAPH_KMER_SIZE=auto\n")
        config_content.append(f"USE_LINKING_MATES={use_linking_mates}\n")
        config_content.append(f"CLOSE_GAPS={close_gaps}\n")
        config_content.append("MEGA_READS_ONE_PASS=0\n")
        config_content.append("LIMIT_JUMP_COVERAGE=300\n")
        config_content.append("CA_PARAMETERS=cgwErrorRate=0.15\n")
        config_content.append(f"NUM_THREADS={cpu_threads}\n")
        config_content.append(f"JF_SIZE={jf_size}\n")
        config_content.append(f"SOAP_ASSEMBLY={soap_assembly}\n")
        config_content.append(f"FLYE_ASSEMBLY={flye_assembly}\n")
        config_content.append("END\n")

        # Debug: Log configuration file contents
        print(f"DEBUG: Config content:\n{''.join(config_content)}")

        # Write the configuration file
        config_path = os.path.join(current_work_dir, "masurca_config_file.txt")
        with open(config_path, "w") as file:
            for entry in config_content:
                file.write(entry)

        # Run the MaSuRCA configuration command
        masurca_config_cmd = ["masurca", "masurca_config_file.txt"]
        _ = run_subprocess_cmd(masurca_config_cmd, False)

        # Modify the assemble.sh to skip gap closing
        assemble_sh_path = os.path.join(current_work_dir, "assemble.sh")
        modified_assemble_sh_path = skip_gap_closing_section(assemble_sh_path)

        # Run the modified assembly script
        masurca_assemble_cmd = ["bash", modified_assemble_sh_path]
        _ = run_subprocess_cmd(masurca_assemble_cmd, False)

        # Refresh assembly file paths post run
        data_output_folder = find_ca_folder(current_work_dir)
        primary_genome_scf = os.path.join(data_output_folder, "primary.genome.scf.fasta")  # (NEW: Fallback for interrupted assemblies)
        terminator_genome_scf = os.path.join(data_output_folder, "9-terminator", "genome.scf.fasta")  # (NEW: Check for successful assembly output)

        # Handle assembly output (NEW: Check for successful or interrupted assembly)
        if os.path.exists(terminator_genome_scf):
            print(f"NOTE:\tSaving successful assembly to {egap_masurca_assembly_path}")
            shutil.copy(terminator_genome_scf, egap_masurca_assembly_path)
        elif os.path.exists(primary_genome_scf):
            print(f"NOTE:\tSaving interrupted assembly to {egap_masurca_assembly_path}")
            shutil.copy(primary_genome_scf, egap_masurca_assembly_path)
        else:
            print(f"ERROR:\tNo assembly file found at {terminator_genome_scf} or {primary_genome_scf}")
            return None

    print(f"DEBUG: Final assembly path: {egap_masurca_assembly_path}")
    egap_masurca_assembly_path, masurca_stats_list, _ = qc_assessment("masurca", input_csv, sample_id, output_dir, cpu_threads, ram_gb)

    # Return paths in work directory (NEW: Rely on publishDir for final output)
    return egap_masurca_assembly_path


if __name__ == "__main__":
    # Provide a usage message if the correct number of arguments is not supplied.
    if len(sys.argv) != 6:
        print("Usage: python3 assemble_masurca.py <sample_id> <input_csv> "
            "<output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)
    
    egap_masurca_assembly_path = masurca_config_gen(sys.argv[1],       # sample_id
                                                    sys.argv[2],       # input_csv
                                                    sys.argv[3],       # output_dir
                                                    str(sys.argv[4]),  # cpu_threads
                                                    str(sys.argv[5]))  # ram_gb
