#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare_assemblies.py

Updated on Sat Mar 29 2025

This script compares assemblies from MaSuRCA, SPAdes, Flye, and Hifiasm, selecting the best based on Compleasm completeness, contig count, and N50.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import os, sys, shutil
from collections import Counter
import pandas as pd
from utilities import run_subprocess_cmd, get_current_row_data
from collections import Counter


# --------------------------------------------------------------
# Extract Compleasm completeness score
# --------------------------------------------------------------
def get_compleasm_score(assembly, db, cpu_threads):
    """Run Compleasm and extract the completeness score.

    Executes Compleasm on the assembly with the specified database and retrieves
    the completeness percentage.

    Args:
        assembly (str): Path to the assembly FASTA file.
        db (str): Compleasm database identifier (e.g., 'fungi_odb10').
        cpu_threads (int): Number of CPU threads to use.

    Returns:
        float or None: Completeness score if successful, else None.
    """
    if assembly == "None" or not os.path.exists(assembly):
        return None
    out_dir = f"{os.path.basename(assembly)}_{db}_compleasm"
    os.makedirs(out_dir, exist_ok=True)
    summary = f"{out_dir}/summary.txt"
    if not os.path.exists(summary):
        run_subprocess_cmd(["compleasm", "run", "-a", assembly,
                            "-o", out_dir, "--lineage", f"{db}_odb12",
                            "-t", str(cpu_threads)], shell_check=False)
    with open(summary, "r") as f:
        for line in f:
            if "Complete" in line:  # Assuming format: "Complete: XX.X%"
                return float(line.split()[1])
    return None


# --------------------------------------------------------------
# Extract BUSCO completeness score
# --------------------------------------------------------------
def get_busco_score(assembly, db, cpu_threads):
    """Run BUSCO and extract the completeness score.

    Executes BUSCO on the assembly with the specified database and retrieves
    the completeness percentage (Single + Duplicated).

    Args:
        assembly (str): Path to the assembly FASTA file.
        db (str): BUSCO database identifier (e.g., 'fungi_odb10').
        cpu_threads (int): Number of CPU threads to use.

    Returns:
        float or None: Completeness score if successful, else None.
    """
    if assembly == "None" or not os.path.exists(assembly):
        return None
    assembly_short = os.path.basename(assembly)
    summary_short = assembly_short.replace(".fasta", f"_{db}_busco")
    out_dir = assembly.replace(".fasta", f"_{db}_busco")
    os.makedirs(out_dir, exist_ok=True)
    summary = os.path.join(out_dir, f"short_summary.specific.{db}_odb12.{summary_short}.txt")
    if not os.path.exists(summary):
        run_subprocess_cmd(["busco", "-m", "genome",
                            "-i", assembly, "-f",
                            "-l", db,
                            "-c", str(cpu_threads),
                            "-o", summary_short,
                            "--out_path", out_dir], shell_check=False)
    with open(summary, "r") as f:
        for line in f:
            if "C:" in line:  # Assuming format: "C:XX.X%[S:YY.Y%,D:ZZ.Z%]"
                return float(line.split("[")[0].split(":")[1].replace("%",""))
    return None


# --------------------------------------------------------------
# Extract Quast contig count and N50
# --------------------------------------------------------------
def get_quast_stats(assembly, cpu_threads):
    """Run Quast and extract contig count and N50 statistics.

    Executes Quast on the assembly and retrieves the total number of contigs and N50 value.

    Args:
        assembly (str): Path to the assembly FASTA file.
        cpu_threads (int): Number of CPU threads to use.

    Returns:
        tuple: (contig count, N50) as integers, or (None, None) if unsuccessful.
    """
    if assembly == "None" or not os.path.exists(assembly):
        return None, None
    out_dir = f"{os.path.basename(assembly)}_quast"
    os.makedirs(out_dir, exist_ok=True)
    report = f"{out_dir}/report.tsv"
    if not os.path.exists(report):
        run_subprocess_cmd(["quast.py", "-o", out_dir, "-t", str(cpu_threads), assembly], False)
    with open(report, "r") as f:
        lines = f.readlines()
        contig_count = None
        n50 = None
        for line in lines:
            if "N50" in line:
                n50 = int(line.split()[1])
            if "Total number of contigs" in line:
                contig_count = int(line.split()[1])
        return contig_count, n50
    return None, None


# --------------------------------------------------------------
# Compare and select best assembly
# --------------------------------------------------------------
def compare_assemblies(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """Compare assemblies and select the best based on BUSCO and Quast metrics.

    Evaluates MaSuRCA, SPAdes, Flye, and Hifiasm assemblies using BUSCO completeness,
    contig count, and N50, then selects the best assembly.

    Args:
        sample_id (str): Sample identifier.
        input_csv (str): Path to metadata CSV file.
        output_dir (str): Directory for output files.
        cpu_threads (int): Number of CPU threads to use.
        ram_gb (int): Available RAM in GB.

    Returns:
        str: Path to the selected best assembly FASTA (gzipped).
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
    busco_1 = current_series["BUSCO_1"]
    busco_2 = current_series["BUSCO_2"]
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
    print(f"DEBUG - species_dir - {species_dir}")
    print(f"DEBUG - est_size - {est_size}")
    print(f"DEBUG - busco_1 - {busco_1}")
    print(f"DEBUG - busco_2 - {busco_2}")

    all_methods = ["MaSuRCA", "SPAdes", "Flye","Hifiasm"]
    assemblies = {}
    
    for method in all_methods:
         temp_assembly_path = os.path.join(sample_dir, f"{method.lower()}_assembly", f"{sample_id}_{method.lower()}.fasta")
         if os.path.exists(temp_assembly_path):
             assemblies[method] = temp_assembly_path
             print(f"DEBUG - {method}_assembly - {temp_assembly_path}")

    # Check if only reference is provided and skip
    if (pd.isna(illumina_f_raw_reads) or pd.isna(illumina_sra)) and (pd.isna(illumina_r_raw_reads) or pd.isna(illumina_sra)) and (pd.isna(ont_raw_reads) or pd.isna(ont_sra)) and (pd.isna(pacbio_raw_reads) or pd.isna(pacbio_sra)):
        print("SKIP:\tNo valid reads provided, required for assembly comparison.")
        return None

    # Gather stats
    stats_dict = {}
    for method, asm in assemblies.items():
        # remove the old 'assembly != "None"' check since asm is either a real path
        busco1   = get_busco_score(asm, busco_1, cpu_threads)
        busco2   = get_busco_score(asm, busco_2, cpu_threads)
        # busco1   = get_compleasm_score(asm, busco_1, cpu_threads)
        # busco2   = get_compleasm_score(asm, busco_2, cpu_threads)
        cnt, n50 = get_quast_stats(asm, cpu_threads)
        stats_dict[method] = [busco1, busco2, cnt, n50]

    # Dynamic list of methods
    methods = list(assemblies.keys())

    # Compare each metric
    labels = [
        "First BUSCO Completeness (Single + Duplicated)",
        "Second BUSCO Completeness (Single + Duplicated)",
        "Assembly Contig Count",
        "Assembly N50"
    ]
    custom_stats = []
    for idx in range(len(labels)):
        # collect only non‐None
        candidates = [(stats_dict[m][idx], m) for m in methods
                      if stats_dict[m][idx] is not None]
        if not candidates:
            custom_stats.append("Unknown")
            continue

        # minimize contig count (idx==2), else maximize
        if idx == 2:
            best = min(candidates, key=lambda x: x[0])[1]
        else:
            best = max(candidates, key=lambda x: x[0])[1]
        custom_stats.append(best)

    print(f"DEBUG - len(assemblies) - {len(assemblies)}")

    # Pick the overall winner
    if len(assemblies) == 1:
        most_represented_method = list(assemblies.keys())[0]
    elif len(assemblies) > 0:
        counts = Counter(m for m in custom_stats if m != "Unknown")
        most_represented_method = counts.most_common(1)[0][0] if counts else "Unknown"
    if most_represented_method == "Unknown":
        print("ERROR:\tNo valid assemblies to compare")
        sys.exit(1)

    print(f"Best Initial Assembly: {most_represented_method}")
    for lbl, val in zip(labels, stats_dict[most_represented_method]):
        print(f"{lbl}: {val}")

    # Copy into sample_dir
    orig = assemblies[most_represented_method]
    dest_fasta     = os.path.join(sample_dir, f"{species_id}_best_assembly.fasta")
    os.makedirs(sample_dir, exist_ok=True)
    shutil.copy(orig, dest_fasta)
    best_assembly = dest_fasta

    return best_assembly


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python3 compare_assemblies.py <sample_id> <input_csv> "
            "<output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)
        
    best_assembly = compare_assemblies(sys.argv[1],       # sample_id
                                       sys.argv[2],       # input_csv
                                       sys.argv[3],       # output_dir
                                       str(sys.argv[4]),  # cpu_threads
                                       str(sys.argv[5]))  # ram_gb
