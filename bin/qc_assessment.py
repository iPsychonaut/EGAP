#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qc_assessment.py

This script performs final assembly assessment with BUSCO and QUAST,
analyzes and classifies the assembly, and finalizes renaming and compression.

Created on Wed Aug 16 2023

Updated on Wed Sept 3 2025

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os
import sys
import shutil
import pandas as pd
from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt
from utilities import run_subprocess_cmd, pigz_compress, pigz_decompress, get_current_row_data, analyze_nanostats

# === Global BUSCO/Compleasm toggle ===
# If True, prefer Compleasm over BUSCO for completeness analysis.
USE_COMPLEASM = True

# If True, automatically fall back to BUSCO when Compleasm fails
# or produces no metrics (recommended so downstream code always has data).
FALLBACK_TO_BUSCO_ON_FAIL = True

def _normalize_busco_keycase(sample_stats_dict: dict, busco_count: str) -> None:
    """
    Ensure keys are in the UPPERCASE style expected elsewhere, e.g.:
    FIRST_BUSCO_S, FIRST_BUSCO_D, FIRST_BUSCO_F, FIRST_BUSCO_M, FIRST_BUSCO_C
    """
    up = busco_count.upper()
    lo = busco_count.lower()
    for tag in ("S", "D", "F", "M", "C"):
        low_key = f"{lo}_BUSCO_{tag}"
        up_key  = f"{up}_BUSCO_{tag}"
        if low_key in sample_stats_dict and up_key not in sample_stats_dict:
            sample_stats_dict[up_key] = sample_stats_dict[low_key]


def run_lineage_eval(assembly_path: str,
                     sample_id: str,
                     sample_stats_dict: dict,
                     busco_count: str,       # "first" or "second"
                     busco_odb: str,
                     assembly_type: str,
                     cpu_threads) -> str:
    """
    Single entry point for BUSCO-like evaluation.
    Respects USE_COMPLEASM and will optionally fall back to BUSCO.
    """
    if USE_COMPLEASM:
        print("INFO:\tUsing Compleasm for BUSCO-like analysis.")
        try:
            compleasm_assembly(assembly_path, sample_id, sample_stats_dict,
                               busco_count, busco_odb, assembly_type, cpu_threads)
            _normalize_busco_keycase(sample_stats_dict, busco_count)

            # Did we actually get metrics? Accept C or (S and D)
            up = busco_count.upper()
            got_any = (
                f"{up}_BUSCO_C" in sample_stats_dict or
                (f"{up}_BUSCO_S" in sample_stats_dict and f"{up}_BUSCO_D" in sample_stats_dict)
            )
            if not got_any and FALLBACK_TO_BUSCO_ON_FAIL:
                print("WARN:\tCompleasm produced no metrics; falling back to BUSCO.")
                busco_assembly(assembly_path, sample_id, sample_stats_dict,
                               busco_count, busco_odb, assembly_type, cpu_threads)

        except Exception as e:
            print(f"WARN:\tCompleasm raised an exception: {e}")
            if FALLBACK_TO_BUSCO_ON_FAIL:
                print("WARN:\tFalling back to BUSCO.")
                busco_assembly(assembly_path, sample_id, sample_stats_dict,
                               busco_count, busco_odb, assembly_type, cpu_threads)
    else:
        print("INFO:\tUsing BUSCO for completeness analysis.")
        busco_assembly(assembly_path, sample_id, sample_stats_dict,
                       busco_count, busco_odb, assembly_type, cpu_threads)

    return assembly_path

# --------------------------------------------------------------
# Validate FASTA file
# --------------------------------------------------------------
def validate_fasta(file_path):
    """Validate that a FASTA file exists, is non-empty, and contains valid nucleotide sequences.

    Args:
        file_path (str): Path to the FASTA file.

    Returns:
        bool: True if valid, False otherwise.
    """
    if not os.path.exists(file_path):
        print(f"ERROR:\tFASTA file not found: {file_path}")
        return False
    if os.path.getsize(file_path) < 100:
        print(f"ERROR:\tFASTA file is suspiciously small: {file_path}")
        return False
    try:
        with open(file_path, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                if not record.seq:
                    print(f"ERROR:\tFASTA file contains empty sequences: {file_path}")
                    return False
                if not all(c.upper() in "ATCGN" for c in record.seq):
                    print(f"ERROR:\tFASTA file contains non-nucleotide sequences: {file_path}")
                    return False
                return True
    except Exception as e:
        print(f"ERROR:\tInvalid FASTA format in {file_path}: {str(e)}")
        return False
    return False


# --------------------------------------------------------------
# Perform quality control on reads using NanoPlot
# --------------------------------------------------------------
def nanoplot_qc_reads(INPUT_READS, READS_ORIGIN, CPU_THREADS, sample_stats_dict):
    """Run NanoPlot to perform quality control on sequencing reads.

    Executes NanoPlot to generate quality metrics and updates the sample statistics
    dictionary with the results.

    Args:
        INPUT_READS (str): Path to the input reads file.
        READS_ORIGIN (str): Read type and stage (e.g., 'Raw_ONT_', 'Filt_ONT_').
        CPU_THREADS (int): Number of CPU threads to use.
        sample_stats_dict (dict): Dictionary to store QC metrics.

    Returns:
        dict: Updated sample statistics dictionary with NanoPlot metrics.
    """
    print(f"NanoPlotting reads: {INPUT_READS}...")
    output_dir = os.path.join(os.path.dirname(INPUT_READS), f"{READS_ORIGIN}nanoplot_analysis")
    print(f"DEBUG - output_dir - {output_dir}")
    os.makedirs(output_dir, exist_ok=True)
    nanoplot_out_file = os.path.join(output_dir, f"{READS_ORIGIN}NanoStats.txt")
    if os.path.exists(nanoplot_out_file):
        print(f"SKIP:\tNanoPlot output already exists: {nanoplot_out_file}.")
    else:
        raw_nanoplot_cmd = ["NanoPlot", "--fastq", INPUT_READS, "-t", str(CPU_THREADS),
                            "-o", output_dir, "--plots", "kde", "dot", "--loglength",
                            "--N50", "--title", f"{READS_ORIGIN} Reads: Preliminary Data",
                            "--prefix", READS_ORIGIN, "--verbose"]
        result = run_subprocess_cmd(raw_nanoplot_cmd, shell_check=False)
        if result != 0:
            print(f"WARN:\tNanoPlot failed with return code {result}")

    sample_stats_dict = analyze_nanostats(READS_ORIGIN, nanoplot_out_file, sample_stats_dict)
    return sample_stats_dict


# --------------------------------------------------------------
# Parse NanoPlot summary for mean quality
# --------------------------------------------------------------
def get_mean_quality(nanoplot_dir, prefix):
    """Extract mean quality from a NanoPlot summary file.

    Reads the NanoPlot summary file to retrieve the mean quality score.

    Args:
        nanoplot_dir (str): Directory containing the NanoPlot summary file.
        prefix (str): Prefix for the NanoPlot summary file (e.g., 'Corrected_ONT_').

    Returns:
        float or None: Mean quality value if found, else None.
    """
    summary_file = os.path.join(nanoplot_dir, f"{prefix}NanoStats.txt")
    if not os.path.exists(summary_file):
        return None
    with open(summary_file, "r") as f:
        for line in f:
            if "Mean quality:" in line:
                return float(line.split(":")[1].strip())
    return None


# --------------------------------------------------------------
# Classify assembly quality metrics
# --------------------------------------------------------------
def classify_assembly(sample_stats):
    """Classify assembly quality based on BUSCO, contig count, and N50 metrics.

    Assigns ratings ('POOR', 'OK', 'GREAT', 'AMAZING') to BUSCO completeness,
    contig count, and N50, and computes an overall rating.

    Args:
        sample_stats (dict): Dictionary containing assembly metrics (e.g., BUSCO scores).

    Returns:
        dict: Ratings for each metric and an overall rating.
    """
    results = {}
    ranking = ["POOR", "OK", "GREAT", "AMAZING"]

    # Initialize BUSCO metrics if not present
    if "FIRST_BUSCO_C" not in sample_stats:
        sample_stats["FIRST_BUSCO_C"] = 0.0
    if "SECOND_BUSCO_C" not in sample_stats:
        sample_stats["SECOND_BUSCO_C"] = 0.0
    if "ASSEMBLY_CONTIGS" not in sample_stats:
        sample_stats["ASSEMBLY_CONTIGS"] = float("inf")
    if "ASSEMBLY_N50" not in sample_stats:
        sample_stats["ASSEMBLY_N50"] = 0.0

    # Classify BUSCO completeness for the first lineage
    if sample_stats["FIRST_BUSCO_C"] >= 98.5:
        results["FIRST_BUSCO_C"] = "AMAZING"
    elif sample_stats["FIRST_BUSCO_C"] > 90.0:
        results["FIRST_BUSCO_C"] = "GREAT"
    elif sample_stats["FIRST_BUSCO_C"] >= 75.0:
        results["FIRST_BUSCO_C"] = "OK"
    else:
        results["FIRST_BUSCO_C"] = "POOR"

    # Classify BUSCO completeness for the second lineage
    if sample_stats["SECOND_BUSCO_C"] >= 98.5:
        results["SECOND_BUSCO_C"] = "AMAZING"
    elif sample_stats["SECOND_BUSCO_C"] > 90.0:
        results["SECOND_BUSCO_C"] = "GREAT"
    elif sample_stats["SECOND_BUSCO_C"] >= 75.0:
        results["SECOND_BUSCO_C"] = "OK"
    else:
        results["SECOND_BUSCO_C"] = "POOR"

    # Classify assembly contig count
    if sample_stats["ASSEMBLY_CONTIGS"] <= 100:
        results["ASSEMBLY_CONTIGS"] = "AMAZING"
    elif sample_stats["ASSEMBLY_CONTIGS"] <= 1000:
        results["ASSEMBLY_CONTIGS"] = "GREAT"
    elif sample_stats["ASSEMBLY_CONTIGS"] <= 10000:
        results["ASSEMBLY_CONTIGS"] = "OK"
    else:
        results["ASSEMBLY_CONTIGS"] = "POOR"

    # Classify assembly N50
    if sample_stats["ASSEMBLY_N50"] <= 100:
        results["ASSEMBLY_N50"] = "POOR"
    elif sample_stats["ASSEMBLY_N50"] <= 1000:
        results["ASSEMBLY_N50"] = "OK"
    elif sample_stats["ASSEMBLY_N50"] <= 10000:
        results["ASSEMBLY_N50"] = "GREAT"
    else:
        results["ASSEMBLY_N50"] = "AMAZING"

    # Compute overall rating
    count = Counter(results.values())
    max_count = max(count.values())
    most_common = [key for key, value in count.items() if value == max_count]
    most_common.sort(key=lambda x: ranking.index(x))
    results["OVERALL"] = "/".join(most_common)

    return results


# --------------------------------------------------------------
# Generate BUSCO status plots
# --------------------------------------------------------------
def plot_busco(sample_id, busco_type, busco_odb, input_busco_tsv, input_fasta, assembly_type):
    """Create stacked bar plots of BUSCO statuses from TSV results.

    Generates SVG and PNG plots showing BUSCO statuses (Single, Duplicated,
    Incomplete, Fragmented) per sequence.

    Args:
        sample_id (str): Sample identifier for output naming.
        busco_type (str): Type of BUSCO run ('busco' or 'compleasm').
        busco_odb (str): BUSCO lineage identifier (e.g., 'fungi_odb10').
        input_busco_tsv (str): Path to the BUSCO TSV file.
        input_fasta (str): Path to the input FASTA file for naming.
        assembly_type (str): Descriptor for the assembly type.
    """
    print(f"Generating BUSCO plot for {input_busco_tsv}...")

    # Check if TSV exists
    if not os.path.exists(input_busco_tsv):
        print(f"ERROR:\tBUSCO TSV not found: {input_busco_tsv}. Skipping plot generation.")
        return

    # Load BUSCO TSV with appropriate headers
    try:
        if busco_type == "busco":
            busco_df = pd.read_csv(input_busco_tsv, sep="\t", skiprows=2, dtype=str)
        elif busco_type == "compleasm":
            busco_df = pd.read_csv(input_busco_tsv, sep="\t", header=0, dtype=str)
    except Exception as e:
        print(f"ERROR:\tFailed to read BUSCO TSV {input_busco_tsv}: {str(e)}")
        return

    # Handle empty data
    if busco_df.empty:
        print("WARNING: BUSCO input file is empty. Generating placeholder plot.")
        plt.figure(figsize=(12, 8))
        plt.text(0.5, 0.5, "No valid BUSCO data to plot", fontsize=14,
                 ha='center', va='center')
        plt.xticks([])
        plt.yticks([])
        plt.title("BUSCO Status Plot - No Data Available")
        output_busco_svg = os.path.join(os.path.dirname(input_fasta),
                                        f"{sample_id}_{assembly_type}_{busco_odb}_busco.svg")
        output_busco_png = os.path.join(os.path.dirname(input_fasta),
                                        f"{sample_id}_{assembly_type}_{busco_odb}_busco.png")
        plt.savefig(output_busco_svg, format="svg")
        plt.savefig(output_busco_png, format="png")
        plt.close()
        return

    # Prepare data for stacked bar plot
    busco_genes = len(busco_df)
    busco_df['Status'] = busco_df['Status'].replace("Complete", "Single")
    status_counts = busco_df.pivot_table(index='Sequence', columns='Status',
                                        aggfunc='size', fill_value=0)
    desired_order = ['Single', 'Duplicated', 'Incomplete', 'Fragmented']
    status_counts = status_counts.reindex(columns=desired_order, fill_value=0)
    total_sequences = len(status_counts)

    excluded_sequences = len(status_counts.loc[status_counts.drop(columns='Duplicated', errors='ignore').sum(axis=1) == 0])
    included_sequences = total_sequences - excluded_sequences

    filtered_status_counts = status_counts.loc[status_counts.drop(columns='Duplicated', errors='ignore').sum(axis=1) > 0]
    if filtered_status_counts.empty:
        print("WARNING: No valid BUSCO data available for plotting.")
        plt.figure(figsize=(12, 8))
        plt.text(0.5, 0.5, "No valid BUSCO data to plot", fontsize=14,
                 ha='center', va='center')
        plt.xticks([])
        plt.yticks([])
        plt.title("BUSCO Status Plot - No Data Available")
        output_busco_svg = os.path.join(os.path.dirname(input_fasta),
                                        f"{sample_id}_{assembly_type}_{busco_odb}_busco.svg")
        output_busco_png = os.path.join(os.path.dirname(input_fasta),
                                        f"{sample_id}_{assembly_type}_{busco_odb}_busco.png")
        plt.savefig(output_busco_svg, format="svg")
        plt.savefig(output_busco_png, format="png")
        plt.close()
        return

    # Reorder rows by total BUSCO matches
    filtered_status_counts = filtered_status_counts.loc[filtered_status_counts.sum(axis=1).sort_values(ascending=False).index]

    # Compute counts and plot
    status_totals = busco_df['Status'].value_counts()
    colors = {'Single': '#619B8AFF',
              'Duplicated': '#A1C181FF',
              'Incomplete': '#FE7F2DFF',
              'Fragmented': '#FCCA46FF'}
    ax = filtered_status_counts.plot(kind='bar', stacked=True, figsize=(12, 8),
                                     color=[colors[col] for col in filtered_status_counts.columns])
    legend_labels = [f"{status} ({round((status_totals.get(status, 0)/busco_genes)*100, 2)}%)"
                     for status in filtered_status_counts.columns]
    completeness_values = [round((status_totals.get(status, 0)/busco_genes)*100, 2)
                           for status in filtered_status_counts.columns]
    completeness_calc = round(completeness_values[0] + completeness_values[1], 2)

    plt.title(f"Distribution of {busco_odb} BUSCO Status per Sequence\n"
              f"Completeness: {completeness_calc}%")
    plt.xlabel(f"Sequences (Contig/Scaffold/Chromosome)\n"
               f"Included={included_sequences}, Excluded={excluded_sequences}")
    plt.ylabel(f"Number of BUSCO Matches (out of {busco_genes})")
    plt.xticks(rotation=45, ha='right')
    ax.legend(legend_labels, title="BUSCO Status", loc='upper right')
    plt.tight_layout()

    # Save plots (SVG, PNG)
    output_busco_svg = os.path.join(os.path.dirname(input_fasta),
                                    f"{sample_id}_{assembly_type}_{busco_odb}_busco.svg")
    output_busco_png = os.path.join(os.path.dirname(input_fasta),
                                    f"{sample_id}_{assembly_type}_{busco_odb}_busco.png")
    plt.savefig(output_busco_svg, format="svg")
    plt.savefig(output_busco_png, format="png")
    print(f"PASS:\tBUSCO {busco_odb} plot saved: {output_busco_svg} & {output_busco_png}")
    plt.close()


# --------------------------------------------------------------
# Run BUSCO on an assembly
# --------------------------------------------------------------
def busco_assembly(assembly_path, sample_id, sample_stats_dict,
                   busco_count, busco_odb, assembly_type, cpu_threads):
    """
    Evaluate assembly completeness using BUSCO and update statistics.

    Runs BUSCO in genome mode, parses the summary, and stores metrics in the
    sample statistics dictionary.

    Args:
        assembly_path (str): Path to the assembly FASTA file.
        sample_id (str): Sample identifier for output naming (not used in path).
        sample_stats_dict (dict): Dictionary to store BUSCO metrics.
        busco_count (str): Label for metrics ('first' or 'second').
        busco_odb (str): BUSCO lineage dataset (e.g., 'fungi_odb10').
        assembly_type (str): Descriptor for the assembly type.
        cpu_threads (int or str): Number of CPU threads to use.

    Returns:
        str: Original assembly path.
    """
    # Determine the directory containing the assembly
    assembly_dir = os.path.dirname(assembly_path)

    # Derive a clean base name (no extension) from the FASTA filename
    base = os.path.splitext(os.path.basename(assembly_path))[0]

    # Build BUSCO output directory alongside the FASTA
    busco_db_version = "odb12"
    busco_dir = os.path.join(assembly_dir, f"{base}_{busco_odb}_busco")
    os.makedirs(busco_dir, exist_ok=True)

    # Validate the input FASTA
    if not validate_fasta(assembly_path):
        print(f"ERROR:\tInvalid assembly for BUSCO: {assembly_path}. Skipping BUSCO.")
        return assembly_path

    # Path to the BUSCO summary file that BUSCO will generate
    busco_summary = os.path.join(
        busco_dir,
        f"short_summary.specific.{busco_odb}_{busco_db_version}."
        f"{base}_{busco_odb}_busco.txt"
    )

    # Run BUSCO if the summary does not already exist
    if os.path.exists(busco_summary):
        print(f"SKIP:\tBUSCO Summary already exists: {busco_summary}.")
    else:
        busco_cmd = [
            "busco", "-m", "genome",
            "-i", assembly_path, "-f",
            "-l", busco_odb,
            "-c", str(cpu_threads),
            "-o", f"{base}_{busco_odb}_busco",
            "--out_path", assembly_dir
        ]
        print(f"DEBUG - Running BUSCO: {' '.join(busco_cmd)}")
        rc = run_subprocess_cmd(busco_cmd, shell_check=False)
        if rc != 0:
            print(f"WARN:\tBUSCO failed with return code {rc}. Skipping BUSCO metrics.")
            return assembly_path

    # Generate BUSCO plot if needed
    comp_busco_svg = os.path.join(
        assembly_dir,
        f"{base}_{busco_odb}_busco.svg"
    )
    busco_tsv = os.path.join(
        busco_dir,
        f"run_{busco_odb}_{busco_db_version}",
        "full_table.tsv"
    )
    if not os.path.exists(comp_busco_svg):
        plot_busco(base, "busco", busco_odb, busco_tsv, assembly_path, assembly_type)

    # Parse BUSCO summary to update statistics
    try:
        with open(busco_summary, "r") as busco_file:
            for line in busco_file:
                if "[" in line:
                    # Split on '%' and look for C, S, D, F, M
                    parts = line.split("%")
                    for item in parts:
                        if "C" in item:
                            sample_stats_dict[f"{busco_count.upper()}_BUSCO_C"] = float(item.split(":")[1])
                        elif "S" in item:
                            sample_stats_dict[f"{busco_count.upper()}_BUSCO_S"] = float(item.split(":")[1])
                        elif "D" in item:
                            sample_stats_dict[f"{busco_count.upper()}_BUSCO_D"] = float(item.split(":")[1])
                        elif "F" in item:
                            sample_stats_dict[f"{busco_count.upper()}_BUSCO_F"] = float(item.split(":")[1])
                        elif "M" in item:
                            sample_stats_dict[f"{busco_count.upper()}_BUSCO_M"] = float(item.split(":")[1])
    except FileNotFoundError:
        print(f"ERROR:\tBUSCO summary not found: {busco_summary}")

    return assembly_path


# --------------------------------------------------------------
# Run Compleasm on an assembly
# --------------------------------------------------------------
def compleasm_assembly(assembly_path, sample_id, sample_stats_dict, busco_count, busco_odb, assembly_type, cpu_threads):
    """Evaluate assembly completeness using Compleasm and update statistics.

    Runs Compleasm, parses BUSCO-like metrics, generates a plot, and updates the
    sample statistics dictionary.

    Args:
        assembly_path (str): Path to the assembly FASTA file.
        sample_id (str): Sample identifier for output naming.
        sample_stats_dict (dict): Dictionary to store Compleasm metrics.
        busco_count (str): Label for metrics ('first' or 'second').
        busco_odb (str): BUSCO lineage dataset (e.g., 'fungi_odb10').
        assembly_type (str): Descriptor for the assembly type.
        cpu_threads (int or str): Number of threads to use.

    Returns:
        str: Original assembly path.
    """
    sample_dir = os.path.dirname(assembly_path)
    compleasm_dir = os.path.join(sample_dir, f"{os.path.basename(assembly_path).replace('.fasta', '')}_{busco_odb}_busco")
    print(f"DEBUG - busco_odb - {busco_odb}")
    print(f"DEBUG - compleasm_dir - {compleasm_dir}")

    if not os.path.exists(compleasm_dir):
        os.makedirs(compleasm_dir)

    # Validate assembly
    if not validate_fasta(assembly_path):
        print(f"ERROR:\tInvalid assembly for Compleasm: {assembly_path}. Skipping Compleasm.")
        return assembly_path

    # Check if a Compleasm summary already exists
    compleasm_summary = os.path.join(compleasm_dir, "summary.txt")
    if os.path.exists(compleasm_summary):
        print(f"SKIP:\tCompleasm Summary already exists: {compleasm_summary}.")
    else:
        compleasm_cmd = ["compleasm", "run",
                         "-a", assembly_path,
                         "-o", compleasm_dir,
                         "--lineage", busco_odb,
                         "-t", str(cpu_threads)]
        print(f"DEBUG - Running Compleasm: {' '.join(compleasm_cmd)}")
        result = run_subprocess_cmd(compleasm_cmd, shell_check=False)
        if result != 0:
            print(f"WARN:\tCompleasm failed with return code {result}. Skipping Compleasm metrics.")
            return assembly_path

    # Generate a BUSCO-like plot if not already present
    comp_busco_svg = os.path.join(sample_dir, f"{os.path.basename(assembly_path).replace('.fasta', '')}_{busco_odb}_busco.svg")
    compleasm_tsv = os.path.join(compleasm_dir, f"{busco_odb}_odb12", "full_table_busco_format.tsv")
    if not os.path.exists(comp_busco_svg):
        plot_busco(sample_id, "compleasm", busco_odb, compleasm_tsv, assembly_path, assembly_type)

    # Parse Compleasm's summary.txt to extract coverage metrics
    try:
        with open(compleasm_summary, "r") as compleasm_file:
            for line in compleasm_file:
                if "S:" in line:
                    sample_stats_dict[f"{busco_count}_BUSCO_S"] = float(line.split("S:")[-1].split(", ")[0].replace("\n", "").replace("%", ""))
                elif "D:" in line:
                    sample_stats_dict[f"{busco_count}_BUSCO_D"] = float(line.split("D:")[-1].split(", ")[0].replace("\n", "").replace("%", ""))
                elif "F:" in line:
                    sample_stats_dict[f"{busco_count}_BUSCO_F"] = float(line.split("F:")[-1].split(", ")[0].replace("\n", "").replace("%", ""))
                elif "M:" in line:
                    sample_stats_dict[f"{busco_count}_BUSCO_M"] = float(line.split("M:")[-1].split(", ")[0].replace("\n", "").replace("%", ""))
    except FileNotFoundError:
        print(f"ERROR:\tCompleasm summary not found: {compleasm_summary}")
        return assembly_path

    # Compute combined completeness: S + D
    sample_stats_dict[f"{busco_count}_BUSCO_C"] = (sample_stats_dict[f"{busco_count}_BUSCO_S"] +
                                                   sample_stats_dict[f"{busco_count}_BUSCO_D"])

    return assembly_path


# --------------------------------------------------------------
# Perform quality control assessment
# --------------------------------------------------------------
def qc_assessment(assembly_type, input_csv, sample_id, output_dir, cpu_threads, ram_gb):
    """Perform quality control on an assembly using BUSCO, QUAST, and coverage.

    Executes BUSCO and QUAST, calculates coverage, classifies assembly quality,
    and saves statistics.

    Args:
        assembly_type (str): Descriptor for the assembly (e.g., 'final').
        input_csv (str): Path to metadata CSV file.
        sample_id (str): Sample identifier.
        output_dir (str): Directory for output files.
        cpu_threads (int or str): Number of CPU threads to use.
        ram_gb (int or str): Available RAM in GB.

    Returns:
        tuple: (path to compressed assembly, list of key statistics, updated stats dictionary).
    """
    # Parse the CSV and retrieve relevant row data
    input_df = pd.read_csv(input_csv)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]  # Convert to Series (single row)

    # Identify read paths, reference, and BUSCO lineage info from CSV
    ont_raw_reads = current_series["ONT_RAW_READS"]
    illu_raw_f_reads = current_series["ILLUMINA_RAW_F_READS"]
    illu_raw_r_reads = current_series["ILLUMINA_RAW_R_READS"]
    pacbio_raw_reads = current_series["PACBIO_RAW_READS"]
    ref_seq = current_series["REF_SEQ"]
    ref_seq_gca = current_series["REF_SEQ_GCA"]
    first_busco_odb = current_series["BUSCO_1"]
    second_busco_odb = current_series["BUSCO_2"]
    kingdom_id = current_series["ORGANISM_KINGDOM"]
    karyote_id = current_series["ORGANISM_KARYOTE"]
    species_id = current_series["SPECIES_ID"]

    species_dir = os.path.join(output_dir, species_id)
    sample_dir = os.path.join(species_dir, sample_id)
    assembly_path = os.path.join(sample_dir, f"{assembly_type}_assembly", f"{sample_id}_{assembly_type}.fasta")

    if pd.notna(ref_seq_gca) and pd.isna(ref_seq):
        ref_seq = os.path.join(species_dir, "RefSeq", f"{species_id}_{ref_seq_gca}_RefSeq.fasta")

    print(f"Parsing assembly for index {current_index} from {input_csv}:\n{current_row}")

    # Ensure working directory is sample_dir
    os.makedirs(sample_dir, exist_ok=True)
    os.chdir(sample_dir)
    print(f"DEBUG - Set working directory to: {os.getcwd()}")

    # Validate assembly path
    if not validate_fasta(assembly_path):
        print(f"ERROR:\tInvalid or missing assembly: {assembly_path}")
        return None, None, sample_stats_dict

    # --------------------------------------------------------------
    # Compleasm/BUSCO QC: Run on two different lineages
    # --------------------------------------------------------------
    run_lineage_eval(assembly_path, sample_id, sample_stats_dict,
                     "first", first_busco_odb, assembly_type, cpu_threads)
    run_lineage_eval(assembly_path, sample_id, sample_stats_dict,
                     "second", second_busco_odb, assembly_type, cpu_threads)

    # --------------------------------------------------------------
    # QUAST QC
    # --------------------------------------------------------------
    quast_dir = os.path.join(sample_dir, f"{assembly_type}_assembly", f"{os.path.basename(assembly_path)}_quast")
    os.makedirs(quast_dir, exist_ok=True)
    quast_report_tsv = os.path.join(quast_dir, "report.tsv")

    # Run QUAST only if no existing report
    if os.path.exists(quast_report_tsv):
        print(f"SKIP:\tQUAST Report already exists: {quast_report_tsv}.")
    else:
        quast_cmd = ["quast", "--threads", str(cpu_threads)]
        if karyote_id == "eukaryote":
            quast_cmd.append("--eukaryote")
        if kingdom_id == "Funga":
            quast_cmd.append("--fungus")
        if pd.notna(ref_seq) and os.path.exists(ref_seq):
            quast_cmd.extend(["-r", ref_seq])
        quast_cmd.extend(["-o", quast_dir, assembly_path])
        print(f"DEBUG - Running QUAST: {' '.join(quast_cmd)}")
        result = run_subprocess_cmd(quast_cmd, shell_check=False)
        if result != 0:
            print(f"WARN:\tQUAST failed with return code {result}. Skipping QUAST metrics.")
            return assembly_path, None, sample_stats_dict

    # Parse QUAST report to populate sample_stats_dict
    try:
        with open(quast_report_tsv, "r") as quast_file:
            for line in quast_file:
                if "Total length (>= 0 bp)" in line:
                    sample_stats_dict["GENOME_SIZE"] = float(line.split("\t")[-1].strip())
                elif "# contigs" in line:
                    sample_stats_dict["ASSEMBLY_CONTIGS"] = float(line.split("\t")[-1].strip())
                elif "N50" in line:
                    sample_stats_dict["ASSEMBLY_N50"] = float(line.split("\t")[-1].strip())
                elif "L50" in line:
                    sample_stats_dict["ASSEMBLY_L50"] = float(line.split("\t")[-1].strip())
                elif "GC (%)" in line:
                    sample_stats_dict["ASSEMBLY_GC"] = float(line.split("\t")[-1].strip())
                if pd.notna(ref_seq) and os.path.exists(ref_seq):
                    if "# misassemblies" in line:
                        sample_stats_dict["MISASSEMBLIES"] = float(line.split("\t")[-1].strip())
                    elif "# N's per 100 kbp" in line:
                        sample_stats_dict["N_PER_100KBP"] = float(line.split("\t")[-1].strip())
                    elif "# mismatches per 100 kbp" in line:
                        sample_stats_dict["MIS_PER_100KBP"] = float(line.split("\t")[-1].strip())
                    elif "# indels per 100 kbp" in line:
                        sample_stats_dict["INDELS_PER_100KPB"] = float(line.split("\t")[-1].strip())
    except FileNotFoundError:
        print(f"ERROR:\tQUAST report not found: {quast_report_tsv}")
        sample_stats_dict["GENOME_SIZE"] = None
        sample_stats_dict["ASSEMBLY_CONTIGS"] = None
        sample_stats_dict["ASSEMBLY_N50"] = None
        sample_stats_dict["ASSEMBLY_L50"] = None
        sample_stats_dict["ASSEMBLY_GC"] = None

    # --------------------------------------------------------------
    # Compute coverage if reference length is known
    # --------------------------------------------------------------
    try:
        if ref_seq and os.path.exists(ref_seq):
            ref_total_bases = 0
            for record in SeqIO.parse(ref_seq, "fasta"):
                ref_total_bases += len(record.seq)
    
            if not pd.isna(illu_raw_f_reads) and not pd.isna(illu_raw_r_reads):
                sample_stats_dict["RAW_ILLU_COVERAGE"] = round(sample_stats_dict["RAW_ILLU_TOTAL_BASES"] / ref_total_bases, 2)
                sample_stats_dict["TRIMMED_ILLU_COVERAGE"] = round(sample_stats_dict["TRIMMED_ILLU_TOTAL_BASES"] / ref_total_bases, 2)
                sample_stats_dict["DEDUPED_ILLU_COVERAGE"] = round(sample_stats_dict["DEDUPED_ILLU_TOTAL_BASES"] / ref_total_bases, 2)
            if not pd.isna(ont_raw_reads):
                sample_stats_dict["RAW_ONT_COVERAGE"] = round(sample_stats_dict["RAW_ONT_TOTAL_BASES"] / ref_total_bases, 2)
                sample_stats_dict["FILT_ONT_COVERAGE"] = round(sample_stats_dict["FILT_ONT_TOTAL_BASES"] / ref_total_bases, 2)
                sample_stats_dict["CORRECT_ONT_COVERAGE"] = round(sample_stats_dict["CORRECT_ONT_TOTAL_BASES"] / ref_total_bases, 2)
            if not pd.isna(pacbio_raw_reads):
                sample_stats_dict["RAW_PACBIO_COVERAGE"] = round(sample_stats_dict["RAW_PACBIO_TOTAL_BASES"] / ref_total_bases, 2)
                sample_stats_dict["HIFI_PACBIO_COVERAGE"] = round(sample_stats_dict["HIFI_PACBIO_TOTAL_BASES"] / ref_total_bases, 2) if "HIFI_PACBIO_TOTAL_BASES" in sample_stats_dict else None
                sample_stats_dict["FILT_PACBIO_COVERAGE"] = round(sample_stats_dict["FILT_PACBIO_TOTAL_BASES"] / ref_total_bases, 2)
        else:
            # If no reference, approximate coverage with QUAST's genome size
            if not pd.isna(illu_raw_f_reads) and not pd.isna(illu_raw_r_reads):
                if pd.notna(sample_stats_dict["GENOME_SIZE"]):
                    sample_stats_dict["RAW_ILLU_COVERAGE"] = round(sample_stats_dict["RAW_ILLU_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
                    sample_stats_dict["TRIMMED_ILLU_COVERAGE"] = round(sample_stats_dict["TRIMMED_ILLU_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
                    sample_stats_dict["DEDUPED_ILLU_COVERAGE"] = round(sample_stats_dict["DEDUPED_ILLU_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
            if not pd.isna(ont_raw_reads):
                if pd.notna(sample_stats_dict["GENOME_SIZE"]):
                    sample_stats_dict["RAW_ONT_COVERAGE"] = round(sample_stats_dict["RAW_ONT_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
                    sample_stats_dict["FILT_ONT_COVERAGE"] = round(sample_stats_dict["FILT_ONT_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
                    sample_stats_dict["CORRECT_ONT_COVERAGE"] = round(sample_stats_dict["CORRECT_ONT_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
            if not pd.isna(pacbio_raw_reads):
                if pd.notna(sample_stats_dict["GENOME_SIZE"]):
                    sample_stats_dict["RAW_PACBIO_COVERAGE"] = round(sample_stats_dict["RAW_PACBIO_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
                    sample_stats_dict["FILT_PACBIO_COVERAGE"] = round(sample_stats_dict["FILT_PACBIO_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
    except TypeError:
        print("SKIP:\tNot updating coverage as raw reads were not analyzed and updated.")

    # Annotate sample metadata
    sample_stats_dict["SPECIES_ID"] = sample_id.split("-")[0]
    sample_stats_dict["SAMPLE_ID"] = sample_id

    # Classify assembly quality
    results = classify_assembly(sample_stats_dict)

    # Prepare simpler summary
    first_busco_c = sample_stats_dict.get("FIRST_BUSCO_S", 0.0) + sample_stats_dict.get("FIRST_BUSCO_D", 0.0)
    second_busco_c = sample_stats_dict.get("SECOND_BUSCO_S", 0.0) + sample_stats_dict.get("SECOND_BUSCO_D", 0.0)
    n50 = sample_stats_dict.get("ASSEMBLY_N50", None)
    contig_count = sample_stats_dict.get("ASSEMBLY_CONTIGS", None)
    sample_stats_list = [first_busco_c, second_busco_c, n50, contig_count]

    stats_filepath = os.path.join(sample_dir, f"{os.path.basename(assembly_path).replace('.fasta', '')}_stats.txt")

    # Save a plain text stats file
    with open(stats_filepath, "w") as stats_file:
        stats_file.write("Assembly Statistics:\n")
        stats_file.write(f"First BUSCO Combined (S+D): {first_busco_c}%\n")
        stats_file.write(f"Second BUSCO Combined (S+D): {second_busco_c}%\n")
        stats_file.write(f"N50: {n50}\n")
        stats_file.write(f"Contig Count: {contig_count}\n")
        stats_file.write("\n")
        stats_file.write(str(results))  # results is a dict

    return assembly_path, sample_stats_list, sample_stats_dict


# --------------------------------------------------------------
# Perform final assembly assessment
# --------------------------------------------------------------
def final_assessment(assembly_type, input_csv, sample_id, output_dir, cpu_threads, ram_gb):
    """Perform final quality control and assessment of an assembly.

    Decompresses the assembly, runs BUSCO and QUAST, calculates coverage,
    classifies quality, and saves compressed assembly and statistics.

    Args:
        assembly_type (str): Descriptor for the assembly (e.g., 'final').
        input_csv (str): Path to metadata CSV file.
        sample_id (str): Sample identifier.
        output_dir (str): Directory for output files.
        cpu_threads (int or str): Number of CPU threads to use.
        ram_gb (int or str): Available RAM in GB.

    Returns:
        tuple: (path to compressed assembly, path to final stats CSV).
    """
    # Read the CSV file and filter to the row corresponding to the sample of interest
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
    first_busco_odb = current_series["BUSCO_1"]
    second_busco_odb = current_series["BUSCO_2"]
    kingdom_id = current_series["ORGANISM_KINGDOM"]
    karyote_id = current_series["ORGANISM_KARYOTE"]
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
        ref_seq_gz = ref_seq + ".gz"
        if os.path.exists(ref_seq_gz) and not os.path.exists(ref_seq):
            _ = pigz_decompress(ref_seq_gz, cpu_threads)

    # Set Illumina deduplicated read paths only if Illumina reads are present
    illu_dedup_f_reads = None
    illu_dedup_r_reads = None
    if pd.notna(illumina_f_raw_reads) and pd.notna(illumina_r_raw_reads):
        illu_dedup_f_reads = os.path.join(species_dir, "Illumina", f"{species_id}_illu_forward_dedup.fastq")
        illu_dedup_f_reads_gz = illu_dedup_f_reads + ".gz"
        illu_dedup_r_reads = os.path.join(species_dir, "Illumina", f"{species_id}_illu_reverse_dedup.fastq")
        illu_dedup_r_reads_gz = illu_dedup_r_reads + ".gz"
        if os.path.exists(illu_dedup_f_reads_gz) and not os.path.exists(illu_dedup_f_reads):
            _ = pigz_decompress(illu_dedup_f_reads_gz, cpu_threads)
        if os.path.exists(illu_dedup_r_reads_gz) and not os.path.exists(illu_dedup_r_reads):
            _ = pigz_decompress(illu_dedup_r_reads_gz, cpu_threads)
        if not os.path.exists(illu_dedup_f_reads) or not os.path.exists(illu_dedup_r_reads):
            print(f"ERROR:\tIllumina deduplicated reads not found: {illu_dedup_f_reads}, {illu_dedup_r_reads}")
            return None, None

    labeled_assembly = None
    if pd.isna(est_size):
        print(f"Processing assembly for quality control only: {ref_seq}")
        assembly_path = ref_seq
        sample_dir = os.path.dirname(assembly_path)
        labeled_assembly = os.path.join(sample_dir, f"{sample_id}_EGAP_assembly.fasta")
    else:
        assembly_path = os.path.join(sample_dir, f"{sample_id}_final_curated.fasta")
        sample_dir = os.path.dirname(assembly_path)
        labeled_assembly = os.path.join(sample_dir, f"{sample_id}_final_EGAP_assembly.fasta")

    # Ensure working directory is sample_dir
    os.makedirs(sample_dir, exist_ok=True)
    os.chdir(sample_dir)
    print(f"DEBUG - Set working directory to: {os.getcwd()}")
    
    # ---- NEW: use the final (homogeneous) filename BEFORE any QC ----
    def _promote_to_final(src_path: str, dst_path: str, copy_only: bool) -> str:
        """Ensure QC runs on the final-labeled path (rename or copy)."""
        if src_path == dst_path and os.path.exists(dst_path):
            return dst_path
        if not os.path.exists(src_path) and os.path.exists(src_path + ".gz"):
            print(f"INFO:\tDecompressing {src_path}.gz ...")
            _ = pigz_decompress(src_path + ".gz", cpu_threads)
        if not os.path.exists(src_path):
            raise FileNotFoundError(f"Assembly path not found: {src_path}")
        os.makedirs(os.path.dirname(dst_path), exist_ok=True)
        if copy_only:
            print(f"INFO:\tCopying assembly -> {dst_path}")
            shutil.copy2(src_path, dst_path)
        else:
            try:
                print(f"INFO:\tRenaming assembly -> {dst_path}")
                os.replace(src_path, dst_path)
            except OSError:
                print("WARN:\tRename failed; copying instead.")
                shutil.copy2(src_path, dst_path)
        return dst_path
    
    # QC-only mode (no curated final): copy the ref to a stable EGAP name
    if pd.isna(est_size):
        if not ref_seq or not os.path.exists(ref_seq):
            print(f"ERROR:\tReference path not found: {ref_seq}")
            return None, None
        if not os.path.exists(labeled_assembly):
            assembly_path = _promote_to_final(ref_seq, labeled_assembly, copy_only=True)
        else:
            assembly_path = labeled_assembly
    else:
        # True final mode: rename curated (or polished) to the final EGAP name
        if not os.path.exists(assembly_path) and os.path.exists(assembly_path + ".gz"):
            assembly_path = pigz_decompress(assembly_path + ".gz", cpu_threads)
        if not os.path.exists(assembly_path):
            # try polished fallback before giving up
            polished_assembly = os.path.join(sample_dir, f"{sample_id}_final_polish_assembly.fasta")
            if os.path.exists(polished_assembly) and validate_fasta(polished_assembly):
                assembly_path = polished_assembly
            else:
                print(f"ERROR:\tAssembly path not found: {assembly_path}")
                return None, None
        assembly_path = _promote_to_final(assembly_path, labeled_assembly, copy_only=False)
    
    # Validate final-named assembly
    if not validate_fasta(assembly_path):
        print(f"ERROR:\tInvalid assembly after promotion: {assembly_path}")
        return None, None
    
    print(f"Parsing final assembly for index {current_index} from {input_csv}:\n{current_row}")
    
    # -------- Compleasm/BUSCO (now runs on the final-labeled path) --------
    run_lineage_eval(assembly_path, sample_id, sample_stats_dict,
                     "first", first_busco_odb, assembly_type, cpu_threads)
    run_lineage_eval(assembly_path, sample_id, sample_stats_dict,
                     "second", second_busco_odb, assembly_type, cpu_threads)
    
    # -------- QUAST (outputs now match the final-labeled basename) --------
    quast_dir = os.path.join(sample_dir, f"{os.path.basename(assembly_path).replace('.fasta', '')}_quast")
    os.makedirs(quast_dir, exist_ok=True)
    quast_report_tsv = os.path.join(quast_dir, "report.tsv")
    
    if os.path.exists(quast_report_tsv):
        print(f"SKIP:\tQUAST Report already exists: {quast_report_tsv}.")
    else:
        quast_cmd = ["quast", "--threads", str(cpu_threads)]
        if karyote_id == "eukaryote":
            quast_cmd.append("--eukaryote")
        if kingdom_id == "Funga":
            quast_cmd.append("--fungus")
        if pd.notna(ref_seq) and os.path.exists(ref_seq):
            quast_cmd.extend(["-r", ref_seq])
        quast_cmd.extend(["-o", quast_dir, assembly_path])
        print(f"DEBUG - Running QUAST: {' '.join(quast_cmd)}")
        result = run_subprocess_cmd(quast_cmd, shell_check=False)
        if result != 0:
            print(f"WARN:\tQUAST failed with return code {result}. Skipping QUAST metrics.")
            return None, None, sample_stats_dict
    
    # Parse QUAST report to populate sample_stats_dict
    try:
        with open(quast_report_tsv, "r") as quast_file:
            for line in quast_file:
                if "Total length (>= 0 bp)" in line:
                    sample_stats_dict["GENOME_SIZE"] = float(line.split("\t")[-1].strip())
                elif "# contigs" in line:
                    sample_stats_dict["ASSEMBLY_CONTIGS"] = float(line.split("\t")[-1].strip())
                elif "N50" in line:
                    sample_stats_dict["ASSEMBLY_N50"] = float(line.split("\t")[-1].strip())
                elif "L50" in line:
                    sample_stats_dict["ASSEMBLY_L50"] = float(line.split("\t")[-1].strip())
                elif "GC (%)" in line:
                    sample_stats_dict["ASSEMBLY_GC"] = float(line.split("\t")[-1].strip())
                if pd.notna(ref_seq) and os.path.exists(ref_seq):
                    if "# misassemblies" in line:
                        sample_stats_dict["MISASSEMBLIES"] = float(line.split("\t")[-1].strip())
                    elif "# N's per 100 kbp" in line:
                        sample_stats_dict["N_PER_100KBP"] = float(line.split("\t")[-1].strip())
                    elif "# mismatches per 100 kbp" in line:
                        sample_stats_dict["MIS_PER_100KBP"] = float(line.split("\t")[-1].strip())
                    elif "# indels per 100 kbp" in line:
                        sample_stats_dict["INDELS_PER_100KPB"] = float(line.split("\t")[-1].strip())
    except FileNotFoundError:
        print(f"ERROR:\tQUAST report not found: {quast_report_tsv}")
        sample_stats_dict["GENOME_SIZE"] = None
        sample_stats_dict["ASSEMBLY_CONTIGS"] = None
        sample_stats_dict["ASSEMBLY_N50"] = None
        sample_stats_dict["ASSEMBLY_L50"] = None
        sample_stats_dict["ASSEMBLY_GC"] = None

    # --------------------------------------------------------------
    # Compute coverage if reference length is known
    # --------------------------------------------------------------
    try:
        if ref_seq and os.path.exists(ref_seq):
            ref_total_bases = 0
            for record in SeqIO.parse(ref_seq, "fasta"):
                ref_total_bases += len(record.seq)

            if not pd.isna(illumina_f_raw_reads) and not pd.isna(illumina_r_raw_reads):
                sample_stats_dict["RAW_ILLU_COVERAGE"] = round(sample_stats_dict["RAW_ILLU_TOTAL_BASES"] / ref_total_bases, 2)
                sample_stats_dict["TRIMMED_ILLU_COVERAGE"] = round(sample_stats_dict["TRIMMED_ILLU_TOTAL_BASES"] / ref_total_bases, 2)
                sample_stats_dict["DEDUPED_ILLU_COVERAGE"] = round(sample_stats_dict["DEDUPED_ILLU_TOTAL_BASES"] / ref_total_bases, 2)
            if not pd.isna(ont_raw_reads):
                sample_stats_dict["RAW_ONT_COVERAGE"] = round(sample_stats_dict["RAW_ONT_TOTAL_BASES"] / ref_total_bases, 2)
                sample_stats_dict["FILT_ONT_COVERAGE"] = round(sample_stats_dict["FILT_ONT_TOTAL_BASES"] / ref_total_bases, 2)
                sample_stats_dict["CORRECT_ONT_COVERAGE"] = round(sample_stats_dict["CORRECT_ONT_TOTAL_BASES"] / ref_total_bases, 2)
            if not pd.isna(pacbio_raw_reads):
                sample_stats_dict["RAW_PACBIO_COVERAGE"] = round(sample_stats_dict["RAW_PACBIO_TOTAL_BASES"] / ref_total_bases, 2)
                sample_stats_dict["HIFI_PACBIO_COVERAGE"] = round(sample_stats_dict["HIFI_PACBIO_TOTAL_BASES"] / ref_total_bases, 2) if "HIFI_PACBIO_TOTAL_BASES" in sample_stats_dict else None
                sample_stats_dict["FILT_PACBIO_COVERAGE"] = round(sample_stats_dict["FILT_PACBIO_TOTAL_BASES"] / ref_total_bases, 2)
        else:
            # If no reference, approximate coverage with QUAST's genome size
            if not pd.isna(illumina_f_raw_reads) and not pd.isna(illumina_r_raw_reads):
                if pd.notna(sample_stats_dict["GENOME_SIZE"]):
                    sample_stats_dict["RAW_ILLU_COVERAGE"] = round(sample_stats_dict["RAW_ILLU_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
                    sample_stats_dict["TRIMMED_ILLU_COVERAGE"] = round(sample_stats_dict["TRIMMED_ILLU_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
                    sample_stats_dict["DEDUPED_ILLU_COVERAGE"] = round(sample_stats_dict["DEDUPED_ILLU_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
            if not pd.isna(ont_raw_reads):
                if pd.notna(sample_stats_dict["GENOME_SIZE"]):
                    sample_stats_dict["RAW_ONT_COVERAGE"] = round(sample_stats_dict["RAW_ONT_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
                    sample_stats_dict["FILT_ONT_COVERAGE"] = round(sample_stats_dict["FILT_ONT_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
                    sample_stats_dict["CORRECT_ONT_COVERAGE"] = round(sample_stats_dict["CORRECT_ONT_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
            if not pd.isna(pacbio_raw_reads):
                if pd.notna(sample_stats_dict["GENOME_SIZE"]):
                    sample_stats_dict["RAW_PACBIO_COVERAGE"] = round(sample_stats_dict["RAW_PACBIO_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
                    sample_stats_dict["FILT_PACBIO_COVERAGE"] = round(sample_stats_dict["FILT_PACBIO_TOTAL_BASES"] / sample_stats_dict["GENOME_SIZE"], 2)
    except TypeError:
        print("SKIP:\tNot updating coverage as raw reads were not analyzed and updated.")

    # Annotate sample metadata
    sample_stats_dict["SPECIES_ID"] = sample_id.split("-")[0]
    sample_stats_dict["SAMPLE_ID"] = sample_id

    # Classify assembly quality
    results = classify_assembly(sample_stats_dict)

    # Prepare simpler summary
    first_busco_c = sample_stats_dict.get("FIRST_BUSCO_S", 0.0) + sample_stats_dict.get("FIRST_BUSCO_D", 0.0)
    second_busco_c = sample_stats_dict.get("SECOND_BUSCO_S", 0.0) + sample_stats_dict.get("SECOND_BUSCO_D", 0.0)
    n50 = sample_stats_dict.get("ASSEMBLY_N50", None)
    contig_count = sample_stats_dict.get("ASSEMBLY_CONTIGS", None)
    sample_stats_list = [first_busco_c, second_busco_c, n50, contig_count]

    stats_filepath = os.path.join(sample_dir, f"{os.path.basename(labeled_assembly).replace('.fasta', '')}_stats.txt")

    if not os.path.exists(ref_seq) and pd.isna(est_size):
        shutil.copy(labeled_assembly, ref_seq)

    # Save a plain text stats file
    with open(stats_filepath, "w") as stats_file:
        stats_file.write("Assembly Statistics:\n")
        stats_file.write(f"First BUSCO Combined (S+D): {first_busco_c}%\n")
        stats_file.write(f"Second BUSCO Combined (S+D): {second_busco_c}%\n")
        stats_file.write(f"N50: {n50}\n")
        stats_file.write(f"Contig Count: {contig_count}\n")
        stats_file.write("\n")
        stats_file.write(str(results))  # results is a dict

    # Save sample_stats_dict to a CSV for final reference
    final_stats_csv = labeled_assembly.replace(".fasta", "_final_stats.csv")
    pd.DataFrame([sample_stats_dict]).to_csv(final_stats_csv, index=False)
    print(f"PASS: Full stats CSV saved: {final_stats_csv}")

    print(f"PASS:\tAssembly Stats for {labeled_assembly}:\n{sample_stats_list}")

    # # Walk through directory and subdirectories and multi-thread compress ALL FASTA or FASTQ files
    # for root, dirs, files in os.walk(sample_dir):
    #     for file in files:
    #         if file.endswith(('.fasta', '.fastq')):
    #             full_path = os.path.join(root, file)
    #             print(f"Compressing: {full_path}")
    #             _ = pigz_compress(full_path, cpu_threads)

    return labeled_assembly, final_stats_csv


if __name__ == "__main__":
    # Handle command-line arguments
    if len(sys.argv) != 7:
        print("Usage: python3 qc_assessment.py <assembly_type> <input_csv> "
              "<sample_id> <output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)

    if sys.argv[1] == "final":
        labeled_assembly, final_stats_csv = final_assessment(sys.argv[1],       # assembly_type
                                                            sys.argv[2],       # input_csv
                                                            sys.argv[3],       # sample_id
                                                            sys.argv[4],       # output_dir
                                                            str(sys.argv[5]),  # cpu_threads
                                                            str(sys.argv[6]))  # ram_gb
    else:
        assembly_path, sample_stats_list, sample_stats_dict = qc_assessment(sys.argv[1],       # assembly_type
                                                                            sys.argv[2],       # input_csv
                                                                            sys.argv[3],       # sample_id
                                                                            sys.argv[4],       # output_dir
                                                                            str(sys.argv[5]),  # cpu_threads
                                                                            str(sys.argv[6]))  # ram_gb
