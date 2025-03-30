#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare_assemblies.py

Updated on Sat Mar 29 2025

This script compares assemblies from MaSuRCA, SPAdes, Flye, and Hifiasm, selecting the best based on Compleasm completeness, contig count, and N50.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import os
import sys
import subprocess
import shutil
from collections import Counter

def run_subprocess_cmd(cmd_list, shell_check):
    """Run a subprocess command and print its output."""
    cmd_str = ' '.join(cmd_list) if isinstance(cmd_list, list) else cmd_list
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

def get_compleasm_score(assembly, db, cpu_threads):
    """Run Compleasm and extract completeness score."""
    if assembly == "None" or not os.path.exists(assembly):
        return None
    out_dir = f"{os.path.basename(assembly)}_{db}_compleasm"
    os.makedirs(out_dir, exist_ok=True)
    summary = f"{out_dir}/summary.txt"
    if not os.path.exists(summary):
        run_subprocess_cmd(["compleasm", "run", "-a", assembly, "-o", out_dir, "--lineage", f"{db}_odb10", "-t", str(cpu_threads)], False)
    with open(summary, "r") as f:
        for line in f:
            if "Complete" in line:  # Assuming format: "Complete: XX.X%"
                return float(line.split()[1])
    return None

def get_quast_stats(assembly, cpu_threads):
    """Run Quast and extract contig count and N50."""
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

def compare_assemblies(masurca, spades, flye, hifiasm, species_id, compleasm_1, compleasm_2, karyote, kingdom, cpu_threads):
    """Compare assemblies and select the best based on EGAP.py logic."""
    assemblies = {"MaSuRCA": masurca, "SPAdes": spades, "Flye": flye, "Hifiasm": hifiasm}
    stats_dict = {}

    # Gather stats for each assembly
    for method, assembly in assemblies.items():
        if assembly != "None" and os.path.exists(assembly):
            comp1 = get_compleasm_score(assembly, compleasm_1, cpu_threads)
            comp2 = get_compleasm_score(assembly, compleasm_2, cpu_threads)
            contigs, n50 = get_quast_stats(assembly, cpu_threads)
            stats_dict[method] = [comp1, comp2, contigs, n50] if all(v is not None for v in [comp1, comp2, contigs, n50]) else [None, None, None, None]
        else:
            stats_dict[method] = [None, None, None, None]

    # Transpose stats for comparison
    stats_combined = list(zip(stats_dict["MaSuRCA"], stats_dict["Flye"], stats_dict["SPAdes"], stats_dict["Hifiasm"]))
    custom_stats = []
    methods = ["MaSuRCA", "Flye", "SPAdes", "Hifiasm"]

    # Compare each metric
    for index, values in enumerate(stats_combined):
        valid_values = [(val, methods[i]) for i, val in enumerate(values) if val is not None]
        if not valid_values:
            custom_stats.append("Unknown")
            continue
        if index == 2:  # Contig count (minimize)
            _, best_method = min(valid_values, key=lambda x: x[0])
        else:  # Compleasm completeness and N50 (maximize)
            _, best_method = max(valid_values, key=lambda x: x[0])
        custom_stats.append(best_method)

    # Determine most frequent method
    method_counts = Counter(m for m in custom_stats if m != "Unknown")
    most_represented_method = max(method_counts, key=method_counts.get) if method_counts else "Unknown"

    if most_represented_method == "Unknown":
        print("ERROR:\tNo valid assemblies to compare")
        sys.exit(1)

    # Output results
    best_stats = stats_dict[most_represented_method]
    print(f"Best Initial Assembly: {most_represented_method}.")
    labels = ["First Compleasm Completeness (Single + Duplicated)", "Second Compleasm Completeness (Single + Duplicated)",
              "Assembly Contig Count", "Assembly N50"]
    for index, item in enumerate(best_stats):
        print(f"{labels[index]}: {item}")

    # Copy best assembly
    final_assembly = f"{species_id}_best_assembly.fasta"
    shutil.copy(assemblies[most_represented_method], final_assembly)

if __name__ == "__main__":
    if len(sys.argv) != 11:
        print("Usage: python compare_assemblies.py <masurca> <spades> <flye> <hifiasm> <species_id> <compleasm_1> <compleasm_2> <karyote> <kingdom> <cpu_threads>", file=sys.stderr)
        sys.exit(1)
    compare_assemblies(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], int(sys.argv[10]))