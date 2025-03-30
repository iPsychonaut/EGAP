#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
assemble_masurca.py

Updated on Sat Mar 29 2025

This script runs MaSuRCA assembly with Illumina and optional long reads.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import os
import sys
import subprocess
import shutil

def run_subprocess_cmd(cmd_list, shell_check):
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

def assemble_masurca(illumina_f, illumina_r, long_reads, species_id, est_size, ref_seq, cpu_threads, ram_gb):
    if illumina_f == "None" or illumina_r == "None":
        print("SKIP:\tMaSuRCA assembly; no Illumina reads provided")
        return

    out_dir = "masurca_assembly"
    final_assembly = f"{species_id}_masurca.fasta"
    if os.path.exists(final_assembly):
        print(f"SKIP:\tFinal MaSuRCA Assembly already exists: {final_assembly}")
        return

    os.makedirs(out_dir, exist_ok=True)
    os.chdir(out_dir)

    # Parse genome size
    size_str = est_size.lower().strip() if est_size != "None" else "25m"
    multipliers = {'m': 10**6, 'g': 10**9}
    est_size_bp = int(float(size_str[:-1]) * multipliers[size_str[-1]]) if size_str[-1] in multipliers else 25000000

    # Simplified config (full logic from EGAP.py can be expanded)
    config_content = ["DATA\n"]
    config_content.append(f"PE= pe 251 30 {illumina_f} {illumina_r}\n")
    if long_reads != "None":
        config_content.append(f"NANOPORE={long_reads}\n" if "ont" in long_reads else f"PACBIO={long_reads}\n")
    if ref_seq != "None":
        config_content.append(f"REFERENCE={ref_seq}\n")
    config_content.append("END\nPARAMETERS\n")
    config_content.append(f"NUM_THREADS={cpu_threads}\n")
    config_content.append(f"JF_SIZE={est_size_bp}\n")
    config_content.append("SOAP_ASSEMBLY=0\nFLYE_ASSEMBLY=0\n")
    config_content.append("END\n")

    with open("masurca_config.txt", "w") as f:
        f.writelines(config_content)

    run_subprocess_cmd(["masurca", "masurca_config.txt"], False)
    run_subprocess_cmd(["bash", "assemble.sh"], False)

    default_assembly = "CA/9-terminator/genome.scf.fasta"
    if os.path.exists(default_assembly):
        shutil.move(default_assembly, f"../{final_assembly}")
    os.chdir("..")

if __name__ == "__main__":
    if len(sys.argv) != 9:
        print("Usage: python assemble_masurca.py <illumina_f> <illumina_r> <long_reads> <species_id> <est_size> <ref_seq> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)
    assemble_masurca(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], int(sys.argv[7]), int(sys.argv[8]))