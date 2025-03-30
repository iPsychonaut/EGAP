#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
assemble_flye.py

Updated on Sat Mar 29 2025

This script runs Flye assembly with long reads.

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

def assemble_flye(long_reads, species_id, est_size, cpu_threads):
    if long_reads == "None":
        print("SKIP:\tFlye assembly; no long reads provided")
        return

    out_dir = "flye_assembly"
    final_assembly = f"{species_id}_flye.fasta"
    if os.path.exists(final_assembly):
        print(f"SKIP:\tFinal Flye Assembly already exists: {final_assembly}")
        return

    os.makedirs(out_dir, exist_ok=True)
    size_str = est_size.lower().strip() if est_size != "None" else "25m"
    multipliers = {'m': 10**6, 'g': 10**9}
    est_size_bp = int(float(size_str[:-1]) * multipliers[size_str[-1]]) if size_str[-1] in multipliers else 25000000

    flye_cmd = ["flye", "--out-dir", out_dir, "--genome-size", str(est_size_bp), "--threads", str(cpu_threads),
                "--iterations", "3", "--keep-haplotypes"]
    flye_cmd += ["--nano-corr" if "ont" in long_reads else "--pacbio-corr", long_reads]

    run_subprocess_cmd(flye_cmd, False)
    if os.path.exists(f"{out_dir}/assembly.fasta"):
        shutil.move(f"{out_dir}/assembly.fasta", final_assembly)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python assemble_flye.py <long_reads> <species_id> <est_size> <cpu_threads>", file=sys.stderr)
        sys.exit(1)
    assemble_flye(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))