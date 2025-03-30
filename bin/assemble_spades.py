#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
assemble_spades.py

Updated on Sat Mar 29 2025

This script runs SPAdes assembly with Illumina and optional long reads.

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

def assemble_spades(illumina_f, illumina_r, long_reads, species_id, ref_seq, cpu_threads, ram_gb):
    if illumina_f == "None" or illumina_r == "None":
        print("SKIP:\tSPAdes assembly; no Illumina reads provided")
        return

    out_dir = "spades_assembly"
    final_assembly = f"{species_id}_spades.fasta"
    if os.path.exists(final_assembly):
        print(f"SKIP:\tFinal SPAdes Assembly already exists: {final_assembly}")
        return

    os.makedirs(out_dir, exist_ok=True)
    spades_cmd = ["spades.py", "--careful", "-1", illumina_f, "-2", illumina_r,
                  "-o", out_dir, "-t", str(cpu_threads), "-m", str(ram_gb), "--cov-cutoff", "auto"]
    if long_reads != "None":
        spades_cmd += ["--nanopore" if "ont" in long_reads else "--pacbio", long_reads]
    if ref_seq != "None":
        spades_cmd += ["--trusted-contigs", ref_seq]
    spades_cmd += ["-k", "21,33,55,77,99,127"]

    run_subprocess_cmd(spades_cmd, False)
    if os.path.exists(f"{out_dir}/scaffolds.fasta"):
        shutil.move(f"{out_dir}/scaffolds.fasta", final_assembly)

if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("Usage: python assemble_spades.py <illumina_f> <illumina_r> <long_reads> <species_id> <ref_seq> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)
    assemble_spades(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], int(sys.argv[6]), int(sys.argv[7]))