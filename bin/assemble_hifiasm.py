#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
assemble_hifiasm.py

Updated on Sat Mar 29 2025

This script runs Hifiasm assembly with PacBio reads.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import os
import sys
import subprocess

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

def assemble_hifiasm(long_reads, species_id, cpu_threads):
    if long_reads == "None" or "pacbio" not in long_reads.lower():
        print("SKIP:\tHifiasm assembly; no PacBio reads provided")
        return

    final_assembly = f"{species_id}_hifiasm.fasta"
    if os.path.exists(final_assembly):
        print(f"SKIP:\tFinal Hifiasm Assembly already exists: {final_assembly}")
        return

    prefix = "hifiasm_assembly"
    hifiasm_cmd = ["hifiasm", "-o", prefix, "-t", str(cpu_threads), long_reads]
    run_subprocess_cmd(hifiasm_cmd, False)

    gfa_file = f"{prefix}.asm.bp.p_ctg.gfa"
    if os.path.exists(gfa_file):
        run_subprocess_cmd(["gfatools", "gfa2fa", gfa_file, ">", final_assembly], True)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python assemble_hifiasm.py <long_reads> <species_id> <cpu_threads>", file=sys.stderr)
        sys.exit(1)
    assemble_hifiasm(sys.argv[1], sys.argv[2], int(sys.argv[3]))