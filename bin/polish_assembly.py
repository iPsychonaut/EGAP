#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polish_assembly.py

This module processes genomic assembly data by applying polishing steps using Racon
for long reads (ONT or PacBio) and Polypolish + Pypolca for Illumina reads. It handles input validation,
subprocess execution, and file management for assembly refinement.

Stage:
    Polishing (minimap2 + Racon, for long reads)
    Polishing (bwa + Polypolish + Pypolca, for Illumina)

Created on Wed Aug 16 2023

Updated on Wed Sept 3 2025

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os
import sys
import shutil
from typing import Optional

import pandas as pd
from pathlib import Path  # <<< NEW
from utilities import run_subprocess_cmd, initialize_logging_environment, load_sample_context
from file_manager import remove_file, remove_dir, remove_glob


# --------------------------------------------------------------
# Polish assembly with Racon using long reads
# --------------------------------------------------------------
def racon_polish_assembly(input_assembly, long_reads, racon_out_dir, sample_id, cpu_threads, iteration_count):
    """Polish a genomic assembly with Racon using long reads."""
    if not os.path.exists(input_assembly):
        print(f"ERROR:\tInput assembly not found: {input_assembly}")
        return input_assembly
    if not os.path.exists(long_reads):
        print(f"ERROR:\tLong reads not found: {long_reads}")
        return input_assembly

    racon_paf = os.path.join(racon_out_dir, f"racon_round{iteration_count}.paf")
    if os.path.exists(racon_paf):
        print(f"SKIP:\tRacon PAF {iteration_count} already exists: {racon_paf}.")
    else:
        minimap2_cmd = f"minimap2 -t {cpu_threads} -x map-ont {input_assembly} {long_reads} > {racon_paf}"
        print(f"DEBUG - Running minimap2: {minimap2_cmd}")
        result = run_subprocess_cmd(minimap2_cmd, shell_check=True)
        if result != 0 or not os.path.exists(racon_paf):
            print(f"WARN:\tminimap2 failed for iteration {iteration_count}. Skipping Racon.")
            return input_assembly

    racon_assembly = os.path.join(racon_out_dir, f"{sample_id}_racon_polish_{iteration_count}.fasta")
    if os.path.exists(racon_assembly):
        print(f"SKIP:\tRacon Assembly {iteration_count} already exists: {racon_assembly}.")
    else:
        racon_cmd = f"racon -t {cpu_threads} {long_reads} {racon_paf} {input_assembly} > {racon_assembly}"
        print(f"DEBUG - Running Racon: {racon_cmd}")
        result = run_subprocess_cmd(racon_cmd, shell_check=True)
        if result != 0 or not os.path.exists(racon_assembly):
            print(f"WARN:\tRacon failed for iteration {iteration_count}. Returning input assembly.")
            return input_assembly

    return racon_assembly


# --------------------------------------------------------------
# Polish assembly with Illumina reads: Polypolish, then Pypolca
# --------------------------------------------------------------
def short_read_polish(input_assembly, illu_f_dedup, illu_r_dedup, work_dir, sample_id, cpu_threads):
    """Polish with Illumina reads: Polypolish, then Pypolca --careful.

    Replaced Pilon after the speed/memory A/B (K. pneumoniae KPC2_202310,
    graded against its same-strain reference): Racon x2 + Polypolish + Pypolca
    left ~6 errors against ~12 for Racon x2 + Pilon, with a faster short-read
    step and ~30% less peak memory. Polypolish needs every alignment of each
    read (``bwa mem -a``), with the two read files aligned separately.

    Returns the Pypolca assembly, the Polypolish assembly if Pypolca fails, or
    the input assembly if Polypolish fails.
    """
    for path in (input_assembly, illu_f_dedup, illu_r_dedup):
        if not path or not os.path.exists(path):
            print(f"ERROR:\tShort-read polishing input not found: {path}")
            return input_assembly
    os.makedirs(work_dir, exist_ok=True)
    pypolca_dir = os.path.join(work_dir, "pypolca")
    pypolca_out = os.path.join(pypolca_dir, "pypolca_corrected.fasta")
    if os.path.exists(pypolca_out):
        print(f"SKIP:\tPypolca assembly already exists: {pypolca_out}.")
        return pypolca_out

    asm = os.path.join(work_dir, f"{sample_id}_short_polish_input.fasta")
    shutil.copy(input_assembly, asm)
    sam1, sam2 = os.path.join(work_dir, "r1.sam"), os.path.join(work_dir, "r2.sam")
    filt1, filt2 = os.path.join(work_dir, "r1.filtered.sam"), os.path.join(work_dir, "r2.filtered.sam")
    polypolish_out = os.path.join(work_dir, f"{sample_id}_polypolish.fasta")
    for cmd in (f"bwa index {asm}",
                f"bwa mem -t {cpu_threads} -a {asm} {illu_f_dedup} > {sam1}",
                f"bwa mem -t {cpu_threads} -a {asm} {illu_r_dedup} > {sam2}",
                f"polypolish filter --in1 {sam1} --in2 {sam2} --out1 {filt1} --out2 {filt2}",
                f"polypolish polish {asm} {filt1} {filt2} > {polypolish_out}"):
        print(f"DEBUG - Running: {cmd}")
        if run_subprocess_cmd(cmd, shell_check=True) != 0:
            print(f"WARN:\tPolypolish step failed: {cmd.split()[0]} {cmd.split()[1]}. Returning input assembly.")
            return input_assembly
    if not os.path.exists(polypolish_out) or os.path.getsize(polypolish_out) < 1000:
        print("WARN:\tPolypolish output missing or suspiciously small. Returning input assembly.")
        return input_assembly
    # The SAM files can be several GB; Polypolish is done with them.
    for p in (sam1, sam2, filt1, filt2):
        remove_file(p)

    pypolca_cmd = (f"pypolca run -a {polypolish_out} -1 {illu_f_dedup} -2 {illu_r_dedup} "
                   f"-t {cpu_threads} -o {pypolca_dir} --careful --force")
    print(f"DEBUG - Running Pypolca: {pypolca_cmd}")
    if run_subprocess_cmd(pypolca_cmd, shell_check=True) != 0 or not os.path.exists(pypolca_out):
        print("WARN:\tPypolca failed. Returning the Polypolish assembly.")
        return polypolish_out
    return pypolca_out


# --------------------------------------------------------------
# Perform full assembly polishing
# --------------------------------------------------------------
def polish_assembly(
    sample_id: str,
    input_tsv: str,
    output_dir: str,
    cpu_threads: int,
    ram_gb: int,
) -> Optional[str]:
    """Polish an assembly: Racon (two rounds, long reads), then Polypolish + Pypolca (Illumina)."""
    ctx = load_sample_context(sample_id, input_tsv, output_dir, cpu_threads, ram_gb)
    current_series = ctx.current_series

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

    print(f"DEBUG - ctx.output_dir      - {ctx.output_dir}")

    species_dir = os.path.join(ctx.output_dir, species_id)
    sample_dir = os.path.join(species_dir, sample_id)

    # ---------- FAST SKIP if final polished output exists ----------
    final_polished_assembly = os.path.join(sample_dir, f"{sample_id}_final_polish_assembly.fasta")
    if os.path.exists(final_polished_assembly) and os.path.getsize(final_polished_assembly) > 0:
        print(f"SKIP:\tFinal polished assembly already present: {final_polished_assembly}")
        return final_polished_assembly

    if pd.notna(ont_sra) and pd.isna(ont_raw_reads):
        ont_raw_reads = os.path.join(species_dir, "ONT", f"{ont_sra}.fastq")
    if pd.notna(illumina_sra) and pd.isna(illumina_f_raw_reads) and pd.isna(illumina_r_raw_reads):
        illumina_f_raw_reads = os.path.join(species_dir, "Illumina", f"{illumina_sra}_1.fastq")
        illumina_r_raw_reads = os.path.join(species_dir, "Illumina", f"{illumina_sra}_2.fastq")
    if pd.notna(pacbio_sra) and pd.isna(pacbio_raw_reads):
        pacbio_raw_reads = os.path.join(species_dir, "PacBio", f"{pacbio_sra}.fastq")
    if pd.notna(ref_seq_gca) and pd.isna(ref_seq):
        ref_seq = os.path.join(species_dir, "RefSeq", f"{species_id}_{ref_seq_gca}_RefSeq.fasta")

    print(f"DEBUG - illumina_sra        - {illumina_sra}")
    print(f"DEBUG - illumina_f_raw_reads- {illumina_f_raw_reads}")
    print(f"DEBUG - illumina_r_raw_reads- {illumina_r_raw_reads}")
    print(f"DEBUG - ont_sra             - {ont_sra}")
    print(f"DEBUG - ont_raw_reads       - {ont_raw_reads}")
    print(f"DEBUG - pacbio_sra          - {pacbio_sra}")
    print(f"DEBUG - pacbio_raw_reads    - {pacbio_raw_reads}")
    print(f"DEBUG - ref_seq_gca         - {ref_seq_gca}")
    print(f"DEBUG - ref_seq             - {ref_seq}")
    print(f"DEBUG - species_id          - {species_id}")
    print(f"DEBUG - est_size            - {est_size}")

    # Skip if no reads at all
    if (pd.isna(illumina_f_raw_reads) and pd.isna(illumina_sra)) and \
       (pd.isna(illumina_r_raw_reads) and pd.isna(illumina_sra)) and \
       (pd.isna(ont_raw_reads) and pd.isna(ont_sra)) and \
       (pd.isna(pacbio_raw_reads) and pd.isna(pacbio_sra)):
        print("SKIP:\tNo valid reads provided, required for assembly comparison.")
        return None

    # Create output directory for polishing (absolute)
    os.makedirs(sample_dir, exist_ok=True)
    polish_out_dir = os.path.join(sample_dir, "polished_assembly")
    os.makedirs(polish_out_dir, exist_ok=True)

    # Best assembly must exist
    best_assembly = os.path.join(sample_dir, f"{species_id}_best_assembly.fasta")
    print(f"DEBUG - best_assembly - {best_assembly}")
    if not os.path.exists(best_assembly):
        print(f"ERROR:\tBest assembly not found: {best_assembly}")
        return None

    # Illumina dedup paths only if present
    illu_dedup_f_reads = None
    illu_dedup_r_reads = None
    if pd.notna(illumina_sra) or (pd.notna(illumina_f_raw_reads) and pd.notna(illumina_r_raw_reads)):
        illu_dedup_f_reads = os.path.join(species_dir, "Illumina", f"{species_id}_illu_forward_dedup.fastq")
        illu_dedup_r_reads = os.path.join(species_dir, "Illumina", f"{species_id}_illu_reverse_dedup.fastq")

    print(f"DEBUG - illu_dedup_f_reads - {illu_dedup_f_reads}")
    print(f"DEBUG - illu_dedup_r_reads - {illu_dedup_r_reads}")

    # Long-read paths (prefer prefiltered)
    highest_mean_qual_long_reads = None
    if pd.notna(ont_raw_reads):
        print("DEBUG - ONT RAW READS EXIST!")
        candidate = os.path.join(species_dir, "ONT", f"{species_id}_ONT_highest_mean_qual_long_reads.fastq")
        highest_mean_qual_long_reads = candidate if os.path.exists(candidate) else ont_raw_reads
    elif pd.notna(pacbio_raw_reads):
        print("DEBUG - PACBIO RAW READS EXIST!")
        candidate = os.path.join(species_dir, "PacBio", f"{species_id}_PacBio_highest_mean_qual_long_reads.fastq")
        highest_mean_qual_long_reads = candidate if os.path.exists(candidate) else pacbio_raw_reads

    print(f"DEBUG - highest_mean_qual_long_reads - {highest_mean_qual_long_reads}")

    # Work inside the absolute polish_out_dir
    os.chdir(polish_out_dir)
    print(f"DEBUG - CWD (polish_out_dir) - {os.getcwd()}")

    # ---- Step 1: Two rounds of Racon polishing ----
    racon_work_dir = os.path.join(polish_out_dir, "racon_polish")
    os.makedirs(racon_work_dir, exist_ok=True)
    racon_final = os.path.join(polish_out_dir, f"{sample_id}_racon.fasta")

    if pd.notna(ont_raw_reads) or pd.notna(pacbio_raw_reads):
        print("2x Racon Polishing Long (ONT or PacBio) Reads...")
        initial_racon_assembly = racon_polish_assembly(best_assembly,
                                                      highest_mean_qual_long_reads,
                                                      racon_work_dir,
                                                      sample_id,
                                                      cpu_threads,
                                                      1)
        racon_assembly = racon_polish_assembly(initial_racon_assembly,
                                               highest_mean_qual_long_reads,
                                               racon_work_dir,
                                               sample_id,
                                               cpu_threads,
                                               2)
    else:
        print("SKIP:\t2x-Racon Polish; No Long Reads Provided.")
        racon_assembly = best_assembly

    if not os.path.exists(racon_assembly):
        print(f"ERROR:\tRacon output not found: {racon_assembly}. Falling back to best assembly.")
        racon_assembly = best_assembly
    if os.path.exists(racon_final):
        print(f"SKIP:\tRacon final already exists: {racon_final}")
    else:
        print(f"DEBUG - Copying Racon assembly: {racon_assembly} to {racon_final}")
        shutil.copy(racon_assembly, racon_final)
        if not os.path.exists(racon_final):
            print(f"ERROR:\tFailed to copy Racon assembly to {racon_final}")
            return None

    # ---- Step 2: Short-read polishing (if Illumina): Polypolish, then Pypolca ----
    short_work_dir = os.path.join(polish_out_dir, "short_read_polish")
    short_final = os.path.join(polish_out_dir, f"{sample_id}_short_polish.fasta")

    if pd.notna(illumina_f_raw_reads) and pd.notna(illumina_r_raw_reads):
        print("Polypolish + Pypolca Polishing with Illumina Reads...")
        polished_assembly = short_read_polish(racon_final, illu_dedup_f_reads, illu_dedup_r_reads,
                                              short_work_dir, sample_id, cpu_threads)
    else:
        print("SKIP:\tShort-read Polish; No Illumina Reads Provided.")
        polished_assembly = racon_final

    print(f"DEBUG - polished_assembly - {polished_assembly}")

    if not os.path.exists(polished_assembly):
        print(f"WARN:\tPolished assembly not found: {polished_assembly}. Falling back to Racon assembly.")
        polished_assembly = racon_final
        if not os.path.exists(polished_assembly):
            print(f"ERROR:\tRacon assembly not found: {polished_assembly}")
            return None

    if os.path.getsize(polished_assembly) < 1000:
        print(f"WARN:\tPolished assembly is suspiciously small: {polished_assembly}. Falling back to Racon assembly.")
        polished_assembly = racon_final

    if os.path.exists(short_final):
        print(f"SKIP:\tShort-read polished assembly already exists: {short_final}")
    else:
        print(f"DEBUG - Copying polished assembly: {polished_assembly} to {short_final}")
        shutil.copy(polished_assembly, short_final)
        if not os.path.exists(short_final):
            print(f"ERROR:\tFailed to copy polished assembly to {short_final}")
            return None

    print(f"DEBUG - short_final - {short_final}")

    # Prepare final polished assembly path (in the main assembly directory)
    final_polished_assembly = os.path.join(sample_dir, f"{sample_id}_final_polish_assembly.fasta")
    print(f"DEBUG - final_polished_assembly - {final_polished_assembly}")

    if os.path.exists(final_polished_assembly):
        print(f"SKIP:\tFinal polished assembly already exists: {final_polished_assembly}")
    else:
        if os.path.exists(short_final):
            print(f"DEBUG - Copying short-read polished assembly: {short_final} to {final_polished_assembly}")
            shutil.copy(short_final, final_polished_assembly)
            if not os.path.exists(final_polished_assembly):
                print(f"ERROR:\tFailed to copy short-read polished assembly to {final_polished_assembly}")
                return None
        else:
            print(f"WARN:\tShort-read polished assembly not found: {short_final}. Falling back to Racon assembly.")
            if os.path.exists(racon_final):
                print(f"DEBUG - Copying Racon final: {racon_final} to {final_polished_assembly}")
                shutil.copy(racon_final, final_polished_assembly)
                if not os.path.exists(final_polished_assembly):
                    print(f"ERROR:\tFailed to copy Racon final to {final_polished_assembly}")
                    return None
            else:
                print(f"ERROR:\tRacon final not found: {racon_final}")
                return None

    print(f"PASS:\tPolishing complete for {best_assembly} -> {final_polished_assembly}")

    # --- Cleanup intermediate files now that final polished assembly is confirmed ---
    if os.path.exists(final_polished_assembly) and os.path.getsize(final_polished_assembly) > 0:
        # Large minimap2 PAF alignments used by Racon
        remove_glob(os.path.join(racon_work_dir, "*.paf"))
        # Per-round Racon FASTA outputs (final is kept as racon_final)
        remove_file(os.path.join(racon_work_dir, f"{sample_id}_racon_polish_1.fasta"))
        remove_file(os.path.join(racon_work_dir, f"{sample_id}_racon_polish_2.fasta"))
        # Short-read polishing work dir (bwa index, Polypolish/Pypolca intermediates)
        remove_dir(short_work_dir)

    return final_polished_assembly


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python3 polish_assembly.py <sample_id> <input_tsv> "
              "<output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)

    initialize_logging_environment(sys.argv[3], sys.argv[1])

    final_polished_assembly = polish_assembly(
        sys.argv[1],       # sample_id
        sys.argv[2],       # input_tsv
        sys.argv[3],       # output_dir
        str(sys.argv[4]),  # cpu_threads
        str(sys.argv[5])   # ram_gb
    )
