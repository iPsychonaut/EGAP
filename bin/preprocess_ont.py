#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
preprocess_ont.py

This script preprocesses ONT reads with NanoPlot, Filtlong, and Ratatosk.
Handles raw directory concatenation and SRA downloads.

Stage:
    Filtering (Filtlong - ONT)
    Error Correction (Ratatosk)

Created on Wed Aug 16 2023

Updated on 2026-04-16

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os
import shutil
import sys
import subprocess
import glob
import re
from typing import Optional

import pandas as pd
from utilities import run_subprocess_cmd, select_long_reads, initialize_logging_environment, load_sample_context
from qc_assessment import nanoplot_qc_reads
from record_provenance import record_file

# String values written by older pipeline runs (na_rep="None") that must be
# treated as missing/null rather than real file paths or accessions.
_NULLS = {"", "none", "nan", "null", "na"}


def null(val):
    """Return None if *val* is NaN or a null-sentinel string, else return *val*."""
    if val is None:
        return None
    try:
        if pd.isna(val):
            return None
    except (TypeError, ValueError):
        pass
    if isinstance(val, str) and val.strip().lower() in _NULLS:
        return None
    return val


# --------------------------------------------------------------
# Preprocess ONT sequencing reads
# --------------------------------------------------------------
def preprocess_ont(
    sample_id: str,
    input_tsv: str,
    output_dir: str,
    cpu_threads: int,
    ram_gb: int,
) -> Optional[str]:
    """Preprocess ONT reads with NanoPlot, Filtlong, and Ratatosk.

    Handles concatenation of raw directory files, SRA downloads, quality control with
    NanoPlot, filtering with Filtlong, and correction with Ratatosk (if Illumina reads
    are available), then selects the highest quality reads.

    Args:
        sample_id (str): Sample identifier.
        input_tsv (str): Path to metadata TSV file.
        output_dir (str): Directory for output files.
        cpu_threads (int or str): Number of CPU threads to use.
        ram_gb (int or str): Available RAM in GB.

    Returns:
        str or None: Path to the highest quality ONT reads file, or None if no reads are available.
    """
    print(f"Preprocessing ONT reads for {sample_id}...")

    # Always anchor everything to absolute output_dir to avoid nested paths after chdir
    ctx = load_sample_context(sample_id, input_tsv, output_dir, cpu_threads, ram_gb)
    current_series = ctx.current_series
    # IMPORTANT: also pull sample_stats_dict out of the context, otherwise the
    # later ``sample_stats_dict = nanoplot_qc_reads(..., sample_stats_dict)``
    # call sites would trip Python's "local variable referenced before
    # assignment" rule (the LHS assignment marks the name local for the whole
    # function, hiding any outer ctx attribute access on the RHS).
    sample_stats_dict = ctx.sample_stats_dict

    # Identify read paths, reference, and BUSCO lineage info from TSV
    ont_raw_reads = null(current_series["ONT_RAW_READS"])
    ont_raw_dir   = null(current_series["ONT_RAW_DIR"])
    ont_sra       = null(current_series["ONT_SRA"])
    species_id    = current_series["SPECIES_ID"]
    est_size      = current_series["EST_SIZE"]

    # Recompute species_dir_abs now that species_id is known
    species_dir_abs = os.path.join(ctx.output_dir, species_id)
    ont_dir_abs = os.path.join(species_dir_abs, "ONT")
    os.makedirs(ont_dir_abs, exist_ok=True)

    # Input source is resolved in priority order (file > directory > SRA)
    # after the chdir below.

    # Illumina dedup (absolute paths so downstream tools can find them regardless of cwd)
    illu_dedup_f_reads = os.path.join(species_dir_abs, "Illumina", f"{species_id}_illu_forward_dedup.fastq")
    illu_dedup_r_reads = os.path.join(species_dir_abs, "Illumina", f"{species_id}_illu_reverse_dedup.fastq")

    print(f"DEBUG - ont_sra - {ont_sra}")
    print(f"DEBUG - ont_raw_reads - {ont_raw_reads}")
    print(f"DEBUG - ont_raw_dir - {ont_raw_dir}")
    print(f"DEBUG - illu_dedup_f_reads - {illu_dedup_f_reads}")
    print(f"DEBUG - illu_dedup_r_reads - {illu_dedup_r_reads}")

    if ont_sra is None and ont_raw_dir is None and ont_raw_reads is None:
        print("SKIP:\tONT preprocessing; no reads provided")
        return None

    # Work within the ONT directory
    prev_cwd = os.getcwd()
    os.chdir(ont_dir_abs)
    
    try:
        # Resolve the ONT input source in priority order:
        #   1) explicit ONT_RAW_READS file, if it is real and non-empty
        #   2) else concatenate a directory of *.fastq (ONT_RAW_DIR)
        #   3) else download ONT_SRA
        # Use the first source that resolves; fall through to the next otherwise.
        resolved_reads = None

        # 1) Explicit file
        if ont_raw_reads is not None:
            candidate = ont_raw_reads if os.path.isabs(ont_raw_reads) \
                else os.path.abspath(os.path.join(ctx.output_dir, ont_raw_reads))
            if os.path.exists(candidate) and os.path.getsize(candidate) > 0:
                resolved_reads = candidate
                print(f"NOTE:\tUsing provided ONT reads file: {resolved_reads}")
            else:
                print(f"WARN:\tONT_RAW_READS provided but not found/empty ({candidate}); trying directory, then SRA.")

        # 2) Directory of fastqs to concatenate
        if resolved_reads is None and ont_raw_dir is not None:
            ont_raw_dir_abs = ont_raw_dir if os.path.isabs(ont_raw_dir) else os.path.join(ctx.output_dir, ont_raw_dir)
            ont_files = sorted(glob.glob(os.path.join(ont_raw_dir_abs, "*.fastq")))
            if ont_files:
                print(f"NOTE:\tConcatenating {len(ont_files)} ONT file(s) from {ont_raw_dir_abs}")
                combined = os.path.join(ont_dir_abs, f"{species_id}_ont_combined.fastq")
                with open(combined, "wb") as w:
                    for f in ont_files:
                        with open(f, "rb") as r:
                            shutil.copyfileobj(r, w)
                resolved_reads = combined
            else:
                print(f"WARN:\tONT_RAW_DIR provided but no .fastq files in {ont_raw_dir_abs}; trying SRA.")

        # 3) SRA download
        if resolved_reads is None and ont_sra is not None:
            dumped = os.path.join(ont_dir_abs, f"{ont_sra}.fastq")
            if not os.path.exists(dumped) or os.path.getsize(dumped) == 0:
                print(f"Downloading SRA {ont_sra} from GenBank...")
                _ = run_subprocess_cmd(["prefetch", "--force", "yes", ont_sra], False)
                _ = run_subprocess_cmd(["fasterq-dump", "-e", str(cpu_threads), "-O", ont_dir_abs, ont_sra], False)
            if not os.path.exists(dumped) or os.path.getsize(dumped) == 0:
                print(f"ERROR:\tExpected FASTQ not found or empty after fasterq-dump: {dumped}")
                return None
            resolved_reads = dumped
            print(f"PASS:\tSRA converted to FASTQ: {resolved_reads}")

        if resolved_reads is None:
            print("SKIP:\tONT preprocessing; no usable reads resolved from file, directory, or SRA")
            return None

        ont_raw_reads = resolved_reads
    
        # Parse estimated genome size
        m = re.match(r"^(\d+(?:\.\d+)?)(\D+)$", str(est_size))
        if m:
            est_size_numb, est_size_mult = m.group(1), m.group(2)
            multipliers = {'m': 10**6, 'g': 10**9}
            est_size_bp = int(float(est_size_numb) * multipliers.get(est_size_mult.lower(), 25_000_000))
        else:
            print(f"NOTE:\tUnable to parse input estimated size {est_size}, using default: 25000000")
            est_size_bp = 25_000_000
    
        if not (isinstance(ont_raw_reads, str) and os.path.exists(ont_raw_reads)):
            print(f"ERROR:\tONT reads not found: {ont_raw_reads}")
            return None
    
        # NanoPlot Raw Reads (soft-guard)
        try:
            sample_stats_dict = nanoplot_qc_reads(ont_raw_reads, "Raw_ONT_", cpu_threads, sample_stats_dict, out_dir=ont_dir_abs)
        except Exception as e:
            print(f"WARN:\tNanoPlot failed on raw ONT reads ({e}); continuing without NanoStats.")
    
        # Filtlong (only include Illumina if those files actually exist).
        # Canonical filtered name regardless of which input source resolved.
        filtered_ont = os.path.join(ont_dir_abs, f"{species_id}_ont_filtered.fastq")
        coverage = 75
        target_bases = est_size_bp * coverage
        use_illumina = os.path.exists(illu_dedup_f_reads) and os.path.exists(illu_dedup_r_reads)
        illumina_opt = f"-1 {illu_dedup_f_reads} -2 {illu_dedup_r_reads}" if use_illumina else ""
    
        if not os.path.exists(filtered_ont):
            filtlong_cmd = f"filtlong {illumina_opt} --trim --min_length 1000 --min_mean_q 8 --keep_percent 90 --target_bases {target_bases} {ont_raw_reads} > {filtered_ont}"
            _ = run_subprocess_cmd(filtlong_cmd, True)
            if (not os.path.exists(filtered_ont)) or os.path.getsize(filtered_ont) == 0:
                print(f"ERROR:\tFiltlong did not produce reads at {filtered_ont}")
                return None
        else:
            print(f"SKIP\tFiltlong filtered reads exist: {filtered_ont}.")
    
        # NanoPlot Filtered Reads (soft-guard)
        try:
            sample_stats_dict = nanoplot_qc_reads(filtered_ont, "Filt_ONT_", cpu_threads, sample_stats_dict, out_dir=ont_dir_abs)
        except Exception as e:
            print(f"WARN:\tNanoPlot failed on filtered ONT reads ({e}); continuing.")
    
        # Ratatosk (only if Illumina reads available)
        corrected_out = os.path.join(ont_dir_abs, "ratatosk_corrected")
        if (not use_illumina):
            print("SKIP:\tRatatosk correction; no usable Illumina reads provided")
            final_corrected_ont = filtered_ont
            safe_alias = os.path.join(ont_dir_abs, f"{species_id}_ont_corrected.fastq")
            if not os.path.exists(safe_alias):
                subprocess.run(["ln", "-sf", filtered_ont, safe_alias], check=True)
        else:
            final_corrected_ont = os.path.join(ont_dir_abs, f"{species_id}_ont_corrected.fastq")
            if os.path.exists(corrected_out + ".fastq"):
                shutil.move(corrected_out + ".fastq", final_corrected_ont)
            if not os.path.exists(final_corrected_ont):
                ratatosk_cmd = ["Ratatosk", "correct", "-s", illu_dedup_f_reads, "-s", illu_dedup_r_reads,
                                "-l", filtered_ont, "-o", corrected_out, "-c", str(cpu_threads), "-v"]
                _ = run_subprocess_cmd(ratatosk_cmd, False)
                if os.path.exists(corrected_out + ".fastq"):
                    shutil.move(corrected_out + ".fastq", final_corrected_ont)
                else:
                    print(f"WARN:\tRatatosk did not produce {corrected_out}.fastq; using filtered reads.")
                    final_corrected_ont = filtered_ont
    
        # NanoPlot Corrected Reads (soft-guard)
        try:
            sample_stats_dict = nanoplot_qc_reads(final_corrected_ont, "Corr_ONT_", cpu_threads, sample_stats_dict, out_dir=ont_dir_abs)
        except Exception as e:
            print(f"WARN:\tNanoPlot failed on corrected ONT reads ({e}); continuing.")
    
        # Select best long reads (your helper uses output_dir/input_tsv paths, unchanged)
        highest_mean_qual_long_reads = select_long_reads(ctx.output_dir, ctx.input_tsv, sample_id, cpu_threads)
        highest = select_long_reads(ctx.output_dir, ctx.input_tsv, sample_id, cpu_threads)
        if not highest:
            highest = final_corrected_ont
    
    finally:
        os.chdir(prev_cwd)

    print(f"PASS:\tPreprocessed Raw ONT Reads for {sample_id}: {highest_mean_qual_long_reads}.")

    record_file("ONT highest-quality long reads", highest_mean_qual_long_reads)
    return highest_mean_qual_long_reads


if __name__ == "__main__":
    # Log raw sys.argv immediately
    print(f"DEBUG: Raw sys.argv = {sys.argv}")
    print(f"DEBUG: Length of sys.argv = {len(sys.argv)}")
    
    # Check argument count
    if len(sys.argv) != 6:
        print(f"ERROR: Expected 5 arguments (plus script name), got {len(sys.argv)-1}: {sys.argv[1:]}", 
              file=sys.stderr)
        print("Usage: python3 preprocess_ont.py <sample_id> <input_tsv> <output_dir> <cpu_threads> <ram_gb>", 
              file=sys.stderr)
        sys.exit(1)

    initialize_logging_environment(sys.argv[3], sys.argv[1])

    # Log each argument
    for i, arg in enumerate(sys.argv):
        print(f"DEBUG: sys.argv[{i}] = '{arg}'")
    
    sample_id = sys.argv[1]
    input_tsv = sys.argv[2]
    output_dir = sys.argv[3]
    cpu_threads = sys.argv[4]
    ram_gb = int(sys.argv[5]) if sys.argv[5] != " " else 8
    
    print(f"DEBUG: Parsed sample_id = '{sample_id}'")
    print(f"DEBUG: Parsed input_tsv = '{input_tsv}'")
    print(f"DEBUG: Parsed output_dir = '{output_dir}'")
    print(f"DEBUG: Parsed cpu_threads = '{sys.argv[4]}' (converted to {cpu_threads})")
    print(f"DEBUG: Parsed ram_gb = '{sys.argv[5]}' (converted to {ram_gb})")
    
    highest_mean_qual_long_reads = preprocess_ont(sample_id, input_tsv, output_dir, cpu_threads, ram_gb)
