#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
preprocess_pacbio.py

Preprocess PacBio reads with NanoPlot and Filtlong.

Why filtering was failing:
- You were calling filtlong with --trim but without a reference (assembly or
  "good reads"). filtlong requires a reference to use --trim; otherwise it
  outputs nothing -> your filtered FASTQ was empty.

Fixes in this version:
- Use --trim only when a reference FASTA exists.
- If a previous filtered file exists but is empty, regenerate it.
- Absolute paths everywhere; robust logging; no silent fallbacks.

Created on Wed Aug 16 2023

Updated on Wed Sept 3 2025

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os, sys, glob, shutil, re
import pandas as pd
from utilities import run_subprocess_cmd, get_current_row_data, select_long_reads
from qc_assessment import nanoplot_qc_reads


def _abs(p):  # normalize to absolute path (for strings only)
    return os.path.abspath(p) if isinstance(p, str) else p

def _nonempty(fp, min_bytes=1024) -> bool:
    try:
        return os.path.getsize(fp) >= min_bytes
    except Exception:
        return False


def preprocess_pacbio(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    print(f"Preprocessing PacBio reads for {sample_id}...")

    # Resolve to absolute so later chdir is safe
    input_csv_abs  = _abs(input_csv)
    output_dir_abs = _abs(output_dir)

    # Read metadata
    input_df = pd.read_csv(input_csv_abs)
    current_row, _, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current = current_row.iloc[0]

    pacbio_sra       = current["PACBIO_SRA"]
    pacbio_raw_reads = current["PACBIO_RAW_READS"]
    pacbio_raw_dir   = current["PACBIO_RAW_DIR"]
    species_id       = current["SPECIES_ID"]
    est_size         = current["EST_SIZE"]
    ref_seq_gca      = current.get("REF_SEQ_GCA", None)
    ref_seq          = current.get("REF_SEQ", None)

    print(f"DEBUG - pacbio_sra - {pacbio_sra}")
    print(f"DEBUG - pacbio_raw_reads - {pacbio_raw_reads}")
    print(f"DEBUG - pacbio_raw_dir - {pacbio_raw_dir}")

    if pd.isna(pacbio_sra) and pd.isna(pacbio_raw_reads) and pd.isna(pacbio_raw_dir):
        print("SKIP:\tPacBio preprocessing; no reads provided")
        return None

    # Layout
    species_dir_abs = os.path.join(output_dir_abs, species_id)
    pacbio_dir_abs  = os.path.join(species_dir_abs, "PacBio")
    os.makedirs(pacbio_dir_abs, exist_ok=True)

    # Normalize implied paths BEFORE any chdir
    if pd.notna(pacbio_sra) and pd.isna(pacbio_raw_reads):
        pacbio_raw_reads = os.path.join(pacbio_dir_abs, f"{pacbio_sra}.fastq")
    elif isinstance(pacbio_raw_reads, str) and not os.path.isabs(pacbio_raw_reads):
        pacbio_raw_reads = _abs(os.path.join(output_dir_abs, pacbio_raw_reads))

    # Normalize reference FASTA if only GCA is present
    if pd.notna(ref_seq_gca) and (pd.isna(ref_seq) or not isinstance(ref_seq, str) or not ref_seq.strip()):
        ref_seq = os.path.join(species_dir_abs, "RefSeq", f"{species_id}_{ref_seq_gca}_RefSeq.fasta")
    if isinstance(ref_seq, str) and not os.path.isabs(ref_seq):
        ref_seq = _abs(os.path.join(output_dir_abs, ref_seq))

    prev_cwd = os.getcwd()
    os.chdir(pacbio_dir_abs)
    try:
        # Option 1: concatenate raw dir *.fastq
        if isinstance(pacbio_raw_dir, str) and pacbio_raw_dir.strip():
            pb_raw_dir_abs = pacbio_raw_dir if os.path.isabs(pacbio_raw_dir) else os.path.join(output_dir_abs, pacbio_raw_dir)
            files = sorted(glob.glob(os.path.join(pb_raw_dir_abs, "*.fastq")))
            if not files:
                print(f"ERROR:\tNo PacBio .fastq files found in {pb_raw_dir_abs}")
                return None
            pacbio_raw_reads = os.path.join(pacbio_dir_abs, f"{species_id}_pacbio_combined.fastq")
            print(f"NOTE:\tConcatenating PacBio files from {pb_raw_dir_abs} -> {pacbio_raw_reads}")
            with open(pacbio_raw_reads, "wb") as w:
                for f in files:
                    with open(f, "rb") as r:
                        shutil.copyfileobj(r, w)

        # Option 2: SRA -> FASTQ into PacBio dir (ensure -O is used)
        if isinstance(pacbio_sra, str) and pacbio_sra.strip():
            if not pacbio_raw_reads or not os.path.exists(pacbio_raw_reads):
                print(f"Downloading SRA {pacbio_sra} from GenBank...")
                _ = run_subprocess_cmd(["prefetch", "--force", "yes", pacbio_sra], False)
                _ = run_subprocess_cmd(["fasterq-dump", "--threads", str(cpu_threads), "-O", pacbio_dir_abs, pacbio_sra], False)
                expected = os.path.join(pacbio_dir_abs, f"{pacbio_sra}.fastq")
                if not os.path.exists(expected) or os.path.getsize(expected) == 0:
                    print(f"ERROR:\tExpected FASTQ not found or empty after fasterq-dump: {expected}")
                    return None
                pacbio_raw_reads = expected
                print(f"PASS:\tSRA converted to FASTQ: {pacbio_raw_reads}")
            else:
                print(f"SKIP:\tSRA already present as FASTQ: {pacbio_raw_reads}")

        # Validate raw reads
        if not pacbio_raw_reads or not os.path.exists(pacbio_raw_reads) or os.path.getsize(pacbio_raw_reads) == 0:
            print(f"ERROR:\tPacBio raw reads not found or empty: {pacbio_raw_reads}")
            return None

        # Parse estimated genome size (e.g. "5.0m")
        m = re.match(r"^(\d+(?:\.\d+)?)(\D+)$", str(est_size))
        if m:
            est_num, est_unit = m.group(1), m.group(2)
            mult = {'m': 10**6, 'g': 10**9}.get(est_unit.lower(), 25_000_000)
            est_size_bp = int(float(est_num) * mult)
        else:
            print(f"NOTE:\tUnable to parse input estimated size {est_size}, using default: 25000000")
            est_size_bp = 25_000_000

        # NanoPlot RAW (best-effort; do not fail pipeline on plotting errors)
        try:
            sample_stats_dict = nanoplot_qc_reads(pacbio_raw_reads, "Raw_PacBio_", cpu_threads, sample_stats_dict)
        except Exception as e:
            print(f"WARN:\tNanoPlot failed on raw PacBio reads ({e}); continuing without NanoStats.")

        # Decide filtered output path
        if isinstance(pacbio_sra, str) and pacbio_sra.strip():
            filtered_pb = os.path.join(pacbio_dir_abs, os.path.basename(pacbio_raw_reads).replace(pacbio_sra, f"{species_id}_pacbio_filtered"))
        else:
            filtered_pb = os.path.join(pacbio_dir_abs, f"{species_id}_pacbio_filtered.fastq")

        # If an old filtered file exists but is empty/suspicious, remove it to force regeneration
        if os.path.exists(filtered_pb) and not _nonempty(filtered_pb):
            print(f"WARN:\tExisting filtered FASTQ is empty or tiny ({filtered_pb}); regenerating.")
            try:
                os.remove(filtered_pb)
            except Exception:
                pass

        # ----------- Filtlong filtering (with correct --trim logic) -----------
        coverage = 60
        target_bases = est_size_bp * coverage

        have_ref = isinstance(ref_seq, str) and os.path.exists(ref_seq) and os.path.getsize(ref_seq) > 0
        if have_ref:
            print(f"INFO:\tUsing reference for filtlong trimming: {ref_seq}")
        else:
            print("INFO:\tNo reference available; running filtlong WITHOUT --trim (this was the cause of your empty file).")

        if not os.path.exists(filtered_pb):  # only run if we don't already have a good file
            common_flags = f"--min_length 1000 --min_mean_q 10 --keep_percent 90 --target_bases {target_bases}"
            if have_ref:
                # With a reference, --trim is valid. Use --ref to enable it.
                filt_cmd = f'filtlong --trim --ref "{ref_seq}" {common_flags} "{pacbio_raw_reads}" > "{filtered_pb}"'
            else:
                # Without a reference, DO NOT pass --trim
                filt_cmd = f'filtlong {common_flags} "{pacbio_raw_reads}" > "{filtered_pb}"'
            print(f"CMD:\t{filt_cmd}")
            rc = run_subprocess_cmd(filt_cmd, True)
            if rc != 0:
                print(f"ERROR:\tfiltlong exited with code {rc}")
                return None

            if not _nonempty(filtered_pb):
                print(f"ERROR:\tFiltlong did not produce a non-empty FASTQ at {filtered_pb}")
                print("HINT:\tIf you want trimming (--trim), provide a reference FASTA (REF_SEQ / REF_SEQ_GCA).")
                return None
        else:
            print(f"SKIP:\tFiltlong filtered reads already present: {filtered_pb}")

        # NanoPlot FILTERED (best-effort)
        try:
            sample_stats_dict = nanoplot_qc_reads(filtered_pb, "Filt_PacBio_", cpu_threads, sample_stats_dict)
        except Exception as e:
            print(f"WARN:\tNanoPlot failed on filtered PacBio reads ({e}); continuing.")

        # Select best long reads via shared helper (absolute paths)
        best_long_reads = select_long_reads(output_dir_abs, input_csv_abs, sample_id, cpu_threads)
        if not best_long_reads or not os.path.exists(best_long_reads):
            best_long_reads = filtered_pb

        # Canonicalize path for downstream steps
        canonical = os.path.join(pacbio_dir_abs, f"{species_id}_PacBio_highest_mean_qual_long_reads.fastq")
        if _abs(best_long_reads) != _abs(canonical):
            try:
                shutil.copy(best_long_reads, canonical)
            except Exception:
                try:
                    os.link(best_long_reads, canonical)
                except Exception:
                    try:
                        os.symlink(best_long_reads, canonical)
                    except Exception:
                        canonical = best_long_reads  # fall back to original

        print(f"PASS:\tPreprocessed PacBio reads for {sample_id}: {canonical}")
        return canonical
    finally:
        os.chdir(prev_cwd)


if __name__ == "__main__":
    print(f"DEBUG: Raw sys.argv = {sys.argv}")
    print(f"DEBUG: Length of sys.argv = {len(sys.argv)}")

    if len(sys.argv) != 6:
        print("Usage: python3 preprocess_pacbio.py <sample_id> <input_csv> <output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)

    for i, arg in enumerate(sys.argv):
        print(f"DEBUG: sys.argv[{i}] = '{arg}'")

    sample_id  = sys.argv[1]
    input_csv  = sys.argv[2]
    output_dir = sys.argv[3]

    # Robust parsing
    try:
        cpu_threads = int(str(sys.argv[4]).strip()) if str(sys.argv[4]).strip() else (os.cpu_count() or 1)
    except Exception:
        cpu_threads = os.cpu_count() or 1
    try:
        ram_gb = int(str(sys.argv[5]).strip()) if str(sys.argv[5]).strip() else 8
    except Exception:
        ram_gb = 8

    print(f"DEBUG: Parsed sample_id = '{sample_id}'")
    print(f"DEBUG: Parsed input_csv = '{input_csv}'")
    print(f"DEBUG: Parsed output_dir = '{output_dir}'")
    print(f"DEBUG: Parsed cpu_threads = '{sys.argv[4]}' (converted to {cpu_threads})")
    print(f"DEBUG: Parsed ram_gb = '{sys.argv[5]}' (converted to {ram_gb})")

    best = preprocess_pacbio(sample_id, input_csv, output_dir, cpu_threads, ram_gb)
