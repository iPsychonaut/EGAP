#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
assemble_hifiasm.py

Run hifiasm on PacBio HiFi reads (FASTQ/FASTA, gz ok). CWD-safe and resilient:
- Uses absolute paths for CSV/output.
- Prefers preprocessed highest-mean-quality PacBio reads when present.
- Stable output prefix (<sample_id>) so downstream paths are predictable.
- Tolerant to hifiasm output filename variants.
- GFA->FASTA via gfatools (fallback to AWK if needed).

Created on Wed Aug 16 2023
Updated on Wed Sept 3 2025

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os, sys, shutil
import pandas as pd
from utilities import run_subprocess_cmd, get_current_row_data
from qc_assessment import qc_assessment


def _abs(p):  # helper for absolute paths
    return os.path.abspath(p) if isinstance(p, str) else p


def assemble_hifiasm(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """Assemble with hifiasm using PacBio reads and QC the result.

    Returns:
        str or None: Path to <sample_id>_hifiasm.fasta or None on failure/skip.
    """
    # --- Resolve critical paths to absolute early ---
    input_csv_abs = _abs(input_csv)
    output_dir_abs = _abs(output_dir)

    # Read metadata safely (independent of CWD changes)
    input_df = pd.read_csv(input_csv_abs)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]

    # CSV fields
    pacbio_sra = current_series["PACBIO_SRA"]
    pacbio_raw_reads = current_series["PACBIO_RAW_READS"]
    ref_seq_gca = current_series["REF_SEQ_GCA"]
    ref_seq = current_series["REF_SEQ"]
    species_id = current_series["SPECIES_ID"]

    species_dir = os.path.join(output_dir_abs, species_id)

    # Normalize implied PacBio path if only SRA is provided
    if (pd.notna(pacbio_sra) and pd.isna(pacbio_raw_reads)):
        pacbio_raw_reads = os.path.join(species_dir, "PacBio", f"{pacbio_sra}.fastq")
    # If CSV had a relative path, anchor it under the project
    if isinstance(pacbio_raw_reads, str) and not os.path.isabs(pacbio_raw_reads):
        pacbio_raw_reads = _abs(os.path.join(output_dir_abs, pacbio_raw_reads))

    print(f"DEBUG - pacbio_sra - {pacbio_sra}")
    print(f"DEBUG - pacbio_raw_reads - {pacbio_raw_reads}")
    print(f"DEBUG - ref_seq_gca - {ref_seq_gca}")
    print(f"DEBUG - ref_seq - {ref_seq}")
    print(f"DEBUG - species_id - {species_id}")

    # Require at least some PacBio path or SRA-derived path
    if pd.isna(pacbio_sra) and (not isinstance(pacbio_raw_reads, str) or pacbio_raw_reads.strip() == ""):
        print("SKIP:\tNo PacBio reads files provided.")
        return None

    # Prefer preprocessed highest-quality reads if present
    highest = os.path.join(species_dir, "PacBio", f"{species_id}_PacBio_highest_mean_qual_long_reads.fastq")
    if not os.path.exists(highest):
        highest = pacbio_raw_reads

    print(f"DEBUG - highest_mean_qual_long_reads    - {highest}")

    # Validate input reads
    if not highest or not os.path.exists(highest) or os.path.getsize(highest) == 0:
        print(f"ERROR:\tPacBio reads not found or empty: {highest}")
        print("HINT:\tRun preprocess_pacbio first, or check your CSV paths.")
        return None

    # Output dirs/files
    sample_dir = os.path.join(species_dir, sample_id)
    hifiasm_out_dir = os.path.join(sample_dir, "hifiasm_assembly")
    os.makedirs(hifiasm_out_dir, exist_ok=True)

    # Use a stable prefix (sample_id) for predictable outputs
    prefix_base = os.path.join(hifiasm_out_dir, sample_id)
    egap_hifiasm_assembly_path = os.path.join(hifiasm_out_dir, f"{sample_id}_hifiasm.fasta")

    # ---------- FAST SKIP if final output exists (re-run QC) ----------
    if os.path.exists(egap_hifiasm_assembly_path) and os.path.getsize(egap_hifiasm_assembly_path) > 0:
        print(f"SKIP:\tHiFi assembly already present: {egap_hifiasm_assembly_path}")
        egap_hifiasm_assembly_path, hifiasm_stats_list, _ = qc_assessment(
            "hifiasm", input_csv_abs, sample_id, output_dir_abs, cpu_threads, ram_gb
        )
        return egap_hifiasm_assembly_path

    # Common hifiasm GFA names across versions
    gfa_candidates = [
        f"{prefix_base}.bp.p_ctg.gfa",
        f"{prefix_base}.p_ctg.gfa",
        f"{prefix_base}.asm.bp.p_ctg.gfa",   # some older scripts used ".asm" in prefix
        f"{prefix_base}.asm.p_ctg.gfa",
    ]

    prev_cwd = os.getcwd()
    try:
        # Work in the output dir (safer for hifiasm temp files)
        os.chdir(hifiasm_out_dir)

        # If a GFA already exists, skip running hifiasm (we'll still do conversion + QC below)
        already_has_gfa = any(os.path.exists(g) for g in gfa_candidates)
        if already_has_gfa:
            print("SKIP:\tHiFi GFA already present; skipping hifiasm run.")
        else:
            # hifiasm command
            hifiasm_cmd = [
                "hifiasm",
                "-o", prefix_base,
                "-t", str(cpu_threads),
                highest
            ]
            print(f"CMD:\t{' '.join(hifiasm_cmd)}")
            rc = run_subprocess_cmd(hifiasm_cmd, shell_check=False)
            if rc != 0:
                print(f"WARN:\thifiasm exited with code {rc}")

        # Pick the first existing GFA
        gfa_path = next((g for g in gfa_candidates if os.path.exists(g)), None)
        if not gfa_path:
            print("ERROR:\tNo hifiasm primary contig GFA found after run.")
            for g in gfa_candidates: print(f"DEBUG: missing -> {g}")
            return None

        # Convert GFA -> FASTA (only if final FASTA not already present—handled by fast-skip earlier)
        if not os.path.exists(egap_hifiasm_assembly_path):
            if shutil.which("gfatools"):
                gfa_cmd = f'gfatools gfa2fa "{gfa_path}" > "{egap_hifiasm_assembly_path}"'
                print(f"CMD:\t{gfa_cmd}")
                rc = run_subprocess_cmd(gfa_cmd, shell_check=True)
                if rc != 0:
                    print("WARN:\tgfatools failed; attempting AWK fallback")
            if not os.path.exists(egap_hifiasm_assembly_path):
                # AWK fallback: write segments (S) as FASTA
                awk_cmd = (
                    f"awk 'BEGIN{{OFS=\"\\t\"}} /^S\\t/ "
                    f"{{print \">\"$2\"\\n\"$3}}' '{gfa_path}' > '{egap_hifiasm_assembly_path}'"
                )
                print(f"CMD:\t{awk_cmd}")
                rc = run_subprocess_cmd(awk_cmd, shell_check=True)
                if rc != 0 or (not os.path.exists(egap_hifiasm_assembly_path)):
                    print("ERROR:\tFailed to convert GFA to FASTA.")
                    return None

        # QC (absolute paths so CWD is irrelevant)
        egap_hifiasm_assembly_path, hifiasm_stats_list, _ = qc_assessment(
            "hifiasm", input_csv_abs, sample_id, output_dir_abs, cpu_threads, ram_gb
        )
        return egap_hifiasm_assembly_path
    finally:
        os.chdir(prev_cwd)


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python3 assemble_hifiasm.py <sample_id> <input_csv> <output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)

    egap_hifiasm_assembly_path = assemble_hifiasm(
        sys.argv[1],              # sample_id
        sys.argv[2],              # input_csv
        sys.argv[3],              # output_dir
        str(sys.argv[4]),         # cpu_threads
        str(sys.argv[5])          # ram_gb
    )
