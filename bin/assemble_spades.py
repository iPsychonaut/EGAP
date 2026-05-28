#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
assemble_spades.py

This script runs SPAdes assembly with Illumina and optional long reads.

Stage:
    Short-read Assembly (SPAdes)

Created on Wed Aug 16 2023

Updated on 2026-04-16

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os
import sys
import shutil
from typing import Optional

import pandas as pd
from utilities import run_subprocess_cmd, log_print, initialize_logging_environment, load_sample_context
from qc_assessment import qc_assessment
from file_manager import remove_file, remove_dir


# --------------------------------------------------------------
# Run SPAdes assembly with Illumina and optional long reads
# --------------------------------------------------------------
def assemble_spades(
    sample_id: str,
    input_csv: str,
    output_dir: str,
    cpu_threads: int,
    ram_gb: int,
) -> Optional[str]:
    """Assemble genomic data using SPAdes with Illumina and optional ONT reads.

    Reads metadata from *input_csv*, resolves read paths, runs SPAdes in
    ``--isolate`` mode (with optional ``--nanopore`` supplemental reads),
    and performs quality control via ``qc_assessment``.  Skips the assembly
    step if the final output file already exists.

    Parameters
    ----------
    sample_id : str
        Sample identifier used to look up the row in *input_csv*.
    input_csv : str
        Path to the metadata CSV file.
    output_dir : str
        Root output directory; per-species subdirectories are created here.
    cpu_threads : int or str
        Number of CPU threads to pass to SPAdes.
    ram_gb : int or str
        Maximum RAM in GB to pass to SPAdes (``-m`` flag).

    Returns
    -------
    str or None
        Absolute path to the final SPAdes assembly FASTA, or ``None``
        if the assembly could not be completed or no suitable reads are
        available.

    Raises
    ------
    SystemExit
        Propagated from ``run_subprocess_cmd`` if the SPAdes binary is not
        found on ``PATH``.
    """
    ctx = load_sample_context(sample_id, input_csv, output_dir, cpu_threads, ram_gb)
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

    species_dir = os.path.join(ctx.output_dir, species_id)

    if pd.notna(ont_sra) and pd.isna(ont_raw_reads):
        ont_raw_reads = os.path.join(species_dir, "ONT", f"{ont_sra}.fastq")
    if pd.notna(illumina_sra) and pd.isna(illumina_f_raw_reads) and pd.isna(illumina_r_raw_reads):
        illumina_f_raw_reads = os.path.join(species_dir, "Illumina", f"{illumina_sra}_1.fastq")
        illumina_r_raw_reads = os.path.join(species_dir, "Illumina", f"{illumina_sra}_2.fastq")
    if pd.notna(pacbio_sra) and pd.isna(pacbio_raw_reads):
        pacbio_raw_reads = os.path.join(species_dir, "PacBio", f"{pacbio_sra}.fastq")
    if pd.notna(ref_seq_gca) and pd.isna(ref_seq):
        ref_seq = os.path.join(species_dir, "RefSeq", f"{species_id}_{ref_seq_gca}_RefSeq.fasta")

    print(f"DEBUG - illumina_sra - {illumina_sra}")
    print(f"DEBUG - illumina_f_raw_reads - {illumina_f_raw_reads}")
    print(f"DEBUG - illumina_r_raw_reads - {illumina_r_raw_reads}")
    print(f"DEBUG - ont_sra - {ont_sra}")
    print(f"DEBUG - ont_raw_reads - {ont_raw_reads}")
    print(f"DEBUG - pacbio_sra - {pacbio_sra}")
    print(f"DEBUG - pacbio_raw_reads - {pacbio_raw_reads}")
    print(f"DEBUG - ref_seq_gca - {ref_seq_gca}")
    print(f"DEBUG - ref_seq - {ref_seq}")
    print(f"DEBUG - species_id - {species_id}")
    print(f"DEBUG - est_size - {est_size}")

    # Set Illumina deduplicated read paths only if Illumina reads are present
    illu_dedup_f_reads = None
    illu_dedup_r_reads = None
    if pd.notna(illumina_f_raw_reads) and pd.notna(illumina_r_raw_reads):
        illu_dedup_f_reads = os.path.join(species_dir, "Illumina", f"{species_id}_illu_forward_dedup.fastq")
        illu_dedup_r_reads = os.path.join(species_dir, "Illumina", f"{species_id}_illu_reverse_dedup.fastq")

    print(f"DEBUG - illu_dedup_f_reads - {illu_dedup_f_reads}")
    print(f"DEBUG - illu_dedup_r_reads - {illu_dedup_r_reads}")

    # Prepare output locations early so we can fast-skip (and still QC) even if reads are missing
    sample_dir = os.path.join(species_dir, sample_id)
    spades_out_dir = os.path.join(sample_dir, "spades_assembly")
    os.makedirs(spades_out_dir, exist_ok=True)
    egap_spades_assembly_path = os.path.join(spades_out_dir, f"{sample_id}_spades.fasta")

    # ---------- FAST SKIP if final output exists (re-run QC) ----------
    if os.path.exists(egap_spades_assembly_path) and os.path.getsize(egap_spades_assembly_path) > 0:
        log_print(f"SKIP:\tSPAdes assembly already present: {egap_spades_assembly_path}")
        egap_spades_assembly_path, spades_stats_list, _ = qc_assessment(
            "spades", ctx.input_csv, sample_id, ctx.output_dir, cpu_threads, ram_gb
        )
        return egap_spades_assembly_path

    # Set long-read paths (ONT only for SPAdes), prefer prefiltered, fallback to raw
    highest_mean_qual_long_reads = None
    if pd.notna(ont_raw_reads):
        log_print("DEBUG - ONT RAW READS EXIST!")
        candidate = os.path.join(species_dir, "ONT", f"{species_id}_ONT_highest_mean_qual_long_reads.fastq")
        highest_mean_qual_long_reads = candidate if os.path.exists(candidate) else ont_raw_reads
    elif pd.notna(pacbio_raw_reads):
        log_print("SKIP:\tSPAdes cannot be used to assemble PacBio reads...")
        return None

    print(f"DEBUG - highest_mean_qual_long_reads    - {highest_mean_qual_long_reads}")

    # If no usable reads at all, bail
    if pd.isna(ont_raw_reads) and pd.isna(illumina_f_raw_reads) and pd.isna(illumina_r_raw_reads) and pd.isna(pacbio_raw_reads):
        log_print("SKIP:\tNo reads available for processing")
        return None

    # -------- Resolve EVERYTHING to absolute paths (no chdir, no double nesting) --------
    start_dir = os.getcwd()

    def abs_safe(p):
        if p is None or pd.isna(p):
            return None
        p = str(p)
        return p if os.path.isabs(p) else os.path.realpath(os.path.join(start_dir, p))

    # Make all potentially used files absolute
    illu_dedup_f_reads = abs_safe(illu_dedup_f_reads)
    illu_dedup_r_reads = abs_safe(illu_dedup_r_reads)
    ref_seq = abs_safe(ref_seq)
    highest_mean_qual_long_reads = abs_safe(highest_mean_qual_long_reads)
    illumina_f_raw_reads = abs_safe(illumina_f_raw_reads)
    illumina_r_raw_reads = abs_safe(illumina_r_raw_reads)

    # Only use dedup FASTQs if they truly exist and are non-empty; otherwise fall back to raw
    use_dedup = all([
        illu_dedup_f_reads, illu_dedup_r_reads,
        os.path.exists(illu_dedup_f_reads) if illu_dedup_f_reads else False,
        os.path.exists(illu_dedup_r_reads) if illu_dedup_r_reads else False,
        os.path.getsize(illu_dedup_f_reads) > 0 if illu_dedup_f_reads else False,
        os.path.getsize(illu_dedup_r_reads) > 0 if illu_dedup_r_reads else False,
    ])
    if not use_dedup:
        # Fall back to raw (already made absolute)
        illu_dedup_f_reads = illumina_f_raw_reads
        illu_dedup_r_reads = illumina_r_raw_reads

    print(f"DEBUG - spades_out_dir - {spades_out_dir}")

    # Establish command (no chdir; output is spades_out_dir)
    kmer_list = ["21", "33", "55", "77", "99"]
    spades_path = os.path.join(spades_out_dir, "scaffolds.fasta")

    spades_cmd = [
        "spades.py",
        "--isolate",
        "-t", str(cpu_threads),
        "-m", str(ram_gb),
        "--cov-cutoff", "auto",
    ]
    if illu_dedup_f_reads and illu_dedup_r_reads:
        spades_cmd += ["-1", illu_dedup_f_reads, "-2", illu_dedup_r_reads]
    if highest_mean_qual_long_reads and pd.notna(ont_raw_reads):
        spades_cmd += ["--nanopore", highest_mean_qual_long_reads]
    if pd.notna(ref_seq):
        spades_cmd += ["--trusted-contigs", ref_seq]
    spades_cmd += ["-o", spades_out_dir, "-k", ",".join(kmer_list)]

    print(f"DEBUG - spades_path - {spades_path}")
    print(f"DEBUG - spades_cmd - {spades_cmd}")

    _ = run_subprocess_cmd(spades_cmd, shell_check=False)

    if not os.path.exists(spades_path):
        log_print("ERROR:\tSPAdes finished but scaffolds.fasta not found; check spades.log and params.txt.")
        return None

    shutil.move(spades_path, egap_spades_assembly_path)

    # QC using absolute paths
    egap_spades_assembly_path, spades_stats_list, _ = qc_assessment(
        "spades", ctx.input_csv, sample_id, ctx.output_dir, cpu_threads, ram_gb
    )

    # --- Cleanup SPAdes intermediates once final assembly is confirmed ---
    if egap_spades_assembly_path and os.path.exists(str(egap_spades_assembly_path)):
        # Per-kmer working directories (K21 … K99)
        for _kmer in ["K21", "K33", "K55", "K77", "K99"]:
            remove_dir(os.path.join(spades_out_dir, _kmer))
        # Intermediate FASTA outputs superseded by the final scaffold
        remove_file(os.path.join(spades_out_dir, "contigs.fasta"))
        remove_file(os.path.join(spades_out_dir, "before_rr.fasta"))
        # Scratch directories
        remove_dir(os.path.join(spades_out_dir, "misc"))
        remove_dir(os.path.join(spades_out_dir, "tmp"))

    return egap_spades_assembly_path


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python3 assemble_spades.py <sample_id> <input_csv> "
              "<output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)

    initialize_logging_environment(sys.argv[3], sys.argv[1])

    egap_spades_assembly_path = assemble_spades(
        sys.argv[1],       # sample_id
        sys.argv[2],       # input_csv
        sys.argv[3],       # output_dir
        str(sys.argv[4]),  # cpu_threads
        str(sys.argv[5])   # ram_gb
    )
