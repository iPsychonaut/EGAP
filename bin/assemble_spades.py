#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
assemble_spades.py

This script runs SPAdes assembly with Illumina and optional long reads.

Created on Wed Aug 16 2023

Updated on Wed Sept 3 2025

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import os, sys, shutil
import pandas as pd
from utilities import run_subprocess_cmd, get_current_row_data
from qc_assessment import qc_assessment


# --------------------------------------------------------------
# Run SPAdes assembly with Illumina and optional long reads
# --------------------------------------------------------------
def assemble_spades(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """Assemble genomic data using SPAdes with Illumina and optional ONT reads.

    Executes SPAdes assembly, processes input reads from a CSV, and performs
    quality control on the resulting assembly.

    Args are identical to your original docstring.
    """
    # Keep absolute paths so later code never loses them
    input_csv_abs  = os.path.abspath(input_csv)
    output_dir_abs = os.path.abspath(output_dir)

    input_df = pd.read_csv(input_csv_abs)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]

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

    species_dir = os.path.join(output_dir_abs, species_id)

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
        print(f"SKIP:\tSPAdes assembly already present: {egap_spades_assembly_path}")
        egap_spades_assembly_path, spades_stats_list, _ = qc_assessment(
            "spades", input_csv_abs, sample_id, output_dir_abs, cpu_threads, ram_gb
        )
        return egap_spades_assembly_path

    # Set long-read paths (ONT only for SPAdes), prefer prefiltered, fallback to raw
    highest_mean_qual_long_reads = None
    if pd.notna(ont_raw_reads):
        print("DEBUG - ONT RAW READS EXIST!")
        candidate = os.path.join(species_dir, "ONT", f"{species_id}_ONT_highest_mean_qual_long_reads.fastq")
        highest_mean_qual_long_reads = candidate if os.path.exists(candidate) else ont_raw_reads
    elif pd.notna(pacbio_raw_reads):
        print("SKIP:\tSPAdes cannot be used to assemble PacBio reads...")
        return None

    print(f"DEBUG - highest_mean_qual_long_reads    - {highest_mean_qual_long_reads}")

    # If no usable reads at all, bail
    if pd.isna(ont_raw_reads) and pd.isna(illumina_f_raw_reads) and pd.isna(illumina_r_raw_reads) and pd.isna(pacbio_raw_reads):
        print("SKIP:\tNo reads available for processing")
        return None

    # -------- Resolve EVERYTHING to absolute paths (no chdir, no double nesting) --------
    start_dir = os.getcwd()

    def _abs_safe(p):
        if p is None or (isinstance(p, float) and pd.isna(p)):
            return p
        return p if os.path.isabs(p) else os.path.realpath(os.path.join(start_dir, p))

    # Make all potentially used files absolute
    illu_dedup_f_reads = _abs_safe(illu_dedup_f_reads)
    illu_dedup_r_reads = _abs_safe(illu_dedup_r_reads)
    ref_seq = _abs_safe(ref_seq)
    highest_mean_qual_long_reads = _abs_safe(highest_mean_qual_long_reads)
    illumina_f_raw_reads = _abs_safe(illumina_f_raw_reads)
    illumina_r_raw_reads = _abs_safe(illumina_r_raw_reads)

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
        print("ERROR:\tSPAdes finished but scaffolds.fasta not found; check spades.log and params.txt.")
        return None

    shutil.move(spades_path, egap_spades_assembly_path)

    # QC using absolute paths
    egap_spades_assembly_path, spades_stats_list, _ = qc_assessment(
        "spades", input_csv_abs, sample_id, output_dir_abs, cpu_threads, ram_gb
    )
    return egap_spades_assembly_path


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python3 assemble_spades.py <sample_id> <input_csv> "
              "<output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)

    egap_spades_assembly_path = assemble_spades(
        sys.argv[1],       # sample_id
        sys.argv[2],       # input_csv
        sys.argv[3],       # output_dir
        str(sys.argv[4]),  # cpu_threads
        str(sys.argv[5])   # ram_gb
    )
