#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
decontaminate_assembly.py

Remove contaminating sequences (bacteria, archaea, prokarya) from the
selected best assembly using Tiara taxonomic classification. Runs between
compare_assemblies and polish_assembly in the EGAP pipeline.

Classification logic:
  KEEP:   eukarya, organelle, unknown
           - eukarya  : target genome sequences
           - organelle: mitochondria / plastid (useful downstream)
           - unknown  : ambiguous sequences kept conservatively to avoid
                        losing real eukaryotic sequence
  REMOVE: bacteria, archaea, prokarya
           - high-confidence contamination

Sequences below Tiara's min_len threshold are not classified and are kept
by default (too short to reliably classify as contamination).

min_len is automatically adjusted per sequencing strategy:
  - Illumina-only : 1100 bp  (shorter contigs from short-read assembly)
  - Hybrid / ONT / PacBio: 3000 bp (default, longer contigs)

Outputs:
  {sample_dir}/tiara_decontamination/
      {species_id}_tiara_classification.tsv  -- raw Tiara output
      tiara_normalized.tsv                   -- normalized column names
      summary_*.tsv                          -- per-class statistics
      tiara_report.html                      -- interactive Plotly report
  {sample_dir}/{species_id}_best_assembly_pretiara.fasta  -- backup original
  {sample_dir}/{species_id}_best_assembly.fasta           -- decontaminated

Created on Mon Mar 16 2026

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""

import os
import sys
import shutil
import pandas as pd
from pathlib import Path
from Bio import SeqIO

from utilities import get_current_row_data
from tiara_classify import run_tiara, normalize_tiara_columns


# Classes that indicate contamination — remove these
REMOVE_FIRST_STAGE = {"bacteria", "archaea", "prokarya"}

# Classes to keep (eukarya = target genome, organelle = mtDNA/plastid,
# unknown = conservative — may be real sequence too short/unusual to classify)
KEEP_FIRST_STAGE = {"eukarya", "organelle", "unknown"}


def decontaminate_assembly(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """
    Run Tiara on the best assembly and remove contaminating sequences.

    Reads the best assembly from {sample_dir}/{species_id}_best_assembly.fasta,
    classifies all sequences above min_len with Tiara, removes bacterial/archaeal/
    prokaryal contigs, backs up the original, and overwrites the best assembly with
    the decontaminated version.

    Args:
        sample_id (str): Sample identifier.
        input_csv (str): Path to metadata CSV file.
        output_dir (str): Directory for pipeline output files.
        cpu_threads (int or str): Number of CPU threads to use.
        ram_gb (int or str): Available RAM in GB (passed for consistency, not used).

    Returns:
        str: Path to the (decontaminated) best assembly FASTA, or None on error.
    """
    cpu_threads = int(cpu_threads)

    # ------------------------------------------------------------------
    # 1. Load CSV and extract sample metadata
    # ------------------------------------------------------------------
    input_df = pd.read_csv(input_csv)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]

    species_id        = current_series["SPECIES_ID"]
    ont_sra           = current_series["ONT_SRA"]
    ont_raw_reads     = current_series["ONT_RAW_READS"]
    illumina_sra      = current_series["ILLUMINA_SRA"]
    illumina_f_reads  = current_series["ILLUMINA_RAW_F_READS"]
    illumina_r_reads  = current_series["ILLUMINA_RAW_R_READS"]
    pacbio_sra        = current_series["PACBIO_SRA"]
    pacbio_raw_reads  = current_series["PACBIO_RAW_READS"]

    # ------------------------------------------------------------------
    # 2. Build paths
    # ------------------------------------------------------------------
    output_dir_abs = str(Path(output_dir).resolve())
    species_dir    = os.path.join(output_dir_abs, species_id)
    sample_dir     = os.path.join(species_dir, sample_id)

    best_assembly  = os.path.join(sample_dir, f"{species_id}_best_assembly.fasta")
    tiara_dir      = os.path.join(sample_dir, "tiara_decontamination")
    output_tsv     = os.path.join(tiara_dir, f"{species_id}_tiara_classification.tsv")
    done_marker    = os.path.join(tiara_dir, "decontamination.done")
    backup_fasta   = os.path.join(sample_dir, f"{species_id}_best_assembly_pretiara.fasta")

    print(f"DEBUG - species_id    : {species_id}")
    print(f"DEBUG - sample_dir    : {sample_dir}")
    print(f"DEBUG - best_assembly : {best_assembly}")
    print(f"DEBUG - tiara_dir     : {tiara_dir}")

    # ------------------------------------------------------------------
    # 3. Guard: best assembly must exist
    # ------------------------------------------------------------------
    if not os.path.exists(best_assembly):
        print(f"ERROR:\tBest assembly not found: {best_assembly}")
        print("HINT:\tRun compare_assemblies before decontaminate_assembly.")
        return None

    # ------------------------------------------------------------------
    # 4. Skip if already done
    # ------------------------------------------------------------------
    if os.path.exists(done_marker):
        print(f"SKIP:\tTiara decontamination already completed: {done_marker}")
        return best_assembly

    os.makedirs(tiara_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # 5. Detect sequencing strategy → choose appropriate min_len
    # ------------------------------------------------------------------
    def _has_data(sra_col, reads_col):
        """Return True if either the SRA accession or a file path is non-null."""
        sra_ok   = pd.notna(sra_col)   and str(sra_col)   not in ("None", "nan", "")
        reads_ok = pd.notna(reads_col) and str(reads_col) not in ("None", "nan", "")
        return sra_ok or reads_ok

    has_ont    = _has_data(ont_sra, ont_raw_reads)
    has_pacbio = _has_data(pacbio_sra, pacbio_raw_reads)
    is_illumina_only = not has_ont and not has_pacbio

    min_length_bp = 1100 if is_illumina_only else 3000
    strategy_label = "Illumina-only" if is_illumina_only else "Hybrid/Long-read"

    print(f"DEBUG - Sequencing strategy : {strategy_label}")
    print(f"DEBUG - Tiara min_len       : {min_length_bp} bp")

    # ------------------------------------------------------------------
    # 6. Run Tiara classification + visualization
    # ------------------------------------------------------------------
    if os.path.exists(output_tsv):
        print(f"SKIP:\tTiara TSV already exists: {output_tsv}")
    else:
        print(f"\n[Tiara] Classifying {best_assembly} ...")
        try:
            run_tiara(
                input_assembly=best_assembly,
                output_tsv=output_tsv,
                min_length_bp=min_length_bp,
                probability_cutoff="0.65,0.65",
                to_fasta="all",
                cpu_threads=cpu_threads,
                first_stage_kmer=6,
                second_stage_kmer=7,
                run_visualization=True,
                visualization_output_dir=tiara_dir,
                visualization_output_html=os.path.join(tiara_dir, "tiara_report.html"),
            )
        except FileNotFoundError as exc:
            print(f"WARN:\tTiara executable not found: {exc}")
            print("WARN:\tSkipping decontamination — ensure tiara is installed in PATH.")
            return best_assembly
        except RuntimeError as exc:
            print(f"WARN:\tTiara subprocess failed: {exc}")
            print("WARN:\tSkipping decontamination.")
            return best_assembly

    if not os.path.exists(output_tsv):
        print(f"WARN:\tTiara TSV not generated: {output_tsv}")
        print("WARN:\tSkipping decontamination.")
        return best_assembly

    # ------------------------------------------------------------------
    # 7. Parse Tiara TSV → identify contaminating sequence IDs
    # ------------------------------------------------------------------
    try:
        tiara_df = pd.read_csv(output_tsv, sep="\t")
        tiara_df = normalize_tiara_columns(tiara_df)
    except Exception as exc:
        print(f"WARN:\tFailed to parse Tiara TSV ({exc}). Skipping decontamination.")
        return best_assembly

    if tiara_df.empty:
        print("WARN:\tTiara TSV is empty — no sequences were classified.")
        print("WARN:\tThis can happen if all contigs are below min_len. Skipping decontamination.")
        _write_done_marker(done_marker, removed_count=0, extra="TSV empty — all below min_len")
        return best_assembly

    classified_ids  = set(tiara_df["sequence_id"].astype(str).tolist())
    remove_ids      = set(
        tiara_df.loc[
            tiara_df["first_stage_class"].isin(REMOVE_FIRST_STAGE),
            "sequence_id"
        ].astype(str).tolist()
    )

    # Summary of what Tiara found
    class_counts = tiara_df["first_stage_class"].value_counts().to_dict()
    print("\n[Tiara] Classification summary:")
    for cls, cnt in sorted(class_counts.items()):
        tag = "REMOVE" if cls in REMOVE_FIRST_STAGE else "KEEP"
        print(f"  {cls:<12} : {cnt:>6} sequences  [{tag}]")
    print(f"  {'(unclassified)':<12} : short contigs below {min_length_bp} bp  [KEEP]")

    # ------------------------------------------------------------------
    # 8. Filter the assembly FASTA
    # ------------------------------------------------------------------
    all_records = list(SeqIO.parse(best_assembly, "fasta"))
    total_count = len(all_records)
    total_bp    = sum(len(r.seq) for r in all_records)

    filtered_records = [
        r for r in all_records
        if str(r.id) not in remove_ids
    ]

    removed_records = [r for r in all_records if str(r.id) in remove_ids]
    removed_count   = len(removed_records)
    removed_bp      = sum(len(r.seq) for r in removed_records)

    print(f"\n[Tiara] Decontamination results:")
    print(f"  Total sequences  : {total_count:>8,}")
    print(f"  Removed sequences: {removed_count:>8,}  ({removed_count / total_count * 100:.1f}%)")
    print(f"  Removed bases    : {removed_bp:>8,} bp  ({removed_bp / total_bp * 100:.1f}%)")
    print(f"  Kept sequences   : {len(filtered_records):>8,}")
    print(f"  Kept bases       : {total_bp - removed_bp:>8,} bp")

    # ------------------------------------------------------------------
    # 9. If nothing to remove, mark done and return
    # ------------------------------------------------------------------
    if removed_count == 0:
        print("\nPASS:\tNo contaminating sequences detected by Tiara.")
        _write_done_marker(done_marker, removed_count=0)
        return best_assembly

    # ------------------------------------------------------------------
    # 10. Safety check — warn if decontamination is very aggressive
    # ------------------------------------------------------------------
    pct_seqs_removed = removed_count / total_count * 100
    pct_bp_removed   = removed_bp   / total_bp    * 100

    if pct_seqs_removed > 50 or pct_bp_removed > 50:
        print(
            f"\nWARN:\tTiara removed >50% of sequences ({pct_seqs_removed:.1f}% seqs, "
            f"{pct_bp_removed:.1f}% bp). Proceeding, but review the tiara_report.html "
            f"to confirm this is genuine contamination."
        )

    # ------------------------------------------------------------------
    # 11. Backup original → overwrite best_assembly with filtered FASTA
    # ------------------------------------------------------------------
    if not os.path.exists(backup_fasta):
        shutil.copy(best_assembly, backup_fasta)
        print(f"\nDEBUG - Original assembly backed up: {backup_fasta}")
    else:
        print(f"SKIP:\tBackup already exists: {backup_fasta}")

    SeqIO.write(filtered_records, best_assembly, "fasta")
    print(f"PASS:\tDecontaminated assembly written: {best_assembly}")

    # ------------------------------------------------------------------
    # 12. Write done marker
    # ------------------------------------------------------------------
    _write_done_marker(done_marker, removed_count=removed_count,
                       removed_ids=[r.id for r in removed_records])

    print(f"PASS:\tTiara decontamination complete for {sample_id}")
    return best_assembly


def _write_done_marker(marker_path, removed_count, removed_ids=None, extra=None):
    """Write a decontamination.done marker file with a brief log."""
    with open(marker_path, "w") as fh:
        fh.write(f"Tiara decontamination complete.\n")
        fh.write(f"Sequences removed: {removed_count}\n")
        if extra:
            fh.write(f"Note: {extra}\n")
        if removed_ids:
            fh.write("Removed sequence IDs:\n")
            for rid in removed_ids:
                fh.write(f"  {rid}\n")


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print(
            "Usage: python3 decontaminate_assembly.py <sample_id> <input_csv> "
            "<output_dir> <cpu_threads> <ram_gb>",
            file=sys.stderr,
        )
        sys.exit(1)

    result = decontaminate_assembly(
        sys.argv[1],        # sample_id
        sys.argv[2],        # input_csv
        sys.argv[3],        # output_dir
        sys.argv[4],        # cpu_threads
        sys.argv[5],        # ram_gb
    )

    if result is None:
        sys.exit(1)
