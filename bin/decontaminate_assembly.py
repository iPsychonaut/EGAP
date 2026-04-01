#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
decontaminate_assembly.py

Removes contaminant sequences from a genome assembly using Tiara, a deep-learning
classifier that labels each contig as one of:

    eukarya | bacteria | archaea | prokarya | organelle | unknown

The keep/remove split is driven by the organism's kingdom, allowing the same
script to decontaminate both eukaryotic and prokaryotic assemblies correctly:

    Kingdom            Keep                           Remove
    -----------------  -----------------------------  ----------------------------
    bacteria           bacteria, prokarya, unknown    eukarya, archaea, organelle
    archaea            archaea, prokarya, unknown     eukarya, bacteria, organelle
    flora/funga/fauna  eukarya, organelle, unknown    bacteria, archaea, prokarya
    (unrecognised)     eukarya, organelle, unknown    bacteria, archaea, prokarya

Rationale for 'organelle':
  - Eukaryotes: kept because mitochondria and plastid sequences belong to the
    organism and would otherwise be lost.
  - Prokaryotes: removed because organelle-classified contigs are most likely
    eukaryotic host contamination, not genuine prokaryotic sequence.

'unknown' is always kept to avoid discarding genuine but low-complexity sequence.

Created on Tue Apr 01 2026

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""

import os
import sys
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from utilities import run_subprocess_cmd, get_current_row_data, initialize_logging_environment, log_print


# --------------------------------------------------------------
# All labels produced by Tiara
# --------------------------------------------------------------
ALL_TIARA_CLASSES = {"eukarya", "bacteria", "archaea", "prokarya", "organelle", "unknown"}


# --------------------------------------------------------------
# Resolve keep/remove sets from organism kingdom
# --------------------------------------------------------------
def get_decontamination_classes(kingdom_id, karyote_id):
    """Return (keep_set, remove_set) for Tiara labels based on organism kingdom.

    Args:
        kingdom_id (str or None): Value of the ORGANISM_KINGDOM CSV column
            (e.g. "funga", "bacteria", "archaea").  Case-insensitive.
        karyote_id (str or None): Value of the ORGANISM_KARYOTE CSV column.
            Currently informational only; kept for future use.

    Returns:
        tuple[set, set]: (keep_classes, remove_classes) where both sets are
            non-overlapping subsets of ALL_TIARA_CLASSES whose union equals
            ALL_TIARA_CLASSES.
    """
    euk_keep    = {"eukarya", "organelle", "unknown"}
    euk_remove  = {"bacteria", "archaea", "prokarya"}
    bac_keep    = {"bacteria", "prokarya", "unknown"}
    bac_remove  = {"eukarya", "archaea", "organelle"}
    arc_keep    = {"archaea", "prokarya", "unknown"}
    arc_remove  = {"eukarya", "bacteria", "organelle"}

    # Normalise kingdom_id
    if pd.isna(kingdom_id) if not isinstance(kingdom_id, str) else False:
        normalised = None
    elif not isinstance(kingdom_id, str) or not kingdom_id.strip():
        normalised = None
    else:
        normalised = kingdom_id.strip().lower()

    if normalised == "bacteria":
        keep_classes, remove_classes = bac_keep, bac_remove
    elif normalised == "archaea":
        keep_classes, remove_classes = arc_keep, arc_remove
    elif normalised in ("flora", "funga", "fauna"):
        keep_classes, remove_classes = euk_keep, euk_remove
    else:
        log_print(f"WARN:\tUnrecognised or missing kingdom_id '{kingdom_id}'. "
                  f"Falling back to eukaryote decontamination profile.")
        keep_classes, remove_classes = euk_keep, euk_remove

    assert keep_classes | remove_classes == ALL_TIARA_CLASSES, \
        "keep_classes and remove_classes must cover ALL_TIARA_CLASSES exactly."
    assert keep_classes & remove_classes == set(), \
        "keep_classes and remove_classes must not overlap."

    return keep_classes, remove_classes


# --------------------------------------------------------------
# Write a completion marker for skip-detection on re-runs
# --------------------------------------------------------------
def _write_done_marker(marker_path, input_assembly, kept_count, removed_count,
                       kingdom=None, karyote=None,
                       keep_classes=None, remove_classes=None):
    """Write a plain-text done marker so subsequent runs can skip this step.

    Args:
        marker_path (str): Path to write the marker file.
        input_assembly (str): Path to the assembly that was decontaminated.
        kept_count (int): Number of sequences retained.
        removed_count (int): Number of sequences removed.
        kingdom (str, optional): Organism kingdom used for classification.
        karyote (str, optional): Organism karyote value.
        keep_classes (set, optional): Tiara classes that were kept.
        remove_classes (set, optional): Tiara classes that were removed.
    """
    with open(marker_path, "w") as fh:
        fh.write(f"Input assembly   : {input_assembly}\n")
        if kingdom is not None:
            fh.write(f"Organism kingdom : {kingdom}\n")
        if karyote is not None:
            fh.write(f"Organism karyote : {karyote}\n")
        if keep_classes is not None:
            fh.write(f"Kept classes     : {', '.join(sorted(keep_classes))}\n")
        if remove_classes is not None:
            fh.write(f"Removed classes  : {', '.join(sorted(remove_classes))}\n")
        fh.write(f"Sequences kept   : {kept_count}\n")
        fh.write(f"Sequences removed: {removed_count}\n")


# --------------------------------------------------------------
# Run Tiara and return path to its classification output file
# --------------------------------------------------------------
def _run_tiara(input_assembly, decontam_dir, cpu_threads):
    """Run Tiara on *input_assembly* and return the path to the output TSV.

    Args:
        input_assembly (str): Path to input FASTA.
        decontam_dir (str): Directory to write Tiara output into.
        cpu_threads (int): Number of threads.

    Returns:
        str or None: Path to Tiara output TSV, or None on failure.
    """
    tiara_out = os.path.join(decontam_dir, "tiara_output.txt")
    if os.path.exists(tiara_out):
        log_print(f"SKIP:\tTiara output already exists: {tiara_out}")
        return tiara_out

    tiara_cmd = (
        f"tiara -i {input_assembly} "
        f"-o {tiara_out} "
        f"--threads {cpu_threads} "
        f"--probabilities"
    )
    rc = run_subprocess_cmd(tiara_cmd, shell_check=True)
    if rc != 0 or not os.path.exists(tiara_out):
        log_print(f"ERROR:\tTiara failed (rc={rc}). Output not found: {tiara_out}")
        return None
    return tiara_out


# --------------------------------------------------------------
# Parse Tiara TSV into {sequence_id: class_label}
# --------------------------------------------------------------
def _parse_tiara(tiara_out):
    """Parse a Tiara output file into a sequence-to-class mapping.

    Tiara writes tab-separated lines: sequence_id <TAB> class <TAB> ...

    Args:
        tiara_out (str): Path to Tiara output TSV.

    Returns:
        dict[str, str]: Mapping of sequence ID to lowercase class label.
    """
    classifications = {}
    with open(tiara_out, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            seq_id = parts[0].strip()
            label  = parts[1].strip().lower()
            classifications[seq_id] = label
    return classifications


# --------------------------------------------------------------
# Main decontamination function
# --------------------------------------------------------------
def decontaminate_assembly(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """Decontaminate an assembly by removing non-target Tiara-classified sequences.

    Args:
        sample_id (str): Sample identifier.
        input_csv (str): Path to metadata CSV.
        output_dir (str): Root output directory.
        cpu_threads (int): CPU threads available.
        ram_gb (int): RAM in GB available (reserved for future use).

    Returns:
        str or None: Path to the decontaminated assembly FASTA, or None on failure.
    """
    cpu_threads = int(cpu_threads)
    ram_gb      = int(ram_gb)

    # ----------------------------------------------------------
    # Section 1: Read CSV metadata
    # ----------------------------------------------------------
    input_df = pd.read_csv(input_csv)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]

    species_id          = current_series["SPECIES_ID"]
    kingdom_id          = current_series["ORGANISM_KINGDOM"]
    karyote_id          = current_series["ORGANISM_KARYOTE"]
    ont_raw_reads       = current_series["ONT_RAW_READS"]
    illumina_f_raw      = current_series["ILLUMINA_RAW_F_READS"]
    illumina_r_raw      = current_series["ILLUMINA_RAW_R_READS"]
    pacbio_raw_reads    = current_series["PACBIO_RAW_READS"]

    species_dir = os.path.join(output_dir, species_id)
    sample_dir  = os.path.join(species_dir, sample_id)
    os.makedirs(sample_dir, exist_ok=True)

    log_print(f"NOTE:\tDecontaminating assembly for sample: {sample_id}")
    log_print(f"NOTE:\tSpecies directory: {species_dir}")

    # ----------------------------------------------------------
    # Section 2: Resolve Tiara decontamination profile
    # ----------------------------------------------------------
    keep_classes, remove_classes = get_decontamination_classes(kingdom_id, karyote_id)
    log_print(f"\n[Tiara] Decontamination profile for '{kingdom_id}' ({karyote_id}):")
    log_print(f"  KEEP   : {', '.join(sorted(keep_classes))}")
    log_print(f"  REMOVE : {', '.join(sorted(remove_classes))}")

    # ----------------------------------------------------------
    # Section 3: Check for existing done marker
    # ----------------------------------------------------------
    decontam_dir  = os.path.join(sample_dir, "decontamination")
    done_marker   = os.path.join(decontam_dir, "decontamination_done.txt")
    final_out     = os.path.join(sample_dir, f"{sample_id}_decontaminated.fasta")

    if os.path.exists(done_marker) and os.path.exists(final_out):
        log_print(f"SKIP:\tDecontamination already complete. Done marker: {done_marker}")
        return final_out

    os.makedirs(decontam_dir, exist_ok=True)

    # ----------------------------------------------------------
    # Section 4: Locate and validate input assembly
    # ----------------------------------------------------------
    curated_assembly = os.path.join(sample_dir, f"{sample_id}_final_curated.fasta")
    if not os.path.exists(curated_assembly):
        # Fall back to polished assembly if curation was skipped
        curated_assembly = os.path.join(sample_dir, f"{sample_id}_final_polish_assembly.fasta")

    if not os.path.exists(curated_assembly):
        log_print(f"ERROR:\tNo input assembly found for decontamination: {curated_assembly}")
        return None

    if os.path.getsize(curated_assembly) < 100:
        log_print(f"ERROR:\tInput assembly is suspiciously small: {curated_assembly}")
        return None

    log_print(f"NOTE:\tInput assembly: {curated_assembly}")

    # ----------------------------------------------------------
    # Section 5: Count input sequences
    # ----------------------------------------------------------
    input_records = list(SeqIO.parse(curated_assembly, "fasta"))
    if not input_records:
        log_print(f"ERROR:\tNo sequences parsed from input assembly: {curated_assembly}")
        return None
    log_print(f"NOTE:\tInput sequences: {len(input_records)}")

    # ----------------------------------------------------------
    # Section 6: Run Tiara
    # ----------------------------------------------------------
    tiara_out = _run_tiara(curated_assembly, decontam_dir, cpu_threads)
    if tiara_out is None:
        log_print("ERROR:\tTiara classification failed. Aborting decontamination.")
        return None

    # ----------------------------------------------------------
    # Section 7: Parse Tiara output
    # ----------------------------------------------------------
    classifications = _parse_tiara(tiara_out)
    if not classifications:
        log_print("ERROR:\tTiara output parsed as empty. Aborting decontamination.")
        return None
    log_print(f"NOTE:\tTiara classified {len(classifications)} sequences.")

    # ----------------------------------------------------------
    # Section 8: Classification summary
    # ----------------------------------------------------------
    class_counts = {cls: 0 for cls in ALL_TIARA_CLASSES}
    class_counts["other"] = 0

    for seq_id, label in classifications.items():
        if label in class_counts:
            class_counts[label] += 1
        else:
            class_counts["other"] += 1

    log_print("\n[Tiara] Classification summary:")
    for cls in sorted(ALL_TIARA_CLASSES):
        count  = class_counts[cls]
        action = "KEEP" if cls in keep_classes else "REMOVE"
        log_print(f"  {cls:<12}: {count:>6}  [{action}]")
    if class_counts["other"]:
        log_print(f"  {'other':<12}: {class_counts['other']:>6}  [REMOVE]")

    # ----------------------------------------------------------
    # Section 9: Write decontaminated (kept) FASTA
    # ----------------------------------------------------------
    kept_records    = []
    removed_records = []

    for record in input_records:
        seq_id = record.id.split()[0]
        label  = classifications.get(seq_id, "unknown")
        if label in keep_classes:
            kept_records.append(record)
        else:
            removed_records.append(record)

    decontam_fasta = os.path.join(decontam_dir, f"{sample_id}_tiara_kept.fasta")
    with open(decontam_fasta, "w") as fh:
        SeqIO.write(kept_records, fh, "fasta")
    log_print(f"NOTE:\tWrote {len(kept_records)} kept sequences to: {decontam_fasta}")

    # ----------------------------------------------------------
    # Section 10: Write removed sequences FASTA
    # ----------------------------------------------------------
    removed_fasta = os.path.join(decontam_dir, f"{sample_id}_tiara_removed.fasta")
    with open(removed_fasta, "w") as fh:
        SeqIO.write(removed_records, fh, "fasta")
    log_print(f"NOTE:\tWrote {len(removed_records)} removed sequences to: {removed_fasta}")

    # ----------------------------------------------------------
    # Section 11: Validate decontaminated assembly
    # ----------------------------------------------------------
    if not os.path.exists(decontam_fasta) or os.path.getsize(decontam_fasta) < 10:
        log_print(f"ERROR:\tDecontaminated FASTA missing or empty: {decontam_fasta}")
        return None

    if len(kept_records) == 0:
        log_print("ERROR:\tDecontamination removed ALL sequences. Aborting.")
        return None

    pct_removed = 100.0 * len(removed_records) / max(len(input_records), 1)
    if pct_removed > 50.0:
        log_print(f"WARN:\tMore than 50% of sequences were removed ({pct_removed:.1f}%). "
                  f"Check Tiara output and organism kingdom assignment.")

    # ----------------------------------------------------------
    # Section 12: Copy decontaminated assembly to final path
    # ----------------------------------------------------------
    import shutil
    if os.path.exists(final_out):
        log_print(f"SKIP:\tFinal decontaminated assembly already exists: {final_out}. Overwriting.")
        os.remove(final_out)

    shutil.copy(decontam_fasta, final_out)
    if not os.path.exists(final_out):
        log_print(f"ERROR:\tFailed to copy decontaminated assembly to: {final_out}")
        return None

    log_print(f"PASS:\tDecontaminated assembly written to: {final_out}")

    # ----------------------------------------------------------
    # Section 13: Write done marker
    # ----------------------------------------------------------
    _write_done_marker(
        done_marker,
        input_assembly=curated_assembly,
        kept_count=len(kept_records),
        removed_count=len(removed_records),
        kingdom=kingdom_id,
        karyote=karyote_id,
        keep_classes=keep_classes,
        remove_classes=remove_classes,
    )
    log_print(f"PASS:\tDone marker written: {done_marker}")
    log_print(f"PASS:\tDecontamination complete for {sample_id}: "
              f"{len(kept_records)} kept, {len(removed_records)} removed.")

    return final_out


# --------------------------------------------------------------
# Entry point
# --------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python3 decontaminate_assembly.py <sample_id> <input_csv> "
              "<output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)

    initialize_logging_environment(sys.argv[3])

    result = decontaminate_assembly(
        sys.argv[1],   # sample_id
        sys.argv[2],   # input_csv
        sys.argv[3],   # output_dir
        sys.argv[4],   # cpu_threads
        sys.argv[5],   # ram_gb
    )

    if result is None:
        log_print("FAIL:\tDecontamination did not produce an output assembly.")
        sys.exit(1)

    log_print(f"PASS:\tDecontamination finished: {result}")
    sys.exit(0)
