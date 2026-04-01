#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
decontaminate_reads.py

Removes contaminant reads from ONT and PacBio long-read FASTQ files using
Kraken2 classification before genome assembly.

Kraken2 assigns each read a taxon ID from the NCBI taxonomy.  This script
parses the Kraken2 report to build a taxon-ID-to-domain map, then keeps
only reads whose domain matches the target organism's kingdom:

    Kingdom            Domains kept
    -----------------  ----------------------------------
    bacteria           bacteria, unclassified
    archaea            archaea,  unclassified
    flora/funga/fauna  eukarya,  unclassified
    (unrecognised)     eukarya,  unclassified  (with WARN)

'unclassified' reads are always kept because they may represent genuine
target sequence that is simply absent from the Kraken2 database.

The Kraken2 database path is resolved in this order:
  1. KRAKEN2_DB environment variable
  2. A 'KRAKEN2_DB' column in the input CSV (if present)
  If neither is available the step is skipped with a WARN rather than
  aborting the pipeline.

Decontaminated reads overwrite the highest-mean-quality long-read file
that the assemblers already expect, preserving the original as
*_pre_decontam.fastq so it can be recovered if needed.

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
# Domain labels produced by Kraken2 (mapped from NCBI rank "D")
# --------------------------------------------------------------
ALL_KRAKEN_DOMAINS = {"bacteria", "archaea", "eukarya", "viruses", "unclassified", "other"}

# Domain taxon IDs at rank D in the NCBI taxonomy (used as fallback labels)
_DOMAIN_TAXIDS = {
    2:     "bacteria",
    2157:  "archaea",
    2759:  "eukarya",
    10239: "viruses",
    0:     "unclassified",
}


# --------------------------------------------------------------
# Resolve which Kraken2 domains to keep for a given kingdom
# --------------------------------------------------------------
def get_kraken_keep_domains(kingdom_id):
    """Return the set of Kraken2 domain labels to KEEP for *kingdom_id*.

    Args:
        kingdom_id (str or None): Value of ORGANISM_KINGDOM from the CSV.
            Case-insensitive.

    Returns:
        set[str]: Domain labels (subset of ALL_KRAKEN_DOMAINS) whose reads
            should be retained.  Always includes 'unclassified'.
    """
    if not isinstance(kingdom_id, str) or not kingdom_id.strip():
        log_print(f"WARN:\tMissing or blank kingdom_id. "
                  f"Falling back to eukaryote Kraken2 keep-domain profile.")
        return {"eukarya", "unclassified"}

    normalised = kingdom_id.strip().lower()
    if normalised == "bacteria":
        return {"bacteria", "unclassified"}
    if normalised == "archaea":
        return {"archaea", "unclassified"}
    if normalised in ("flora", "funga", "fauna"):
        return {"eukarya", "unclassified"}

    log_print(f"WARN:\tUnrecognised kingdom_id '{kingdom_id}'. "
              f"Falling back to eukaryote Kraken2 keep-domain profile.")
    return {"eukarya", "unclassified"}


# --------------------------------------------------------------
# Locate the Kraken2 database
# --------------------------------------------------------------
def _get_kraken2_db(current_series):
    """Return the Kraken2 database path or None if not configured.

    Checks (in order):
      1. KRAKEN2_DB environment variable
      2. 'KRAKEN2_DB' column in the CSV row (if it exists)

    Args:
        current_series (pandas.Series): One row from the metadata CSV.

    Returns:
        str or None: Absolute path to the Kraken2 database directory.
    """
    db = os.environ.get("KRAKEN2_DB", "").strip()
    if not db and "KRAKEN2_DB" in current_series.index:
        val = current_series["KRAKEN2_DB"]
        if isinstance(val, str) and val.strip():
            db = val.strip()
    if not db:
        return None
    if not os.path.isdir(db):
        log_print(f"WARN:\tKraken2 database path does not exist: {db}")
        return None
    return db


# --------------------------------------------------------------
# Build a taxid -> domain label map from a Kraken2 report file
# --------------------------------------------------------------
def _build_taxid_domain_map(report_path):
    """Parse a kreport2 file and return a {taxid: domain_label} dict.

    The kreport2 format is a tab-separated file:
        pct  cum_reads  direct_reads  rank  taxid  name

    The 'name' column is indented with two spaces per taxonomic level.
    This function tracks the current domain label as it descends the tree.

    Args:
        report_path (str): Path to the Kraken2 --report output file.

    Returns:
        dict[int, str]: taxid -> domain label string.
    """
    taxid_domain = {0: "unclassified"}
    domain_stack = []   # list of (indent_depth, domain_label)

    with open(report_path, "r") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            rank  = parts[3].strip()
            try:
                taxid = int(parts[4].strip())
            except ValueError:
                continue
            name_field   = parts[5]                          # may be indented
            name_clean   = name_field.lstrip()
            indent        = len(name_field) - len(name_clean)
            name_lower   = name_clean.lower()

            # Pop stack entries at same or greater indentation
            while domain_stack and domain_stack[-1][0] >= indent:
                domain_stack.pop()

            if rank == "D":
                # Assign a canonical domain label from the name
                if "bacteria" in name_lower:
                    label = "bacteria"
                elif "archaea" in name_lower:
                    label = "archaea"
                elif "eukaryota" in name_lower or "eukarya" in name_lower:
                    label = "eukarya"
                elif "virus" in name_lower or "viroid" in name_lower:
                    label = "viruses"
                else:
                    label = "other"
                domain_stack.append((indent, label))
                taxid_domain[taxid] = label
            elif domain_stack:
                taxid_domain[taxid] = domain_stack[-1][1]
            else:
                taxid_domain[taxid] = "other"

    return taxid_domain


# --------------------------------------------------------------
# Parse Kraken2 per-read output -> {read_id: taxid}
# --------------------------------------------------------------
def _parse_kraken_reads(kraken_output_path):
    """Parse a Kraken2 per-read output file into a read-ID-to-taxid dict.

    Kraken2 --output format (one line per read):
        C/U  read_id  taxid  length  kmer_mapping

    Args:
        kraken_output_path (str): Path to the Kraken2 --output file.

    Returns:
        dict[str, int]: read_id -> taxid (0 for unclassified).
    """
    read_taxid = {}
    with open(kraken_output_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            read_id = parts[1].strip()
            try:
                taxid = int(parts[2].strip())
            except ValueError:
                taxid = 0
            read_taxid[read_id] = taxid
    return read_taxid


# --------------------------------------------------------------
# Run Kraken2 on a FASTQ file
# --------------------------------------------------------------
def _run_kraken2(input_reads, kraken_dir, kraken2_db, cpu_threads, label):
    """Run Kraken2 classification on *input_reads*.

    Args:
        input_reads (str): Path to input FASTQ (may be .gz).
        kraken_dir (str): Working directory for Kraken2 output files.
        kraken2_db (str): Path to the Kraken2 database directory.
        cpu_threads (int): Number of threads.
        label (str): Short label used in output filenames (e.g. "ONT").

    Returns:
        tuple[str, str] or (None, None): (kraken_output_path, kraken_report_path)
            or (None, None) on failure.
    """
    os.makedirs(kraken_dir, exist_ok=True)
    kraken_out    = os.path.join(kraken_dir, f"{label}_kraken2.out")
    kraken_report = os.path.join(kraken_dir, f"{label}_kraken2_report.txt")

    if os.path.exists(kraken_out) and os.path.exists(kraken_report):
        log_print(f"SKIP:\tKraken2 output already exists for {label}: {kraken_out}")
        return kraken_out, kraken_report

    kraken_cmd = (
        f"kraken2 "
        f"--db {kraken2_db} "
        f"--threads {cpu_threads} "
        f"--output {kraken_out} "
        f"--report {kraken_report} "
        f"--report-zero-counts "
        f"{input_reads}"
    )
    rc = run_subprocess_cmd(kraken_cmd, shell_check=True)
    if rc != 0 or not os.path.exists(kraken_out) or not os.path.exists(kraken_report):
        log_print(f"ERROR:\tKraken2 failed for {label} (rc={rc}).")
        return None, None

    return kraken_out, kraken_report


# --------------------------------------------------------------
# Filter a FASTQ by a set of read IDs to keep
# --------------------------------------------------------------
def _filter_fastq(input_fastq, keep_read_ids, output_fastq):
    """Write records from *input_fastq* whose ID is in *keep_read_ids*.

    Args:
        input_fastq (str): Path to input FASTQ (plain or .gz).
        keep_read_ids (set[str]): Read IDs to retain.
        output_fastq (str): Path to write filtered FASTQ.

    Returns:
        int: Number of records written.
    """
    fmt = "fastq"
    kept = 0
    with open(output_fastq, "w") as out_fh:
        with open(input_fastq, "rt") as in_fh:
            for record in SeqIO.parse(in_fh, fmt):
                # Kraken2 uses the portion before the first space as the ID
                rid = record.id.split()[0]
                if rid in keep_read_ids:
                    SeqIO.write(record, out_fh, fmt)
                    kept += 1
    return kept


# --------------------------------------------------------------
# Decontaminate one long-read FASTQ using Kraken2
# --------------------------------------------------------------
def _decontaminate_one(reads_path, reads_type, kraken_dir, kraken2_db,
                       keep_domains, cpu_threads):
    """Run Kraken2 on *reads_path* and return the filtered FASTQ path.

    Keeps reads whose Kraken2-assigned domain is in *keep_domains*.
    Saves the original reads as *_pre_decontam.fastq before overwriting.

    Args:
        reads_path (str): Path to the input FASTQ to decontaminate.
        reads_type (str): 'ONT' or 'PacBio' (used in filenames and logs).
        kraken_dir (str): Directory for Kraken2 intermediates.
        kraken2_db (str): Kraken2 database path.
        keep_domains (set[str]): Domain labels to keep.
        cpu_threads (int): Thread count.

    Returns:
        str or None: Path to the decontaminated FASTQ, or None on failure.
    """
    label = reads_type

    # Run Kraken2
    kraken_out, kraken_report = _run_kraken2(
        reads_path, kraken_dir, kraken2_db, cpu_threads, label
    )
    if kraken_out is None:
        return None

    # Build taxid -> domain map from the report
    taxid_domain = _build_taxid_domain_map(kraken_report)

    # Parse per-read classifications
    read_taxid = _parse_kraken_reads(kraken_out)

    # Determine which read IDs to keep
    keep_read_ids = set()
    domain_counts = {d: 0 for d in ALL_KRAKEN_DOMAINS}
    for read_id, taxid in read_taxid.items():
        domain = taxid_domain.get(taxid, "other")
        domain_counts[domain] = domain_counts.get(domain, 0) + 1
        if domain in keep_domains:
            keep_read_ids.add(read_id)

    # Classification summary
    total = max(len(read_taxid), 1)
    log_print(f"\n[Kraken2] {reads_type} classification summary:")
    for domain in sorted(ALL_KRAKEN_DOMAINS):
        count  = domain_counts.get(domain, 0)
        action = "KEEP" if domain in keep_domains else "REMOVE"
        log_print(f"  {domain:<14}: {count:>8}  ({100*count/total:5.1f}%)  [{action}]")

    pct_kept = 100.0 * len(keep_read_ids) / total
    log_print(f"  Total kept : {len(keep_read_ids):>8}  ({pct_kept:.1f}%)")

    if pct_kept < 50.0:
        log_print(f"WARN:\tFewer than 50% of {reads_type} reads were kept "
                  f"({pct_kept:.1f}%). Check kingdom assignment and Kraken2 database.")

    # Back up the original reads then write decontaminated file in its place
    pre_decontam = reads_path.replace(".fastq", "_pre_decontam.fastq")
    if not os.path.exists(pre_decontam):
        import shutil
        shutil.copy(reads_path, pre_decontam)
        log_print(f"NOTE:\tOriginal reads preserved as: {pre_decontam}")

    decontam_path = reads_path  # overwrite the file the assemblers expect
    kept_count = _filter_fastq(pre_decontam, keep_read_ids, decontam_path)

    if kept_count == 0:
        log_print(f"ERROR:\tKraken2 decontamination removed ALL {reads_type} reads. "
                  f"Restoring original.")
        import shutil
        shutil.copy(pre_decontam, reads_path)
        return reads_path  # return original so pipeline can continue

    log_print(f"PASS:\t{reads_type} decontamination complete: "
              f"{kept_count} reads kept -> {decontam_path}")
    return decontam_path


# --------------------------------------------------------------
# Main entry point
# --------------------------------------------------------------
def decontaminate_reads(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """Kraken2 read-level decontamination for ONT and PacBio long reads.

    Runs after preprocessing (filtering/correction) and before assembly.
    Decontaminated reads overwrite the *_highest_mean_qual_long_reads.fastq
    files that the assemblers look for, so no assembler changes are needed.

    Args:
        sample_id (str): Sample identifier.
        input_csv (str): Path to metadata CSV.
        output_dir (str): Root output directory.
        cpu_threads (int): CPU threads available.
        ram_gb (int): RAM in GB (reserved for future use).

    Returns:
        bool: True if at least one read set was processed or skipped cleanly,
              False on unrecoverable failure.
    """
    cpu_threads = int(cpu_threads)

    # ----------------------------------------------------------
    # Section 1: Read CSV metadata
    # ----------------------------------------------------------
    input_df = pd.read_csv(os.path.abspath(input_csv))
    current_row, _, _ = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]

    species_id       = current_series["SPECIES_ID"]
    kingdom_id       = current_series["ORGANISM_KINGDOM"]
    karyote_id       = current_series["ORGANISM_KARYOTE"]
    ont_raw_reads    = current_series["ONT_RAW_READS"]
    pacbio_raw_reads = current_series["PACBIO_RAW_READS"]
    ont_sra          = current_series["ONT_SRA"]
    pacbio_sra       = current_series["PACBIO_SRA"]

    species_dir  = os.path.join(os.path.abspath(output_dir), species_id)
    ont_dir      = os.path.join(species_dir, "ONT")
    pacbio_dir   = os.path.join(species_dir, "PacBio")

    log_print(f"NOTE:\tKraken2 read decontamination for sample: {sample_id}")

    # ----------------------------------------------------------
    # Section 2: Resolve Kraken2 database
    # ----------------------------------------------------------
    kraken2_db = _get_kraken2_db(current_series)
    if kraken2_db is None:
        log_print("WARN:\tNo Kraken2 database found (set KRAKEN2_DB env var or "
                  "add a KRAKEN2_DB column to the CSV). Skipping read decontamination.")
        return True  # non-fatal skip

    log_print(f"NOTE:\tKraken2 database: {kraken2_db}")

    # ----------------------------------------------------------
    # Section 3: Resolve keep domains from kingdom
    # ----------------------------------------------------------
    keep_domains = get_kraken_keep_domains(kingdom_id)
    log_print(f"[Kraken2] Decontamination profile for '{kingdom_id}' ({karyote_id}):")
    log_print(f"  KEEP domains: {', '.join(sorted(keep_domains))}")

    # ----------------------------------------------------------
    # Section 4: Locate highest-quality long-read files
    # ----------------------------------------------------------
    # These are the canonical paths that assemblers look for.
    ont_hq   = os.path.join(ont_dir,    f"{species_id}_ONT_highest_mean_qual_long_reads.fastq")
    pacbio_hq = os.path.join(pacbio_dir, f"{species_id}_PacBio_highest_mean_qual_long_reads.fastq")

    has_ont    = (pd.notna(ont_sra) or pd.notna(ont_raw_reads)) and os.path.exists(ont_hq)
    has_pacbio = (pd.notna(pacbio_sra) or pd.notna(pacbio_raw_reads)) and os.path.exists(pacbio_hq)

    if not has_ont and not has_pacbio:
        log_print("SKIP:\tNo preprocessed long-read files found for Kraken2 decontamination.")
        return True

    # ----------------------------------------------------------
    # Section 5: Check done marker
    # ----------------------------------------------------------
    kraken_dir  = os.path.join(species_dir, "kraken2_reads")
    done_marker = os.path.join(kraken_dir, "decontaminate_reads_done.txt")

    if os.path.exists(done_marker):
        log_print(f"SKIP:\tKraken2 read decontamination already complete: {done_marker}")
        return True

    os.makedirs(kraken_dir, exist_ok=True)

    # ----------------------------------------------------------
    # Section 6: Decontaminate ONT reads
    # ----------------------------------------------------------
    if has_ont:
        log_print(f"NOTE:\tDecontaminating ONT reads: {ont_hq}")
        result = _decontaminate_one(
            ont_hq, "ONT", kraken_dir, kraken2_db, keep_domains, cpu_threads
        )
        if result is None:
            log_print("ERROR:\tONT Kraken2 decontamination failed.")
            return False
    else:
        log_print("SKIP:\tNo ONT highest-quality reads found for Kraken2.")

    # ----------------------------------------------------------
    # Section 7: Decontaminate PacBio reads
    # ----------------------------------------------------------
    if has_pacbio:
        log_print(f"NOTE:\tDecontaminating PacBio reads: {pacbio_hq}")
        result = _decontaminate_one(
            pacbio_hq, "PacBio", kraken_dir, kraken2_db, keep_domains, cpu_threads
        )
        if result is None:
            log_print("ERROR:\tPacBio Kraken2 decontamination failed.")
            return False
    else:
        log_print("SKIP:\tNo PacBio highest-quality reads found for Kraken2.")

    # ----------------------------------------------------------
    # Section 8: Write done marker
    # ----------------------------------------------------------
    with open(done_marker, "w") as fh:
        fh.write(f"Sample           : {sample_id}\n")
        fh.write(f"Organism kingdom : {kingdom_id}\n")
        fh.write(f"Organism karyote : {karyote_id}\n")
        fh.write(f"Kept domains     : {', '.join(sorted(keep_domains))}\n")
        fh.write(f"Kraken2 database : {kraken2_db}\n")
        fh.write(f"ONT processed    : {has_ont}\n")
        fh.write(f"PacBio processed : {has_pacbio}\n")

    log_print(f"PASS:\tKraken2 read decontamination complete for {sample_id}.")
    return True


# --------------------------------------------------------------
# Entry point
# --------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python3 decontaminate_reads.py <sample_id> <input_csv> "
              "<output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)

    initialize_logging_environment(sys.argv[3])

    success = decontaminate_reads(
        sys.argv[1],   # sample_id
        sys.argv[2],   # input_csv
        sys.argv[3],   # output_dir
        sys.argv[4],   # cpu_threads
        sys.argv[5],   # ram_gb
    )

    if not success:
        log_print("FAIL:\tKraken2 read decontamination did not complete successfully.")
        sys.exit(1)

    sys.exit(0)
