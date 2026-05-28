#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
file_operations.py

Filesystem, path, and FASTA/FASTQ I/O helpers for EGAP modules.

Extracted from :mod:`utilities` in v3.4.1.  Hosts path-normalisation,
FASTA validation, parallel gzip wrappers (``pigz``), MD5 verification,
base-count tallies, coverage calculation, and a small move helper.

Depends on :mod:`subprocess_runner` (``run_subprocess_cmd`` for the ``pigz`` wrappers).
Callers that still ``from utilities import to_abs`` continue to work via
the re-export shim in :mod:`utilities`.

Stage:
    Cross-cutting (filesystem layer used by every pipeline stage)

Created on Wed Aug 16 2023

Updated on 2026-05-24

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import hashlib
import os
import shutil
import tempfile
from pathlib import Path
from typing import List

from Bio import SeqIO

from subprocess_runner import run_subprocess_cmd


# --------------------------------------------------------------
# Normalize a value to an absolute filesystem path
# --------------------------------------------------------------
def to_abs(path):
    """Return an absolute filesystem path, leaving non-strings unchanged.

    Many pipeline call sites pass values that may be ``str`` paths or
    ``None``/``pandas.NA``/``float('nan')`` (when the corresponding CSV cell
    is empty).  Wrapping the value with this helper keeps the
    ``isinstance(p, str)`` guard centralized so callers can do
    ``to_abs(maybe_path)`` without first checking the type themselves.

    Parameters
    ----------
    path : str or other
        A string path, or any non-string value (typically ``None`` or
        ``pandas.NA``) which should pass through untouched.

    Returns
    -------
    str or original
        ``os.path.abspath(path)`` when *path* is a string, otherwise the
        input is returned untouched so that downstream
        ``pd.notna(...)`` / ``pd.isna(...)`` guards still work on
        sentinel values.
    """
    return os.path.abspath(path) if isinstance(path, str) else path


# --------------------------------------------------------------
# Validate that a FASTA file exists and looks well-formed
# --------------------------------------------------------------
# IUPAC nucleotide alphabet -- the set of characters that may legitimately
# appear in a polished assembly.  Pilon emits ambiguity codes (K, R, Y, S,
# W, M, B, D, H, V) at positions where read evidence is mixed; rejecting
# them would mean rejecting every Pilon-polished assembly that has any
# residual uncertainty after polishing.  T and U are both allowed so the
# same validator works for RNA-derived assemblies too.
_VALID_NUCLEOTIDES = set("ACGTUNRYSWKMBDHV")


def validate_fasta(file_path):
    """Validate that a FASTA file exists, is non-empty, and parses as nucleotide.

    The canonical FASTA sanity check used across the pipeline: rejects
    empty/missing paths, suspiciously small files (<100 bytes), records
    with empty sequences, and records that contain characters outside
    the IUPAC nucleotide alphabet (case-insensitive).

    Notes
    -----
    The accepted character set is the full IUPAC nucleotide alphabet
    (``A C G T U N R Y S W K M B D H V``).  Earlier versions of this
    validator accepted only ``ATCGN`` and so rejected every Pilon-
    polished hybrid assembly (Pilon emits ``K``/``R`` ambiguity codes
    at mixed-evidence positions).  Only the first record is inspected
    once it parses successfully -- matching the historical fast-path
    behaviour of the per-module copies this function replaces.

    Parameters
    ----------
    file_path : str
        Filesystem path to a ``.fasta`` / ``.fa`` file.  ``None`` or an
        empty string is treated as a validation failure.

    Returns
    -------
    bool
        ``True`` when the file exists, is >=100 bytes, and the first
        record's sequence is non-empty and contains only IUPAC
        nucleotide characters.  ``False`` in every failure case
        (missing file, parse error, etc.); diagnostics are printed to
        stdout, including the offending character set when validation
        fails on alphabet grounds so callers can see *which* code was
        unexpected.
    """
    if not file_path:
        print("ERROR:\tFASTA path is empty/None")
        return False
    if not os.path.exists(file_path):
        print(f"ERROR:\tFASTA file not found: {file_path}")
        return False
    if os.path.getsize(file_path) < 100:
        print(f"ERROR:\tFASTA file is suspiciously small: {file_path}")
        return False
    try:
        with open(file_path, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                if not record.seq:
                    print(f"ERROR:\tFASTA file contains empty sequences: {file_path}")
                    return False
                seq_chars = set(str(record.seq).upper())
                bad = seq_chars - _VALID_NUCLEOTIDES
                if bad:
                    print(
                        f"ERROR:\tFASTA file contains non-nucleotide sequences: {file_path} "
                        f"(unexpected chars: {sorted(bad)})"
                    )
                    return False
                return True
    except Exception as e:
        print(f"ERROR:\tInvalid FASTA format in {file_path}: {str(e)}")
        return False
    return False


# --------------------------------------------------------------
# Compress a file using pigz
# --------------------------------------------------------------
def pigz_compress(input_file, cpu_threads):
    """Compress a file using pigz with multiple threads.

    Parameters
    ----------
    input_file : str
        Path to the file to compress.
    cpu_threads : int
        Number of threads to use for compression.

    Returns
    -------
    str
        Path to the compressed ``.gz`` file.
    """
    pigz_cmd = f"pigz -p {cpu_threads} {input_file}"
    _ = run_subprocess_cmd(pigz_cmd, shell_check=True)
    gzip_file = input_file + ".gz"
    return gzip_file


# --------------------------------------------------------------
# Decompress a file using pigz
# --------------------------------------------------------------
def pigz_decompress(input_file, cpu_threads):
    """Decompress a file using pigz with multiple threads.

    Parameters
    ----------
    input_file : str
        Path to the ``.gz`` file to decompress.
    cpu_threads : int
        Number of threads to use for decompression.

    Returns
    -------
    str
        Path to the decompressed file (with ``.gz`` extension removed).
    """
    pigz_cmd = f"pigz -p {cpu_threads} -d -f {input_file}"
    _ = run_subprocess_cmd(pigz_cmd, shell_check=True)
    unzip_file = input_file.replace(".gz", "")
    return unzip_file


# --------------------------------------------------------------
# Catches and unzips compressed files for FASTQ
# --------------------------------------------------------------
def sum_fastq_bases_with_pigz_safe(fq_path: Path, cpu_threads: int) -> int:
    """
    Safely sum read lengths from a FASTQ that may be .gz by:
    - copying the .gz to a temp dir,
    - pigz_decompress() on the *copy*,
    - parsing the decompressed temp file,
    - removing the temp dir.
    """
    total = 0
    if str(fq_path).endswith(".gz"):
        with tempfile.TemporaryDirectory(prefix="egap_pigz_") as tdir:
            tmp_gz = Path(tdir) / fq_path.name
            shutil.copy2(str(fq_path), str(tmp_gz))                 # safe copy
            tmp_fastq = pigz_decompress(str(tmp_gz), cpu_threads)   # your function returns str path
            with open(tmp_fastq, "rt") as handle:
                for rec in SeqIO.parse(handle, "fastq"):
                    total += len(rec.seq)
            # temp dir cleanup happens automatically
    else:
        with open(fq_path, "rt") as handle:
            for rec in SeqIO.parse(handle, "fastq"):
                total += len(rec.seq)
    return total


# --------------------------------------------------------------
# Catches and unzips compressed files for FASTA
# --------------------------------------------------------------
def sum_fasta_bases_with_pigz_safe(fa_path: Path, cpu_threads: int) -> int:
    """
    Same safety pattern for FASTA/FA files that may be .gz.
    """
    total = 0
    if str(fa_path).endswith(".gz"):
        with tempfile.TemporaryDirectory(prefix="egap_pigz_") as tdir:
            tmp_gz = Path(tdir) / fa_path.name
            shutil.copy2(str(fa_path), str(tmp_gz))
            tmp_fa = pigz_decompress(str(tmp_gz), cpu_threads)      # returns str
            with open(tmp_fa, "rt") as handle:
                for rec in SeqIO.parse(handle, "fasta"):
                    total += len(rec.seq)
    else:
        with open(fa_path, "rt") as handle:
            for rec in SeqIO.parse(handle, "fasta"):
                total += len(rec.seq)
    return total


# --------------------------------------------------------------
# Calculates Coverage for a given Genome
# --------------------------------------------------------------
def calculate_genome_coverage(read_fastqs: List[str], assembly_fasta: str, cpu_threads: int) -> float:
    """
    Compute genome coverage = (total bases in all reads) / (total bases in the assembly).
    Uses your pigz_{compress,decompress} functions safely (no mutation of inputs).
    Supports .fastq/.fq(.gz) for reads and .fa/.fasta(.gz) for assembly.
    """
    total_bases = 0
    for fq in read_fastqs:
        fq_path = Path(fq)
        total_bases += sum_fastq_bases_with_pigz_safe(fq_path, cpu_threads)

    assembly_bases = sum_fasta_bases_with_pigz_safe(Path(assembly_fasta), cpu_threads)
    if assembly_bases == 0:
        raise ValueError(f"No contigs found in {assembly_fasta}")

    return total_bases / assembly_bases


# --------------------------------------------------------------
# Verify MD5 checksums for Illumina files
# --------------------------------------------------------------
def md5_check(folder_name, illumina_df):
    """Verify MD5 checksums for Illumina files in the specified folder.

    Compares computed MD5 checksums against those listed in ``MD5.txt``.
    Prints a PASS or ERROR message for each file.

    Parameters
    ----------
    folder_name : str
        Directory containing Illumina FASTQ files and ``MD5.txt``.
    illumina_df : pandas.DataFrame
        DataFrame used to accumulate MD5 and filename entries.
    """
    md5_file = os.path.join(folder_name, "MD5.txt")
    if not os.path.exists(md5_file):
        print(f"WARNING: MD5.txt not found in {folder_name}. Skipping MD5 check.")
        return
    with open(md5_file, "r") as f:
        for line in f:
            md5, filename = line.strip().split()
            illumina_df = illumina_df.append({"MD5": md5, "Filename": filename}, ignore_index=True)

    for index, row in illumina_df.iterrows():
        file_path = os.path.join(folder_name, row["Filename"])
        if os.path.exists(file_path):
            with open(file_path, "rb") as f:
                file_hash = hashlib.md5(f.read()).hexdigest()
            if file_hash == row["MD5"]:
                print(f"PASS: MD5 check passed for {row['Filename']}")
            else:
                print(f"ERROR: MD5 check failed for {row['Filename']}. Expected {row['MD5']}, got {file_hash}")
        else:
            print(f"ERROR: File not found for MD5 check: {file_path}")


# --------------------------------------------------------------
# Move a file up the directory tree
# --------------------------------------------------------------
def move_file_up(input_file, up_count):
    """Move a file up the directory hierarchy by a specified number of levels.

    Parameters
    ----------
    input_file : str
        Path to the file to move.
    up_count : int
        Number of directory levels to ascend.

    Returns
    -------
    str
        New path to the moved file, or the original *input_file* path if
        the file does not exist.
    """
    if os.path.exists(input_file):
        move_dir = "/".join(os.path.dirname(input_file).split("/")[:-int(up_count)])
        os.makedirs(move_dir, exist_ok=True)
        new_path = os.path.join(move_dir, os.path.basename(input_file))
        shutil.move(input_file, new_path)
        return new_path
    else:
        print(f"ERROR:\tCannot move a non-existing file: {input_file}")
        return input_file
