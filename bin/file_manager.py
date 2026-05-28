#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
file_manager.py

Centralised file and directory management for the EGAP pipeline.
Handles safe removal of intermediate files with size logging and
optional dry-run mode.

Set the environment variable EGAP_DRY_RUN=1 before launching the
pipeline to log what would be removed without actually deleting
anything — useful for auditing before committing to cleanup.

Created on Mon Apr 07 2026

Updated on 2026-04-16

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""

import os
import glob as glob_module
import shutil
from pathlib import Path
from utilities import log_print

def _is_dry_run() -> bool:
    """Return True if dry-run mode is active.

    Checked lazily at each call so that EGAP.py (or EGAP_TUI.py) can set
    ``os.environ["EGAP_DRY_RUN"] = "1"`` *after* this module is imported
    (e.g. after argparse resolves ``--dry_run``) and have the setting take
    effect immediately.  Equivalent to exporting ``EGAP_DRY_RUN=1`` from
    the shell before launching the pipeline.
    """
    return os.environ.get("EGAP_DRY_RUN", "0").strip() == "1"


# --------------------------------------------------------------
# Human-readable byte formatter
# --------------------------------------------------------------
def format_size(n_bytes: int) -> str:
    """Return *n_bytes* as a concise human-readable string (e.g. '1.4 GB')."""
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n_bytes < 1024.0:
            return f"{n_bytes:.1f} {unit}"
        n_bytes /= 1024.0
    return f"{n_bytes:.1f} PB"


# --------------------------------------------------------------
# Measure size of a file or directory tree
# --------------------------------------------------------------
def get_size(path: str) -> int:
    """Return the total byte size of *path* (file or full directory tree).

    Returns 0 if *path* does not exist or cannot be measured.
    """
    try:
        p = Path(path)
        if p.is_file():
            return p.stat().st_size
        if p.is_dir():
            return sum(f.stat().st_size for f in p.rglob("*") if f.is_file())
    except OSError:
        pass
    return 0


# --------------------------------------------------------------
# Remove a single file
# --------------------------------------------------------------
def remove_file(path: str, log_file=None) -> bool:
    """Remove a single file, logging the space freed.

    Parameters
    ----------
    path : str
        Absolute path to the file to remove.
    log_file : str, optional
        Path to the sample log file.  Falls back to the module-level
        ``DEFAULT_LOG_FILE`` set by ``initialize_logging_environment``.

    Returns
    -------
    bool
        ``True`` if the file was removed (or would be removed in dry-run
        mode), ``False`` if the file did not exist or removal failed.
    """
    if not os.path.isfile(path):
        return False
    size_str = format_size(get_size(path))
    if _is_dry_run():
        log_print(f"DRY_RUN:\tWould remove file ({size_str}): {path}", log_file)
        return True
    try:
        os.remove(path)
        log_print(f"NOTE:\tRemoved intermediate file ({size_str} freed): {path}", log_file)
        return True
    except OSError as exc:
        log_print(f"WARN:\tCould not remove file: {path} — {exc}", log_file)
        return False


# --------------------------------------------------------------
# Remove a directory tree
# --------------------------------------------------------------
def remove_dir(path: str, log_file=None) -> bool:
    """Recursively remove a directory, logging the total space freed.

    Parameters
    ----------
    path : str
        Absolute path to the directory to remove.
    log_file : str, optional
        Path to the sample log file.

    Returns
    -------
    bool
        ``True`` if the directory was removed (or would be in dry-run),
        ``False`` if it did not exist or removal failed.
    """
    if not os.path.isdir(path):
        return False
    size_str = format_size(get_size(path))
    if _is_dry_run():
        log_print(f"DRY_RUN:\tWould remove directory ({size_str}): {path}", log_file)
        return True
    try:
        shutil.rmtree(path)
        log_print(f"NOTE:\tRemoved intermediate directory ({size_str} freed): {path}", log_file)
        return True
    except OSError as exc:
        log_print(f"WARN:\tCould not remove directory: {path} — {exc}", log_file)
        return False


# --------------------------------------------------------------
# Remove all paths matching a glob pattern
# --------------------------------------------------------------
def remove_glob(pattern: str, log_file=None) -> int:
    """Remove every file or directory matching *pattern*.

    Parameters
    ----------
    pattern : str
        Shell-style glob pattern (forwarded to ``glob.glob``).
    log_file : str, optional
        Path to the sample log file.

    Returns
    -------
    int
        Number of paths successfully removed.
    """
    matches = glob_module.glob(pattern)
    removed = 0
    for path in matches:
        if os.path.isfile(path):
            if remove_file(path, log_file):
                removed += 1
        elif os.path.isdir(path):
            if remove_dir(path, log_file):
                removed += 1
    return removed


# --------------------------------------------------------------
# Clean up Compleasm / BUSCO lineage download archives
# --------------------------------------------------------------
def clean_lineage_downloads(mb_downloads_dir: str, log_file=None) -> None:
    """Remove lineage ``.tar.gz`` archives once they have been extracted.

    Compleasm writes downloaded lineage archives into ``mb_downloads/``
    alongside a zero-byte ``<lineage>.done`` marker and an extracted
    subdirectory named ``<lineage>/``.  Once those two artefacts exist
    the compressed archive is redundant.

    Parameters
    ----------
    mb_downloads_dir : str
        Path to the ``mb_downloads/`` directory created by Compleasm.
    log_file : str, optional
        Path to the sample log file.
    """
    if not os.path.isdir(mb_downloads_dir):
        return
    for tar_path in glob_module.glob(os.path.join(mb_downloads_dir, "*.tar.gz")):
        base = os.path.basename(tar_path)
        # Archive names look like: agaricales_odb12.2025-07-01.tar.gz
        # The lineage key is everything before the first '.'
        lineage = base.split(".")[0]
        done_marker  = os.path.join(mb_downloads_dir, f"{lineage}.done")
        extracted_dir = os.path.join(mb_downloads_dir, lineage)
        if os.path.exists(done_marker) and os.path.isdir(extracted_dir):
            remove_file(tar_path, log_file)


# --------------------------------------------------------------
# Clean up BWA-mem2 index files alongside a FASTA
# --------------------------------------------------------------
def remove_bwa_indices(fasta_path: str, log_file=None) -> None:
    """Remove all BWA-mem2 index files generated for *fasta_path*.

    ``bwa-mem2 index`` creates five sidecar files next to the FASTA:
    ``.0123``, ``.bwt.2bit.64``, ``.pac``, ``.amb``, ``.ann``.

    Parameters
    ----------
    fasta_path : str
        Absolute path to the FASTA that was indexed.
    log_file : str, optional
        Path to the sample log file.
    """
    for suffix in (".0123", ".bwt.2bit.64", ".pac", ".amb", ".ann"):
        remove_file(fasta_path + suffix, log_file)
