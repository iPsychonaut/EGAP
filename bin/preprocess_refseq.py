#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
preprocess_refseq.py

Prepare a standardized Reference Sequence FASTA for a given sample so downstream
assemblies (polish/compare/QC) can consume a predictable file path and name.

Stage:
    Reference Sequence Acquisition (RefSeq / NCBI)

Created on Wed Aug 16 2023

Updated on 2026-04-16

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os
import sys
import glob
import shutil
import zipfile
from typing import Optional

import pandas as pd
from utilities import run_subprocess_cmd, initialize_logging_environment, load_sample_context


# --------------------------------------------------------------
# Preprocess reference sequence data (ABSOLUTE-PATH SAFE)
# --------------------------------------------------------------
def _resolve_cell(series, col):
    """Return a cleaned string for *col*, or None for a blank/placeholder cell."""
    _NULLS = {"", "none", "nan", "null", "na"}
    tmp = series.get(col)
    if tmp is not None and pd.notna(tmp) and str(tmp).strip().lower() not in _NULLS:
        return str(tmp).strip()
    return None


def place_assembly(gca, path, species_id, output_dir_abs, label, cpu_threads):
    """Fetch or normalize an existing assembly into ``<output>/<species>/<label>/``.

    Shared by both supply modes EGAP recognizes for an already-built assembly:

    - *path* given and present  -> copy/normalize into place.
    - *gca* given (or *path* fell through) -> download via the NCBI ``datasets``
      CLI, unzip, locate ``*_genomic.fna`` (or ``.fna.gz``), gunzip if needed.

    *label* tags both the sub-directory and the filename: when *gca* is present
    the stem is ``<species>_<gca>_<label>``, otherwise ``<species>_<label>``.
    ``qc_assessment.final_assessment`` reconstructs the same path, so the two
    must stay in sync. Never uses ``os.chdir``; subprocess calls run via
    ``cd '<dir>' && ...``. Returns the absolute finalized FASTA path, or None.
    """
    import gzip
    from pathlib import Path

    output_dir_abs = Path(output_dir_abs)
    species_dir = output_dir_abs / species_id
    dest_dir = species_dir / label
    dest_dir.mkdir(parents=True, exist_ok=True)

    # Compute target filename
    if gca:
        if "." not in gca:
            print(f"WARN:\t{label} GCA has no version number: {gca}")
        target_fasta = dest_dir / f"{species_id}_{gca}_{label}.fasta"
    else:
        target_fasta = dest_dir / f"{species_id}_{label}.fasta"

    print(f"DEBUG - {label} gca  - {gca}")
    print(f"DEBUG - {label} path - {path}")
    print(f"DEBUG - target_fasta - {target_fasta}")

    # If neither provided, nothing to do
    if not gca and not path:
        return None

    # Case A: explicit path provided — normalize into place
    if path:
        # Anchor relative paths under the project output directory
        src = Path(path)
        if not src.is_absolute():
            src = (output_dir_abs / path)
        src = src.expanduser()

        if src.exists():
            try:
                # If already in place, short-circuit
                if target_fasta.exists():
                    try:
                        if src.resolve().samefile(target_fasta.resolve()):
                            print(f"PASS:\t{label} already in place: {target_fasta.resolve()}")
                            return str(target_fasta.resolve())
                    except Exception:
                        # samefile may fail across filesystems; compare sizes as a weak check
                        if src.stat().st_size == target_fasta.stat().st_size:
                            print(f"PASS:\t{label} appears already normalized: {target_fasta.resolve()}")
                            return str(target_fasta.resolve())
                # Copy into expected location/name
                shutil.copy2(str(src), str(target_fasta))
                print(f"PASS:\tCopied {label} into expected location: {target_fasta.resolve()}")
                return str(target_fasta.resolve())
            except Exception as e:
                print(f"WARN:\tFailed to normalize provided {label} '{src}': {e}. Will try GCA if available.")
        else:
            print(f"ERROR:\tProvided {label} does not exist on disk: {src}. Will try GCA if available.")

    # Case B: Need (or fall back) to download via GCA
    if gca:
        if target_fasta.exists():
            print(f"SKIP:\t{label} GCA already exists: {target_fasta.resolve()}")
            return str(target_fasta.resolve())

        pkg_zip = dest_dir / "ncbi_dataset.zip"
        pkg_dir = dest_dir / "ncbi_dataset"
        gca_dir = pkg_dir / "data" / gca

        # Clean prior partials
        try:
            if pkg_zip.exists():
                pkg_zip.unlink()
        except Exception:
            pass
        if pkg_dir.exists():
            try:
                shutil.rmtree(pkg_dir)
            except Exception:
                pass

        print(f"Downloading {label} Assembly for: {gca}")

        # 1) Download with datasets (no chdir; use cd '<dir>' && ...)
        dl_cmd = (
            f"cd '{dest_dir}' && "
            f"datasets download genome accession {gca} --include genome --filename ncbi_dataset.zip"
        )
        _ = run_subprocess_cmd(dl_cmd, shell_check=True)

        if not pkg_zip.exists():
            raise FileNotFoundError(f"NCBI download failed: {pkg_zip} not found.")

        # 2) Unzip in-place under dest_dir using Python's stdlib (no system
        #    ``unzip`` dep -- some minimal Linux installs and most WSL Ubuntu
        #    images ship without it, and adding a hard system dependency just
        #    for NCBI datasets archives is unnecessary when zipfile is built
        #    into the interpreter).
        try:
            with zipfile.ZipFile(str(pkg_zip)) as zf:
                zf.extractall(str(dest_dir))
            print(f"PASS:\tExtracted {pkg_zip.name} via zipfile to {dest_dir}")
        except zipfile.BadZipFile as exc:
            raise FileNotFoundError(
                f"NCBI dataset archive is corrupt or not a valid zip: {pkg_zip} ({exc})"
            )

        if not gca_dir.exists():
            raise FileNotFoundError(f"Expected unpack dir not found: {gca_dir}")

        # 3) Locate *_genomic.fna or *_genomic.fna.gz
        fna_candidates = list(gca_dir.glob("*_genomic.fna")) + list(gca_dir.glob("*_genomic.fna.gz"))
        if not fna_candidates:
            raise FileNotFoundError(f"No '*_genomic.fna(.gz)' found in {gca_dir}")
        src_fna = fna_candidates[0]

        # 4) Gunzip if needed, then normalize name
        if src_fna.suffix == ".gz":
            tmp_out = dest_dir / f"{species_id}_{gca}_{label}.tmp.fasta"
            with gzip.open(src_fna, "rb") as fin, open(tmp_out, "wb") as fout:
                shutil.copyfileobj(fin, fout)
            tmp_out.replace(target_fasta)
        else:
            shutil.copy2(str(src_fna), str(target_fasta))

        print(f"PASS:\tSuccessfully placed {label}: {target_fasta.resolve()}")
        return str(target_fasta.resolve())

    return None


def preprocess_refseq(
    sample_id: str,
    input_tsv: str,
    output_dir: str,
    cpu_threads: int,
    ram_gb: int,
) -> Optional[str]:
    """Preprocess an existing assembly (REF_SEQ or NT_ASSEMBLY) for the pipeline.

    Behavior
    --------
    - When REF_SEQ / REF_SEQ_GCA is given, fetch/normalize it into
      <output_dir>/<SPECIES_ID>/RefSeq/ (EGAP's reference/QC target).
    - Otherwise, when NT_ASSEMBLY_PATH / NT_ASSEMBLY_GCA is given (the
      annotation-pipeline assembly carried in the shared table), fetch/normalize
      it into <output_dir>/<SPECIES_ID>/NtAssembly/ so the same final QC can run
      on it.
    - Returns the absolute path to the finalized FASTA, or None if neither is
      available.

    See :func:`place_assembly` for the per-source fetch/normalize logic.
    """
    from pathlib import Path

    print(f"Preprocessing existing assembly for {sample_id.split('-')[0]}...")

    ctx = load_sample_context(sample_id, input_tsv, output_dir, cpu_threads, ram_gb)
    output_dir_abs = Path(ctx.output_dir)
    print(f"DEBUG - input_tsv  - {ctx.input_tsv}")
    print(f"DEBUG - output_dir - {output_dir_abs}")
    current_series = ctx.current_series

    species_id = str(current_series["SPECIES_ID"]).strip()

    ref_seq_gca = _resolve_cell(current_series, "REF_SEQ_GCA")
    ref_seq = _resolve_cell(current_series, "REF_SEQ")
    nt_assembly_gca = _resolve_cell(current_series, "NT_ASSEMBLY_GCA")
    nt_assembly_path = _resolve_cell(current_series, "NT_ASSEMBLY_PATH")

    # Priority: EGAP's REF_SEQ first, then the annotation-pipeline NT_ASSEMBLY.
    if ref_seq_gca or ref_seq:
        return place_assembly(ref_seq_gca, ref_seq, species_id,
                              str(output_dir_abs), "RefSeq", cpu_threads)
    if nt_assembly_gca or nt_assembly_path:
        print("NOTE:\tNo REF_SEQ given; fetching NT_ASSEMBLY for final QC.")
        return place_assembly(nt_assembly_gca, nt_assembly_path, species_id,
                              str(output_dir_abs), "NtAssembly", cpu_threads)

    print(f"SKIP:\tSample includes no REF_SEQ or NT_ASSEMBLY: {sample_id}.")
    return None


if __name__ == "__main__":
    # Log raw sys.argv immediately
    print(f"DEBUG: Raw sys.argv = {sys.argv}")
    print(f"DEBUG: Length of sys.argv = {len(sys.argv)}")

    # Check argument count
    if len(sys.argv) != 6:
        print(f"ERROR: Expected 5 arguments (plus script name), got {len(sys.argv)-1}: {sys.argv[1:]}", 
              file=sys.stderr)
        print("Usage: python3 preprocess_refseq.py <sample_id> <input_tsv> <output_dir> <cpu_threads> <ram_gb>", 
              file=sys.stderr)
        sys.exit(1)

    initialize_logging_environment(sys.argv[3], sys.argv[1])

    # Log each argument
    for i, arg in enumerate(sys.argv):
        print(f"DEBUG: sys.argv[{i}] = '{arg}'")

    sample_id   = sys.argv[1]
    input_tsv   = sys.argv[2]
    output_dir  = sys.argv[3]
    cpu_threads = int(sys.argv[4])
    ram_gb      = int(sys.argv[5]) if sys.argv[5].strip() else 8

    print(f"DEBUG: Parsed sample_id = '{sample_id}' {type(sample_id)}")
    print(f"DEBUG: Parsed input_tsv = '{input_tsv}' {type(input_tsv)}")
    print(f"DEBUG: Parsed output_dir = '{output_dir}' {type(output_dir)}")
    print(f"DEBUG: Parsed cpu_threads = '{sys.argv[4]}' {sys.argv[4]} (converted to {cpu_threads}) {type(cpu_threads)}")
    print(f"DEBUG: Parsed ram_gb = '{sys.argv[5]}' {sys.argv[5]} (converted to {ram_gb}) {type(ram_gb)}")

    renamed_gca = preprocess_refseq(sample_id, input_tsv, output_dir, cpu_threads, ram_gb)
