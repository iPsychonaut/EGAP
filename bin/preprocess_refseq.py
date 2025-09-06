#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
preprocess_refseq.py

Prepare a standardized Reference Sequence FASTA for a given sample so downstream
assemblies (polish/compare/QC) can consume a predictable file path and name.

Created on Wed Aug 16 2023

Updated on Wed Sept 3 2025

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import os, sys, glob, shutil
import pandas as pd
from utilities import run_subprocess_cmd, get_current_row_data


# --------------------------------------------------------------
# Preprocess reference sequence data (ABSOLUTE-PATH SAFE)
# --------------------------------------------------------------
def preprocess_refseq(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """Preprocess reference sequence data for the assembly pipeline.

    Behavior
    --------
    - Creates <output_dir>/<SPECIES_ID>/RefSeq/ if missing.
    - If REF_SEQ (path) is provided and exists, copy/normalize it into RefSeq/
      as SPECIES_GCA_RefSeq.fasta (or SPECIES_RefSeq.fasta if no GCA).
    - If REF_SEQ_GCA is provided but the file is not present, download it with
      NCBI 'datasets' CLI, unzip, locate '*_genomic.fna' (or .fna.gz), gunzip
      if needed, and place as SPECIES_GCA_RefSeq.fasta.
    - Never uses os.chdir; subprocess calls run via "cd '<dir>' && ...".
    - Returns absolute path to finalized FASTA, or None if unavailable.
    """
    import gzip
    from pathlib import Path

    print(f"Preprocessing Reference Sequence assembly for {sample_id.split('-')[0]}...")

    # --- Resolve critical paths to absolute early ---
    input_csv_abs  = Path(input_csv).expanduser().resolve()
    output_dir_abs = Path(output_dir).expanduser().resolve()
    print(f"DEBUG - input_csv_abs  - {input_csv_abs}")
    print(f"DEBUG - output_dir_abs - {output_dir_abs}")

    # Load CSV & select current row
    input_df = pd.read_csv(str(input_csv_abs))
    print(f"DEBUG - input_df - {input_df}")

    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]

    # Pull fields (normalize NaNs → None)
    ref_seq_gca = None
    tmp = current_series.get("REF_SEQ_GCA")
    if pd.notna(tmp):
        ref_seq_gca = str(tmp).strip()

    ref_seq = None
    tmp = current_series.get("REF_SEQ")
    if pd.notna(tmp):
        ref_seq = str(tmp).strip()

    species_id = str(current_series["SPECIES_ID"]).strip()

    # Prepare directories (absolute)
    species_dir = output_dir_abs / species_id
    refseq_dir  = species_dir / "RefSeq"
    refseq_dir.mkdir(parents=True, exist_ok=True)

    # Compute target filename
    if ref_seq_gca and ref_seq_gca.lower() != "nan":
        if "." not in ref_seq_gca:
            print(f"ERROR:\tReference Sequence GCA requires version number: {ref_seq_gca} has no '.#'")
        target_fasta = refseq_dir / f"{species_id}_{ref_seq_gca}_RefSeq.fasta"
    else:
        target_fasta = refseq_dir / f"{species_id}_RefSeq.fasta"

    print(f"DEBUG - ref_seq_gca - {ref_seq_gca}")
    print(f"DEBUG - ref_seq - {ref_seq}")
    print(f"DEBUG - target_fasta - {target_fasta}")

    # If neither provided, nothing to do
    if (not ref_seq_gca or ref_seq_gca.lower() == "nan") and (not ref_seq or ref_seq.lower() == "nan"):
        print(f"SKIP:\tSample does not include Reference Sequence: {sample_id}.")
        return None

    # Case A: REF_SEQ path provided — normalize into place
    if ref_seq and ref_seq.lower() != "nan":
        # Anchor relative REF_SEQ paths under the project output directory
        src = Path(ref_seq)
        if not src.is_absolute():
            src = (output_dir_abs / ref_seq)
        src = src.expanduser()

        if src.exists():
            try:
                # If already in place, short-circuit
                if target_fasta.exists():
                    try:
                        if src.resolve().samefile(target_fasta.resolve()):
                            print(f"PASS:\tReference sequence already in place: {target_fasta.resolve()}")
                            return str(target_fasta.resolve())
                    except Exception:
                        # samefile may fail across filesystems; compare sizes as a weak check
                        if src.stat().st_size == target_fasta.stat().st_size:
                            print(f"PASS:\tReference sequence appears already normalized: {target_fasta.resolve()}")
                            return str(target_fasta.resolve())
                # Copy into expected location/name
                shutil.copy2(str(src), str(target_fasta))
                print(f"PASS:\tCopied REF_SEQ into expected location: {target_fasta.resolve()}")
                return str(target_fasta.resolve())
            except Exception as e:
                print(f"WARN:\tFailed to normalize provided REF_SEQ '{src}': {e}. Will try GCA if available.")
        else:
            print(f"ERROR:\tProvided REF_SEQ does not exist on disk: {src}. Will try GCA if available.")

    # Case B: Need (or fall back) to download via REF_SEQ_GCA
    if ref_seq_gca and ref_seq_gca.lower() != "nan":
        if target_fasta.exists():
            print(f"SKIP:\tREF_SEQ GCA already exists: {target_fasta.resolve()}")
            return str(target_fasta.resolve())

        pkg_zip = refseq_dir / "ncbi_dataset.zip"
        pkg_dir = refseq_dir / "ncbi_dataset"
        gca_dir = pkg_dir / "data" / ref_seq_gca

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

        print(f"Downloading Reference Sequence Assembly for: {ref_seq_gca}")

        # 1) Download with datasets (no chdir; use cd '<dir>' && ...)
        dl_cmd = (
            f"cd '{refseq_dir}' && "
            f"datasets download genome accession {ref_seq_gca} --include genome --filename ncbi_dataset.zip"
        )
        _ = run_subprocess_cmd(dl_cmd, shell_check=True)

        if not pkg_zip.exists():
            raise FileNotFoundError(f"NCBI download failed: {pkg_zip} not found.")

        # 2) Unzip in-place under refseq_dir
        unzip_cmd = f"cd '{refseq_dir}' && unzip -o ncbi_dataset.zip -d ."
        _ = run_subprocess_cmd(unzip_cmd, shell_check=True)

        if not gca_dir.exists():
            raise FileNotFoundError(f"Expected unpack dir not found: {gca_dir}")

        # 3) Locate *_genomic.fna or *_genomic.fna.gz
        fna_candidates = list(gca_dir.glob("*_genomic.fna")) + list(gca_dir.glob("*_genomic.fna.gz"))
        if not fna_candidates:
            raise FileNotFoundError(f"No '*_genomic.fna(.gz)' found in {gca_dir}")
        src_fna = fna_candidates[0]

        # 4) Gunzip if needed, then normalize name
        if src_fna.suffix == ".gz":
            tmp_out = refseq_dir / f"{species_id}_{ref_seq_gca}_RefSeq.tmp.fasta"
            with gzip.open(src_fna, "rb") as fin, open(tmp_out, "wb") as fout:
                shutil.copyfileobj(fin, fout)
            tmp_out.replace(target_fasta)
        else:
            shutil.copy2(str(src_fna), str(target_fasta))

        print(f"PASS:\tSuccessfully placed Reference Sequence: {target_fasta.resolve()}")
        return str(target_fasta.resolve())

    # Fallback
    print(f"ERROR:\tNo reference sequence could be resolved for {sample_id}")
    return None


if __name__ == "__main__":
    # Log raw sys.argv immediately
    print(f"DEBUG: Raw sys.argv = {sys.argv}")
    print(f"DEBUG: Length of sys.argv = {len(sys.argv)}")

    # Check argument count
    if len(sys.argv) != 6:
        print(f"ERROR: Expected 5 arguments (plus script name), got {len(sys.argv)-1}: {sys.argv[1:]}", 
              file=sys.stderr)
        print("Usage: python3 preprocess_refseq.py <sample_id> <input_csv> <output_dir> <cpu_threads> <ram_gb>", 
              file=sys.stderr)
        sys.exit(1)

    # Log each argument
    for i, arg in enumerate(sys.argv):
        print(f"DEBUG: sys.argv[{i}] = '{arg}'")

    sample_id   = sys.argv[1]
    input_csv   = sys.argv[2]
    output_dir  = sys.argv[3]
    cpu_threads = int(sys.argv[4])
    ram_gb      = int(sys.argv[5]) if sys.argv[5].strip() else 8

    print(f"DEBUG: Parsed sample_id = '{sample_id}' {type(sample_id)}")
    print(f"DEBUG: Parsed input_csv = '{input_csv}' {type(input_csv)}")
    print(f"DEBUG: Parsed output_dir = '{output_dir}' {type(output_dir)}")
    print(f"DEBUG: Parsed cpu_threads = '{sys.argv[4]}' {sys.argv[4]} (converted to {cpu_threads}) {type(cpu_threads)}")
    print(f"DEBUG: Parsed ram_gb = '{sys.argv[5]}' {sys.argv[5]} (converted to {ram_gb}) {type(ram_gb)}")

    renamed_gca = preprocess_refseq(sample_id, input_csv, output_dir, cpu_threads, ram_gb)
