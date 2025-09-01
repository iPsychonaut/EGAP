#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
preprocess_refseq.py

Updated on Sat Mar 29 2025


@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import os, sys, glob, shutil
import pandas as pd
from utilities import run_subprocess_cmd, get_current_row_data

# --------------------------------------------------------------
# Preprocess reference sequence data
# --------------------------------------------------------------
def preprocess_refseq(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """Preprocess reference sequence data for the assembly pipeline.

    Behavior:
      - Creates <output_dir>/<SPECIES_ID>/RefSeq/ if missing.
      - If REF_SEQ (path) exists, copy/move it into RefSeq/ as SPECIES_GCA_RefSeq.fasta (or SPECIES_RefSeq.fasta if no GCA).
      - If REF_SEQ_GCA is provided but file not found, download via NCBI datasets, unzip, find *_genomic.fna(.gz),
        gunzip if necessary, and place as SPECIES_GCA_RefSeq.fasta.

    Returns:
        str or None: Absolute path to the finalized reference FASTA, or None if unavailable.
    """
    from pathlib import Path
    import gzip

    print(f"Preprocessing Reference Sequence assembly for {sample_id.split('-')[0]}...")

    # Load CSV & select current row
    input_df = pd.read_csv(input_csv)
    print(f"DEBUG - input_df - {input_df}")

    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]

    # Pull fields
    ref_seq_gca = (str(current_series["REF_SEQ_GCA"]).strip()
                   if pd.notna(current_series.get("REF_SEQ_GCA")) else None)
    ref_seq = (str(current_series["REF_SEQ"]).strip()
               if pd.notna(current_series.get("REF_SEQ")) else None)
    species_id = str(current_series["SPECIES_ID"]).strip()

    # Prepare directories (no chdir)
    species_dir = Path(output_dir).resolve() / species_id
    refseq_dir = species_dir / "RefSeq"
    refseq_dir.mkdir(parents=True, exist_ok=True)

    # If only GCA was supplied, compute the expected REF_SEQ target path
    # (this is the *final* filename we want to create)
    if ref_seq_gca and ref_seq_gca.lower() != "nan":
        # basic sanity: requires version (e.g., ".1")
        if "." not in ref_seq_gca:
            print(f"ERROR:\tReference Sequence GCA requires version number: {ref_seq_gca} has no '.#'")
        target_fasta = refseq_dir / f"{species_id}_{ref_seq_gca}_RefSeq.fasta"
    else:
        target_fasta = refseq_dir / f"{species_id}_RefSeq.fasta"

    print(f"DEBUG - ref_seq_gca - {ref_seq_gca}")
    print(f"DEBUG - ref_seq - {ref_seq}")

    # Nothing to do if neither a GCA nor a REF_SEQ path is provided
    if (not ref_seq_gca or ref_seq_gca.lower() == "nan") and (not ref_seq or ref_seq.lower() == "nan"):
        print(f"SKIP:\tSample does not include Reference Sequence: {sample_id}.")
        return None

    # Case A: REF_SEQ path provided by CSV (user-supplied file)
    # Normalize and copy/move into expected location/name
    if ref_seq and ref_seq.lower() != "nan":
        src = Path(ref_seq).expanduser().resolve()
        if not src.exists():
            print(f"ERROR:\tProvided REF_SEQ does not exist on disk: {src}")
            # fallthrough to GCA download if available
        else:
            # If user provided a file, normalize it into our expected target filename
            if src.samefile(target_fasta):
                print(f"PASS:\tReference sequence already in place: {target_fasta}")
                return str(target_fasta)
            else:
                # Copy/rename into place
                shutil.copy2(str(src), str(target_fasta))
                print(f"PASS:\tCopied REF_SEQ into expected location: {target_fasta}")
                return str(target_fasta)

    # Case B: Need to download via REF_SEQ_GCA (NCBI datasets)
    if ref_seq_gca and ref_seq_gca.lower() != "nan":
        if target_fasta.exists():
            print(f"SKIP:\tREF_SEQ GCA already exists: {target_fasta}")
            return str(target_fasta)

        # Where datasets will unpack files:
        # datasets creates ./ncbi_dataset.zip and ./ncbi_dataset/...
        # We keep its outputs inside refseq_dir
        pkg_zip = refseq_dir / "ncbi_dataset.zip"
        pkg_dir = refseq_dir / "ncbi_dataset"
        gca_dir = pkg_dir / "data" / ref_seq_gca

        # Clean any previous partials
        if pkg_zip.exists():
            try:
                pkg_zip.unlink()
            except Exception:
                pass
        if pkg_dir.exists():
            try:
                shutil.rmtree(pkg_dir)
            except Exception:
                pass

        print(f"Downloading Reference Sequence Assembly for: {ref_seq_gca}")
        # Use NCBI datasets CLI to fetch genome FASTA
        # (no chdir; run in refseq_dir via cwd=)
        dl_cmd = (
            f"datasets download genome accession {ref_seq_gca} "
            f"--include genome --filename {pkg_zip.name}"
        )
        _ = run_subprocess_cmd(dl_cmd, shell_check=True, cwd=str(refseq_dir))

        # Unzip
        if not pkg_zip.exists():
            raise FileNotFoundError(f"NCBI download failed: {pkg_zip} not found.")
        unzip_cmd = f"unzip -o {pkg_zip.name} -d {refseq_dir.name}"
        # Since -d path is relative in unzip, we run from parent of refseq_dir to place under refseq_dir cleanly
        _ = run_subprocess_cmd(unzip_cmd, shell_check=True, cwd=str(refseq_dir.parent))

        # Locate *_genomic.fna or *_genomic.fna.gz
        if not gca_dir.exists():
            raise FileNotFoundError(f"Expected unpack dir not found: {gca_dir}")

        fna_candidates = list(gca_dir.glob("*_genomic.fna")) + list(gca_dir.glob("*_genomic.fna.gz"))
        if not fna_candidates:
            raise FileNotFoundError(f"No ‘*_genomic.fna(.gz)’ found in {gca_dir}")

        src_fna = fna_candidates[0]

        # If gzipped, gunzip to a temp path within refseq_dir, then move to final target
        if src_fna.suffix == ".gz":
            temp_out = refseq_dir / f"{species_id}_{ref_seq_gca}_RefSeq.tmp.fasta"
            with gzip.open(src_fna, "rb") as fin, open(temp_out, "wb") as fout:
                shutil.copyfileobj(fin, fout)
            # Move into place
            temp_out.replace(target_fasta)
        else:
            shutil.copy2(str(src_fna), str(target_fasta))

        print(f"PASS:\tSuccessfully placed Reference Sequence: {target_fasta}")
        return str(target_fasta)

    # Fallback: nothing obtained
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
        print("Usage: python3 preprocess_illumina.py <sample_id> <input_csv> <output_dir> <cpu_threads> <ram_gb>", 
              file=sys.stderr)
        sys.exit(1)
    
    # Log each argument
    for i, arg in enumerate(sys.argv):
        print(f"DEBUG: sys.argv[{i}] = '{arg}'")
    
    sample_id = sys.argv[1]
    input_csv = sys.argv[2]
    output_dir = sys.argv[3]
    cpu_threads = int(sys.argv[4])
    ram_gb = int(sys.argv[5]) if sys.argv[5] != " " else 8
    
    print(f"DEBUG: Parsed sample_id = '{sample_id}' {type(sample_id)}")
    print(f"DEBUG: Parsed input_csv = '{input_csv}' {type(input_csv)}")
    print(f"DEBUG: Parsed output_dir = '{output_dir}' {type(output_dir)}")
    print(f"DEBUG: Parsed cpu_threads = '{sys.argv[4]}' {sys.argv[4]} (converted to {cpu_threads}) {type(cpu_threads)}")
    print(f"DEBUG: Parsed ram_gb = '{sys.argv[5]}' {sys.argv[5]} (converted to {ram_gb}) {type(ram_gb)}")
    
    renamed_gca = preprocess_refseq(sample_id, input_csv, output_dir, cpu_threads, ram_gb)

