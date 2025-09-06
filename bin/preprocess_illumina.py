#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
preprocess_illumina.py

This script preprocesses Illumina reads with FastQC, Trimmomatic, BBDuk, and Clumpify.
Reads input from a CSV file, processes all Illumina datasets, and updates the CSV
with declumpified FASTQ file paths in ${params.output_dir}/${sample_prefix}/Illumina/.

Created on Wed Aug 16 2023

Updated on Wed Sept 3 2025

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import os, sys, shutil
import pandas as pd
from utilities import run_subprocess_cmd, get_current_row_data, md5_check


# --------------------------------------------------------------
# Combine and verify Illumina reads
# --------------------------------------------------------------
def illumina_extract_and_check(folder_name, SAMPLE_ID):
    """Combine paired-end Illumina reads after MD5 verification.

    Verifies MD5 checksums, concatenates forward and reverse FASTQ files.

    Args:
        folder_name (str): Directory containing Illumina FASTQ files and MD5.txt.
        SAMPLE_ID (str): Sample identifier for naming output files.

    Returns:
        list or None: Paths to combined forward and reverse files, or None if failed.
    """
    print(f"Running MD5 Checksum Analysis on Raw Illumina FASTQ files in {folder_name}...")
    illumina_df = pd.DataFrame(columns=["MD5", "Filename"])
    base_folder = SAMPLE_ID.split("-")[0]  # e.g., "Es_coli"
    illumina_dir = os.path.join(os.getcwd(), base_folder, "Illumina")  # Use current dir, not hardcoded EGAP_Test_Data
    os.makedirs(illumina_dir, exist_ok=True)
    combined_1_file = os.path.join(illumina_dir, f"{SAMPLE_ID}_combined_1.fq")
    combined_2_file = os.path.join(illumina_dir, f"{SAMPLE_ID}_combined_2.fq")
    combined_list = [combined_1_file, combined_2_file]
    
    if not os.path.isfile(combined_list[0]) or not os.path.isfile(combined_list[1]):
        if not os.path.isfile(combined_1_file) or not os.path.isfile(combined_2_file):
            md5_check(folder_name, illumina_df)
            raw_1_list = []
            raw_2_list = []
            for filename in os.listdir(folder_name):
                if "_1.fastq" in filename:
                    raw_1_list.append(os.path.join(folder_name, filename))
                elif "_2.fastq" in filename:
                    raw_2_list.append(os.path.join(folder_name, filename))
            if not raw_1_list or not raw_2_list:
                print(f"ERROR:\tNo paired Illumina files found in {folder_name}")
                return None
            fwd_cat_cmd = f"cat {' '.join(raw_1_list)} > {combined_1_file}"
            _ = run_subprocess_cmd(fwd_cat_cmd, shell_check=True)
            rev_cat_cmd = f"cat {' '.join(raw_2_list)} > {combined_2_file}"
            _ = run_subprocess_cmd(rev_cat_cmd, shell_check=True)
        else:
            print(f"SKIP:\tCombined FASTQ files already exist: {combined_list[0]}; {combined_list[1]}.")
    else:
        print(f"SKIP:\tGzipped Combined FASTQ files already exist: {combined_list[0]}; {combined_list[1]}.")

    return combined_list if os.path.exists(combined_list[0]) and os.path.exists(combined_list[1]) else None


# --------------------------------------------------------------
# Preprocess Illumina sequencing reads
# --------------------------------------------------------------
def preprocess_illumina(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """Preprocess Illumina reads with FastQC, Trimmomatic, BBDuk, and Clumpify.

    Behavior
    --------
    - Creates <output_dir>/<SPECIES_ID>/Illumina/ if missing.
    - If ILLUMINA_SRA is provided (and no raw read paths), ensure {SRA}_1.fastq and {SRA}_2.fastq exist
      in Illumina/. If missing, run:
        prefetch --force yes <SRA>
        fasterq-dump --split-files -e <threads> -O <Illumina/> <SRA>
      Fallback: fastq-dump --split-files -O <Illumina/> <SRA>
    - If ILLUMINA_RAW_DIR is provided, combine verified pairs via illumina_extract_and_check(...).
    - Then run FastQC → Trimmomatic → BBDuk → Clumpify.
    - Returns (illu_dedup_f_reads, illu_dedup_r_reads) or (None, None) if unavailable/failed.
    """
    from pathlib import Path

    print(f"Preprocessing Illumina reads for {sample_id.split('-')[0]}...")

    # Parse the CSV and retrieve relevant row data
    input_df = pd.read_csv(input_csv)
    print(f"DEBUG - input_df - {input_df}")

    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]  # Convert to Series (single row)

    print(f"DEBUG - current_series - {current_series}")

    # Identify read paths/reference info from CSV
    illu_sra = current_series["ILLUMINA_SRA"]
    illu_raw_f_reads = current_series["ILLUMINA_RAW_F_READS"]
    illu_raw_r_reads = current_series["ILLUMINA_RAW_R_READS"]
    illu_raw_dir = current_series["ILLUMINA_RAW_DIR"]
    species_id = current_series["SPECIES_ID"]

    print(f"DEBUG - illu_sra - {illu_sra}")
    print(f"DEBUG - illu_raw_f_reads - {illu_raw_f_reads}")
    print(f"DEBUG - illu_raw_r_reads - {illu_raw_r_reads}")
    print(f"DEBUG - illu_raw_dir - {illu_raw_dir}")

    # Early skip if nothing Illumina-like is provided
    if pd.isna(illu_sra) and pd.isna(illu_raw_f_reads) and pd.isna(illu_raw_r_reads) and pd.isna(illu_raw_dir):
        print(f"SKIP:\tSample does not include Illumina Reads: {sample_id}.")
        return None, None

    # Prepare directories
    species_dir = Path(output_dir).resolve() / str(species_id)
    illumina_dir = species_dir / "Illumina"
    illumina_dir.mkdir(parents=True, exist_ok=True)

    # ---- CASE A: SRA accession provided; ensure R1/R2 exist or download ----
    if pd.notna(illu_sra) and pd.isna(illu_raw_f_reads) and pd.isna(illu_raw_r_reads):
        sra_id = str(illu_sra).strip()
        r1 = illumina_dir / f"{sra_id}_1.fastq"
        r2 = illumina_dir / f"{sra_id}_2.fastq"

        if not (r1.exists() and r2.exists()):
            print(f"Downloading Illumina SRA: {sra_id}...")

            # 1) prefetch directly into the Illumina directory (avoid cwd clutter)
            #    This writes SRR*/SRR*.sra under <Illumina/>.
            prefetch_cmd = [
                "prefetch",
                "--force", "yes",
                "--output-directory", str(illumina_dir),
                sra_id,
            ]
            _ = run_subprocess_cmd(prefetch_cmd, False)

            # 2) fasterq-dump with explicit output dir AND temp dir in <Illumina/>.
            #    Many environments honor TMPDIR for large temps; we set it only for this call.
            #    We also set VDB_CONFIG to avoid surprises from user/global configs.
            env_prefix = (
                f"TMPDIR='{illumina_dir}' "
                f"VDB_CONFIG='{illumina_dir}' "
            )
            fq_cmd = (
                env_prefix +
                "fasterq-dump "
                f"--split-files -e {int(cpu_threads)} "
                f"-O '{illumina_dir}' "
                f"'{illumina_dir}/{sra_id}'"
            )
            ret = run_subprocess_cmd(fq_cmd, True)

            # Fallback to fastq-dump if fasterq-dump failed or outputs missing
            if ret != 0 or not (r1.exists() and r2.exists()):
                print("WARN:\tfasterq-dump failed or files not found; trying fastq-dump fallback...")
                fq2_cmd = (
                    env_prefix +
                    "fastq-dump "
                    "--split-files "
                    f"-O '{illumina_dir}' "
                    f"'{illumina_dir}/{sra_id}'"
                )
                ret2 = run_subprocess_cmd(fq2_cmd, True)

                if ret2 != 0 or not (r1.exists() and r2.exists()):
                    print(f"ERROR:\tFailed to produce FASTQ for {sra_id}.")
                    return None, None

            print(f"PASS:\tIllumina SRA processed: {r1}, {r2}")

        illu_raw_f_reads = str(r1)
        illu_raw_r_reads = str(r2)

    # ---- CASE B: Raw directory provided; combine/verify there ----
    elif pd.notna(illu_raw_dir) and (pd.isna(illu_raw_f_reads) or pd.isna(illu_raw_r_reads)):
        print(f"Process Illumina Raw Directory: {illu_raw_dir}...")
        combined_list = illumina_extract_and_check(str(illu_raw_dir), sample_id)
        if combined_list:
            print("PASS:\tSuccessfully processed Illumina Raw Directory.")
            illu_raw_f_reads, illu_raw_r_reads = combined_list[0], combined_list[1]
            # Update the 'SRA' field logically for downstream logging (won't rewrite CSV here)
        else:
            print(f"ERROR:\tFailed to process Illumina Raw Directory for {sample_id}")
            return None, None

    # If both explicit files provided in CSV, just trust them
    else:
        # Normalize strings if they exist
        if pd.notna(illu_raw_f_reads):
            illu_raw_f_reads = str(illu_raw_f_reads).strip()
        if pd.notna(illu_raw_r_reads):
            illu_raw_r_reads = str(illu_raw_r_reads).strip()

    # Final guard: we must have both R1 and R2 paths now
    if not illu_raw_f_reads or not illu_raw_r_reads or not os.path.exists(illu_raw_f_reads) or not os.path.exists(illu_raw_r_reads):
        print("ERROR:\tIllumina paired-end files are missing after preprocessing:")
        print(f"      R1: {illu_raw_f_reads}")
        print(f"      R2: {illu_raw_r_reads}")
        return None, None

    # Output filenames
    illu_dedup_f_reads = str(illumina_dir / f"{species_id}_illu_forward_dedup.fastq")
    illu_dedup_r_reads = str(illumina_dir / f"{species_id}_illu_reverse_dedup.fastq")

    # Short-circuit if already done
    if os.path.exists(illu_dedup_f_reads) and os.path.exists(illu_dedup_r_reads):
        print(f"SKIP:\tIllumina preprocessing already completed: {illu_dedup_f_reads}, {illu_dedup_r_reads}.")
        return illu_dedup_f_reads, illu_dedup_r_reads

    # ---------- FastQC on raw ----------
    fastqc_results_dir = str(illumina_dir / "fastqc_results")
    os.makedirs(fastqc_results_dir, exist_ok=True)
    run_subprocess_cmd(["fastqc", "-t", str(cpu_threads), "-o", fastqc_results_dir, illu_raw_f_reads, illu_raw_r_reads], False)

    # ---------- Trimmomatic ----------
    # Derive paired/unpaired output names from input
    def _swap_suffix(pth, old, new):
        if pth.endswith(old + ".gz"):
            return pth.replace(old + ".gz", new + ".gz")
        return pth.replace(old, new)

    def _nonempty(p):
        try:
            return os.path.isfile(p) and os.path.getsize(p) > 0
        except OSError:
            return False

    trimmo_f_pair    = _swap_suffix(illu_raw_f_reads, "_1.fastq", "_forward_paired.fastq")
    trimmo_r_pair    = _swap_suffix(illu_raw_r_reads, "_2.fastq", "_reverse_paired.fastq")
    trimmo_f_unpair  = _swap_suffix(illu_raw_f_reads, "_1.fastq", "_forward_unpaired.fastq")
    trimmo_r_unpair  = _swap_suffix(illu_raw_r_reads, "_2.fastq", "_reverse_unpaired.fastq")


    def _run_trimmomatic(t_threads, use_heap=False):
        truseq3_path = shutil.which("TruSeq3-PE.fa") or shutil.which("TruSeq3-PE")
        if not truseq3_path:
            print("WARN:\tTruSeq3-PE adapters file not found on PATH; Trimmomatic may fail.")
        base_cmd = [
            "trimmomatic", "PE",
            "-threads", str(t_threads),
            "-phred33",
            illu_raw_f_reads,           # Input R1
            illu_raw_r_reads,           # Input R2
            trimmo_f_pair,              # R1 paired output
            trimmo_f_unpair,            # R1 unpaired output
            trimmo_r_pair,              # R2 paired output
            trimmo_r_unpair,            # R2 unpaired output
            f"ILLUMINACLIP:{truseq3_path}:2:30:10:11",
            "HEADCROP:10",
            "CROP:145",
            "SLIDINGWINDOW:50:25",
            "MINLEN:125"
        ]
        if use_heap:
            # The conda 'trimmomatic' wrapper honors TRIMMOMATIC_OPTS for JVM flags.
            cmd = "TRIMMOMATIC_OPTS='-Xmx8g' " + " ".join(base_cmd)
            return run_subprocess_cmd(cmd, True)
        else:
            return run_subprocess_cmd(base_cmd, False)

    if os.path.exists(trimmo_f_pair) and os.path.exists(trimmo_r_pair) and _nonempty(trimmo_f_pair) and _nonempty(trimmo_r_pair):
        print(f"SKIP:\tTrimmomatic files exist: {trimmo_f_pair} & {trimmo_r_pair}.")
    else:
        # 1st attempt: as-is
        rc = _run_trimmomatic(cpu_threads, use_heap=False)

        # Check success: return code + non-empty outputs
        if rc != 0 or not (_nonempty(trimmo_f_pair) and _nonempty(trimmo_r_pair)):
            print("WARN:\tTrimmomatic failed or produced empty outputs; retrying with fewer threads and larger Java heap...")
            # Clean any empty outputs to avoid confusion
            for p in [trimmo_f_pair, trimmo_r_pair, trimmo_f_unpair, trimmo_r_unpair]:
                try:
                    if os.path.isfile(p) and os.path.getsize(p) == 0:
                        os.remove(p)
                except OSError:
                    pass

            # 2nd attempt: safer settings
            safer_threads = max(1, min(4, int(cpu_threads)))
            rc2 = _run_trimmomatic(safer_threads, use_heap=True)

            if rc2 != 0 or not (_nonempty(trimmo_f_pair) and _nonempty(trimmo_r_pair)):
                raise RuntimeError(
                    "Trimmomatic failed to generate non-empty paired outputs after retry.\n"
                    f"  R1 paired: {trimmo_f_pair} (exists={os.path.exists(trimmo_f_pair)}, size={os.path.getsize(trimmo_f_pair) if os.path.exists(trimmo_f_pair) else 'NA'})\n"
                    f"  R2 paired: {trimmo_r_pair} (exists={os.path.exists(trimmo_r_pair)}, size={os.path.getsize(trimmo_r_pair) if os.path.exists(trimmo_r_pair) else 'NA'})\n"
                    "Check Java memory/threads or adapter path."
                )

    # ---------- BBDuk ----------
    bbduk_f_map = str(Path(trimmo_f_pair).with_name(Path(trimmo_f_pair).name.replace("_forward_paired", "_forward_mapped")))
    bbduk_r_map = str(Path(trimmo_r_pair).with_name(Path(trimmo_r_pair).name.replace("_reverse_paired", "_reverse_mapped")))

    if os.path.exists(bbduk_f_map) and os.path.exists(bbduk_r_map):
        print(f"SKIP:\tbbduk mapped files exist: {bbduk_f_map} & {bbduk_r_map}.")
    else:
        run_subprocess_cmd([
            "bbduk.sh",
            f"in1={trimmo_f_pair}", f"in2={trimmo_r_pair}",
            f"out1={bbduk_f_map}",  f"out2={bbduk_r_map}",
            "ktrim=r", "k=23", "mink=11", "hdist=1",
            "tpe", "tbo", "qtrim=rl", "trimq=20"
        ], False)

    # ---------- Clumpify (dedupe) ----------
    if os.path.exists(illu_dedup_f_reads) and os.path.exists(illu_dedup_r_reads):
        print(f"SKIP:\tClumpify deduplicated files exist: {illu_dedup_f_reads} & {illu_dedup_r_reads}.")
    else:
        run_subprocess_cmd([
            "clumpify.sh",
            f"in={bbduk_f_map}", f"in2={bbduk_r_map}",
            f"out={illu_dedup_f_reads}", f"out2={illu_dedup_r_reads}",
            "dedupe"
        ], False)

    print(f"PASS:\tPreprocessed Raw Illumina Reads for {species_id}: {illu_dedup_f_reads}, {illu_dedup_r_reads}.")

    # ---------- FastQC on dedup ----------
    dedup_fastqc_results_dir = str(illumina_dir / "dedup_fastqc_results")
    os.makedirs(dedup_fastqc_results_dir, exist_ok=True)
    run_subprocess_cmd(["fastqc", "-t", str(cpu_threads), "-o", dedup_fastqc_results_dir, illu_dedup_f_reads, illu_dedup_r_reads], False)

    return illu_dedup_f_reads, illu_dedup_r_reads


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
    
    illu_dedup_f_reads, illu_dedup_r_reads = preprocess_illumina(sample_id, input_csv, output_dir, cpu_threads, ram_gb)
