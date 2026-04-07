#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
preprocess_illumina.py

This script preprocesses Illumina reads with FastQC, Trimmomatic, BBDuk, and Clumpify,
then optionally decontaminates with Kraken2.
Reads input from a CSV file, processes all Illumina datasets, and updates the CSV
with declumpified FASTQ file paths in ${params.output_dir}/${sample_prefix}/Illumina/.

Created on Wed Aug 16 2023

Updated on Fri Apr 04 2026

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import os
import sys
import shutil
import glob
import pandas as pd
from utilities import run_subprocess_cmd, get_current_row_data, md5_check, initialize_logging_environment, log_print
from decontaminate_reads import (
    get_kraken2_db,
    get_kraken_keep_domains,
    build_taxid_domain_map,
    parse_kraken_reads,
    ALL_KRAKEN_DOMAINS,
)


# --------------------------------------------------------------
# Helper: treat pd.NA / NaN / "None" string as missing
# --------------------------------------------------------------
def clean(val):
    """Return *val* or ``None`` if it represents a missing value.

    Parameters
    ----------
    val : object
        Value to inspect; may be any type including ``pd.NA``.

    Returns
    -------
    object or None
        *val* unchanged, or ``None`` when *val* is NA/NaN or the
        literal string ``"none"`` (case-insensitive).
    """
    if pd.isna(val) or (isinstance(val, str) and val.strip().lower() == "none"):
        return None
    return val


# --------------------------------------------------------------
# Helper: swap a filename suffix
# --------------------------------------------------------------
def swap_suffix(pth, old, new):
    """Replace *old* suffix with *new* in *pth*, preserving optional .gz.

    Parameters
    ----------
    pth : str
        File path whose suffix should be swapped.
    old : str
        Suffix to replace (without .gz).
    new : str
        Replacement suffix (without .gz).

    Returns
    -------
    str
        Path with the suffix replaced.
    """
    if pth.endswith(old + ".gz"):
        return pth.replace(old + ".gz", new + ".gz")
    return pth.replace(old, new)


# --------------------------------------------------------------
# Helper: check a path is a non-empty file
# --------------------------------------------------------------
def nonempty(p):
    """Return True if *p* is an existing, non-empty regular file.

    Parameters
    ----------
    p : str
        Path to check.

    Returns
    -------
    bool
        ``True`` when the file exists and has size > 0, ``False`` otherwise.
    """
    try:
        return os.path.isfile(p) and os.path.getsize(p) > 0
    except OSError:
        return False


# --------------------------------------------------------------
# Helper: locate TruSeq3-PE adapter file
# --------------------------------------------------------------
def find_truseq_adapters():
    """Search common locations for the TruSeq3-PE.fa adapter file.

    Searches (in order): the active conda environment, the directory
    adjacent to the trimmomatic binary, and common system paths.

    Returns
    -------
    str or None
        Absolute path to TruSeq3-PE.fa, or ``None`` if not found.
    """
    # 1) Check CONDA_PREFIX (active conda env)
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    if conda_prefix:
        candidates = glob.glob(
            os.path.join(conda_prefix, "share", "trimmomatic*", "adapters", "TruSeq3-PE.fa")
        )
        if candidates:
            return candidates[0]
    # 2) Check next to the trimmomatic binary
    trimmo_bin = shutil.which("trimmomatic")
    if trimmo_bin:
        trimmo_root = os.path.dirname(os.path.dirname(os.path.realpath(trimmo_bin)))
        candidates = glob.glob(
            os.path.join(trimmo_root, "share", "trimmomatic*", "adapters", "TruSeq3-PE.fa")
        )
        if candidates:
            return candidates[0]
    # 3) Fallback: search common paths
    for search_root in ["/usr/share", "/usr/local/share", os.path.expanduser("~/miniconda3")]:
        candidates = glob.glob(os.path.join(search_root, "**", "TruSeq3-PE.fa"), recursive=True)
        if candidates:
            return candidates[0]
    return None


# --------------------------------------------------------------
# Combine and verify Illumina reads
# --------------------------------------------------------------
def illumina_extract_and_check(folder_name, SAMPLE_ID):
    """Combine paired-end Illumina reads after MD5 verification.

    Verifies MD5 checksums, then concatenates all forward and reverse FASTQ
    files found in *folder_name*.

    Parameters
    ----------
    folder_name : str
        Directory containing Illumina FASTQ files and ``MD5.txt``.
    SAMPLE_ID : str
        Sample identifier used to name the combined output files.

    Returns
    -------
    list of str or None
        Two-element list ``[combined_1_file, combined_2_file]``, or ``None``
        if the files could not be produced.
    """
    log_print(f"Running MD5 Checksum Analysis on Raw Illumina FASTQ files in {folder_name}...")
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
                log_print(f"ERROR:\tNo paired Illumina files found in {folder_name}")
                return None
            fwd_cat_cmd = f"cat {' '.join(raw_1_list)} > {combined_1_file}"
            _ = run_subprocess_cmd(fwd_cat_cmd, shell_check=True)
            rev_cat_cmd = f"cat {' '.join(raw_2_list)} > {combined_2_file}"
            _ = run_subprocess_cmd(rev_cat_cmd, shell_check=True)
        else:
            log_print(f"SKIP:\tCombined FASTQ files already exist: {combined_list[0]}; {combined_list[1]}.")
    else:
        log_print(f"SKIP:\tGzipped Combined FASTQ files already exist: {combined_list[0]}; {combined_list[1]}.")

    return combined_list if os.path.exists(combined_list[0]) and os.path.exists(combined_list[1]) else None


# --------------------------------------------------------------
# Kraken2 paired-end Illumina decontamination
# --------------------------------------------------------------
def decontaminate_illumina_paired(f_reads, r_reads, kraken_dir, kraken2_db,
                                   keep_domains, cpu_threads, species_id):
    """Run Kraken2 in paired-end mode and filter Illumina reads by domain.

    Parameters
    ----------
    f_reads : str
        Path to the forward deduped FASTQ file.
    r_reads : str
        Path to the reverse deduped FASTQ file.
    kraken_dir : str
        Directory for Kraken2 intermediate output files.
    kraken2_db : str
        Path to the Kraken2 database directory.
    keep_domains : set of str
        Domain labels (from ``ALL_KRAKEN_DOMAINS``) whose reads should be
        retained.
    cpu_threads : int
        Number of threads to pass to Kraken2.
    species_id : str
        Species identifier used in log messages.

    Returns
    -------
    tuple of (str, str) or (None, None)
        Paths ``(f_reads, r_reads)`` to the decontaminated FASTQ files,
        or ``(None, None)`` on failure.
    """
    os.makedirs(kraken_dir, exist_ok=True)
    label = "Illumina"
    kraken_out    = os.path.join(kraken_dir, f"{label}_kraken2.out")
    kraken_report = os.path.join(kraken_dir, f"{label}_kraken2_report.txt")

    # -- Run Kraken2 in paired-end mode --
    if os.path.exists(kraken_out) and os.path.exists(kraken_report):
        log_print(f"SKIP:\tKraken2 output already exists for {label}: {kraken_out}")
    else:
        kraken_cmd = (
            f"kraken2 "
            f"--db {kraken2_db} "
            f"--threads {cpu_threads} "
            f"--paired "
            f"--output {kraken_out} "
            f"--report {kraken_report} "
            f"--report-zero-counts "
            f"{f_reads} {r_reads}"
        )
        rc = run_subprocess_cmd(kraken_cmd, shell_check=True)
        if rc != 0 or not os.path.exists(kraken_out) or not os.path.exists(kraken_report):
            log_print(f"ERROR:\tKraken2 failed for Illumina paired-end (rc={rc}).")
            return None, None

    # -- Build taxid -> domain map from the report --
    taxid_domain = build_taxid_domain_map(kraken_report)

    # -- Parse per-read classifications --
    read_taxid = parse_kraken_reads(kraken_out)

    # -- Determine which read IDs to keep --
    keep_read_ids = set()
    domain_counts = {d: 0 for d in ALL_KRAKEN_DOMAINS}
    for read_id, taxid in read_taxid.items():
        domain = taxid_domain.get(taxid, "other")
        domain_counts[domain] = domain_counts.get(domain, 0) + 1
        if domain in keep_domains:
            keep_read_ids.add(read_id)

    # -- Classification summary --
    total = max(len(read_taxid), 1)
    log_print(f"\n[Kraken2] Illumina paired-end classification summary:")
    for domain in sorted(ALL_KRAKEN_DOMAINS):
        count  = domain_counts.get(domain, 0)
        action = "KEEP" if domain in keep_domains else "REMOVE"
        log_print(f"  {domain:<14}: {count:>8}  ({100*count/total:5.1f}%)  [{action}]")

    pct_kept = 100.0 * len(keep_read_ids) / total
    log_print(f"  Total kept : {len(keep_read_ids):>8}  ({pct_kept:.1f}%)")

    if pct_kept < 50.0:
        log_print(f"WARN:\tFewer than 50% of Illumina reads were kept "
                  f"({pct_kept:.1f}%). Check kingdom assignment and Kraken2 database.")

    # -- Back up originals, then write filtered reads --
    # os.rename is an atomic move — avoids duplicating multi-GB files.
    f_pre = f_reads.replace(".fastq", "_pre_decontam.fastq")
    r_pre = r_reads.replace(".fastq", "_pre_decontam.fastq")
    if not os.path.exists(f_pre):
        os.rename(f_reads, f_pre)
        log_print(f"NOTE:\tOriginal forward reads preserved as: {f_pre}")
    if not os.path.exists(r_pre):
        os.rename(r_reads, r_pre)
        log_print(f"NOTE:\tOriginal reverse reads preserved as: {r_pre}")

    # Filter both files in sync — read paired FASTQ records together
    f_kept = 0
    with open(f_pre, "r", buffering=1 << 20) as fh_f_in, \
         open(r_pre, "r", buffering=1 << 20) as fh_r_in, \
         open(f_reads, "w", buffering=1 << 20) as fh_f_out, \
         open(r_reads, "w", buffering=1 << 20) as fh_r_out:
        while True:
            # Read one 4-line FASTQ record from each file
            f_lines = [fh_f_in.readline() for _ in range(4)]
            r_lines = [fh_r_in.readline() for _ in range(4)]

            # EOF check
            if not f_lines[0] or not r_lines[0]:
                break

            # Extract read ID (portion before first space, strip leading @)
            f_rid = f_lines[0].split()[0].lstrip("@")

            if f_rid in keep_read_ids:
                fh_f_out.writelines(f_lines)
                fh_r_out.writelines(r_lines)
                f_kept += 1

    if f_kept == 0:
        log_print(f"ERROR:\tKraken2 decontamination removed ALL Illumina reads. Restoring originals.")
        shutil.copy(f_pre, f_reads)
        shutil.copy(r_pre, r_reads)
        return f_reads, r_reads

    log_print(f"PASS:\tIllumina Kraken2 decontamination complete: "
              f"{f_kept} read pairs kept -> {f_reads}, {r_reads}")
    return f_reads, r_reads


# --------------------------------------------------------------
# Preprocess Illumina sequencing reads
# --------------------------------------------------------------
def preprocess_illumina(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """Preprocess Illumina reads with FastQC, Trimmomatic, BBDuk, Clumpify, and Kraken2.

    Creates ``<output_dir>/<SPECIES_ID>/Illumina/`` if absent, then runs the
    full preprocessing chain.  Returns the paths to the final deduplicated
    FASTQ files that downstream assemblers expect.

    Parameters
    ----------
    sample_id : str
        Sample identifier used to look up the row in *input_csv*.
    input_csv : str
        Path to the metadata CSV file.
    output_dir : str
        Root output directory; per-species subdirectories are created here.
    cpu_threads : int
        Number of CPU threads available for tools such as FastQC and
        Trimmomatic.
    ram_gb : int
        RAM (GB) available; passed to tools that accept a memory limit.

    Returns
    -------
    tuple of (str, str) or (None, None)
        ``(illu_dedup_f_reads, illu_dedup_r_reads)`` — absolute paths to the
        deduplicated forward and reverse FASTQ files — or ``(None, None)``
        when no Illumina data are available or preprocessing fails.

    Raises
    ------
    RuntimeError
        If Trimmomatic fails to produce non-empty paired outputs after retry,
        or if ``reformat.sh trd=t`` fails to fix Trimmomatic headers.
    """
    from pathlib import Path

    log_print(f"Preprocessing Illumina reads for {sample_id.split('-')[0]}...")

    # Parse the CSV and retrieve relevant row data
    input_df = pd.read_csv(input_csv)
    log_print(f"DEBUG - input_df - {input_df}")

    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]  # Convert to Series (single row)

    log_print(f"DEBUG - current_series - {current_series}")

    # Identify read paths/reference info from CSV
    illu_sra = clean(current_series["ILLUMINA_SRA"])
    illu_raw_f_reads = clean(current_series["ILLUMINA_RAW_F_READS"])
    illu_raw_r_reads = clean(current_series["ILLUMINA_RAW_R_READS"])
    illu_raw_dir = clean(current_series["ILLUMINA_RAW_DIR"])
    species_id = current_series["SPECIES_ID"]

    log_print(f"DEBUG - illu_sra - {illu_sra}")
    log_print(f"DEBUG - illu_raw_f_reads - {illu_raw_f_reads}")
    log_print(f"DEBUG - illu_raw_r_reads - {illu_raw_r_reads}")
    log_print(f"DEBUG - illu_raw_dir - {illu_raw_dir}")

    # Early skip if nothing Illumina-like is provided
    if pd.isna(illu_sra) and pd.isna(illu_raw_f_reads) and pd.isna(illu_raw_r_reads) and pd.isna(illu_raw_dir):
        log_print(f"SKIP:\tSample does not include Illumina Reads: {sample_id}.")
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
            log_print(f"Downloading Illumina SRA: {sra_id}...")

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
                log_print("WARN:\tfasterq-dump failed or files not found; trying fastq-dump fallback...")
                fq2_cmd = (
                    env_prefix +
                    "fastq-dump "
                    "--split-files "
                    f"-O '{illumina_dir}' "
                    f"'{illumina_dir}/{sra_id}'"
                )
                ret2 = run_subprocess_cmd(fq2_cmd, True)

                if ret2 != 0 or not (r1.exists() and r2.exists()):
                    log_print(f"ERROR:\tFailed to produce FASTQ for {sra_id}.")
                    return None, None

            log_print(f"PASS:\tIllumina SRA processed: {r1}, {r2}")

        illu_raw_f_reads = str(r1)
        illu_raw_r_reads = str(r2)

    # ---- CASE B: Raw directory provided; combine/verify there ----
    elif pd.notna(illu_raw_dir) and (pd.isna(illu_raw_f_reads) or pd.isna(illu_raw_r_reads)):
        log_print(f"Process Illumina Raw Directory: {illu_raw_dir}...")
        combined_list = illumina_extract_and_check(str(illu_raw_dir), sample_id)
        if combined_list:
            log_print("PASS:\tSuccessfully processed Illumina Raw Directory.")
            illu_raw_f_reads, illu_raw_r_reads = combined_list[0], combined_list[1]
            # Update the 'SRA' field logically for downstream logging (won't rewrite CSV here)
        else:
            log_print(f"ERROR:\tFailed to process Illumina Raw Directory for {sample_id}")
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
        log_print("ERROR:\tIllumina paired-end files are missing after preprocessing:")
        log_print(f"      R1: {illu_raw_f_reads}")
        log_print(f"      R2: {illu_raw_r_reads}")
        return None, None

    # Output filenames
    illu_dedup_f_reads = str(illumina_dir / f"{species_id}_illu_forward_dedup.fastq")
    illu_dedup_r_reads = str(illumina_dir / f"{species_id}_illu_reverse_dedup.fastq")

    # Short-circuit if already done
    if os.path.exists(illu_dedup_f_reads) and os.path.exists(illu_dedup_r_reads):
        log_print(f"SKIP:\tIllumina preprocessing already completed: {illu_dedup_f_reads}, {illu_dedup_r_reads}.")
        return illu_dedup_f_reads, illu_dedup_r_reads

    # ---------- FastQC on raw ----------
    fastqc_results_dir = str(illumina_dir / "fastqc_results")
    os.makedirs(fastqc_results_dir, exist_ok=True)
    # Derive expected FastQC output names to check for skip
    _r1_base = os.path.splitext(os.path.basename(illu_raw_f_reads))[0]
    _r2_base = os.path.splitext(os.path.basename(illu_raw_r_reads))[0]
    _fqc_r1 = os.path.join(fastqc_results_dir, f"{_r1_base}_fastqc.html")
    _fqc_r2 = os.path.join(fastqc_results_dir, f"{_r2_base}_fastqc.html")
    if os.path.exists(_fqc_r1) and os.path.exists(_fqc_r2):
        log_print(f"SKIP:\tFastQC results already exist: {_fqc_r1} & {_fqc_r2}")
    else:
        run_subprocess_cmd(["fastqc", "-t", str(cpu_threads), "-o", fastqc_results_dir, illu_raw_f_reads, illu_raw_r_reads], False)

    # ---------- Trimmomatic ----------
    # Derive paired/unpaired output names from input
    trimmo_f_pair    = swap_suffix(illu_raw_f_reads, "_1.fastq", "_forward_paired.fastq")
    trimmo_r_pair    = swap_suffix(illu_raw_r_reads, "_2.fastq", "_reverse_paired.fastq")
    trimmo_f_unpair  = swap_suffix(illu_raw_f_reads, "_1.fastq", "_forward_unpaired.fastq")
    trimmo_r_unpair  = swap_suffix(illu_raw_r_reads, "_2.fastq", "_reverse_unpaired.fastq")


    def run_trimmomatic(t_threads, use_heap=False):
        truseq3_path = find_truseq_adapters()
        if not truseq3_path:
            log_print("WARN:\tTruSeq3-PE adapters file not found; Trimmomatic may fail.")
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

    if os.path.exists(trimmo_f_pair) and os.path.exists(trimmo_r_pair) and nonempty(trimmo_f_pair) and nonempty(trimmo_r_pair):
        log_print(f"SKIP:\tTrimmomatic files exist: {trimmo_f_pair} & {trimmo_r_pair}.")
    else:
        # 1st attempt: as-is
        rc = run_trimmomatic(cpu_threads, use_heap=False)

        # Check success: return code + non-empty outputs
        if rc != 0 or not (nonempty(trimmo_f_pair) and nonempty(trimmo_r_pair)):
            log_print("WARN:\tTrimmomatic failed or produced empty outputs; retrying with fewer threads and larger Java heap...")
            # Clean any empty outputs to avoid confusion
            for p in [trimmo_f_pair, trimmo_r_pair, trimmo_f_unpair, trimmo_r_unpair]:
                try:
                    if os.path.isfile(p) and os.path.getsize(p) == 0:
                        os.remove(p)
                except OSError:
                    pass

            # 2nd attempt: safer settings
            safer_threads = max(1, min(4, int(cpu_threads)))
            rc2 = run_trimmomatic(safer_threads, use_heap=True)

            if rc2 != 0 or not (nonempty(trimmo_f_pair) and nonempty(trimmo_r_pair)):
                raise RuntimeError(
                    "Trimmomatic failed to generate non-empty paired outputs after retry.\n"
                    f"  R1 paired: {trimmo_f_pair} (exists={os.path.exists(trimmo_f_pair)}, size={os.path.getsize(trimmo_f_pair) if os.path.exists(trimmo_f_pair) else 'NA'})\n"
                    f"  R2 paired: {trimmo_r_pair} (exists={os.path.exists(trimmo_r_pair)}, size={os.path.getsize(trimmo_r_pair) if os.path.exists(trimmo_r_pair) else 'NA'})\n"
                    "Check Java memory/threads or adapter path."
                )

    # ---------- Reformat: strip stale length= from Trimmomatic headers ----------
    # Trimmomatic HEADCROP/CROP trims bases but keeps "length=151" in headers.
    # BBDuk's parallel FASTQ reader uses that annotation for byte-offset seeking
    # BEFORE any processing flags like trd take effect, causing "missing plus"
    # assertions. reformat.sh trd=t is BBTools' native fix for this.
    reformat_f = str(Path(trimmo_f_pair).with_name(
        Path(trimmo_f_pair).name.replace("_forward_paired", "_forward_reformat")))
    reformat_r = str(Path(trimmo_r_pair).with_name(
        Path(trimmo_r_pair).name.replace("_reverse_paired", "_reverse_reformat")))

    if nonempty(reformat_f) and nonempty(reformat_r):
        log_print(f"SKIP:\tReformatted FASTQ files already exist: {reformat_f} & {reformat_r}.")
    else:
        log_print("NOTE:\tStripping stale length= annotations from Trimmomatic output via reformat.sh...")
        rc_f = run_subprocess_cmd([
            "reformat.sh",
            f"in={trimmo_f_pair}", f"out={reformat_f}", "trd=t"
        ], False)
        rc_r = run_subprocess_cmd([
            "reformat.sh",
            f"in={trimmo_r_pair}", f"out={reformat_r}", "trd=t"
        ], False)
        if rc_f != 0 or rc_r != 0 or not (nonempty(reformat_f) and nonempty(reformat_r)):
            for p in [reformat_f, reformat_r]:
                try:
                    if os.path.isfile(p) and os.path.getsize(p) == 0:
                        os.remove(p)
                except OSError:
                    pass
            raise RuntimeError(
                "reformat.sh trd=t failed to fix Trimmomatic headers.\n"
                f"  Forward: {reformat_f} (rc={rc_f})\n"
                f"  Reverse: {reformat_r} (rc={rc_r})"
            )

    # ---------- BBDuk ----------
    # ref=adapters : load built-in adapter kmers
    bbduk_f_map = str(Path(trimmo_f_pair).with_name(Path(trimmo_f_pair).name.replace("_forward_paired", "_forward_mapped")))
    bbduk_r_map = str(Path(trimmo_r_pair).with_name(Path(trimmo_r_pair).name.replace("_reverse_paired", "_reverse_mapped")))

    if nonempty(bbduk_f_map) and nonempty(bbduk_r_map):
        log_print(f"SKIP:\tbbduk mapped files exist: {bbduk_f_map} & {bbduk_r_map}.")
    else:
        run_subprocess_cmd([
            "bbduk.sh",
            f"in1={reformat_f}", f"in2={reformat_r}",
            f"out1={bbduk_f_map}",  f"out2={bbduk_r_map}",
            "ref=adapters",
            "ktrim=r", "k=23", "mink=11", "hdist=1",
            "tpe", "tbo", "qtrim=rl", "trimq=20"
        ], False)

    # ---------- Clumpify (dedupe) ----------
    if os.path.exists(illu_dedup_f_reads) and os.path.exists(illu_dedup_r_reads):
        log_print(f"SKIP:\tClumpify deduplicated files exist: {illu_dedup_f_reads} & {illu_dedup_r_reads}.")
    else:
        run_subprocess_cmd([
            "clumpify.sh",
            f"in={bbduk_f_map}", f"in2={bbduk_r_map}",
            f"out={illu_dedup_f_reads}", f"out2={illu_dedup_r_reads}",
            "dedupe"
        ], False)

    log_print(f"PASS:\tPreprocessed Raw Illumina Reads for {species_id}: {illu_dedup_f_reads}, {illu_dedup_r_reads}.")

    # ---------- FastQC on dedup ----------
    dedup_fastqc_results_dir = str(illumina_dir / "dedup_fastqc_results")
    os.makedirs(dedup_fastqc_results_dir, exist_ok=True)
    _dd_f_base = os.path.splitext(os.path.basename(illu_dedup_f_reads))[0]
    _dd_r_base = os.path.splitext(os.path.basename(illu_dedup_r_reads))[0]
    _dd_fqc_f = os.path.join(dedup_fastqc_results_dir, f"{_dd_f_base}_fastqc.html")
    _dd_fqc_r = os.path.join(dedup_fastqc_results_dir, f"{_dd_r_base}_fastqc.html")
    if os.path.exists(_dd_fqc_f) and os.path.exists(_dd_fqc_r):
        log_print(f"SKIP:\tDedup FastQC results already exist: {_dd_fqc_f} & {_dd_fqc_r}")
    else:
        run_subprocess_cmd(["fastqc", "-t", str(cpu_threads), "-o", dedup_fastqc_results_dir, illu_dedup_f_reads, illu_dedup_r_reads], False)

    # ---------- Kraken2 Decontamination ----------
    # Resolve Kraken2 DB — if not configured, skip gracefully (non-fatal)
    kraken2_db = get_kraken2_db(current_series)
    if kraken2_db is None:
        log_print("WARN:\tNo Kraken2 database found (set KRAKEN2_DB env var or "
                  "add a KRAKEN2_DB column to the CSV). Skipping Illumina decontamination.")
    elif shutil.which("kraken2") is None:
        log_print("WARN:\tkraken2 is not on PATH. Skipping Illumina decontamination.")
    else:
        kingdom_id = current_series.get("ORGANISM_KINGDOM", None)
        keep_domains = get_kraken_keep_domains(kingdom_id)
        log_print(f"[Kraken2] Illumina decontamination profile for '{kingdom_id}':")
        log_print(f"  KEEP domains: {', '.join(sorted(keep_domains))}")

        kraken_dir = str(species_dir / "kraken2_reads")
        decontam_f, decontam_r = decontaminate_illumina_paired(
            illu_dedup_f_reads, illu_dedup_r_reads,
            kraken_dir, kraken2_db, keep_domains,
            cpu_threads, species_id
        )
        if decontam_f is None:
            log_print("WARN:\tKraken2 Illumina decontamination failed. "
                      "Continuing with unfiltered deduped reads.")

    return illu_dedup_f_reads, illu_dedup_r_reads


if __name__ == "__main__":
    # Check argument count before any log_print calls
    if len(sys.argv) != 6:
        print(f"ERROR: Expected 5 arguments (plus script name), got {len(sys.argv)-1}: {sys.argv[1:]}",
              file=sys.stderr)
        print("Usage: python3 preprocess_illumina.py <sample_id> <input_csv> <output_dir> <cpu_threads> <ram_gb>",
              file=sys.stderr)
        sys.exit(1)

    initialize_logging_environment(sys.argv[3], sys.argv[1])

    # Log raw sys.argv immediately
    log_print(f"DEBUG: Raw sys.argv = {sys.argv}")
    log_print(f"DEBUG: Length of sys.argv = {len(sys.argv)}")

    # Log each argument
    for i, arg in enumerate(sys.argv):
        log_print(f"DEBUG: sys.argv[{i}] = '{arg}'")

    sample_id = sys.argv[1]
    input_csv = sys.argv[2]
    output_dir = sys.argv[3]
    cpu_threads = int(sys.argv[4])
    ram_gb = int(sys.argv[5]) if sys.argv[5] != " " else 8

    log_print(f"DEBUG: Parsed sample_id = '{sample_id}' {type(sample_id)}")
    log_print(f"DEBUG: Parsed input_csv = '{input_csv}' {type(input_csv)}")
    log_print(f"DEBUG: Parsed output_dir = '{output_dir}' {type(output_dir)}")
    log_print(f"DEBUG: Parsed cpu_threads = '{sys.argv[4]}' {sys.argv[4]} (converted to {cpu_threads}) {type(cpu_threads)}")
    log_print(f"DEBUG: Parsed ram_gb = '{sys.argv[5]}' {sys.argv[5]} (converted to {ram_gb}) {type(ram_gb)}")

    illu_dedup_f_reads, illu_dedup_r_reads = preprocess_illumina(sample_id, input_csv, output_dir, cpu_threads, ram_gb)
