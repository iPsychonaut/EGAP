#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
assemble_masurca.py

This script runs MaSuRCA assembly with Illumina and optional long reads.

Stage:
    Hybrid Assembly (MaSuRCA)

Created on Wed Aug 16 2023

Updated on 2026-04-16

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os
import sys
import shutil
import re
import stat
from typing import Optional

import pandas as pd
from pathlib import Path
from utilities import run_subprocess_cmd, log_print, initialize_logging_environment, load_sample_context
from qc_assessment import qc_assessment
from file_manager import remove_file, remove_dir
from monitor_assembly import AssemblyMonitor
from estimate_runtime import log_estimate_for
from record_provenance import record_file

# --------------------------------------------------------------
# Locate MaSuRCA CA folder
# --------------------------------------------------------------
def find_ca_folder(current_work_dir):
    """Return the MaSuRCA CA folder path under *current_work_dir*.

    Prefers the most-recently modified directory whose name starts with
    ``'CA'``; falls back to ``<current_work_dir>/CA`` if none is found.
    Creates *current_work_dir* if it does not exist.

    Parameters
    ----------
    current_work_dir : str
        Directory to scan for CA output subdirectories.

    Returns
    -------
    str
        Absolute path to the CA folder.
    """
    os.makedirs(current_work_dir, exist_ok=True)
    candidates = []
    with os.scandir(current_work_dir) as it:
        for e in it:
            if e.is_dir() and os.path.basename(e.path).startswith("CA"):
                try:
                    mtime = os.path.getmtime(e.path)
                except OSError:
                    mtime = 0
                candidates.append((mtime, e.path))
    if candidates:
        candidates.sort(reverse=True)
        return candidates[0][1]
    return os.path.join(current_work_dir, "CA")


# --------------------------------------------------------------
# Parse BBMerge insert size statistics
# --------------------------------------------------------------
def parse_bbmerge_output(insert_size_histogram_txt):
    """Extract insert size statistics from a BBMerge histogram file.

    Parameters
    ----------
    insert_size_histogram_txt : str
        Path to the BBMerge insert-size histogram text file.

    Returns
    -------
    tuple of (float, float)
        ``(avg_insert, std_dev)`` — average insert size and standard
        deviation, both rounded to the nearest integer.

    Raises
    ------
    ValueError
        If ``#Mean`` or ``#STDev`` lines cannot be found or parsed.
    """
    log_print(f"NOTE:\tProcessing insert size histogram: {insert_size_histogram_txt}...")
    avg_insert = None
    std_dev = None
    with open(insert_size_histogram_txt, "r") as file:
        for line in file:
            if "#Mean\t" in line:
                avg_insert = round(float(line.replace("#Mean\t", "").replace("\n", "")), 0)
            if "#STDev\t" in line:
                std_dev = round(float(line.replace("#STDev\t", "").replace("\n", "")), 0)
    if avg_insert is None or std_dev is None:
        raise ValueError("Could not find average insert size and/or standard deviation in the output.")
    return avg_insert, std_dev


# --------------------------------------------------------------
# Compute insert size statistics with BBMerge
# --------------------------------------------------------------
def bbmap_stats(input_folder, reads_list, cpu_threads):
    """Compute Illumina insert size statistics using BBMerge.

    If an insert-size histogram already exists in *input_folder* it is
    parsed directly; otherwise BBMerge is run.  Paired reads are
    interleaved via ``reformat.sh`` before calling BBMerge to avoid the
    PairStreamer byte-offset desync bug.

    Parameters
    ----------
    input_folder : str
        Directory used for BBMerge intermediate and output files.
    reads_list : list of str
        File path list where index 1 is the forward FASTQ and index 2 is
        the reverse FASTQ (matches the MaSuRCA reads_list convention).
    cpu_threads : int
        Number of threads available (passed to ``reformat.sh``).

    Returns
    -------
    tuple of (float, float)
        ``(avg_insert, std_dev)`` — average insert size and standard
        deviation.
    """
    bbmap_out_path = f"{input_folder}/bbmap_data.fq"
    insert_size_histogram_txt = f"{input_folder}/insert_size_histogram.txt"
    log_print(f"DEBUG - bbmap_out_path - {bbmap_out_path}")
    avg_insert = 251
    std_dev = 30
    if os.path.isfile(insert_size_histogram_txt):
        log_print(f"SKIP:\tbbmap insert size histogram output already exists: {insert_size_histogram_txt}")
        avg_insert, std_dev = parse_bbmerge_output(insert_size_histogram_txt)
    else:
        log_print("NOTE:\tProcessing fastq files for bbmap stats...")
        # Interleave F/R reads first — BBMerge's parallel PairStreamer
        # splits paired files by byte offset, which desyncs when reads
        # have variable lengths after quality trimming. A single
        # interleaved file avoids this entirely.
        reformat_path = shutil.which("reformat.sh")
        interleaved_path = f"{input_folder}/bbmerge_interleaved.fq"
        if reformat_path and not (os.path.isfile(interleaved_path) and os.path.getsize(interleaved_path) > 0):
            log_print("NOTE:\tInterleaving paired reads for BBMerge...")
            _ = run_subprocess_cmd([
                reformat_path,
                f"in1={reads_list[1]}", f"in2={reads_list[2]}",
                f"out={interleaved_path}"
            ], False)

        bbmerge_path = shutil.which("bbmerge.sh") or shutil.which("bbmerge")
        if not bbmerge_path:
            raise FileNotFoundError("bbmerge not found in PATH")

        # Use interleaved input if available, otherwise fall back to paired
        if os.path.isfile(interleaved_path):
            bbmerge_cmd = [bbmerge_path,
                           f"in={interleaved_path}",
                           f"out={bbmap_out_path}",
                           f"ihist={input_folder}/insert_size_histogram.txt"]
        else:
            bbmerge_cmd = [bbmerge_path,
                           f"in1={reads_list[1]}",
                           f"in2={reads_list[2]}",
                           f"out={bbmap_out_path}",
                           f"ihist={input_folder}/insert_size_histogram.txt"]
        _ = run_subprocess_cmd(bbmerge_cmd, False)
        try:
            avg_insert, std_dev = parse_bbmerge_output(insert_size_histogram_txt)
        except:
            log_print("NOTE:\tUnable to parse avg_insert or std_dev, using default values: 251 and 30 respectively.")
            avg_insert = 251
            std_dev = 30
    return avg_insert, std_dev


# --------------------------------------------------------------
# Modify MaSuRCA script to skip gap closing
# --------------------------------------------------------------
def skip_gap_closing_section(assembly_sh_path):
    """Modify a MaSuRCA ``assemble.sh`` script to bypass gap closing.

    Writes a new ``assemble_skip_gap.sh`` alongside the original that
    skips the gap-closing step and keeps the 9-terminator scaffolds.

    Parameters
    ----------
    assembly_sh_path : str
        Path to the original MaSuRCA ``assemble.sh`` script.

    Returns
    -------
    str
        Path to the modified ``assemble_skip_gap.sh`` script.
    """
    with open(assembly_sh_path, "r") as f_in:
        original_script = f_in.read()
    pattern = re.compile(r"(if \[ -s \$CA_DIR/9-terminator/genome\.scf\.fasta \];then)"
                         r"(.*?)"
                         r"(else\s+fail 'Assembly stopped or failed, see \$CA_DIR\.log'\nfi)",
                         re.DOTALL)
    replacement_snippet = (r"\1\n"
                           r"  # Force the final terminator to remain '9-terminator'\n"
                           r"  log \"Skipping gap closing step; using 9-terminator as final.\"\n"
                           r"  TERMINATOR='9-terminator'\n"
                           r"\3")
    modified_text = re.sub(pattern, replacement_snippet, original_script)
    modified_script_path = os.path.join(os.path.dirname(assembly_sh_path), "assemble_skip_gap.sh")
    with open(modified_script_path, "w") as f_out:
        f_out.write(modified_text)
    return modified_script_path


# --------------------------------------------------------------
# NEW: Relax Flye check to allow PATH-based Flye
# --------------------------------------------------------------
def _masurca_vendor_flye_ensure_shim():
    """
    Make MaSuRCA happy regardless of where Conda installed Flye by ensuring
    $CONDA_PREFIX/Flye/bin/flye exists (the path MaSuRCA probes).
    Creates a symlink to the real flye if possible, else a tiny wrapper script.
    Safe on Linux and macOS. No effect if already present.
    """
    masurca_bin = shutil.which("masurca")
    flye_bin = shutil.which("flye")
    if not masurca_bin or not flye_bin:
        return  # nothing we can do

    # $PREFIX/bin/masurca -> $PREFIX
    env_prefix = os.path.realpath(os.path.join(os.path.dirname(masurca_bin), ".."))
    vendor_dir = os.path.join(env_prefix, "Flye", "bin")
    vendor_exe = os.path.join(vendor_dir, "flye")

    # If already present and executable, we’re done
    if os.path.exists(vendor_exe) and os.access(vendor_exe, os.X_OK):
        return

    os.makedirs(vendor_dir, exist_ok=True)
    try:
        if os.path.exists(vendor_exe):
            os.remove(vendor_exe)
        os.symlink(flye_bin, vendor_exe)
        made = "symlink"
    except Exception:
        # Symlinks may be blocked; fall back to a tiny wrapper
        with open(vendor_exe, "w") as w:
            w.write(f'#!/usr/bin/env bash\nexec "{flye_bin}" "$@"\n')
        os.chmod(vendor_exe, 0o755)
        made = "wrapper"

    print(f"INFO:\tInstalled Flye vendor {made} for MaSuRCA: {vendor_exe} -> {flye_bin}")


def _masurca_vendor_flye_patch_environment(env_sh_path: str) -> str:
    """
    Prepend a small block to environment.sh so it prefers a PATH flye if present,
    provides a python2.7 fallback shim (to python3) if python2.7 is missing,
    and only fails hard when FLYE_ASSEMBLY=1. Leaves other failure paths intact.
    """
    if not os.path.exists(env_sh_path):
        return env_sh_path

    with open(env_sh_path, "r") as f:
        original = f.read()

    marker = "# EGAP_FLYE_PATH_PATCH"
    if marker not in original:
        # NOTE: environment.sh is sourced from masurca_out_dir, so $(pwd) is that dir.
        prepend = (
            f"{marker}\n"
            'if command -v flye >/dev/null 2>&1; then\n'
            '  export FLYE="$(command -v flye)"\n'
            'fi\n'
            'if ! command -v python2.7 >/dev/null 2>&1; then\n'
            '  if command -v python2 >/dev/null 2>&1; then\n'
            '    export PYTHON2="$(command -v python2)"\n'
            '  elif command -v python3 >/dev/null 2>&1; then\n'
            '    export PYTHON2="$(command -v python3)"\n'
            '    if [ -w "$(pwd)" ]; then\n'
            '      printf \'#!/usr/bin/env bash\\nexec "%s" "$@"\\n\' "$PYTHON2" > "$(pwd)/python2.7"\n'
            '      chmod +x "$(pwd)/python2.7"\n'
            '      export PATH="$(pwd):$PATH"\n'
            '    fi\n'
            '  fi\n'
            'fi\n'
            '# Only hard-fail for missing Flye when FLYE_ASSEMBLY=1\n'
            f"{marker}\n"
        )
        patched = prepend + original
    else:
        patched = original

    patched = re.sub(
        r'(flye not found[^\n]*\n)(\s*exit\s+1)',
        r'\1if [ "${FLYE_ASSEMBLY:-0}" = "1" ]; then exit 1; fi',
        patched, count=1, flags=re.IGNORECASE
    )

    with open(env_sh_path, "w") as f:
        f.write(patched)

    st = os.stat(env_sh_path)
    os.chmod(env_sh_path, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return env_sh_path


def _masurca_vendor_flye_patch_check(assembly_sh_path: str) -> str:
    """
    Make assemble.sh accept Flye from PATH if the bundled Flye is absent,
    and only hard-fail when FLYE_ASSEMBLY=1.
    """
    with open(assembly_sh_path, "r") as f:
        txt = f.read()

    pat_block = re.compile(
        r'if\s+\[\s*-x\s+"\$MASURCA/../Flye/bin/flye"\s*\]\s*;?\s*then\s*'
        r'FLYE=.*?\n'
        r'\s*else\s*\n'
        r'\s*echo\s+"flye not found[^"\n]*"\s*\n'
        r'\s*exit\s+1\s*\n'
        r'\s*fi',
        re.DOTALL
    )
    repl_block = (
        'if [ -x "$MASURCA/../Flye/bin/flye" ]; then\n'
        '  FLYE="$MASURCA/../Flye/bin/flye"\n'
        'elif command -v flye >/dev/null 2>&1; then\n'
        '  FLYE="$(command -v flye)"\n'
        'else\n'
        '  if [ "${FLYE_ASSEMBLY:-0}" = "1" ]; then\n'
        '    echo "ERROR: flye required (FLYE_ASSEMBLY=1) but not found on PATH"; exit 1\n'
        '  else\n'
        '    echo "WARN: flye not found; continuing because FLYE_ASSEMBLY=${FLYE_ASSEMBLY:-0}"\n'
        '  fi\n'
        'fi'
    )
    new_txt, n = pat_block.subn(repl_block, txt, count=1)

    if n == 0:
        pat_exit = re.compile(r'echo\s+"flye not found[^\n]*"\s*\n\s*exit\s+1', re.IGNORECASE)
        new_txt, _ = pat_exit.subn(
            'echo "WARN: flye not found; will try PATH and/or continue if FLYE_ASSEMBLY=0"\n'
            'if [ "${FLYE_ASSEMBLY:-0}" = "1" ]; then exit 1; fi',
            new_txt, count=1
        )

    if "EGAP FLYE PATH PATCH" not in new_txt:
        shebang_idx = new_txt.find("\n")
        if shebang_idx != -1:
            inject = (
                '\n# --- EGAP FLYE PATH PATCH ---\n'
                'if [ -z "${FLYE:-}" ]; then\n'
                '  if [ -x "$MASURCA/../Flye/bin/flye" ]; then\n'
                '    FLYE="$MASURCA/../Flye/bin/flye"\n'
                '  elif command -v flye >/dev/null 2>&1; then\n'
                '    FLYE="$(command -v flye)"\n'
                '  fi\n'
                'fi\n'
                '# --- END EGAP FLYE PATH PATCH ---\n'
            )
            new_txt = new_txt[:shebang_idx] + inject + new_txt[shebang_idx:]

    with open(assembly_sh_path, "w") as f:
        f.write(new_txt)
    return assembly_sh_path


def clear_stale_overlap_store(masurca_out_dir) -> None:
    """Remove an orphaned CABOG ``*.ovlStore.BUILDING`` left by an interrupted run.

    CABOG's ``overlapStoreBuild`` writes the overlap store into
    ``<name>.ovlStore.BUILDING`` and only renames it to the final
    ``<name>.ovlStore`` once the build completes. If that build is killed
    mid-flight (for example the controlling terminal sends SIGHUP), the partial
    ``.BUILDING`` directory survives. Every later ``runCA`` then aborts with
    "is a valid overlap store, will not overwrite", so MaSuRCA can never finish
    on a re-run into the same output directory. This preflight clears such an
    orphan (and its stale ``*.ovlStore.list``) when there is no completed sibling
    store, so the next attempt rebuilds the overlap store cleanly. A completed
    ``*.ovlStore`` is never touched.
    """
    base = Path(masurca_out_dir)
    if not base.exists():
        return
    for building in base.rglob("*.ovlStore.BUILDING"):
        if not building.is_dir():
            continue
        completed = building.with_suffix("")  # drop the trailing ".BUILDING"
        if completed.exists():
            # A finished store sits alongside the orphan; leave both alone.
            continue
        try:
            shutil.rmtree(building)
            log_print(f"NOTE:\tCleared stale CABOG overlap store: {building}")
        except OSError as e:
            log_print(f"WARN:\tCould not remove stale overlap store {building}: {e}")
        stale_list = Path(f"{completed}.list")
        if stale_list.exists():
            try:
                stale_list.unlink()
                log_print(f"NOTE:\tCleared stale overlap store list: {stale_list}")
            except OSError as e:
                log_print(f"WARN:\tCould not remove {stale_list}: {e}")


# --------------------------------------------------------------
# Generate and run MaSuRCA assembly
# --------------------------------------------------------------
def masurca_config_gen(
    sample_id: str,
    input_tsv: str,
    output_dir: str,
    cpu_threads: int,
    ram_gb: int,
) -> Optional[str]:
    """Generate and execute MaSuRCA configuration for genome assembly.

    Configures MaSuRCA for CABOG assembly with Illumina paired-end reads
    and optional ONT/PacBio long reads, runs the assembly, and performs
    quality control via ``qc_assessment``.

    Parameters
    ----------
    sample_id : str
        Sample identifier used to look up the row in *input_tsv*.
    input_tsv : str
        Path to the metadata TSV file.
    output_dir : str
        Root output directory; per-species subdirectories are created here.
    cpu_threads : int
        Number of CPU threads to pass to MaSuRCA and BBMerge.
    ram_gb : int
        Available RAM in GB (passed to quality-control tools).

    Returns
    -------
    str or None
        Absolute path to the final MaSuRCA assembly FASTA, or ``None``
        if the assembly could not be completed.
    """
    ctx = load_sample_context(sample_id, input_tsv, output_dir, cpu_threads, ram_gb)
    print(f"DEBUG - input_tsv - {ctx.input_tsv}")
    print(f"DEBUG - output_dir - {ctx.output_dir}")
    current_series = ctx.current_series

    # TSV fields
    illumina_sra = current_series["ILLUMINA_SRA"]
    illumina_f_raw_reads = current_series["ILLUMINA_RAW_F_READS"]
    illumina_r_raw_reads = current_series["ILLUMINA_RAW_R_READS"]
    ont_sra = current_series["ONT_SRA"]
    ont_raw_reads = current_series["ONT_RAW_READS"]
    pacbio_sra = current_series["PACBIO_SRA"]
    pacbio_raw_reads = current_series["PACBIO_RAW_READS"]
    ref_seq_gca = current_series["REF_SEQ_GCA"]
    ref_seq = current_series["REF_SEQ"]
    species_id = current_series["SPECIES_ID"]
    est_size = current_series["EST_SIZE"]

    species_dir = Path(ctx.output_dir).resolve() / str(species_id)
    sample_dir = species_dir / str(sample_id)
    masurca_out_dir = (sample_dir / "masurca_assembly").resolve()
    masurca_out_dir.mkdir(parents=True, exist_ok=True)

    # ---------- FAST SKIP if final output exists (re-run QC) ----------
    expected_final = masurca_out_dir / f"{sample_id}_masurca.fasta"
    if expected_final.exists() and expected_final.stat().st_size > 0:
        log_print(f"SKIP:\tFinal MaSuRCA assembly already present: {expected_final}")
        # Re-run QC on the existing assembly/output structure
        egap_masurca_assembly_path, masurca_stats_list, _ = qc_assessment(
            "masurca", ctx.input_tsv, sample_id, ctx.output_dir, cpu_threads, ram_gb
        )
        return str(egap_masurca_assembly_path)

    # Helper: absolute path for str/Path alike
    def to_abs(p):
        if p is None or pd.isna(p):
            return None
        return str(Path(str(p)).resolve())

    # Fill in implied paths from SRA when raw paths are empty
    if pd.notna(ont_sra) and pd.isna(ont_raw_reads):
        ont_raw_reads = to_abs(species_dir / "ONT" / f"{ont_sra}.fastq")
    if pd.notna(illumina_sra) and pd.isna(illumina_f_raw_reads) and pd.isna(illumina_r_raw_reads):
        illumina_f_raw_reads = to_abs(species_dir / "Illumina" / f"{illumina_sra}_1.fastq")
        illumina_r_raw_reads = to_abs(species_dir / "Illumina" / f"{illumina_sra}_2.fastq")
    if pd.notna(pacbio_sra) and pd.isna(pacbio_raw_reads):
        pacbio_raw_reads = to_abs(species_dir / "PacBio" / f"{pacbio_sra}.fastq")
    if pd.notna(ref_seq_gca) and pd.isna(ref_seq):
        ref_seq = to_abs(species_dir / "RefSeq" / f"{species_id}_{ref_seq_gca}_RefSeq.fasta")

    # Preferred Illumina inputs = your dedup outputs (fall back to raw if missing)
    dedup_f = to_abs(species_dir / "Illumina" / f"{species_id}_illu_forward_dedup.fastq")
    dedup_r = to_abs(species_dir / "Illumina" / f"{species_id}_illu_reverse_dedup.fastq")
    use_dedup = dedup_f and dedup_r and os.path.exists(dedup_f) and os.path.exists(dedup_r) \
                and os.path.getsize(dedup_f) > 0 and os.path.getsize(dedup_r) > 0

    # Long-read selection (prefer prefiltered)
    highest_mean_qual_long_reads = None
    if pd.notna(ont_raw_reads):
        candidate = species_dir / "ONT" / f"{species_id}_ONT_highest_mean_qual_long_reads.fastq"
        highest_mean_qual_long_reads = to_abs(candidate if candidate.exists() else ont_raw_reads)
    elif pd.notna(pacbio_raw_reads):
        candidate = species_dir / "PacBio" / f"{species_id}_PacBio_highest_mean_qual_long_reads.fastq"
        highest_mean_qual_long_reads = to_abs(candidate if candidate.exists() else pacbio_raw_reads)

    print(f"DEBUG - illumina_sra - {illumina_sra}")
    print(f"DEBUG - illumina_f_raw_reads - {illumina_f_raw_reads}")
    print(f"DEBUG - illumina_r_raw_reads - {illumina_r_raw_reads}")
    print(f"DEBUG - ont_sra - {ont_sra}")
    print(f"DEBUG - ont_raw_reads - {ont_raw_reads}")
    print(f"DEBUG - pacbio_sra - {pacbio_sra}")
    print(f"DEBUG - pacbio_raw_reads - {pacbio_raw_reads}")
    print(f"DEBUG - ref_seq_gca - {ref_seq_gca}")
    print(f"DEBUG - ref_seq - {ref_seq}")
    print(f"DEBUG - species_id - {species_id}")
    print(f"DEBUG - est_size - {est_size}")
    print(f"DEBUG - masurca_out_dir - {masurca_out_dir}")
    print(f"DEBUG - dedup_f - {dedup_f}")
    print(f"DEBUG - dedup_r - {dedup_r}")
    print(f"DEBUG - highest_mean_qual_long_reads    - {highest_mean_qual_long_reads}")

    # Require Illumina for MaSuRCA
    have_raw_illumina = (isinstance(illumina_f_raw_reads, str) and os.path.exists(illumina_f_raw_reads)) \
                        and (isinstance(illumina_r_raw_reads, str) and os.path.exists(illumina_r_raw_reads))
    if not (use_dedup or have_raw_illumina):
        log_print("SKIP:\tNo Illumina paired-end reads provided, required for MaSuRCA assembly")
        return None

    # Insert size (prefer computing from available Illumina)
    if use_dedup:
        avg_insert, std_dev = bbmap_stats(str(masurca_out_dir),
                                          [ont_raw_reads, dedup_f, dedup_r, pacbio_raw_reads],
                                          cpu_threads)
    elif have_raw_illumina:
        avg_insert, std_dev = bbmap_stats(str(masurca_out_dir),
                                          [ont_raw_reads, illumina_f_raw_reads, illumina_r_raw_reads, pacbio_raw_reads],
                                          cpu_threads)
    else:
        avg_insert, std_dev = 251, 30

    # Jellyfish hash size
    m = re.match(r"^(\d+(?:\.\d+)?)(\D+)$", str(est_size))
    if m:
        est_size_numb, est_size_mult = m.group(1), m.group(2)
        multipliers = {'m': 10**6, 'g': 10**9}
        jf_size = int(float(est_size_numb) * multipliers.get(est_size_mult.lower(), 25_000_000))
    else:
        log_print(f"NOTE:\tUnable to parse input estimated size {est_size}, using default: 25000000")
        jf_size = 25_000_000

    # ---------- Build MaSuRCA config (absolute paths only) ----------
    config_lines = ["DATA\n"]
    if use_dedup:
        log_print("NOTE:\tUsing deduplicated Illumina reads for MaSuRCA PE.")
        config_lines.append(f"PE= pe {int(avg_insert)} {int(std_dev)} {dedup_f} {dedup_r}\n")
    else:
        log_print("NOTE:\tUsing raw Illumina reads for MaSuRCA PE.")
        config_lines.append(f"PE= pe {int(avg_insert)} {int(std_dev)} {illumina_f_raw_reads} {illumina_r_raw_reads}\n")

    if highest_mean_qual_long_reads:
        if pd.notna(ont_raw_reads):
            config_lines.append(f"NANOPORE={highest_mean_qual_long_reads}\n")
        elif pd.notna(pacbio_raw_reads):
            config_lines.append(f"PACBIO={highest_mean_qual_long_reads}\n")

    if pd.notna(ref_seq):
        config_lines.append(f"REFERENCE={ref_seq}\n")
    config_lines.append("END\n")

    # When long reads are present, MaSuRCA builds "mega-reads" via create_mega_reads.
    # The default two-pass mode (MEGA_READS_ONE_PASS=0) runs an iterative refinement
    # that can hang pathologically on noisy/over-long ONT input (a single core pegged
    # at 100% with no output progress for hours). One-pass mode skips that refinement:
    # it is faster and far less prone to the stall, at the cost of a slightly less
    # contiguous mega-read set. EGAP compares all assemblers and keeps the best, so a
    # completing-but-slightly-rougher MaSuRCA is strictly better than one that hangs.
    has_long_reads = pd.notna(ont_raw_reads) or pd.notna(pacbio_raw_reads)
    config_lines += [
        "PARAMETERS\n",
        "GRAPH_KMER_SIZE=auto\n",
        f"USE_LINKING_MATES={0 if has_long_reads else 1}\n",
        f"CLOSE_GAPS={0 if has_long_reads else 1}\n",
        f"MEGA_READS_ONE_PASS={1 if has_long_reads else 0}\n",
        "LIMIT_JUMP_COVERAGE=300\n",
        "CA_PARAMETERS=cgwErrorRate=0.15\n",
        f"NUM_THREADS={cpu_threads}\n",
        f"JF_SIZE={jf_size}\n",
        "SOAP_ASSEMBLY=0\n",
        "FLYE_ASSEMBLY=0\n",
        "END\n",
    ]

    config_path = masurca_out_dir / "masurca_config_file.txt"
    with open(config_path, "w") as fh:
        fh.writelines(config_lines)
    print(f"DEBUG: Wrote config to {config_path}")

    # Always run MaSuRCA from its output dir (run_subprocess_cmd has no cwd)
    prev_cwd = os.getcwd()
    os.chdir(str(masurca_out_dir))
    try:
        rc = run_subprocess_cmd(["masurca", "masurca_config_file.txt"], False)
        if rc != 0:
            log_print("ERROR:\tmasurca failed to generate assemble.sh")
            return None

        assemble_sh_path = masurca_out_dir / "assemble.sh"
        if not assemble_sh_path.exists():
            log_print("ERROR:\tassemble.sh not found after masurca step")
            return None
        
        # Ensure MaSuRCA sees Flye regardless of install location
        _masurca_vendor_flye_ensure_shim()
        
        # Patch environment.sh to prefer PATH Flye and avoid fatal exit when FLYE_ASSEMBLY=0
        env_sh_path = masurca_out_dir / "environment.sh"
        _masurca_vendor_flye_patch_environment(str(env_sh_path))

        # Patch assemble.sh's Flye check to accept PATH Flye and only hard-fail when FLYE_ASSEMBLY=1.
        # Must run before skip_gap_closing_section, which reads assemble.sh and writes assemble_skip_gap.sh.
        _masurca_vendor_flye_patch_check(str(assemble_sh_path))

        # Patch assemble.sh to skip gap-closing (your existing helper)
        modified_assemble = skip_gap_closing_section(str(assemble_sh_path))
        
        # (Optional) print a one-liner sanity check
        real_flye = shutil.which("flye")
        print(f"DEBUG - flye on PATH - {real_flye}")
        
        # Print a pre-flight runtime estimate for this assembler given the
        # available resources and read volume.
        log_estimate_for("masurca", sample_id, ctx.input_tsv, ctx.output_dir, cpu_threads, ram_gb)

        # Preflight: an interrupted prior run (e.g. SIGHUP during the overlap
        # build) can leave a partial CABOG *.ovlStore.BUILDING that blocks every
        # retry with "will not overwrite". Clear any such orphan with no completed
        # sibling before assemble.sh runs, so the overlap store rebuilds cleanly.
        clear_stale_overlap_store(masurca_out_dir)

        # Run under a non-lethal runtime monitor. It reports elapsed time and
        # output throughput and warns if growth flatlines (the create_mega_reads
        # hang signature: CPU busy, zero progress) but never kills the assembler
        # -- the operator decides whether to keep waiting.
        with AssemblyMonitor(masurca_out_dir, label=f"{sample_id} masurca"):
            rc = run_subprocess_cmd(["bash", os.path.basename(modified_assemble)], False)

        # Try to salvage even if assemble.sh returned non-zero
        ca_dir = Path(find_ca_folder(str(masurca_out_dir))).resolve()
        primary_genome_scf = ca_dir / "primary.genome.scf.fasta"
        terminator_genome_scf = ca_dir / "9-terminator" / "genome.scf.fasta"
        
        if rc != 0:
            if terminator_genome_scf.exists() or primary_genome_scf.exists():
                log_print("WARN:\tassemble.sh exited non-zero, but CA output exists; continuing.")
            else:
                log_print("ERROR:\tAssembly run failed.")
                return None

        egap_masurca_assembly_path = masurca_out_dir / f"{sample_id}_masurca.fasta"
        if terminator_genome_scf.exists():
            shutil.copy(str(terminator_genome_scf), str(egap_masurca_assembly_path))
        elif primary_genome_scf.exists():
            shutil.copy(str(primary_genome_scf), str(egap_masurca_assembly_path))
        else:
            log_print("ERROR:\tNo assembly file found in CA output.")
            return None
        record_file("MaSuRCA assembly", str(egap_masurca_assembly_path))

        print(f"DEBUG: Final assembly path: {egap_masurca_assembly_path}")
        os.chdir(prev_cwd)
        egap_masurca_assembly_path, masurca_stats_list, _ = qc_assessment(
            "masurca", ctx.input_tsv, sample_id, ctx.output_dir, cpu_threads, ram_gb
        )

        # --- Cleanup MaSuRCA intermediates once final assembly is confirmed ---
        if egap_masurca_assembly_path and os.path.exists(str(egap_masurca_assembly_path)):
            # CA/ CABOG tree can reach ~11 GB
            remove_dir(str(ca_dir))
            # work1/ overlap/unitig scratch space
            remove_dir(str(masurca_out_dir / "work1"))
            # Large working FASTQ/FASTA files produced during config generation
            for _wf in [
                "pe.cor.fa",
                "pe.linking.fa",
                "guillaumeKUnitigsAtLeast32bases_all.fasta",
                "guillaumeKUnitigsAtLeast32bases_all.jump.fasta",
                "bbmerge_interleaved.fq",
                "bbmap_data.fq",
            ]:
                remove_file(str(masurca_out_dir / _wf))

        return str(egap_masurca_assembly_path)
    finally:
        os.chdir(prev_cwd)


if __name__ == "__main__":
    # Provide a usage message if the correct number of arguments is not supplied.
    if len(sys.argv) != 6:
        print("Usage: python3 assemble_masurca.py <sample_id> <input_tsv> "
            "<output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)

    initialize_logging_environment(sys.argv[3], sys.argv[1])

    egap_masurca_assembly_path = masurca_config_gen(sys.argv[1],       # sample_id
                                                    sys.argv[2],       # input_tsv
                                                    sys.argv[3],       # output_dir
                                                    str(sys.argv[4]),  # cpu_threads
                                                    str(sys.argv[5]))  # ram_gb
