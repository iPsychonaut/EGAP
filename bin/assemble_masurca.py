#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
assemble_masurca.py

This script runs MaSuRCA assembly with Illumina and optional long reads.

Created on Wed Aug 16 2023

Updated on Wed Sept 3 2025

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import os, sys, shutil, re, stat
import pandas as pd
from utilities import run_subprocess_cmd, get_current_row_data
from qc_assessment import qc_assessment
from pathlib import Path

# --------------------------------------------------------------
# Locate MaSuRCA CA folder
# --------------------------------------------------------------
def find_ca_folder(current_work_dir):
    """
    Return the MaSuRCA CA folder path under current_work_dir.
    Prefer the most-recent directory starting with 'CA', else '<cwd>/CA'.
    Ensures current_work_dir exists.
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
    """Extract insert size statistics from BBMerge histogram file.

    Parses the histogram to retrieve average insert size and standard deviation.

    Args:
        insert_size_histogram_txt (str): Path to the BBMerge histogram file.

    Returns:
        tuple: (average insert size, standard deviation).

    Raises:
        ValueError: If mean or standard deviation cannot be parsed.
    """
    print(f"Processing insert size histogram: {insert_size_histogram_txt}...")
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
    """Compute insert size statistics using BBMerge.

    Runs BBMerge to generate a histogram or parses an existing one to extract
    insert size statistics.

    Args:
        input_folder (str): Directory for BBMerge output files.
        reads_list (list): List of FASTQ file paths (paired-end or with extra read).

    Returns:
        tuple: (average insert size, standard deviation).
    """
    bbmap_out_path = f"{input_folder}/bbmap_data.fq"
    insert_size_histogram_txt = f"{input_folder}/insert_size_histogram.txt"
    print(f"DEBUG - bbmap_out_path - {bbmap_out_path}")
    avg_insert = 251
    std_dev = 30
    if os.path.isfile(insert_size_histogram_txt):
        print(f"SKIP:\tbbmap insert size histogram output already exists: {insert_size_histogram_txt}")
        avg_insert, std_dev = parse_bbmerge_output(insert_size_histogram_txt)
    else:
        print("Processing fastq files for bbmap stats...")
        bbmerge_path = shutil.which("bbmerge.sh") or shutil.which("bbmerge")
        if not bbmerge_path:
            raise FileNotFoundError("bbmerge not found in PATH")
        bbmerge_cmd = [bbmerge_path,
                       f"in1={reads_list[1]}",
                       f"in2={reads_list[2]}",
                       f"out={bbmap_out_path}",
                       f"ihist={input_folder}/insert_size_histogram.txt"]
        _ = run_subprocess_cmd(bbmerge_cmd, False)
        try:
            avg_insert, std_dev = parse_bbmerge_output(insert_size_histogram_txt)
        except:
            print("NOTE:\tUnable to parse avg_insert or std_dev, using default values: 251 and 30 respectively.")
            avg_insert = 251
            std_dev = 30
    return avg_insert, std_dev


# --------------------------------------------------------------
# Modify MaSuRCA script to skip gap closing
# --------------------------------------------------------------
def skip_gap_closing_section(assembly_sh_path):
    """Modify MaSuRCA assemble.sh script to bypass gap closing.

    Edits the script to skip the gap closing step and use the 9-terminator output.

    Args:
        assembly_sh_path (str): Path to the assemble.sh script.

    Returns:
        str: Path to the modified script.
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
def ensure_flye_vendor_shim():
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


def patch_environment_flye(env_sh_path: str) -> str:
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
            # Prefer Flye on PATH if available
            'if command -v flye >/dev/null 2>&1; then\n'
            '  export FLYE="$(command -v flye)"\n'
            'fi\n'
            # Provide PYTHON2 fallback and a local python2.7 shim if python2.7 is absent
            'if ! command -v python2.7 >/dev/null 2>&1; then\n'
            '  if command -v python2 >/dev/null 2>&1; then\n'
            '    export PYTHON2="$(command -v python2)"\n'
            '  elif command -v python3 >/dev/null 2>&1; then\n'
            '    export PYTHON2="$(command -v python3)"\n'
            '    # Drop a local shim named python2.7 that execs $PYTHON2 and put it first on PATH\n'
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

    # Soften exactly one “flye not found … exit 1” occurrence, if any
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


def patch_flye_check(assembly_sh_path: str) -> str:
    """
    Make assemble.sh accept Flye from PATH if the bundled Flye is absent,
    and only hard-fail when FLYE_ASSEMBLY=1.

    This keeps MaSuRCA’s original behavior if its bundled Flye exists,
    but avoids aborting when Flye is installed elsewhere (e.g., conda env).
    """
    with open(assembly_sh_path, "r") as f:
        txt = f.read()

    # 1) Replace the strict "Flye/bin ... exit 1" block with a PATH fallback.
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

    # 2) If that exact block wasn’t present (version drift), soften any generic "flye not found ... exit 1".
    if n == 0:
        pat_exit = re.compile(r'echo\s+"flye not found[^\n]*"\s*\n\s*exit\s+1', re.IGNORECASE)
        new_txt, _ = pat_exit.subn(
            'echo "WARN: flye not found; will try PATH and/or continue if FLYE_ASSEMBLY=0"\n'
            'if [ "${FLYE_ASSEMBLY:-0}" = "1" ]; then exit 1; fi',
            new_txt, count=1
        )

    # 3) Ensure FLYE gets resolved from PATH near the top if unset.
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
            new_txt = new_txt[:shebang_idx] + new_txt[shebang_idx:] + inject

    with open(assembly_sh_path, "w") as f:
        f.write(new_txt)
    return assembly_sh_path


# --------------------------------------------------------------
# Generate and run MaSuRCA assembly
# --------------------------------------------------------------
def masurca_config_gen(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """Generate and execute MaSuRCA configuration for genome assembly.

    Configures MaSuRCA for CABOG assembly with Illumina and optional long reads,
    runs the assembly, and performs quality control.

    Args:
        sample_id (str): Sample identifier.
        input_csv (str): Path to metadata CSV file.
        output_dir (str): Directory for output files.
        cpu_threads (int): Number of CPU threads to use.
        ram_gb (int): Available RAM in GB.

    Returns:
        str or None: Path to the gzipped MaSuRCA assembly FASTA, or None if assembly fails.
    """
    input_csv_abs = os.path.abspath(input_csv)
    print(f"DEBUG - input_csv_abs - {input_csv_abs}")

    output_dir_abs = str(Path(output_dir).resolve())
    print(f"DEBUG - output_dir_abs - {output_dir_abs}")

    input_df = pd.read_csv(input_csv_abs)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]

    # CSV fields
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

    species_dir = Path(output_dir_abs).resolve() / str(species_id)
    sample_dir = species_dir / str(sample_id)
    masurca_out_dir = (sample_dir / "masurca_assembly").resolve()
    masurca_out_dir.mkdir(parents=True, exist_ok=True)

    # Helper: absolute path for str/Path alike
    def _abs(p):
        return str(Path(p).resolve()) if p else None

    # Fill in implied paths from SRA when raw paths are empty
    if pd.notna(ont_sra) and pd.isna(ont_raw_reads):
        ont_raw_reads = _abs(species_dir / "ONT" / f"{ont_sra}.fastq")
    if pd.notna(illumina_sra) and pd.isna(illumina_f_raw_reads) and pd.isna(illumina_r_raw_reads):
        illumina_f_raw_reads = _abs(species_dir / "Illumina" / f"{illumina_sra}_1.fastq")
        illumina_r_raw_reads = _abs(species_dir / "Illumina" / f"{illumina_sra}_2.fastq")
    if pd.notna(pacbio_sra) and pd.isna(pacbio_raw_reads):
        pacbio_raw_reads = _abs(species_dir / "PacBio" / f"{pacbio_sra}.fastq")
    if pd.notna(ref_seq_gca) and pd.isna(ref_seq):
        ref_seq = _abs(species_dir / "RefSeq" / f"{species_id}_{ref_seq_gca}_RefSeq.fasta")

    # Preferred Illumina inputs = your dedup outputs (fall back to raw if missing)
    dedup_f = _abs(species_dir / "Illumina" / f"{species_id}_illu_forward_dedup.fastq")
    dedup_r = _abs(species_dir / "Illumina" / f"{species_id}_illu_reverse_dedup.fastq")
    use_dedup = dedup_f and dedup_r and os.path.exists(dedup_f) and os.path.exists(dedup_r) \
                and os.path.getsize(dedup_f) > 0 and os.path.getsize(dedup_r) > 0

    # Long-read selection (prefer prefiltered)
    highest_mean_qual_long_reads = None
    if pd.notna(ont_raw_reads):
        candidate = species_dir / "ONT" / f"{species_id}_ONT_highest_mean_qual_long_reads.fastq"
        highest_mean_qual_long_reads = _abs(candidate if candidate.exists() else ont_raw_reads)
    elif pd.notna(pacbio_raw_reads):
        candidate = species_dir / "PacBio" / f"{species_id}_PacBio_highest_mean_qual_long_reads.fastq"
        highest_mean_qual_long_reads = _abs(candidate if candidate.exists() else pacbio_raw_reads)

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
        print("SKIP:\tNo Illumina paired-end reads provided, required for MaSuRCA assembly")
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
        print(f"NOTE:\tUnable to parse input estimated size {est_size}, using default: 25000000")
        jf_size = 25_000_000

    # ---------- Build MaSuRCA config (absolute paths only) ----------
    config_lines = ["DATA\n"]
    if use_dedup:
        print("INFO:\tUsing deduplicated Illumina reads for MaSuRCA PE.")
        config_lines.append(f"PE= pe {int(avg_insert)} {int(std_dev)} {dedup_f} {dedup_r}\n")
    else:
        print("INFO:\tUsing raw Illumina reads for MaSuRCA PE.")
        config_lines.append(f"PE= pe {int(avg_insert)} {int(std_dev)} {illumina_f_raw_reads} {illumina_r_raw_reads}\n")

    if highest_mean_qual_long_reads:
        if pd.notna(ont_raw_reads):
            config_lines.append(f"NANOPORE={highest_mean_qual_long_reads}\n")
        elif pd.notna(pacbio_raw_reads):
            config_lines.append(f"PACBIO={highest_mean_qual_long_reads}\n")

    if pd.notna(ref_seq):
        config_lines.append(f"REFERENCE={ref_seq}\n")
    config_lines.append("END\n")

    config_lines += [
        "PARAMETERS\n",
        "GRAPH_KMER_SIZE=auto\n",
        f"USE_LINKING_MATES={0 if (pd.notna(ont_raw_reads) or pd.notna(pacbio_raw_reads)) else 1}\n",
        f"CLOSE_GAPS={0 if (pd.notna(ont_raw_reads) or pd.notna(pacbio_raw_reads)) else 1}\n",
        "MEGA_READS_ONE_PASS=0\n",
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
            print("ERROR:\tmasurca failed to generate assemble.sh")
            return None

        assemble_sh_path = masurca_out_dir / "assemble.sh"
        if not assemble_sh_path.exists():
            print("ERROR:\tassemble.sh not found after masurca step")
            return None
        
        # Ensure MaSuRCA sees Flye regardless of install location
        ensure_flye_vendor_shim()
        
        # Patch environment.sh to prefer PATH Flye and avoid fatal exit when FLYE_ASSEMBLY=0
        env_sh_path = masurca_out_dir / "environment.sh"
        patch_environment_flye(str(env_sh_path))
        
        # Patch assemble.sh to skip gap-closing (your existing helper)
        modified_assemble = skip_gap_closing_section(str(assemble_sh_path))
        
        # (Optional) print a one-liner sanity check
        real_flye = shutil.which("flye")
        print(f"DEBUG - flye on PATH - {real_flye}")
        
        # Run
        rc = run_subprocess_cmd(["bash", os.path.basename(modified_assemble)], False)
        
        # Try to salvage even if assemble.sh returned non-zero (MaSuRCA sometimes fails late
        # while the CA output is already produced, e.g. after synteny step).
        ca_dir = Path(find_ca_folder(str(masurca_out_dir))).resolve()
        primary_genome_scf = ca_dir / "primary.genome.scf.fasta"
        terminator_genome_scf = ca_dir / "9-terminator" / "genome.scf.fasta"
        
        if rc != 0:
            if terminator_genome_scf.exists() or primary_genome_scf.exists():
                print("WARN:\tassemble.sh exited non-zero, but CA output exists; continuing.")
            else:
                print("ERROR:\tAssembly run failed.")
                return None

        egap_masurca_assembly_path = masurca_out_dir / f"{sample_id}_masurca.fasta"
        if terminator_genome_scf.exists():
            shutil.copy(str(terminator_genome_scf), str(egap_masurca_assembly_path))
        elif primary_genome_scf.exists():
            shutil.copy(str(primary_genome_scf), str(egap_masurca_assembly_path))
        else:
            print("ERROR:\tNo assembly file found in CA output.")
            return None

        print(f"DEBUG: Final assembly path: {egap_masurca_assembly_path}")
        os.chdir(prev_cwd)
        egap_masurca_assembly_path, masurca_stats_list, _ = qc_assessment(
            "masurca", input_csv_abs, sample_id, output_dir_abs, cpu_threads, ram_gb
        )
        return str(egap_masurca_assembly_path)
    finally:
        os.chdir(prev_cwd)


if __name__ == "__main__":
    # Provide a usage message if the correct number of arguments is not supplied.
    if len(sys.argv) != 6:
        print("Usage: python3 assemble_masurca.py <sample_id> <input_csv> "
            "<output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)
    
    egap_masurca_assembly_path = masurca_config_gen(sys.argv[1],       # sample_id
                                                    sys.argv[2],       # input_csv
                                                    sys.argv[3],       # output_dir
                                                    str(sys.argv[4]),  # cpu_threads
                                                    str(sys.argv[5]))  # ram_gb
