#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EGAP.py

Orchestrates the Entheome Genome Assembly Pipeline (EGAP) by executing a series of
preprocessing, assembly, comparison, polishing, curation, and quality assessment steps
for genomic data. Processes multiple samples from a CSV input, utilizing specified
CPU threads and RAM resources.

Created on Wed Aug 16 2023

Updated on Wed Feb 25 2026

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""

import argparse
import sys
import time
import subprocess
import re
import os
from pathlib import Path
import pandas as pd
from datetime import datetime


# --------------------------------------------------------------
# Tee: mirror stdout to a log file in real time
# --------------------------------------------------------------
class _Tee:
    """Write to *primary* stream and strip-ANSI copy to *secondary* (log file)."""
    _ANSI_ESCAPE = re.compile(r'\x1b\[[0-9;]*[mK]')

    def __init__(self, primary, secondary):
        self._primary = primary
        self._secondary = secondary

    def write(self, data):
        self._primary.write(data)
        clean = self._ANSI_ESCAPE.sub('', data)
        self._secondary.write(clean)
        self._secondary.flush()

    def flush(self):
        self._primary.flush()
        self._secondary.flush()

    def fileno(self):
        return self._primary.fileno()


# Lines containing any of these substrings are suppressed from console output.
# They are non-fatal noise from numpy C-extension API version mismatches that
# occur when packages in the same env were compiled against different numpy
# minor versions (e.g. numpy=1.19.5 in EGAP_env vs a package built on 1.20+).
_NOISE_SUBSTRINGS = (
    "module compiled against API version",
    "RuntimeError: module compiled against",
)


def run_filtered(cmd: list) -> int:
    """Run *cmd* as a subprocess, stream output to stdout in real time,
    and silently drop lines matching _NOISE_SUBSTRINGS.

    Returns the subprocess exit code.
    """
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    for line in proc.stdout:
        stripped = line.rstrip("\n")
        if not any(ns in stripped for ns in _NOISE_SUBSTRINGS):
            print(stripped, flush=True)
    proc.wait()
    return proc.returncode


# --------------------------------------------------------------
# Establish Global Pipeline Settings (i.e. command variables)
# --------------------------------------------------------------
VERSION = "3.4.2"

TRIMMOMATIC_SETTINGS = {
    "mode": "-PE",
    "phred version": "-phred33",
    "illuminaclip adapter": "/opt/conda/envs/EGAP_env/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa",
    "fastaWithAdaptersEtc": 2,
    "seed mismatches": 30,
    "palindrome clip threshold": 10,
    "simple clip threshold": 11,
    "HEADCROP": 0,
    "CROP": 145,
    "SLIDINGWINDOW": "50:25",
    "MINLEN": 125,
}

BBDUK_SETTINGS = {
    "ktrim": "-r",
    "k": 23,
    "mink": 11,
    "hdist": 1,
    "trimpairsevenly": "-tpe",
    "trimbyoverlap": "-tbo",
    "qtrim": "-rl",
    "trimq": 20,
}

CLUMPIFY_SETTINGS = {
    "remove duplicate reads": "-dedupe (haploid / unspecified PLOIDY only; skipped for diploid+, PLOIDY>1, to preserve allelic coverage)",
}

FILTLONG_SETTINGS = {
    "min_length": 1000,
    "min_mean_q": 8,
    "keep_percent": 90,
    "coverage": 75,
    "target_bases": "estimated size (bp) * 75 (coverage)",
}

MASURCA_SETTINGS = {
    "configuration file sections": "DATA, PARAMETERS",
    "graph kmer size": "auto",
    "use linking mates": "(0 if Hybrid assembly, else 1)",
    "close gaps": "(0 if Hybrid assembly, else 1)",
    "mega reads one pass": 0,
    "limit jump coverage": 300,
    "ca parameters": "cgwErrorRate=0.15",
    "jellyfish hash size": "based on the estimated size of the genome provided",
    "soap assembly": "0 (disabled to force CABOG assembly)",
    "flye assembly": "0 (disabled to force CABOG assembly)",
}

SPADES_SETTINGS = {
    "high-cov. isolate & multi-cell": "--isolate",
    "coverage cutoff ": "auto",
}

FLYE_SETTINGS = {
    "estimated genome size": "est_size",
    "number of polishing iterations": 3,
    "collapse alternative haplotypes": "--keep-haplotypes",
}

HIFIASM_SETTINGS = {
    "settings": "default",
}

RACON_SETTINGS = {
    "settings": "default",
}

PILON_SETTINGS = {
    "change file generation": "-changes",
    "vcf file generation": "-vcf",
    "tracks file generation": "-tracks",
    "largest chunksize limit": 5000000,
    "fix list": "indels, local, snps",
}

RAGTAG_SETTINGS = {
    "concatenate unplaced": "-C",
    "madd suffix to unplaced": "-u",
}

TGSGAPCLOSER_SETTINGS = {
    "do not error correct": "-ne",
}

ABYSSSEALER_SETTINGS = {
    "pseudoreads length used": 400,
    "bloom filter size": "500M",
}

FASTQC_SETTINGS = {
    "settings": "default",
}

NANOPLOT_SETTINGS = {
    "bivariate plots": "kde, dot",
    "show logarithmic lengths scaling": "--loglength",
    "N50 mark in read length histogram": "--N50",
    "log messages to terminal": "--verbose",
}

BUSCO_SETTINGS = {
    "mode": "genome",
    "force overwrite": "-f",
}

COMPLEASM_SETTINGS = {
    "settings": "default",
} 

QUAST_SETTINGS = {
    "settings": "default",
}

KRAKEN2_SETTINGS = {
    "stage": "pre-assembly (reads)",
    "targets": "ONT and PacBio highest-quality long reads",
    "database": "resolved from KRAKEN2_DB env var or CSV column",
    "keep (bacteria)": "bacteria, unclassified",
    "keep (archaea)": "archaea, unclassified",
    "keep (flora/funga/fauna)": "eukarya, unclassified",
    "always kept": "unclassified",
    "missing db behaviour": "WARN and skip (non-fatal)",
}

TIARA_SETTINGS = {
    "stage": "post-assembly (contigs)",
    "targets": "final curated assembly",
    "keep (bacteria)": "bacteria, prokarya, unknown",
    "keep (archaea)": "archaea, prokarya, unknown",
    "keep (flora/funga/fauna)": "eukarya, organelle, unknown",
    "organelle note": "kept for eukaryotes (own mitochondria/plastid); removed for prokaryotes",
    "always kept": "unknown",
    "missing tiara behaviour": "WARN and skip (non-fatal)",
}


def get_pipeline_settings(current_moment, ram_gb, cpu_threads, input_csv, output_dir):
    """Return all pipeline settings as a structured dict for reuse (TUI, logs, etc.)."""
    return {
        "input-setup settings": {
            "started at": str(current_moment),
            "config files": "N/A (running python version)",
            "container": "N/A (running python version)",
            "RAM GB": str(ram_gb),
            "CPU threads": str(cpu_threads),
            "input csv": str(input_csv),
            "output to": str(output_dir),
        },
        "trimmomatic settings": TRIMMOMATIC_SETTINGS,
        "bbduk settings": BBDUK_SETTINGS,
        "clumpify settings": CLUMPIFY_SETTINGS,
        "filtlong settings": FILTLONG_SETTINGS,
        "masurca settings": MASURCA_SETTINGS,
        "spades settings": SPADES_SETTINGS,
        "flye settings": FLYE_SETTINGS,
        "hifiasm settings": HIFIASM_SETTINGS,
        "racon settings": RACON_SETTINGS,
        "pilon settings": PILON_SETTINGS,
        "ragtag settings": RAGTAG_SETTINGS,
        "tgs-gapcloser settings": TGSGAPCLOSER_SETTINGS,
        "abyss-sealer settings": ABYSSSEALER_SETTINGS,
        "fastqc settings": FASTQC_SETTINGS,
        "nanoplot settings": NANOPLOT_SETTINGS,
        "busco settings": BUSCO_SETTINGS,
        "compleasm settings": COMPLEASM_SETTINGS,
        "quast settings": QUAST_SETTINGS,
        "kraken2 settings": KRAKEN2_SETTINGS,
        "tiara settings": TIARA_SETTINGS,
    }


def print_nested_dict(data, indent=0,
                      section_color="\033[92m",
                      key_color="\033[94m",
                      reset="\033[0m"):
    """
    Recursively print nested dictionaries in a readable, indented format.

    - Top-level dict keys (sections) are printed in section_color.
    - Nested key/value lines are printed in key_color.
    """
    pad = " " * indent

    if isinstance(data, dict):
        for k, v in data.items():
            if isinstance(v, dict):
                # Section header
                if indent < 5:
                    print(f"{pad}{section_color}{k}{reset}")
                else:
                    print(f"{pad}{key_color}{k}{reset}")
                print_nested_dict(v, indent=indent + 4,
                                  section_color=section_color,
                                  key_color=key_color,
                                  reset=reset)
            else:
                # Leaf item
                print(f"{pad}{key_color}{str(k):<34}{reset}: {v}")
    else:
        # Non-dict fallback
        print(f"{pad}{data}")


def print_pipeline_settings(current_moment, ram_gb, cpu_threads, input_csv, output_dir,
                            header=True, divider=True, sleep_s=0.0):
    """
    Calls get_pipeline_settings(...) and prints the entire nested structure.
    Returns the settings dict for reuse in logs, TUI, etc.
    """
    settings = get_pipeline_settings(
        current_moment=current_moment,
        ram_gb=ram_gb,
        cpu_threads=cpu_threads,
        input_csv=input_csv,
        output_dir=output_dir,
    )

    if header:
        print(f"\n\033[91m{'=' * 80}\033[0m\n")

    print_nested_dict(settings, indent=4)

    if divider:
        print(f"\n\033[91m{'=' * 80}\033[0m\n")

    if sleep_s and sleep_s > 0:
        time.sleep(sleep_s)

    return settings


# --------------------------------------------------------------
# Updates the input.csv with any .fq -> .fastq
# --------------------------------------------------------------
def preprocess_csv(csv_file_path):
    """
    1) Load the CSV at `csv_file_path` into a pandas DataFrame.
    2) Find any cell whose string ends in .fq or .fq.gz.
    3) Rename that file on disk to use .fastq (e.g. file.fq.gz → file.fastq.gz).
    4) Update the DataFrame so that all old “.fq” or “.fq.gz” entries
       become “.fastq” or “.fastq.gz”.
    5) Save (overwrite) the original CSV with these new paths.
    """
    # (a) Read the CSV
    df = pd.read_csv(csv_file_path)

    # (b) Prepare a regex that matches “something.fq” or “something.fq.gz”
    pattern = re.compile(r'^(.*)\.fq(\.gz)?$')

    # (c) Gather all unique values in the DataFrame (so we don’t rename the same file twice)
    all_vals = pd.unique(df.values.ravel())

    # (d) Build a mapping { old_path_str: new_path_str } for everything that matches
    mapping = {}
    for val in all_vals:
        if isinstance(val, str):
            m = pattern.match(val)
            if m:
                base      = m.group(1)         # “path/to/file” without “.fq”
                gz_suffix = m.group(2) or ''   # either “.gz” or ''
                new_path  = f"{base}.fastq{gz_suffix}"

                # Rename on disk (if possible)
                try:
                    os.rename(val, new_path)
                except FileNotFoundError:
                    print(f"Warning: file not found, skipping rename: {val!r}")
                except FileExistsError:
                    print(f"Warning: target already exists, skipping rename: {new_path!r}")
                # In all cases, update the mapping so the DataFrame can be rewritten
                mapping[val] = new_path

    # (e) Replace all occurrences of old_path → new_path in the DataFrame
    df.replace(mapping, inplace=True)

    # (f) Overwrite the original CSV
    df.to_csv(csv_file_path, index=False)
    
    return df

# --------------------------------------------------------------
# Finds a file in a given folder
# --------------------------------------------------------------
def find_files(search_string, folder=None):
    """
    Find all files in the given directory and its subdirectories containing '_dedup' in their names.
    
    Args:
        directory (str): Path to the directory to search
        
    Returns:
        list: List of full file paths containing '_dedup' in their names
    """
    files = []
    
    # Walk through directory and subdirectories
    for root, _, files in os.walk(folder):
        # Check each file in current directory
        for file in files:
            if search_string in file:
                # Add full path to the list
                files.append(os.path.join(root, file))
    return files

# --------------------------------------------------------------
# Locates the one directory under `search_root` that contains ALL of `required_scripts`
# --------------------------------------------------------------
def locate_bin_dir(required_scripts, search_root):
    """
    Walk `search_root` and its subdirectories. At each directory, check if
    that directory’s file‐list contains every name in `required_scripts`.
    If so, return that Path immediately. If no such directory is found, return None.
    """
    # We only need to look at folders (not individual files) and see if each
    # proc_name+".py" is present.
    for root, subdirs, files in os.walk(search_root):
        # root is a str; convert to Path for convenience
        root_path = Path(root)
        # Check if *all* required script filenames appear in `files`
        if all(f"{proc}.py" in files for proc in required_scripts):
            return root_path
    return None


# --------------------------------------------------------------
# Orchestrate the Entheome Genome Assembly Pipeline
# --------------------------------------------------------------
def _cell_has_value(row, col):
    """True when *row[col]* is present and not a blank/``None`` placeholder."""
    if col not in row:
        return False
    val = row[col]
    return pd.notna(val) and str(val).strip() not in ("", "None", "nan")


def sample_is_qc_only(row):
    """Return ``True`` for a QC-only sample: a reference but no reads/EST_SIZE.

    The documented QC-only mode assesses an existing assembly: provide
    ``REF_SEQ`` / ``REF_SEQ_GCA`` (plus BUSCO dbs) with no reads and no
    ``EST_SIZE``. Such samples only need the reference fetched, then QC and the
    report; the assembly-building steps (preprocess reads, decontaminate reads,
    assemble, compare, polish, curate, decontaminate assembly) do not apply.

    This is mode dispatch, not skip-on-failure: read-based samples still run
    every applicable step and fail loudly on error.
    """
    has_ref = _cell_has_value(row, "REF_SEQ") or _cell_has_value(row, "REF_SEQ_GCA")
    read_cols = ("ILLUMINA_SRA", "ILLUMINA_RAW_DIR", "ILLUMINA_RAW_F_READS",
                 "ILLUMINA_RAW_R_READS", "ONT_SRA", "ONT_RAW_DIR", "ONT_RAW_READS",
                 "PACBIO_SRA", "PACBIO_RAW_DIR", "PACBIO_RAW_READS")
    has_reads = any(_cell_has_value(row, c) for c in read_cols)
    has_est_size = _cell_has_value(row, "EST_SIZE")
    return has_ref and not has_reads and not has_est_size


if __name__ == "__main__":
    # Parse command-line arguments for pipeline configuration
    parser = argparse.ArgumentParser(description=f"Run Entheome Genome Assembly Pipeline (EGAP)\nVersion:{VERSION}")
    parser.add_argument("--input_csv", "-csv",
                        type=str,
                        required=True,
                        help="Path to a CSV containing multiple sample data")
    parser.add_argument("--output_dir", "-o",
                        type=str,
                        required=True,
                        help="Directory for pipeline output files")
    parser.add_argument("--cpu_threads", "-t",
                        type=int,
                        default=1,
                        help="Number of CPU threads to use (default: 1)")
    parser.add_argument("--ram_gb", "-r",
                        type=int,
                        default=8,
                        help="Amount of RAM in GB to allocate (default: 8)")
    parser.add_argument("--dry_run", action="store_true", default=False,
                        help="Log file-management actions (removals/compressions) "
                             "without executing them. Equivalent to EGAP_DRY_RUN=1.")
    parser.add_argument("--tui", action="store_true", default=False,
                        help="Run the pipeline through the interactive TUI instead of "
                             "plain terminal output.")
    # Opt-out flags to skip individual assemblers. EGAP otherwise runs every
    # assembler compatible with the supplied read types and picks the best.
    parser.add_argument("--no-masurca", "-no_m", dest="no_masurca", action="store_true",
                        default=False, help="Skip the MaSuRCA assembler.")
    parser.add_argument("--no-flye", "-no_f", dest="no_flye", action="store_true",
                        default=False, help="Skip the Flye assembler.")
    parser.add_argument("--no-spades", "-no_s", dest="no_spades", action="store_true",
                        default=False, help="Skip the SPAdes assembler.")
    parser.add_argument("--no-hifiasm", "-no_h", dest="no_hifiasm", action="store_true",
                        default=False, help="Skip the hifiasm assembler.")

    args = parser.parse_args()
    input_csv   = args.input_csv
    output_dir  = args.output_dir
    cpu_threads = args.cpu_threads
    ram_gb      = args.ram_gb

    if args.dry_run:
        os.environ["EGAP_DRY_RUN"] = "1"
    dry_run = args.dry_run
    current_moment = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    os.makedirs(output_dir, exist_ok=True)

    # ---- TUI mode: hand off to EGAP_TUI and exit ----
    if args.tui:
        _bin_dir = Path(__file__).resolve().parent / "bin"
        if str(_bin_dir) not in sys.path:
            sys.path.insert(0, str(_bin_dir))
        from EGAP_TUI import ENTHEOME_GENOME_ASSEMBLY_PIPELINE  # noqa: PLC0415
        _tui_app = ENTHEOME_GENOME_ASSEMBLY_PIPELINE()
        # Pass the already-parsed namespace so the TUI skips its own argparse.
        _tui_app._preloaded_args = args
        _tui_app.run()
        sys.exit(0)

    # ---- Session-wide log ----
    # The per-sample logs opened in the loop below only begin once each sample
    # starts, so the startup banner, the full tool/pipeline configuration dump,
    # and the launched commands would otherwise never reach disk. Tee stdout to
    # a session log now so the entire run (banner + settings + commands) is
    # captured. _Tee strips ANSI colour codes from the file copy.
    _session_stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    _session_log_path = os.path.join(output_dir, f"EGAP_session_{_session_stamp}_log.txt")
    _session_log_fh = open(_session_log_path, "a", buffering=1)
    sys.stdout = _Tee(sys.stdout, _session_log_fh)
    print(f"NOTE:\tSession log: {_session_log_path}")

    # Print Banner & Parameters
    print(f"""
\033[91m.---.\033[92m________\\\033[38;5;208m=\033[96m/\033[92m________\033[91m.---.\033[0m 
\033[91m|\033[38;5;208m[\033[93m_\033[38;5;208m]\033[91m|\033[94m--------\033[96m/\033[91m=\033[92m\\\033[94m--------\033[91m|\033[38;5;208m[\033[93m_\033[38;5;208m]\033[91m|   \033[94m.\033[96m---------.  \033[94m.\033[96m------.    \033[94m.\033[96m------.    \033[94m.\033[96m-------.\033[0m
\033[91m`---'\033[96m~~~~~~~(\033[38;5;208m===\033[92m)\033[96m~~~~~~~\033[91m`---'  \033[94m/\033[96m|         |\033[94m/\033[96m/        \\ \033[94m/\033[96m/        \\ \033[94m/\033[96m/         \\\033[0m
 \033[92m|\033[94m|\033[96m| \033[92m.--     \033[96m\\\033[91m=\033[92m/    \033[92m,--. \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  .------\033[94m'\033[96m|   .------\033[94m'\033[96m|   .--\033[94m.   \033[96m. |   .--\033[94m.\033[96m   .\033[0m
 \033[92m|\033[94m|\033[96m| \033[94m|-      \033[92m/\033[38;5;208m=\033[96m\\    \033[94m|  _ \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  |\033[94m----"| \033[96m|   |\033[94m----"| \033[96m|   |\033[94m-'\033[96m|   | |   |\033[94m-'\033[96m|   |\033[0m
 \033[92m|\033[94m|\033[96m| \033[96m`--    \033[92m(\033[91m===\033[96m)   \033[96m`--' \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  `----.\033[94m|\033[96m |   |     \033[94m| \033[96m|   +--+   | |   +--+   |\033[0m
 \033[92m|\033[94m|\033[96m|         \033[92m\\\033[38;5;208m=\033[96m/         \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|       |\033[94m| \033[96m|   | \033[94m.\033[96m----.|          | |          |\033[0m
 \033[92m|\033[94m|\033[96m|         \033[96m/\033[91m=\033[92m\\         \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  .----'\033[94m| \033[96m|   |\033[94m"\033[96m|    ||   +--+   | |   +------'\033[0m
 \033[92m|\033[94m|\033[96m|        \033[96m(\033[38;5;208m===\033[92m)        \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  |\033[94m---" | \033[96m|   |\033[94m"\033[96m`-.  ||   |\033[94m| \033[96m|   | |   |\033[94m-----"\033[0m
 \033[92m|\033[94m|\033[96m|  \033[96m__     \033[96m\\\033[91m=\033[92m/    \033[96m,__. \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  `------.|   `---'  ||   |\033[94m| \033[96m|   | |   |\033[0m
 \033[92m|\033[94m|\033[96m| \033[94m/__\\    \033[92m/\033[38;5;208m=\033[96m\\    \033[94m|__| \033[96m|\033[94m|\033[92m|  \033[94m`.\033[96m|         |`.        \033[96m/\033[94m.\033[96m|   |\033[94m| \033[96m|   |\033[94m.\033[96m|   |\033[0m
 \033[92m|\033[94m|\033[96m| \033[92m|  |   \033[92m(\033[91m===\033[96m)   \033[92m|    \033[96m|\033[94m|\033[92m|    \033[96m`---------'  `------'  `---' \033[94m`\033[96m'---' `---'\033[0m
\033[91m.---.\033[96m_____\033[92m[:\033[94m:\033[96m::::\033[94m:\033[92m]\033[96m_____\033[91m.---.\033[92m    \033[92m╔═══════════════════════════════════════════╗\033[0m
\033[91m|\033[38;5;208m[\033[93m_\033[38;5;208m]\033[91m|\033[94m------\033[92m|:\033[94m:\033[96m::\033[94m:\033[92m|\033[94m------\033[91m|\033[38;5;208m[\033[93m_\033[38;5;208m]\033[91m|    \033[92m║     \033[94mEnthe\033[96mome Ge\033[97mnome Assemb\033[96mly Pip\033[94meline     \033[92m║\033[0m
\033[91m`---'\033[92m~~~~~~\033[92m|::\033[94m:\033[96m:\033[94m:\033[92m|~~~~~~\033[91m`---'    \033[92m╚═══════════════════════════════════════════╝\033[0m

                    Curated & Maintained by Ian M Bollinger
                         (\033[94mian.bollinger@entheome.org\033[0m)

                                   \033[92mEGAP.py\033[0m
                                Version {VERSION}
                                  
   Preprocess \033[94m-\033[92m>\033[0m Assemble \033[94m-\033[92m>\033[0m Compare \033[94m-\033[92m>\033[0m Polish \033[94m-\033[92m>\033[0m Curate \033[94m-\033[92m>\033[0m Assess \033[94m-\033[92m>\033[0m Report
    """)
    time.sleep(1.5)
    print(f"""
\033[91m================================================================================\033[0m

    \033[92minput-setup settings
        \033[94mstarted at                       \033[0m: {current_moment}
        \033[94mconfig files                     \033[0m: N/A (running python version)
        \033[94mcontainer                        \033[0m: N/A (running python version)
        \033[94mRAM GB                           \033[0m: {ram_gb}
        \033[94mCPU threads                      \033[0m: {cpu_threads}
        \033[94mdry run                          \033[0m: {dry_run}
        \033[94minput csv                        \033[0m: {input_csv}
        \033[94moutput to                        \033[0m: {output_dir}

    \033[96m-----------------------------------------------------------------------\033[0m
    """)
    time.sleep(0.25)
    print("""
    \033[92mtrimmomatic settings
        \033[94mmode                             \033[0m: -PE
        \033[94mphred version                    \033[0m: -phred33
        \033[94milluminaclip adapter             \033[0m: /opt/conda/envs/EGAP_env/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa
        \033[94mfastaWithAdaptersEtc             \033[0m: 2
        \033[94mseed mismatches                  \033[0m: 30
        \033[94mpalindrome clip threshold        \033[0m: 10
        \033[94msimple clip threshold            \033[0m: 11
        \033[94mHEADCROP                         \033[0m: 10
        \033[94mCROP                             \033[0m: 145
        \033[94mSLIDINGWINDOW                    \033[0m: 50:25
        \033[94mMINLEN                           \033[0m: 125
        
    \033[92mbbduk settings
        \033[94mktrim                            \033[0m: -r
        \033[94mk                                \033[0m: 23
        \033[94mmink                             \033[0m: 11
        \033[94mhdist                            \033[0m: 1
        \033[94mtrimpairsevenly                  \033[0m: -tpe
        \033[94mtrimbyoverlap                    \033[0m: -tbo
        \033[94mqtrim                            \033[0m: -rl
        \033[94mtrimq                            \033[0m: 20

    \033[92mclumpify settings
        \033[94mremove duplicate reads           \033[0m: -dedupe (haploid/unspecified; skipped for diploid+ PLOIDY>1)

    \033[96m-----------------------------------------------------------------------\033[0m
    """)
    time.sleep(0.25)
    print("""
    \033[92mfiltlong settings
        \033[94mmin_length                       \033[0m: 1000
        \033[94mmin_mean_q                       \033[0m: 8
        \033[94mkeep_percent                     \033[0m: 90 
        \033[94mcoverage                         \033[0m: 75
        \033[94mtarget_bases                     \033[0m: estimated size (bp) * 75 (coverage)
        
    \033[92mratatosk settings
        \033[94mverbose output                   \033[0m: -v

    \033[96m-----------------------------------------------------------------------\033[0m
    """)
    time.sleep(0.25)
    print("""
    \033[92mmasurca settings
        \033[94mconfiguration file sections      \033[0m: DATA, PARAMETERS
        \033[94mgraph kmer size                  \033[0m: auto
        \033[94muse linking mates                \033[0m: (0 if Hybrid assembly, else 1)
        \033[94mclose gaps                       \033[0m: (0 if Hybrid assembly, else 1)
        \033[94mmega reads one pass              \033[0m: (1 if long reads present, else 0)
        \033[94mlimit jump coverage              \033[0m: 300
        \033[94mca parameters                    \033[0m: cgwErrorRate=0.15
        \033[94mjellyfish hash size              \033[0m: based on the estimated size of the genome provided
        \033[94msoap assembly                    \033[0m: 0 (disabled to force CABOG assembly)
        \033[94mflye assembly                    \033[0m: 0 (disabled to force CABOG assembly)
                                    
    \033[92mspades settings
        \033[94mhigh-cov. isolate & multi-cell   \033[0m: --isolate
        \033[94mcoverage cutoff                  \033[0m: auto

    \033[92mflye settings
        \033[94mestimated genome size            \033[0m: est_size
        \033[94mnumber of polishing iterations   \033[0m: 3
        \033[94mcollapse alternative haplotypes  \033[0m: --keep-haplotypes
        
    \033[92mhifasm settings
        \033[0mdefault settings

    \033[96m-----------------------------------------------------------------------\033[0m
    """)
    time.sleep(0.25)
    print("""
    \033[92mracon settings 
        \033[0mdefault settings

    \033[92mpilon settings 
        \033[94mchange file generation           \033[0m: -changes
        \033[94mvcf file generation              \033[0m: -vcf
        \033[94mtracks file generation           \033[0m: -tracks
        \033[94mlargest chunksize limit          \033[0m: 5000000
        \033[94mfix list                         \033[0m: indels, local, snps

    \033[96m-----------------------------------------------------------------------\033[0m
    """)
    time.sleep(0.25)
    print("""
    \033[92mragtag settings 
        \033[0mscaffold
           \033[94mconcatenate unplaced          \033[0m: -C
           \033[94madd suffix to unplaced        \033[0m: -u
        \033[0mcorrect
           \033[94madd suffix to unaltered       \033[0m: -u
        \033[0mpatch
           \033[94madd suffix to unplaced        \033[0m: -u
        
    \033[92mtgs-gapcloser settings
        \033[94mdo not error correct             \033[0m: -ne
            
    \033[92mabyss-sealer settings
        \033[94mpseudoreads length used          \033[0m: 400
        \033[94mbloom filter size                \033[0m: 500M

    \033[96m-----------------------------------------------------------------------\033[0m
    """)
    time.sleep(0.25)
    print("""
    \033[92mfastqc settings
       \033[0mdefault settings
    
    \033[92mnanoplot settings
       \033[94mbivariate plots                   \033[0m: kde, dot
       \033[94mshow logarithmic lengths scaling  \033[0m: --loglength
       \033[94mN50 mark in read length histogram \033[0m: --N50
       \033[94mlog messages to terminal          \033[0m: --verbose

    \033[96m-----------------------------------------------------------------------\033[0m
    """)
    time.sleep(0.25)
    print("""
    \033[92mbusco settings
       \033[94mmode                              \033[0m: genome
       \033[94mforce overwrite                   \033[0m: -f

    \033[92mcompleasm settings
       \033[0mdefault settings

    \033[92mquast settings
       \033[0mdefault settings

    \033[96m-----------------------------------------------------------------------\033[0m
    """)
    time.sleep(0.25)
    print("""
    \033[92mkraken2 settings\033[0m  (decontaminate_reads — pre-assembly)
       \033[94mstage                             \033[0m: pre-assembly (reads)
       \033[94mtargets                           \033[0m: ONT and PacBio highest-quality long reads
       \033[94mdatabase                          \033[0m: resolved from KRAKEN2_DB env var or CSV column
       \033[94mkeep (bacteria)                   \033[0m: bacteria, unclassified
       \033[94mkeep (archaea)                    \033[0m: archaea, unclassified
       \033[94mkeep (flora/funga/fauna)          \033[0m: eukarya, unclassified
       \033[94malways kept                       \033[0m: unclassified
       \033[94mmissing db behaviour              \033[0m: WARN and skip (non-fatal)

    \033[92mtiara settings\033[0m  (decontaminate_assembly — post-assembly)
       \033[94mstage                             \033[0m: post-assembly (contigs)
       \033[94mtargets                           \033[0m: final curated assembly
       \033[94mkeep (bacteria)                   \033[0m: bacteria, prokarya, unknown
       \033[94mkeep (archaea)                    \033[0m: archaea, prokarya, unknown
       \033[94mkeep (flora/funga/fauna)          \033[0m: eukarya, organelle, unknown
       \033[94morganelle note                    \033[0m: kept for eukaryotes; removed for prokaryotes
       \033[94malways kept                       \033[0m: unknown
       \033[94mmissing tiara behaviour           \033[0m: WARN and skip (non-fatal)

\033[91m================================================================================\033[0m
    """)
    print(f"""
\033[91m.---.\033[92m________\\\033[38;5;208m=\033[96m/\033[92m________\033[91m.---.\033[0m 
\033[91m|\033[38;5;208m[\033[93m_\033[38;5;208m]\033[91m|\033[94m--------\033[96m/\033[91m=\033[92m\\\033[94m--------\033[91m|\033[38;5;208m[\033[93m_\033[38;5;208m]\033[91m|   \033[94m.\033[96m---------.  \033[94m.\033[96m------.    \033[94m.\033[96m------.    \033[94m.\033[96m-------.\033[0m
\033[91m`---'\033[96m~~~~~~~(\033[38;5;208m===\033[92m)\033[96m~~~~~~~\033[91m`---'  \033[94m/\033[96m|         |\033[94m/\033[96m/        \\ \033[94m/\033[96m/        \\ \033[94m/\033[96m/         \\\033[0m
 \033[92m|\033[94m|\033[96m| \033[92m.--     \033[96m\\\033[91m=\033[92m/    \033[92m,--. \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  .------\033[94m'\033[96m|   .------\033[94m'\033[96m|   .--\033[94m.   \033[96m. |   .--\033[94m.\033[96m   .\033[0m
 \033[92m|\033[94m|\033[96m| \033[94m|-      \033[92m/\033[38;5;208m=\033[96m\\    \033[94m|  _ \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  |\033[94m----"| \033[96m|   |\033[94m----"| \033[96m|   |\033[94m-'\033[96m|   | |   |\033[94m-'\033[96m|   |\033[0m
 \033[92m|\033[94m|\033[96m| \033[96m`--    \033[92m(\033[91m===\033[96m)   \033[96m`--' \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  `----.\033[94m|\033[96m |   |     \033[94m| \033[96m|   +--+   | |   +--+   |\033[0m
 \033[92m|\033[94m|\033[96m|         \033[92m\\\033[38;5;208m=\033[96m/         \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|       |\033[94m| \033[96m|   | \033[94m.\033[96m----.|          | |          |\033[0m
 \033[92m|\033[94m|\033[96m|         \033[96m/\033[91m=\033[92m\\         \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  .----'\033[94m| \033[96m|   |\033[94m"\033[96m|    ||   +--+   | |   +------'\033[0m
 \033[92m|\033[94m|\033[96m|        \033[96m(\033[38;5;208m===\033[92m)        \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  |\033[94m---" | \033[96m|   |\033[94m"\033[96m`-.  ||   |\033[94m| \033[96m|   | |   |\033[94m-----"\033[0m
 \033[92m|\033[94m|\033[96m|  \033[96m__     \033[96m\\\033[91m=\033[92m/    \033[96m,__. \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  `------.|   `---'  ||   |\033[94m| \033[96m|   | |   |\033[0m
 \033[92m|\033[94m|\033[96m| \033[94m/__\\    \033[92m/\033[38;5;208m=\033[96m\\    \033[94m|__| \033[96m|\033[94m|\033[92m|  \033[94m`.\033[96m|         |`.        \033[96m/\033[94m.\033[96m|   |\033[94m| \033[96m|   |\033[94m.\033[96m|   |\033[0m
 \033[92m|\033[94m|\033[96m| \033[92m|  |   \033[92m(\033[91m===\033[96m)   \033[92m|    \033[96m|\033[94m|\033[92m|    \033[96m`---------'  `------'  `---' \033[94m`\033[96m'---' `---'\033[0m
\033[91m.---.\033[96m_____\033[92m[:\033[94m:\033[96m::::\033[94m:\033[92m]\033[96m_____\033[91m.---.\033[92m    \033[92m╔═══════════════════════════════════════════╗\033[0m
\033[91m|\033[38;5;208m[\033[93m_\033[38;5;208m]\033[91m|\033[94m------\033[92m|:\033[94m:\033[96m::\033[94m:\033[92m|\033[94m------\033[91m|\033[38;5;208m[\033[93m_\033[38;5;208m]\033[91m|    \033[92m║     \033[94mEnthe\033[96mome Ge\033[97mnome Assemb\033[96mly Pip\033[94meline     \033[92m║\033[0m
\033[91m`---'\033[92m~~~~~~\033[92m|::\033[94m:\033[96m:\033[94m:\033[92m|~~~~~~\033[91m`---'    \033[92m╚═══════════════════════════════════════════╝\033[0m

                    Curated & Maintained by Ian M Bollinger
                         (\033[94mian.bollinger@entheome.org\033[0m)

                                   \033[92mEGAP.py\033[0m
                                Version {VERSION}
                                  
   Preprocess \033[94m-\033[92m>\033[0m Assemble \033[94m-\033[92m>\033[0m Compare \033[94m-\033[92m>\033[0m Polish \033[94m-\033[92m>\033[0m Curate \033[94m-\033[92m>\033[0m Assess \033[94m-\033[92m>\033[0m Report
    """)
    time.sleep(1.5)


    # Print settings from the structured dict (single source of truth)
    pipeline_settings = print_pipeline_settings(
        current_moment=current_moment,
        ram_gb=ram_gb,
        cpu_threads=cpu_threads,
        input_csv=input_csv,
        output_dir=output_dir,
        header=True,
        divider=True,
        sleep_s=0.0,
    )
   
    processes = ["preprocess_refseq", "preprocess_illumina", "preprocess_ont",
                 "preprocess_pacbio", "decontaminate_reads",
                 "assemble_masurca", "assemble_flye",
                 "assemble_spades", "assemble_hifiasm", "compare_assemblies",
                 "polish_assembly", "curate_assembly", "decontaminate_assembly"]

    # Honour per-assembler skip flags (--no-masurca/-no_m, etc.). EGAP runs every
    # assembler compatible with the read types unless explicitly disabled here.
    _skip_assemblers = {
        "assemble_masurca": args.no_masurca,
        "assemble_flye":    args.no_flye,
        "assemble_spades":  args.no_spades,
        "assemble_hifiasm": args.no_hifiasm,
    }
    _disabled = [proc for proc, skip in _skip_assemblers.items() if skip]
    if _disabled:
        processes = [p for p in processes if p not in _disabled]
        print(f"NOTE:\tSkipping assembler(s) by request: "
              f"{', '.join(p.replace('assemble_', '') for p in _disabled)}")
    if all(_skip_assemblers.values()):
        print("WARN:\tAll assemblers were disabled; no assembly will be produced.")

    # Determine project root and locate the bin/ directory there
    this_file   = Path(__file__).resolve()
    project_dir = this_file.parent

    # Attempt to find a directory under project_dir that contains all "<proc>.py"
    bin_dir_candidate = locate_bin_dir(processes, project_dir)

    if bin_dir_candidate is not None:
        bin_dir = bin_dir_candidate
    else:
        # Fallback: maybe everything was installed flatly in project_dir/bin
        if (project_dir / "bin").is_dir():
            bin_dir = project_dir / "bin"
        else:
            raise FileNotFoundError("Could not locate a single folder containing all of: "
                                    + ", ".join(f"{p}.py" for p in processes))

    # ---- Fail-fast environment preflight (tools on PATH, RAM sanity) ----
    if str(bin_dir) not in sys.path:
        sys.path.insert(0, str(bin_dir))
    try:
        from preflight_checks import run_preflight
        run_preflight(ram_gb)
    except SystemExit:
        raise  # a failed check aborts the run loudly, as intended
    except Exception as _pf_exc:
        print(f"WARN:\tPre-flight checks could not run ({_pf_exc}); continuing.")

    # Load the Sample Table and correct any ".fq" -> ".fastq"
    input_csv_df = preprocess_csv(input_csv)

    # Resolve QC and reporter scripts once (outside the sample loop)
    qc_script = bin_dir / "qc_assessment.py"
    if not qc_script.exists():
        raise FileNotFoundError(f"Missing QC script: {qc_script}")
    reporter_script = bin_dir / "html_reporter.py"
    if not reporter_script.exists():
        raise FileNotFoundError(f"Missing Reporter script: {reporter_script}")

    # Process each sample fully — pipeline steps, then QC, then HTML — before
    # moving on to the next sample.
    failed_samples = []
    _real_stdout = sys.stdout
    for index, row in input_csv_df.iterrows():
        sample_id = row["SAMPLE_ID"]

        # Open a per-sample log file and tee stdout into it for this sample's run.
        _sample_log_path = os.path.join(output_dir, f"{sample_id}_log.txt")
        _sample_log_fh = open(_sample_log_path, "a", buffering=1)
        sys.stdout = _Tee(_real_stdout, _sample_log_fh)

        print(f"\n{'='*70}")
        print(f"NOTE:\tBeginning pipeline for sample: {sample_id}")
        print(f"NOTE:\tSample log: {_sample_log_path}")
        print(f"{'='*70}\n")

        sample_step_failed = False

        # QC-only samples (reference, no reads, no EST_SIZE) only need the
        # reference fetched; QC + report run afterward. The assembly-building
        # steps do not apply -- mode dispatch, not skip-on-failure. Read-based
        # samples run every step in `processes` and fail loudly on error.
        if sample_is_qc_only(row):
            sample_processes = [p for p in processes if p == "preprocess_refseq"]
            print("NOTE:\tQC-only sample (reference, no reads, no EST_SIZE): "
                  "running reference fetch + QC + report only.")
        else:
            sample_processes = processes

        # ---- Pipeline steps ----
        for proc in sample_processes:
            script = bin_dir / f"{proc}.py"
            if not script.exists():
                print(f"\nERROR:\tMissing script: {script}")
                sample_step_failed = True
                break
            current_cmd = [sys.executable,
                           str(script),
                           sample_id,
                           input_csv,
                           output_dir,
                           str(cpu_threads),
                           str(ram_gb)]
            print(f"\n→ Running {proc}: {' '.join(current_cmd)}\n")
            current_return_code = run_filtered(current_cmd)
            if current_return_code != 0:
                print(f"\nERROR:\t{proc} failed for {sample_id} (rc={current_return_code})")
                sample_step_failed = True
                break

        if sample_step_failed:
            print(f"\nWARN:\tOne or more steps failed for {sample_id}; "
                  f"still running final QC and HTML report.")

        # ---- Final QC assessment (always runs per sample) ----
        qc_cmd = [sys.executable,
                  str(qc_script),
                  "final",
                  input_csv,
                  sample_id,
                  output_dir,
                  str(cpu_threads),
                  str(ram_gb)]
        print(f"\n→ Running qc_assessment: {' '.join(qc_cmd)}\n")
        qc_return_code = run_filtered(qc_cmd)
        if qc_return_code != 0:
            print(f"\nWARN:\tqc_assessment returned non-zero for {sample_id} (rc={qc_return_code})")
            sample_step_failed = True

        # ---- HTML report (always runs per sample) ----
        reporter_cmd = [sys.executable,
                        str(reporter_script),
                        sample_id,
                        input_csv,
                        output_dir,
                        str(cpu_threads),
                        str(ram_gb)]
        print(f"\n→ Running html_reporter: {' '.join(reporter_cmd)}\n")
        reporter_return_code = run_filtered(reporter_cmd)
        if reporter_return_code != 0:
            print(f"\nWARN:\thtml_reporter returned non-zero for {sample_id} (rc={reporter_return_code})")
            sample_step_failed = True

        if sample_step_failed:
            print(f"\nFAIL:\tSample {sample_id} completed with errors. Continuing to next sample.")
            failed_samples.append(sample_id)
        else:
            print(f"\nPASS:\tSample {sample_id} completed successfully.")

        # Close this sample's log and restore stdout before moving to the next sample.
        sys.stdout = _real_stdout
        _sample_log_fh.close()

    if failed_samples:
        print(f"\nWARN:\tAll samples processed. The following had errors: {', '.join(failed_samples)}")
    else:
        print("\nPASS:\tAll samples processed successfully.")
