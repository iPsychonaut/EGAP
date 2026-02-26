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
# Establish Global Pipeline Settings (i.e. command variables)
# --------------------------------------------------------------
VERSION = "3.3.8"

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
    "settings": "default",
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
    df.to_csv(csv_file_path, index=False, na_rep="None")
    
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

    args = parser.parse_args()
    input_csv   = args.input_csv
    output_dir  = args.output_dir
    cpu_threads = args.cpu_threads
    ram_gb      = args.ram_gb
    current_moment = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

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
        \033[94mremove duplicate reads           \033[0m: -dedupe

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
        \033[94mmega reads one pass              \033[0m: 0
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
    
    \033[0m*\033[92mcompleasm settings
       \033[0mdefault settings
    
    \033[92mquast settings
       \033[0mdefault settings

    \033[0m*: Currently disabled since version 3.0.0

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
                 "preprocess_pacbio", "assemble_masurca", "assemble_flye",
                 "assemble_spades", "assemble_hifiasm", "compare_assemblies",
                 "polish_assembly", "curate_assembly"]

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

    # Load the Sample Table and correct any ".fq" -> ".fastq"
    input_csv_df = preprocess_csv(input_csv)

    # Loop through each process for each sample
    for index, row in input_csv_df.iterrows():
        sample_id = row["SAMPLE_ID"]
        for proc in processes:
            script = bin_dir / f"{proc}.py"
            if not script.exists():
                raise FileNotFoundError(f"Missing script: {script}")
            current_cmd = [sys.executable,
                   str(script),
                   sample_id,
                   input_csv,
                   output_dir,
                   str(cpu_threads),
                   str(ram_gb)]
            print(f"\n→ Running {proc}: {' '.join(current_cmd)}\n")
            current_process = subprocess.Popen(current_cmd)
            current_return_code = current_process.wait()
            if current_return_code != 0:
                raise RuntimeError(f"{proc} failed with return code {current_return_code}")

    # final QC
    qc_script = bin_dir / "qc_assessment.py"
    if not qc_script.exists():
        raise FileNotFoundError(f"Missing QC script: {qc_script}")
    for index, row in input_csv_df.iterrows():
        sample_id = row["SAMPLE_ID"]
        qc_cmd = [sys.executable,
                  str(qc_script),
                  "final",
                  input_csv,
                  sample_id,
                  output_dir,
                  str(cpu_threads),
                  str(ram_gb)]
        print(f"\n→ Running qc_assessment: {' '.join(qc_cmd)}\n")
        qc_process = subprocess.Popen(qc_cmd)
        qc_return_code = qc_process.wait()
        if qc_return_code != 0:
            raise RuntimeError(f"qc_assessment failed with return code {qc_return_code}")

    # Report Generation
    reporter_script = bin_dir / "html_reporter.py"
    if not reporter_script.exists():
        raise FileNotFoundError(f"Missing Reporter script: {reporter_script}")
    for index, row in input_csv_df.iterrows():
        sample_id = row["SAMPLE_ID"]
        reporter_cmd = [sys.executable,
                        str(reporter_script),
                        sample_id,
                        input_csv,
                        output_dir,
                        str(cpu_threads),
                        str(ram_gb)]
        print(f"\n→ Running html_reporter: {' '.join(reporter_cmd)}\n")
        reporter_process = subprocess.Popen(reporter_cmd)
        reporter_return_code = reporter_process.wait()
        if reporter_return_code != 0:
            raise RuntimeError(f"html_reporter failed with return code {reporter_return_code}")
    
    print("\nPASS:\tAll samples processed successfully.")
