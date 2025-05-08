#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EGAP.py

Orchestrates the Entheome Genome Assembly Pipeline (EGAP) by executing a series of
preprocessing, assembly, comparison, polishing, curation, and quality assessment steps
for genomic data. Processes multiple samples from a CSV input, utilizing specified
CPU threads and RAM resources.

Created on Fri May  2 21:18:21 2025

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import argparse, sys, time, subprocess
from pathlib import Path
import pandas as pd
from datetime import datetime

version = "3.0.0c"


# --------------------------------------------------------------
# Orchestrate the Entheome Genome Assembly Pipeline
# --------------------------------------------------------------
if __name__ == "__main__":
    # Parse command-line arguments for pipeline configuration
    parser = argparse.ArgumentParser(description="Run Entheome Genome Assembly Pipeline (EGAP)")

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
                         (\033[94mian.bollinger@entheome.org)\033[0m

                              \033[92mdraft_assembly.nf\033[0m
                                version {version}
                                  
 Input-Setup \033[94m-\033[92m>\033[0m Preprocess \033[94m-\033[92m>\033[0m Assemble \033[94m-\033[92m>\033[0m Compare \033[94m-\033[92m>\033[0m Polish \033[94m-\033[92m>\033[0m Curate \033[94m-\033[92m>\033[0m Assess
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
                         (\033[94mian.bollinger@entheome.org)\033[0m

                              \033[92mdraft_assembly.nf\033[0m
                                version {version}
                                  
\033[91m================================================================================\033[0m
    """)

    # Determine project root and locate the bin/ directory there
    this_file = Path(__file__).resolve()
    project_dir = this_file.parent
    bin_dir = project_dir / "bin"
    
    # Load the Sample Table
    input_csv_df = pd.read_csv(input_csv)

    processes = ["preprocess_refseq", "preprocess_ont", "preprocess_illumina",
                 "preprocess_pacbio", "assemble_masurca", "assemble_flye",
                 "assemble_spades", "assemble_hifiasm", "compare_assemblies",
                 "polish_assembly", "curate_assembly"]

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

    print("\nPASS:\tAll samples processed successfully.")
