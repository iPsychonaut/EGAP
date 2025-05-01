#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
assemble_spades.py

Updated on Sat Mar 29 2025

This script runs SPAdes assembly with Illumina and optional long reads.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import os, sys, shutil
import pandas as pd
from utilities import run_subprocess_cmd, get_current_row_data, pigz_compress, pigz_decompress
from qc_assessment import qc_assessment


# --------------------------------------------------------------
# Run SPAdes assembly with Illumina and optional long reads
# --------------------------------------------------------------
def assemble_spades(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """Assemble genomic data using SPAdes with Illumina and optional ONT reads.

    Executes SPAdes assembly, processes input reads from a CSV, and performs
    quality control on the resulting assembly.

    Args:
        sample_id (str): Sample identifier.
        input_csv (str): Path to metadata CSV file.
        output_dir (str): Directory for output files.
        cpu_threads (int): Number of CPU threads to use.
        ram_gb (int): Available RAM in GB.

    Returns:
        str or None: Path to the gzipped SPAdes assembly FASTA, or None if no valid reads are provided.
    """
    input_df = pd.read_csv(input_csv)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]  # Convert to Series (single row)

    # Identify read paths, reference, and BUSCO lineage info from CSV
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

    species_dir = os.path.join(output_dir, species_id)

    if pd.notna(ont_sra) and pd.isna(ont_raw_reads):
        ont_raw_reads = os.path.join(species_dir, "ONT", f"{ont_sra}.fastq.gz")
    if pd.notna(illumina_sra) and pd.isna(illumina_f_raw_reads) and pd.isna(illumina_r_raw_reads):
        illumina_f_raw_reads = os.path.join(species_dir, "Illumina", f"{illumina_sra}_1.fastq.gz")
        illumina_r_raw_reads = os.path.join(species_dir, "Illumina", f"{illumina_sra}_2.fastq.gz")    
    if pd.notna(pacbio_sra) and pd.isna(pacbio_raw_reads):
        pacbio_raw_reads = os.path.join(species_dir, "PacBio", f"{pacbio_sra}.fastq.gz")
    if pd.notna(ref_seq_gca) and pd.isna(ref_seq):
        ref_seq = os.path.join(species_dir, "RefSeq", f"{species_id}_{ref_seq_gca}_RefSeq.fasta.gz")

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
    
    # Set long-read paths (NEW: Simplified to match original, using raw reads if preprocessed unavailable)
    highest_mean_qual_long_reads_gz = None
    highest_mean_qual_long_reads = None
    if pd.notna(ont_raw_reads):
        print("DEBUG - ONT RAW READS EXIST!")
        highest_mean_qual_long_reads_gz = os.path.join(species_dir, "ONT", f"{species_id}_ONT_highest_mean_qual_long_reads.fastq.gz")
        highest_mean_qual_long_reads = highest_mean_qual_long_reads_gz.replace(".gz", "")
    elif pd.notna(pacbio_raw_reads):
        print("SKIP:\tSPAdes will not run on PacBio data only.")
        return None

    # Set Illumina deduplicated read paths only if Illumina reads are present
    illu_dedup_f_reads_gz = None
    illu_dedup_r_reads_gz = None
    illu_dedup_f_reads = None
    illu_dedup_r_reads = None
    if pd.notna(illumina_f_raw_reads) and pd.notna(illumina_r_raw_reads):
        illu_dedup_f_reads_gz = os.path.join(species_dir, "Illumina", f"{species_id}_illu_forward_dedup.fastq.gz")
        illu_dedup_r_reads_gz = os.path.join(species_dir, "Illumina", f"{species_id}_illu_reverse_dedup.fastq.gz")

    # Unzip deduplicated files if they are gzipped
    dedup_unzip_list = []
    if illu_dedup_f_reads_gz and illu_dedup_r_reads_gz:
        dedup_unzip_list = [illu_dedup_f_reads_gz.replace(".gz", ""), illu_dedup_r_reads_gz.replace(".gz", "")]
        for input_file in [illu_dedup_f_reads_gz, illu_dedup_r_reads_gz]:
            if ".gz" in input_file and os.path.exists(input_file):
                _ = pigz_decompress(input_file, cpu_threads)
    if len(dedup_unzip_list) == 2:
        illu_dedup_f_reads = dedup_unzip_list[0]
        illu_dedup_r_reads = dedup_unzip_list[1]

    print(f"DEBUG - illu_dedup_f_reads - {illu_dedup_f_reads}")
    print(f"DEBUG - illu_dedup_r_reads - {illu_dedup_r_reads}")
    print(f"DEBUG - highest_mean_qual_long_reads_gz - {highest_mean_qual_long_reads_gz}")

    # Unzip long-read file if necessary (NEW: Ensure decompression of preprocessed file)
    if pd.notna(pacbio_raw_reads) :
      print("SKIP:\tCannot run SPAdes on PacBio only.")  
      return None
    elif pd.notna(ont_raw_reads) and highest_mean_qual_long_reads_gz:
        if os.path.exists(highest_mean_qual_long_reads_gz):
            if not os.path.exists(highest_mean_qual_long_reads):
                print(f"NOTE:\tUnzipping {highest_mean_qual_long_reads_gz}")
                highest_mean_qual_long_reads = pigz_decompress(highest_mean_qual_long_reads_gz, cpu_threads)
        
    print(f"DEBUG - highest_mean_qual_long_reads - {highest_mean_qual_long_reads}")

    # Unzip reference sequence if necessary (NEW: Corrected to ensure decompressed path)
    if pd.notna(ref_seq):
        ref_seq_unzipped = ref_seq.replace(".gz","")
    if pd.notna(ref_seq) and ref_seq.endswith(".gz") and os.path.exists(ref_seq):
        ref_seq_unzipped = ref_seq.replace(".gz", "")
        if not os.path.exists(ref_seq_unzipped):
            print(f"NOTE:\tUnzipping {ref_seq}")
            ref_seq_unzipped = pigz_decompress(ref_seq, cpu_threads)

    if pd.isna(ont_raw_reads) and pd.isna(illumina_f_raw_reads) and pd.isna(illumina_r_raw_reads) and pd.isna(pacbio_raw_reads):
        print("SKIP:\tNo reads available for processing")
        return None

    sample_dir = os.path.join(species_dir, sample_id)
    os.makedirs(sample_dir, exist_ok=True)    
    spades_out_dir = os.path.join(sample_dir, "spades_assembly")
    os.makedirs(spades_out_dir, exist_ok=True)
    egap_spades_assembly_path = os.path.join(spades_out_dir, f"{sample_id}_spades.fasta")
    egap_spades_assembly_path_gz = egap_spades_assembly_path + ".gz"
    
    print(f"DEBUG - spades_out_dir - {spades_out_dir}")

    # SPAdes Hybrid Assembly (Illumina w/ ONT, or Illumina Only)
    if os.path.exists(egap_spades_assembly_path):
        print(f"SKIP:\tSPAdes Assembly, scaffolded assembly already exists: {egap_spades_assembly_path}.")
        if ".gz" not in egap_spades_assembly_path:
            egap_spades_assembly_path_gz = pigz_compress(egap_spades_assembly_path, cpu_threads)
        else:
            egap_spades_assembly_path_gz = egap_spades_assembly_path
        egap_spades_assembly_path_gz, masurca_stats_list, _ = qc_assessment("spades", input_csv, sample_id, output_dir, cpu_threads, ram_gb)
        return egap_spades_assembly_path_gz
    elif os.path.exists(egap_spades_assembly_path + ".gz"):
        egap_spades_assembly_path_gz = egap_spades_assembly_path + ".gz"
        egap_spades_assembly_path_gz, masurca_stats_list, _ = qc_assessment("spades", input_csv, sample_id, output_dir, cpu_threads, ram_gb)
        print(f"SKIP:\tMaSuRCA Assembly, scaffolded assembly already exists: {egap_spades_assembly_path_gz}.")
        return egap_spades_assembly_path_gz
    else:
        kmer_list = ["21", "33", "55", "77", "99"]        
        spades_work_dir = os.path.join(os.getcwd(),"spades_assembly")
        spades_path = os.path.join(spades_work_dir, "scaffolds.fasta")
        
        print(f"DEBUG - spades_path - {spades_path}")

        spades_base = ["spades.py",
                       "--isolate",
                       "-t", str(cpu_threads),
                       "-m", str(ram_gb),
                       "--cov-cutoff", "auto",]
        # after you’ve un-gzipped as needed:
        spades_cmd = spades_base + []
        
        # 3a) Illumina paired-end?
        if illu_dedup_f_reads and illu_dedup_r_reads:
            spades_cmd += ["-1", illu_dedup_f_reads,
                           "-2", illu_dedup_r_reads]
        
        # 3b) Nanopore hybrid?
        if highest_mean_qual_long_reads and pd.notna(ont_raw_reads):
            spades_cmd += ["--nanopore", highest_mean_qual_long_reads]
                
        # 3c) Reference contigs?
        if pd.notna(ref_seq):
            spades_cmd += ["--trusted-contigs", ref_seq_unzipped]
        
        # 3d) output & k-mer list  
        spades_cmd += ["-o", spades_work_dir,
                       "-k", ",".join(kmer_list)]
        
        print(f"DEBUG - spades_cmd - {spades_cmd}")
        
        _ = run_subprocess_cmd(spades_cmd, shell_check = False)
        
        shutil.move(spades_path, egap_spades_assembly_path)
        
        if not os.path.exists(egap_spades_assembly_path_gz) and os.path.exists(egap_spades_assembly_path):
            egap_spades_assembly_path_gz = pigz_compress(egap_spades_assembly_path, cpu_threads)

        egap_spades_assembly_path_gz, spades_stats_list, _ = qc_assessment("spades", input_csv, sample_id, output_dir, cpu_threads, ram_gb)
    
        return egap_spades_assembly_path_gz


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python3 assemble_spades.py <sample_id> <input_csv> "
            "<output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)
        
    egap_spades_assembly_path_gz = assemble_spades(sys.argv[1],       # sample_id
                                                   sys.argv[2],       # input_csv
                                                   sys.argv[3],       # output_dir
                                                   str(sys.argv[4]),  # cpu_threads
                                                   str(sys.argv[5]))  # ram_gb