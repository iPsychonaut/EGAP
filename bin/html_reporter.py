#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
html_reporter.py

Created on Mon May 12 08:37:29 2025

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""

import os, glob, sys
import pandas as pd
from bs4 import BeautifulSoup
from jinja2 import Environment, FileSystemLoader
from datetime import datetime
from pathlib import Path
from process_metadata import get_inat_obs, reverse_geocode

# Helper: parse a NanoStats.txt file into a pandas.DataFrame
def parse_nanostats(nanostats_path):
    """
    Reads a NanoStats.txt summary file and returns a DataFrame
    with columns ["Metric", "Value"].
    """
    rows = []
    with open(nanostats_path, "r") as fh:
        for line in fh:
            line = line.strip()
            # skip empty or non‐key:value lines
            if not line or ":" not in line:
                continue
            key, val = line.split(":", 1)
            rows.append({"Metric": key.strip(), "Value": val.strip()})
    return pd.DataFrame(rows)


# Helper: parse a FastQC HTML file into a pandas.DataFrame
def parse_fastqc_basic_stats(html_path):
    """
    Parses the FastQC HTML report at `html_path` and returns
    a DataFrame of the “Basic Statistics” table (Metric / Value).
    """
    with open(html_path, "r", encoding="utf-8") as fh:
        soup = BeautifulSoup(fh, "html.parser")

    # Find the <h2> whose text (including children) contains "Basic Statistics"
    header = soup.find(lambda tag: tag.name == "h2" 
                                  and "Basic Statistics" in tag.get_text())
    if not header:
        raise ValueError(f"No Basic Statistics header found in {html_path}")

    # The very next <table> after that header
    table = header.find_next("table")
    if not table:
        raise ValueError(f"No table found after Basic Statistics header in {html_path}")

    # Extract rows
    rows = []
    for tr in table.find_all("tr"):
        tds = tr.find_all("td")
        if len(tds) == 2:
            metric = tds[0].get_text(strip=True)
            value   = tds[1].get_text(strip=True)
            rows.append({"Metric": metric, "Value": value})
    return pd.DataFrame(rows)


# Helper: parse a QUAST TSV file into a pandas.DataFrame and the HTML Link
def process_quast_folder(quast_dir):
    """
    Given a path to a *_quast folder, finds the report.tsv and report.html,
    parses the TSV into a DataFrame, and returns (df, html_path).
    Raises if either file not found.
    """
    # Find report.tsv
    tsv_paths = glob.glob(os.path.join(quast_dir, "report.tsv"))
    if not tsv_paths:
        raise FileNotFoundError(f"No report.tsv in {quast_dir}")
    tsv = tsv_paths[0]
    
    # Read TSV, replacing "#" with "Number" in the first column
    df = pd.read_csv(tsv, sep="\t")
    df.iloc[:, 0] = df.iloc[:, 0].str.replace("#", "Number of", regex=False)
    
    # QUAST report.tsv has columns: “Metric” and “Value”
    df = df.set_index(df.columns[0])
    df.columns = ["Value"]

    # Find HTML report
    html_paths = glob.glob(os.path.join(quast_dir, "report.html"))
    if not html_paths:
        raise FileNotFoundError(f"No report.html in {quast_dir}")
    html = html_paths[0]

    return df, html


# Extract and process sample metadata
def get_current_row_data(input_df, sample_id):
    """Extract row data for a sample ID and generate a stats dictionary.

    Filters a DataFrame for a specific sample ID and creates a statistics dictionary.

    Args:
        input_df (pandas.DataFrame): DataFrame with sample metadata.
        sample_id (str): Sample identifier to filter.

    Returns:
        tuple: (filtered DataFrame row, row index list, sample statistics dictionary).
    """
    # Filter the DataFrame for rows where the "SAMPLE_ID" column equals the provided sample_id
    current_row = input_df[input_df["SAMPLE_ID"] == sample_id]
    sample_stats_dict = gen_sample_stats_dict(current_row)
    current_index = current_row.index.tolist()
    
    return current_row, current_index, sample_stats_dict


# Create a sample statistics dictionary from metadata
def gen_sample_stats_dict(row):
    """Generate a sample statistics dictionary from a metadata row.

    Extracts key fields from a pandas Series and initializes placeholders for metrics.

    Args:
        row (pandas.Series): Metadata row containing sample information.

    Returns:
        dict: Dictionary with initialized sample statistics.
    """
    sample_stats_dict = {"SAMPLE_ID": row["SAMPLE_ID"],
                         "SPECIES_ID": row["SPECIES_ID"],
                         "ONT_SRA": row["ONT_SRA"] if isinstance(row["ONT_SRA"], str) else None,
                         "ONT": os.path.basename(row["ONT_RAW_READS"]) if isinstance(row["ONT_RAW_READS"], str) else None,
                         "ILLU_SRA": row["ILLUMINA_SRA"] if isinstance(row["ILLUMINA_SRA"], str) else None,
                         "ILLU_F": os.path.basename(row["ILLUMINA_RAW_F_READS"]) if isinstance(row["ILLUMINA_RAW_F_READS"], str) else None,
                         "ILLU_R": os.path.basename(row["ILLUMINA_RAW_R_READS"]) if isinstance(row["ILLUMINA_RAW_R_READS"], str) else None,
                         "PACBIO_SRA": row["PACBIO_SRA"] if isinstance(row["PACBIO_SRA"], str) else None,
                         "PACBIO": os.path.basename(row["PACBIO_RAW_READS"]) if isinstance(row["PACBIO_RAW_READS"], str) else None,
                         "REF_SEQ_GCA": row["REF_SEQ_GCA"] if isinstance(row["REF_SEQ_GCA"], str) else None,
                         "REF_SEQ": os.path.basename(row["REF_SEQ"]) if isinstance(row["REF_SEQ"], str) else None,
                         "RAW_ILLU_TOTAL_BASES": None,
                         "RAW_ILLU_COVERAGE": None,
                         "TRIMMED_ILLU_TOTAL_BASES": None,
                         "TRIMMED_ILLU_COVERAGE": None,
                         "DEDUPED_ILLU_TOTAL_BASES": None,
                         "DEDUPED_ILLU_COVERAGE": None,
                         "RAW_ONT_READS": None,
                         "RAW_ONT_MEAN_LENGTH": None,
                         "RAW_ONT_MEAN_QUAL": None,
                         "RAW_ONT_TOTAL_BASES": None,
                         "RAW_ONT_COVERAGE": None,
                         "FILT_ONT_READS": None,
                         "FILT_ONT_MEAN_LENGTH": None,
                         "FILT_ONT_MEAN_QUAL": None,
                         "FILT_ONT_TOTAL_BASES": None,
                         "FILT_ONT_COVERAGE": None,
                         "CORRECT_ONT_READS": None,
                         "CORRECT_ONT_MEAN_LENGTH": None,
                         "CORRECT_ONT_MEAN_QUAL": None,
                         "CORRECT_ONT_TOTAL_BASES": None,
                         "CORRECT_ONT_COVERAGE": None,
                         "KMER_COMPLETENESS": None,
                         "QUAL_VAL": None,
                         "RAW_PACBIO_READS": None,
                         "RAW_PACBIO_MEAN_LENGTH": None,
                         "RAW_PACBIO_MEAN_QUAL": None,
                         "RAW_PACBIO_TOTAL_BASES": None,
                         "RAW_PACBIO_COVERAGE": None,
                         "HIFI_PACBIO_READS": None,
                         "HIFI_PACBIO_MEAN_LENGTH": None,
                         "HIFI_PACBIO_MEAN_QUAL": None,
                         "HIFI_PACBIO_TOTAL_BASES": None,
                         "HIFI_PACBIO_COVERAGE": None,
                         "FILT_PACBIO_READS": None,
                         "FILT_PACBIO_MEAN_LENGTH": None,
                         "FILT_PACBIO_MEAN_QUAL": None,
                         "FILT_PACBIO_TOTAL_BASES": None,
                         "FILT_PACBIO_COVERAGE": None,                        
                         "FIRST_COMPLEASM_S": None,
                         "FIRST_COMPLEASM_D": None,
                         "FIRST_COMPLEASM_F": None,
                         "FIRST_COMPLEASM_M": None,
                         "FIRST_COMPLEASM_C": None,
                         "SECOND_COMPLEASM_S": None,
                         "SECOND_COMPLEASM_D": None,
                         "SECOND_COMPLEASM_F": None,
                         "SECOND_COMPLEASM_M": None,
                         "SECOND_COMPLEASM_C": None,
                         "GENOME_SIZE": None,
                         "ASSEMBLY_READS": None,
                         "ASSEMBLY_CONTIGS": None,
                         "ASSEMBLY_N50": None,
                         "ASSEMBLY_L50": None,
                         "ASSEMBLY_GC": None,
                         "MISASSEMBLIES": None,
                         "N_PER_100KBP": None,
                         "MIS_PER_100KBP": None,
                         "INDELS_PER_100KPB": None,
                         "FINAL_ASSEMBLY": None}
    return sample_stats_dict


def html_reporter(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """
    PEP 8 Documentation
    """
    # Path to this script
    script_file    = Path(__file__).resolve()
    bin_dir        = script_file.parent
    resources_dir  = bin_dir.parent / "resources"
    templates_dir  = resources_dir / "templates"
    print(f"DEBUG - templates_dir - {templates_dir}")

    # point Jinja2 at resources/templates
    env = Environment(loader=FileSystemLoader(str(templates_dir)),
        autoescape=True)

    # load by template name only
    main_tmpl  = env.get_template("EGAP_summary.html")
    busco_tmpl = env.get_template("EGAP_busco.html")
    
    # Format current time as YYYYMMDD-HH:MM:SS
    now = datetime.now().strftime("%Y%m%d-%H:%M:%S")
    
    # Parsing Input CSV for the current row that needs a 
    input_df = pd.read_csv(input_csv)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]
    
    # Identify read paths, reference, and BUSCO lineage info from CSV
    illumina_sra = current_series["ILLUMINA_SRA"]
    illumina_f_raw_reads = current_series["ILLUMINA_RAW_F_READS"]
    illumina_r_raw_reads = current_series["ILLUMINA_RAW_R_READS"]
    ont_sra = current_series["ONT_SRA"]
    ont_raw_reads = current_series["ONT_RAW_READS"]
    pacbio_sra = current_series["PACBIO_SRA"]
    pacbio_raw_reads = current_series["PACBIO_RAW_READS"]
    first_busco_db = current_series["BUSCO_1"]
    second_busco_db = current_series["BUSCO_2"]
    ref_seq_gca = current_series["REF_SEQ_GCA"]
    ref_seq = current_series["REF_SEQ"]
    species_id = current_series["SPECIES_ID"]
    inat_id = current_series["INATRUALIST_ID"]
    est_size = current_series["EST_SIZE"]
    
    # Generate Expected Output Directories
    species_dir = os.path.join(output_dir, species_id)
    sample_dir = os.path.join(species_dir, sample_id)
    
    # Directories that MIGHT EXIST and MIGHT contain data
    ont_dir = os.path.join(species_dir, "ONT")
    illumina_dir = os.path.join(species_dir, "Illumina")
    pacbio_dir = os.path.join(species_dir, "PacBio")
    masurca_dir = os.path.join(sample_dir, "masurca_assembly")
    flye_dir = os.path.join(sample_dir, "flye_assembly")
    spades_dir = os.path.join(sample_dir, "spades_assembly")
    hifiasm_dir = os.path.join(sample_dir, "hifiasm_assembly")
    
    # Generate paths to SRA downloads if necessary
    if pd.notna(ont_sra) and pd.isna(ont_raw_reads):
        ont_raw_reads = os.path.join(species_dir, "ONT", f"{ont_sra}.fastq")
    if pd.notna(illumina_sra) and pd.isna(illumina_f_raw_reads) and pd.isna(illumina_r_raw_reads):
        illumina_f_raw_reads = os.path.join(species_dir, "Illumina", f"{illumina_sra}_1.fastq")
        illumina_r_raw_reads = os.path.join(species_dir, "Illumina", f"{illumina_sra}_2.fastq")    
    if pd.notna(pacbio_sra) and pd.isna(pacbio_raw_reads):
        pacbio_raw_reads = os.path.join(species_dir, "PacBio", f"{pacbio_sra}.fastq")
    if pd.notna(ref_seq_gca) and pd.isna(ref_seq):
        ref_seq = os.path.join(species_dir, "RefSeq", f"{species_id}_{ref_seq_gca}_RefSeq.fasta")
    
    # Only include these specific NanoStats metrics in the summary table:
    inat_metrics = ["Collector",
                    "Observation Date",
                    "Latitude",
                    "Longitude",
                    "Country",
                    "State",
                    "County"]
    ont_metrics = ["Mean read length",
                   "Mean read quality",
                   "Median read length",
                   "Median read quality",
                   "Number of reads",
                   "Read length N50",
                   "STDEV read length",
                   "Total bases",
                   ">Q10",
                   ">Q15",
                   ">Q20",
                   ">Q25",
                   ">Q30"]
    illumina_metrics = ["Filename",
                        "File type",
                        "Encoding",
                        "Total Sequences",
                        "Total Bases",
                        "Sequences flagged as poor quality",
                        "Sequence Length",
                        "%GC"]
    pacbio_metrics = ["Mean read length",
                      "Mean read quality",
                      "Median read length",
                      "Median read quality",
                      "Number of reads",
                      "Read length N50",
                      "STDEV read length",
                      "Total bases",
                      ">Q10",
                      ">Q15",
                      ">Q20",
                      ">Q25",
                      ">Q30"]
    quast_metrics = ["Number of contigs (>= 0 bp)",
                     "Number of contigs (>= 1000 bp)",
                     "Number of contigs (>= 5000 bp)",
                     "Number of contigs (>= 10000 bp)",
                     "Number of contigs (>= 25000 bp)",
                     "Number of contigs (>= 50000 bp)",
                     "Total length (>= 0 bp)",
                     "Total length (>= 1000 bp)",
                     "Total length (>= 5000 bp)",
                     "Total length (>= 10000 bp)",
                     "Total length (>= 25000 bp)",
                     "Total length (>= 50000 bp)",
                     "Number of contigs",
                     "Largest contig",
                     "Total length",
                     "GC (%)",
                     "N50",
                     "N90",
                     "auN",
                     "L50",
                     "L90",
                     "Number of N's per 100 kbp"]
    busco_metrics = ["Single", "Duplicated", "Fragmented", "Missing"]    
    
    # If iNaturalist ID is provided, collect add that data 
    if pd.notna(inat_id):
        inat_html = f"https://www.inaturalist.org/observations/{int(inat_id)}"
        collector, observation_date, latitude, longitude, original_photo = get_inat_obs(inat_id)
        country, state, county = reverse_geocode(latitude, longitude)
    
        # Build the DataFrame with inat_metrics and values
        inat_df = pd.DataFrame({"Metric": inat_metrics,
                                "Value": [collector,
                                          observation_date,
                                          latitude,
                                          longitude,
                                          country,
                                          state,
                                          county]})
        sample_stats_dict["INAT_DF"] = inat_df
        sample_stats_dict["INAT_HTML"] = inat_html   
        sample_stats_dict["INAT_PHOTO"] = original_photo   
    
    else:
        print(f"SKIP:\tNo iNaturalist ID provided {inat_id}.")
    
    
    # Scan for *_ONT_nanoplot_analysis folders, pick up HTML + NanoStats
    if os.path.isdir(ont_dir):
        nanoplot_dirs = glob.glob(os.path.join(ont_dir, "*_ONT_nanoplot_analysis"))
        for ndir in nanoplot_dirs:
            # run_label will be "Raw", "Filtered", or "Corrected" (depending on your folder names)
            run_label = os.path.basename(ndir).split("_ONT")[0]
    
            # Find the NanoPlot HTML report
            html_matches = glob.glob(os.path.join(ndir, "*NanoPlot-report.html"))
            if html_matches:
                sample_stats_dict[f"{run_label.upper()}_ONT_NANOPLOT_HTML"] = html_matches[0]
    
            # Find and parse the NanoStats text summary
            stats_matches = glob.glob(os.path.join(ndir, "*NanoStats.txt"))
            if stats_matches:
                df_stats = parse_nanostats(stats_matches[0])
                # Store the DataFrame for downstream report generation
                sample_stats_dict[f"{run_label.upper()}_ONT_NANOSTATS_DF"] = df_stats
    
                # (Optional) print a markdown‐style table to your console/log
                print(f"\n=== {run_label} ONT NanoStats ===")
                print(df_stats.to_markdown(index=False))
    else:
        print(f"SKIP:\tNo ONT directory found at {ont_dir}.")
    
    
    # Scan for fastqc_results (Raw) and dedup_fastqc_results (Dedup)
    if os.path.isdir(illumina_dir):
        fastqc_dirs = glob.glob(os.path.join(illumina_dir, "*fastqc_results"))
        for ndir in fastqc_dirs:
            base = os.path.basename(ndir)
            if base == "fastqc_results":
                run_label = "RAW"
                f_pattern = "*_1_fastqc.html"
                r_pattern = "*_2_fastqc.html"
            elif base == "dedup_fastqc_results":
                run_label = "DEDUP"
                f_pattern = "*_forward_dedup_fastqc.html"
                r_pattern = "*_reverse_dedup_fastqc.html"
            else:
                # skip any other folder
                continue
    
            # Forward read
            f_html_matches = glob.glob(os.path.join(ndir, f_pattern))
            if f_html_matches:
                f_html = f_html_matches[0]
                # store HTML link
                sample_stats_dict[f"{run_label}_ILLU_F_FASTQC_HTML"] = f_html
                # parse Basic Statistics and store DataFrame
                df_f = parse_fastqc_basic_stats(f_html)
                sample_stats_dict[f"{run_label}_ILLU_F_BASIC_STATS_DF"] = df_f
    
            # Reverse read
            r_html_matches = glob.glob(os.path.join(ndir, r_pattern))
            if r_html_matches:
                r_html = r_html_matches[0]
                sample_stats_dict[f"{run_label}_ILLU_R_FASTQC_HTML"] = r_html
                df_r = parse_fastqc_basic_stats(r_html)
                sample_stats_dict[f"{run_label}_ILLU_R_BASIC_STATS_DF"] = df_r
    else:
        print(f"SKIP:\tNo Illumina directory found at {illumina_dir}.")
    
    
    # Scan for *_PacBio_nanoplot_analysis folders, pick up HTML + NanoStats
    if os.path.isdir(pacbio_dir):
        nanoplot_dirs = glob.glob(os.path.join(pacbio_dir, "*_PacBio_nanoplot_analysis"))
        for ndir in nanoplot_dirs:
            # run_label will be "Raw" or "Filtered" (depending on your folder names)
            run_label = os.path.basename(ndir).split("_PacBio")[0]
    
            # Find the NanoPlot HTML report
            html_matches = glob.glob(os.path.join(ndir, "*NanoPlot-report.html"))
            if html_matches:
                sample_stats_dict[f"{run_label.upper()}_PACBIO_NANOPLOT_HTML"] = html_matches[0]
    
            # Find and parse the NanoStats text summary
            stats_matches = glob.glob(os.path.join(ndir, "*NanoStats.txt"))
            if stats_matches:
                df_stats = parse_nanostats(stats_matches[0])
                # Store the DataFrame for downstream report generation
                sample_stats_dict[f"{run_label.upper()}_PACBIO_NANOSTATS_DF"] = df_stats
    
                # (Optional) print a markdown‐style table to your console/log
                print(f"\n=== {run_label} PacBio NanoStats ===")
                print(df_stats.to_markdown(index=False))
    else:
        print(f"SKIP:\tNo PacBio directory found at {pacbio_dir}.")
    
    
    # Scan for assembly folders, pick up QUAST metrics + link, then BUSCO metrics
    assembly_candidates = [masurca_dir, flye_dir, spades_dir, hifiasm_dir]
    for assembly_dir in assembly_candidates:    
        if os.path.isdir(assembly_dir):
            asm_type = os.path.basename(assembly_dir).split("_")[0].upper()
            assembly_fasta = os.path.join(assembly_dir, f"{sample_id}_{asm_type.lower()}.fasta")
            
            # Process QUAST
            quast_folders = glob.glob(os.path.join(assembly_dir, "*_quast"))
            for quast_dir in quast_folders:
                try:
                    quast_df, quast_html = process_quast_folder(quast_dir)
                    sample_stats_dict[f"{asm_type}_QUAST_DF"] = quast_df
                    sample_stats_dict[f"{asm_type}_QUAST_HTML"] = quast_html
                except FileNotFoundError as e:
                    print(f"QUAST processing failed for {quast_dir}: {e}")
            
            # Process BUSCO
            busco_list = [first_busco_db, second_busco_db]
            for busco_index, busco_db in enumerate(busco_list):
                busco_type = "FIRST" if busco_index == 0 else "SECOND"
                busco_dirs = glob.glob(os.path.join(assembly_dir, f"*_{busco_db}_busco"))
                if not busco_dirs:
                    busco_dirs = glob.glob(os.path.join(assembly_dir, f"*_{busco_db}_compleasm"))
                    print(f"TODO: Process Compleasm files for {busco_db} in {assembly_dir}")
                    continue
                
                for busco_dir in busco_dirs:
                    # Find BUSCO summary
                    summary_files = glob.glob(os.path.join(busco_dir, "short_summary*.txt"))
                    if not summary_files:
                        print(f"No BUSCO summary found in {busco_dir}")
                        continue
                    
                    # Parse summary
                    metrics = {}
                    with open(summary_files[0], "r") as f:
                        for line in f:
                            if "%" in line and "C:" in line:
                                # Example: C:99.0%[S:98.5%,D:0.5%],F:0.5%,M:0.5%,n:1000
                                parts = line.strip().split(",")
                                metrics["Complete"] = parts[0].split("[")[0].replace("C:", "")
                                metrics["Single"] = parts[0].split("[")[1].replace("S:", "")
                                metrics["Duplicated"] = parts[1].replace("D:", "").replace("]", "")
                                metrics["Fragmented"] = parts[2].replace("F:", "")
                                metrics["Missing"] = parts[3].replace("M:", "")
                                metrics["Total"] = parts[4].replace("n:", "")
                                break
    
                    # Convert metrics to DataFrame for sample_stats_dict
                    metrics_df = pd.DataFrame(list(metrics.items()), columns=["Metric", "Value"])
            
                    # Find SVG
                    svg_files = glob.glob(os.path.join(assembly_dir, f"*{busco_db}_busco.svg"))
                    svg_path = svg_files[0] if svg_files else ""
            
                    # Read BUSCO genes table (if exists)
                    genes_file = os.path.join(busco_dir, f"run_{busco_db}_odb12", "full_table.tsv")
                    busco_genes_df = pd.DataFrame()
                    if os.path.exists(genes_file):
                        try:
                            busco_genes_df = pd.read_csv(genes_file, sep="\t", comment="#")
                            # Ensure expected columns are present, set defaults if missing
                            expected_columns = ["Busco id", "Status", "Sequence", "Gene Start", 
                                              "Gene End", "Strand", "Score", "Length"]
                            for col in expected_columns:
                                if col not in busco_genes_df.columns:
                                    busco_genes_df[col] = "-"
                        except Exception as e:
                            print(f"Failed to parse BUSCO genes table {genes_file}: {e}")
            
                    # Generate BUSCO HTML report
                    busco_data = [{"assembly_fasta": assembly_fasta,
                                   "busco_db": busco_db,
                                   "svg_path": svg_path,
                                   "metrics": metrics,  # Dictionary for summary table
                                   "busco_genes_df": busco_genes_df}]  # DataFrame for genes table
            
                    busco_html = busco_tmpl.render(assembly_name = os.path.basename(assembly_fasta),
                                                   generated_time = now,
                                                   busco_data = busco_data)
                    outpath = os.path.join(busco_dir, f"{sample_id}_{busco_db}_EGAP_busco.html")
                    with open(outpath, "w") as fh:
                        fh.write(busco_html)
                    
                    # Store in sample_stats_dict
                    sample_stats_dict[f"{asm_type}_{busco_type}_BUSCO_DF"] = metrics_df
                    sample_stats_dict[f"{asm_type}_{busco_type}_BUSCO_HTML"] = outpath
        else:
            print(f"SKIP:\tNo Assembly directory found at {assembly_dir}.")
    
    # If ONLY Reference sequence is provided (noted by the lack of est_size), use that as assembly
    if pd.notna(ref_seq) and pd.isna(est_size):
        assembly_fasta = ref_seq
    else:
        assembly_fasta = os.path.join(sample_dir, f"{sample_id}_final_EGAP_assembly.fasta")
    asm_type = "FINAL"

    print(f"DEBUG - assembly_fasta - {assembly_fasta}")
    
    # Process QUAST
    final_quast_dir = assembly_fasta.replace(".fasta","_quast")
    
    print(f"DEBUG - final_quast_dir - {final_quast_dir}")
    
    try:
        final_df, final_html = process_quast_folder(final_quast_dir)
        sample_stats_dict["FINAL_QUAST_DF"]   = final_df
        sample_stats_dict["FINAL_QUAST_HTML"] = final_html
    except FileNotFoundError as e:
        print(f"FINAL QUAST processing failed for {quast_dir}: {e}")
    
    # Process BUSCO
    busco_list = [first_busco_db, second_busco_db]
    for busco_index, busco_db in enumerate(busco_list):
        busco_type = "FIRST" if busco_index == 0 else "SECOND"
        busco_dir = assembly_fasta.replace(".fasta",f"_{busco_db}_busco")
        if not os.path.exists(busco_dir):
            busco_dirs = glob.glob(os.path.join(assembly_dir, f"*_{busco_db}_compleasm"))
            print(f"TODO: Process Compleasm files for {busco_db} in {assembly_dir}")
            continue
        
        # Find BUSCO summary
        summary_files = glob.glob(os.path.join(busco_dir, "short_summary*.txt"))
        if not summary_files:
            print(f"No BUSCO summary found in {busco_dir}")
            continue
        
        # Parse summary
        metrics = {}
        with open(summary_files[0], "r") as f:
            for line in f:
                if "%" in line and "C:" in line:
                    # Example: C:99.0%[S:98.5%,D:0.5%],F:0.5%,M:0.5%,n:1000
                    parts = line.strip().split(",")
                    metrics["Complete"] = parts[0].split("[")[0].replace("C:", "")
                    metrics["Single"] = parts[0].split("[")[1].replace("S:", "")
                    metrics["Duplicated"] = parts[1].replace("D:", "").replace("]", "")
                    metrics["Fragmented"] = parts[2].replace("F:", "")
                    metrics["Missing"] = parts[3].replace("M:", "")
                    metrics["Total"] = parts[4].replace("n:", "")
                    break

        # Convert metrics to DataFrame for sample_stats_dict
        metrics_df = pd.DataFrame(list(metrics.items()), columns=["Metric", "Value"])

        # Find SVG
        svg_path = assembly_fasta.replace(".fasta", f"_final_{busco_db}_busco.svg")
        print(f"DEBUG - svg_path - {svg_path}")

        # Read BUSCO genes table (if exists)
        genes_file = os.path.join(busco_dir, f"run_{busco_db}_odb12", "full_table.tsv")
        busco_genes_df = pd.DataFrame()
        if os.path.exists(genes_file):
            try:
                busco_genes_df = pd.read_csv(genes_file, sep="\t", comment="#")

                # Ensure expected columns are present, set defaults if missing
                expected_columns = ["Busco id", "Status", "Sequence", "Gene Start", 
                                  "Gene End", "Strand", "Score", "Length"]
                for col in expected_columns:
                    if col not in busco_genes_df.columns:
                        busco_genes_df[col] = "-"
            except Exception as e:
                print(f"Failed to parse BUSCO genes table {genes_file}: {e}")

        # Generate BUSCO HTML report
        busco_data = [{"assembly_fasta": assembly_fasta,
                       "busco_db": busco_db.capitalize(),
                       "svg_path": svg_path,
                       "metrics": metrics,  # Dictionary for summary table
                       "busco_genes_df": busco_genes_df}]  # DataFrame for genes table
        busco_html = busco_tmpl.render(assembly_name = os.path.basename(assembly_fasta),
                                       generated_time = now,
                                       busco_data = busco_data,
                                       busco_db = busco_db.capitalize())
        outpath = os.path.join(busco_dir, f"{sample_id}_{busco_db}_EGAP_busco.html")
        with open(outpath, "w") as fh:
            fh.write(busco_html)
        
        # Store in sample_stats_dict
        sample_stats_dict[f"FINAL_{busco_type}_BUSCO_DF"] = metrics_df
        sample_stats_dict[f"FINAL_{busco_type}_BUSCO_HTML"] = outpath
    
    
    # Pull out DataFrames and HTML links
    inat_df = sample_stats_dict.get("INAT_DF")
    inat_link = sample_stats_dict.get("INAT_HTML")
    inat_thumbnail = sample_stats_dict.get("INAT_THUMBNAIL")
    inat_photo = sample_stats_dict.get("INAT_PHOTO")
    
    ont_raw_df = sample_stats_dict.get("RAW_ONT_NANOSTATS_DF")
    ont_filt_df = sample_stats_dict.get("FILT_ONT_NANOSTATS_DF")
    ont_corr_df = sample_stats_dict.get("CORR_ONT_NANOSTATS_DF")
    ont_raw_link = sample_stats_dict.get("RAW_ONT_NANOPLOT_HTML", "")
    ont_filt_link = sample_stats_dict.get("FILT_ONT_NANOPLOT_HTML", "")
    ont_corr_link = sample_stats_dict.get("CORR_ONT_NANOPLOT_HTML", "")
    
    illumina_raw_f_df = sample_stats_dict.get("RAW_ILLU_F_BASIC_STATS_DF")
    if illumina_raw_f_df is None:
        raise KeyError("Expected key RAW_ILLU_F_BASIC_STATS_DF not found in sample_stats_dict.\n"
                       f"Available keys: {list(sample_stats_dict.keys())}")
    illumina_raw_r_df = sample_stats_dict.get("RAW_ILLU_R_BASIC_STATS_DF")
    if illumina_raw_r_df is None:
        raise KeyError("Expected key RAW_ILLU_R_BASIC_STATS_DF not found in sample_stats_dict.")
    illumina_raw_f_link = sample_stats_dict.get("RAW_ILLU_F_FASTQC_HTML", "")
    illumina_raw_r_link = sample_stats_dict.get("RAW_ILLU_R_FASTQC_HTML", "")
    
    illumina_dedup_f_df = sample_stats_dict.get("DEDUP_ILLU_F_BASIC_STATS_DF")
    if illumina_dedup_f_df is None:
        raise KeyError("Expected key DEDUP_ILLU_F_BASIC_STATS_DF not found in sample_stats_dict.\n"
                       f"Available keys: {list(sample_stats_dict.keys())}")
    illumina_dedup_r_df = sample_stats_dict.get("DEDUP_ILLU_R_BASIC_STATS_DF")
    if illumina_dedup_r_df is None:
        raise KeyError("Expected key DEDUP_ILLU_R_BASIC_STATS_DF not found in sample_stats_dict.")
    illumina_dedup_f_link = sample_stats_dict.get("DEDUP_ILLU_F_FASTQC_HTML", "")
    illumina_dedup_r_link = sample_stats_dict.get("DEDUP_ILLU_R_FASTQC_HTML", "")
    
    pacbio_raw_df = sample_stats_dict.get("RAW_PACBIO_NANOSTATS_DF")
    pacbio_filt_df = sample_stats_dict.get("FILT_PACBIO_NANOSTATS_DF")
    pacbio_raw_link = sample_stats_dict.get("RAW_PACBIO_NANOPLOT_HTML", "")
    pacbio_filt_link = sample_stats_dict.get("FILT_PACBIO_NANOPLOT_HTML", "")
    
    masurca_quast_df = sample_stats_dict.get("MASURCA_QUAST_DF")
    masurca_first_busco_df = sample_stats_dict.get("MASURCA_FIRST_BUSCO_DF")
    masurca_second_busco_df = sample_stats_dict.get("MASURCA_SECOND_BUSCO_DF")
    masurca_quast_link = sample_stats_dict.get("MASURCA_QUAST_HTML", "")
    masurca_first_busco_link = sample_stats_dict.get("MASURCA_FIRST_BUSCO_HTML", "")
    masurca_second_busco_link = sample_stats_dict.get("MASURCA_SECOND_BUSCO_HTML", "")
    
    flye_quast_df = sample_stats_dict.get("FLYE_QUAST_DF")
    flye_first_busco_df = sample_stats_dict.get("FLYE_FIRST_BUSCO_DF")
    flye_second_busco_df = sample_stats_dict.get("FLYE_SECOND_BUSCO_DF")
    flye_quast_link = sample_stats_dict.get("FLYE_QUAST_HTML", "")
    flye_first_busco_link = sample_stats_dict.get("FLYE_FIRST_BUSCO_HTML", "")
    flye_second_busco_link = sample_stats_dict.get("FLYE_SECOND_BUSCO_HTML", "")
    
    spades_quast_df = sample_stats_dict.get("SPADES_QUAST_DF")
    spades_first_busco_df = sample_stats_dict.get("SPADES_FIRST_BUSCO_DF")
    spades_second_busco_df = sample_stats_dict.get("SPADES_SECOND_BUSCO_DF")
    spades_quast_link = sample_stats_dict.get("SPADES_QUAST_HTML", "")
    spades_first_busco_link = sample_stats_dict.get("SPADES_FIRST_BUSCO_HTML", "")
    spades_second_busco_link = sample_stats_dict.get("SPADES_SECOND_BUSCO_HTML", "")
    
    hifiasm_quast_df = sample_stats_dict.get("HIFIASM_QUAST_DF")
    hifiasm_first_busco_df = sample_stats_dict.get("HIFIASM_FIRST_BUSCO_DF")
    hifiasm_second_busco_df = sample_stats_dict.get("HIFIASM_SECOND_BUSCO_DF")
    hifiasm_quast_link = sample_stats_dict.get("HIFIASM_QUAST_HTML", "")
    hifiasm_first_busco_link = sample_stats_dict.get("HIFIASM_FIRST_BUSCO_HTML", "")
    hifiasm_second_busco_link = sample_stats_dict.get("HIFIASM_SECOND_BUSCO_HTML", "")
    
    final_quast_df = sample_stats_dict.get("FINAL_QUAST_DF")
    final_first_busco_df = sample_stats_dict.get("FINAL_FIRST_BUSCO_DF")
    final_second_busco_df = sample_stats_dict.get("FINAL_SECOND_BUSCO_DF")
    final_quast_link = sample_stats_dict.get("FINAL_QUAST_HTML", "")
    final_first_busco_link = sample_stats_dict.get("FINAL_FIRST_BUSCO_HTML", "")
    final_second_busco_link = sample_stats_dict.get("FINAL_SECOND_BUSCO_HTML", "")
    
    # Turn each DataFrame into a dict Metric→Value
    inat_map = dict(zip(inat_df["Metric"], inat_df["Value"])) if inat_df is not None else {}
    
    ont_raw_map = dict(zip(ont_raw_df["Metric"], ont_raw_df["Value"])) if ont_raw_df is not None else {}
    ont_filt_map = dict(zip(ont_filt_df["Metric"], ont_filt_df["Value"])) if ont_filt_df is not None else {}
    ont_corr_map = dict(zip(ont_corr_df["Metric"], ont_filt_df["Value"])) if ont_corr_df is not None else {}
    
    illumina_raw_f_map = dict(zip(illumina_raw_f_df["Metric"], illumina_raw_f_df["Value"])) if illumina_raw_f_df is not None else {}
    illumina_raw_r_map = dict(zip(illumina_raw_r_df["Metric"], illumina_raw_r_df["Value"])) if illumina_raw_r_df is not None else {}
    illumina_dedup_f_map = dict(zip(illumina_dedup_f_df["Metric"], illumina_dedup_f_df["Value"])) if illumina_dedup_f_df is not None else {}
    illumina_dedup_r_map = dict(zip(illumina_dedup_r_df["Metric"], illumina_dedup_r_df["Value"])) if illumina_dedup_r_df is not None else {}
    
    pacbio_raw_map = dict(zip(pacbio_raw_df["Metric"], pacbio_raw_df["Value"])) if pacbio_raw_df is not None else {}
    pacbio_filt_map = dict(zip(pacbio_filt_df["Metric"], pacbio_filt_df["Value"])) if pacbio_filt_df is not None else {}
    
    masurca_quast_map = dict(zip(masurca_quast_df.index, masurca_quast_df["Value"])) if masurca_quast_df is not None else {}
    masurca_first_busco_map = dict(zip(masurca_first_busco_df["Metric"], masurca_first_busco_df["Value"])) if masurca_first_busco_df is not None else {}
    masurca_second_busco_map = dict(zip(masurca_second_busco_df["Metric"], masurca_second_busco_df["Value"])) if masurca_second_busco_df is not None else {}
    
    flye_quast_map = dict(zip(flye_quast_df.index, flye_quast_df["Value"])) if flye_quast_df is not None else {}
    flye_first_busco_map = dict(zip(flye_first_busco_df["Metric"], flye_first_busco_df["Value"])) if flye_first_busco_df is not None else {}
    flye_second_busco_map = dict(zip(flye_second_busco_df["Metric"], flye_second_busco_df["Value"])) if flye_second_busco_df is not None else {}
    
    spades_quast_map = dict(zip(spades_quast_df.index, spades_quast_df["Value"])) if spades_quast_df is not None else {}
    spades_first_busco_map = dict(zip(spades_first_busco_df["Metric"], spades_first_busco_df["Value"])) if spades_first_busco_df is not None else {}
    spades_second_busco_map = dict(zip(spades_second_busco_df["Metric"], spades_second_busco_df["Value"])) if spades_second_busco_df is not None else {}
    
    hifiasm_quast_map = dict(zip(hifiasm_quast_df.index, hifiasm_quast_df["Value"])) if hifiasm_quast_df is not None else {}
    hifiasm_first_busco_map = dict(zip(hifiasm_first_busco_df["Metric"], hifiasm_first_busco_df["Value"])) if hifiasm_first_busco_df is not None else {}
    hifiasm_second_busco_map = dict(zip(hifiasm_second_busco_df["Metric"], hifiasm_second_busco_df["Value"])) if hifiasm_second_busco_df is not None else {}
    
    final_quast_map = dict(zip(final_quast_df.index, final_quast_df["Value"])) if final_quast_df is not None else {}
    final_first_busco_map = dict(zip(final_first_busco_df["Metric"], final_first_busco_df["Value"])) if final_first_busco_df is not None else {}
    final_second_busco_map = dict(zip(final_second_busco_df["Metric"], final_second_busco_df["Value"])) if final_second_busco_df is not None else {}
        
    # Render HTML file
    render_ctx = {"sample_id": sample_id,
                  "generated_time": now,
                  "inat_metrics": inat_metrics,
                  "inat_map": inat_map,
                  "inat_link": inat_link,
                  "inat_thumbnail": inat_thumbnail,
                  "inat_photo": inat_photo,
                  "ont_metrics": ont_metrics,
                  "ont_raw_map": ont_raw_map,
                  "ont_filt_map": ont_filt_map,
                  "ont_corr_map": ont_corr_map,
                  "ont_raw_link": ont_raw_link,
                  "ont_filtered_link": ont_filt_link,
                  "ont_corrected_link": ont_corr_link,
                  "illumina_metrics": illumina_metrics,
                  "illumina_raw_f_map": illumina_raw_f_map,
                  "illumina_raw_r_map": illumina_raw_r_map,
                  "illumina_dedup_f_map": illumina_dedup_f_map,
                  "illumina_dedup_r_map": illumina_dedup_r_map,
                  "illumina_raw_f_link": illumina_raw_f_link,
                  "illumina_raw_r_link": illumina_raw_r_link,
                  "illumina_dedup_f_link": illumina_dedup_f_link,
                  "illumina_dedup_r_link": illumina_dedup_r_link,
                  "pacbio_metrics": pacbio_metrics,
                  "pacbio_raw_map": pacbio_raw_map,
                  "pacbio_filt_map": pacbio_filt_map,
                  "pacbio_raw_link": pacbio_raw_link,
                  "pacbio_filtered_link": pacbio_filt_link,
                  "quast_metrics": quast_metrics,
                  "busco_metrics": busco_metrics,
                  "masurca_quast_map": masurca_quast_map,
                  "masurca_first_busco_map": masurca_first_busco_map,
                  "masurca_second_busco_map": masurca_second_busco_map,
                  "masurca_quast_link": masurca_quast_link,
                  "masurca_first_busco_link": masurca_first_busco_link,
                  "masurca_second_busco_link": masurca_second_busco_link,
                  "flye_quast_map": flye_quast_map,
                  "flye_first_busco_map": flye_first_busco_map,
                  "flye_second_busco_map": flye_second_busco_map,
                  "flye_quast_link": flye_quast_link,
                  "flye_first_busco_link": flye_first_busco_link,
                  "flye_second_busco_link": flye_second_busco_link,
                  "spades_quast_map": spades_quast_map,
                  "spades_first_busco_map": spades_first_busco_map,
                  "spades_second_busco_map": spades_second_busco_map,
                  "spades_quast_link": spades_quast_link,
                  "spades_first_busco_link": spades_first_busco_link,
                  "spades_second_busco_link": spades_second_busco_link,
                  "hifiasm_quast_map": hifiasm_quast_map,
                  "hifiasm_first_busco_map": hifiasm_first_busco_map,
                  "hifiasm_second_busco_map": hifiasm_second_busco_map,
                  "hifiasm_quast_link": hifiasm_quast_link,
                  "hifiasm_first_busco_link": hifiasm_first_busco_link,
                  "hifiasm_second_busco_link": hifiasm_second_busco_link,
                  "final_quast_map": final_quast_map,
                  "final_first_busco_map": final_first_busco_map,
                  "final_second_busco_map": final_second_busco_map,
                  "final_quast_link": final_quast_link,
                  "final_first_busco_link": final_first_busco_link,
                  "final_second_busco_link": final_second_busco_link}
    html = main_tmpl.render(**render_ctx)
    html_report_outpath = os.path.join(sample_dir, f"{sample_id}_EGAP_summary.html")
    with open(html_report_outpath, "w") as fh:
        fh.write(html)

    print(f"PASS:\tWrote report: {html_report_outpath }")
    
    return html_report_outpath


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python3 assemble_flye.py <sample_id> <input_csv> "
            "<output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)
        
    html_report_outpath = html_reporter(sys.argv[1],       # sample_id
                                        sys.argv[2],       # input_csv
                                        sys.argv[3],       # output_dir
                                        str(sys.argv[4]),  # cpu_threads
                                        str(sys.argv[5]))  # ram_gb