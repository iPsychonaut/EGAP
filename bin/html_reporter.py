#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
html_reporter.py

Generates an EGAP HTML summary using Jinja templates. Robust to missing
artifacts (FastQC / NanoPlot / QUAST / BUSCO), CWD-safe, and tolerant of
metadata quirks.

Created on Wed Aug 16 2023

Updated on Wed Sept 3 2025

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os, glob, sys
from pathlib import Path
from datetime import datetime
import urllib.request
from typing import Optional  # <-- Py3.8/3.9-safe Optional

import pandas as pd
from bs4 import BeautifulSoup
from jinja2 import Environment, FileSystemLoader, select_autoescape

# Optional iNat import (process_metadata may fail due to utilities.calculate_genome_coverage)
_HAS_INAT = True
try:
    from process_metadata import get_inat_obs, reverse_geocode
except Exception as e:
    print(f"WARN:\tprocess_metadata import failed ({e}); iNaturalist fields will be skipped.")
    _HAS_INAT = False
    def get_inat_obs(_): return (None, None, None, None, None)
    def reverse_geocode(_, __): return (None, None, None)

# Use the shared helper from utilities so row extraction stays consistent
try:
    from utilities import get_current_row_data
except Exception as e:
    print(f"ERROR:\tutilities.get_current_row_data not importable: {e}")
    raise


# ------------------------------ Helpers ------------------------------
REQUIRED_TEMPLATES = ("EGAP_summary.html", "EGAP_busco.html")
RAW_URLS = {
    "EGAP_summary.html": "https://raw.githubusercontent.com/iPsychonaut/EGAP/main/resources/templates/EGAP_summary.html",
    "EGAP_busco.html"  : "https://raw.githubusercontent.com/iPsychonaut/EGAP/main/resources/templates/EGAP_busco.html",
}

def _abs(p):
    return os.path.abspath(p) if isinstance(p, str) else p

def _has_both_templates(dir_path: Path) -> bool:
    return all((dir_path / fname).exists() for fname in REQUIRED_TEMPLATES)

def _download_templates_to(dst_dir: Path) -> bool:
    """Try to download both templates into dst_dir. Return True if both exist afterward."""
    dst_dir.mkdir(parents=True, exist_ok=True)
    ok = True
    for name, url in RAW_URLS.items():
        out = dst_dir / name
        try:
            print(f"INFO:\tDownloading template '{name}' -> {out}")
            urllib.request.urlretrieve(url, out.as_posix())
        except Exception as e:
            ok = False
            print(f"WARN:\tFailed to download {name} from {url}: {e}")
    return ok and _has_both_templates(dst_dir)

def _resolve_templates_dir(output_dir_abs: str) -> Optional[Path]:
    """
    Robustly locate a folder containing BOTH 'EGAP_summary.html' and 'EGAP_busco.html'.
    If not found, attempt to download both files into the SAME directory as html_reporter.py.
    If download fails, return None (caller will stub out the report and continue).

    Search order (first folder that contains BOTH wins):
      1) $EGAP_TEMPLATES_DIR
      2) $EGAP_PROJECT_DIR/resources/templates
      3) <script>/../resources/templates
      4) <script>/../../resources/templates
      5) <sys.prefix>/share/egap/templates  (conda share)
      6) <sys.prefix>/share/EGAP/templates
      7) site-packages scans: */EGAP/resources/templates, */egap/resources/templates
      8) <output_dir>/resources/templates
      9) If none above found, download both templates into <script_dir> and use it.
    """
    script_file = Path(__file__).resolve()
    script_dir  = script_file.parent
    sys_prefix  = Path(sys.prefix).resolve()
    output_dir_p = Path(output_dir_abs).resolve()

    tried = []

    def _check_and_return(p: Path):
        tried.append(p)
        if _has_both_templates(p):
            return p
        return None

    # 1) Explicit override
    env_templates = os.environ.get("EGAP_TEMPLATES_DIR")
    if env_templates:
        p = Path(env_templates).expanduser().resolve()
        found = _check_and_return(p)
        if found: return found

    # 2) Project dir provided
    env_project = os.environ.get("EGAP_PROJECT_DIR")
    if env_project:
        p = Path(env_project).expanduser().resolve() / "resources" / "templates"
        found = _check_and_return(p)
        if found: return found

    # 3–4) Repo layouts near the script
    for p in [
        script_dir.parent / "resources" / "templates",
        script_dir.parent.parent / "resources" / "templates",
    ]:
        found = _check_and_return(p)
        if found: return found

    # 5–6) Conda share
    for share_name in ("egap", "EGAP"):
        p = sys_prefix / "share" / share_name / "templates"
        found = _check_and_return(p)
        if found: return found

    # 7) site-packages scan
    for sp in [Path(p) for p in sys.path if isinstance(p, str)]:
        for p in [
            sp / "EGAP" / "resources" / "templates",
            sp / "egap" / "resources" / "templates",
            sp / "resources" / "templates",
        ]:
            found = _check_and_return(p)
            if found: return found

    # 8) project-local
    p = output_dir_p / "resources" / "templates"
    found = _check_and_return(p)
    if found: return found

    # 9) Download into the SAME directory that html_reporter.py is in
    print("INFO:\tTemplates not found in standard locations; attempting download into script directory.")
    if _download_templates_to(script_dir):
        print(f"PASS:\tTemplates downloaded to {script_dir}")
        return script_dir

    # Still nothing: inform and continue (return None; caller will stub).
    msg = "Templates directory with BOTH required files not found.\nTried:\n  " + "\n  ".join(str(x) for x in tried)
    msg += ("\nAlso attempted download into script directory but failed.\n"
            "A report cannot be generated unless 'EGAP_summary.html' and 'EGAP_busco.html' are available.")
    print("WARN:\t" + msg)
    return None


def parse_nanostats(nanostats_path):
    """Read a NanoStats.txt and return DataFrame[Metric, Value]."""
    rows = []
    with open(nanostats_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or ":" not in line:
                continue
            key, val = line.split(":", 1)
            rows.append({"Metric": key.strip(), "Value": val.strip()})
    return pd.DataFrame(rows)

def parse_fastqc_basic_stats(html_path):
    """Parse 'Basic Statistics' table from a FastQC HTML report."""
    with open(html_path, "r", encoding="utf-8") as fh:
        soup = BeautifulSoup(fh, "html.parser")
    header = soup.find(lambda tag: tag.name == "h2" and "Basic Statistics" in tag.get_text())
    if not header:
        raise ValueError(f"No Basic Statistics header found in {html_path}")
    table = header.find_next("table")
    if not table:
        raise ValueError(f"No table after Basic Statistics header in {html_path}")
    rows = []
    for tr in table.find_all("tr"):
        tds = tr.find_all("td")
        if len(tds) == 2:
            rows.append({"Metric": tds[0].get_text(strip=True),
                         "Value":  tds[1].get_text(strip=True)})
    return pd.DataFrame(rows)

def process_quast_folder(quast_dir):
    """Return (df, html_path) from a QUAST <assembly>_quast folder."""
    tsv = os.path.join(quast_dir, "report.tsv")
    html = os.path.join(quast_dir, "report.html")
    if not os.path.exists(tsv):
        raise FileNotFoundError(f"No report.tsv in {quast_dir}")
    if not os.path.exists(html):
        raise FileNotFoundError(f"No report.html in {quast_dir}")
    df = pd.read_csv(tsv, sep="\t")
    df.iloc[:, 0] = df.iloc[:, 0].str.replace("#", "Number of", regex=False)
    df = df.set_index(df.columns[0])
    df.columns = ["Value"]
    return df, html

def _parse_busco_summary_line(line: str) -> dict:
    """
    Parse a BUSCO short_summary line like:
    "C:97.5%[S:97.3%,D:0.2%],F:0.5%,M:2.0%,n:124"
    Returns keys with percent values as strings (no trailing spaces).
    """
    parts = line.strip().split(",")
    m = {}
    try:
        m["Complete"]   = parts[0].split("[")[0].replace("C:", "").strip()
        m["Single"]     = parts[0].split("[")[1].replace("S:", "").strip()
        m["Duplicated"] = parts[1].replace("D:", "").replace("]", "").strip()
        m["Fragmented"] = parts[2].replace("F:", "").strip()
        m["Missing"]    = parts[3].replace("M:", "").strip()
        m["Total"]      = parts[4].replace("n:", "").strip()
    except Exception:
        pass
    return m

def _parse_compleasm_summary(path: str) -> dict:
    """
    Parse Compleasm summary.txt where S:/D:/F:/M: are on separate lines.
    Returns percent values as strings; computes Complete = S + D (rounded 2dp).
    """
    if not os.path.exists(path):
        return {}
    S = D = F = M = None
    with open(path, "r") as fh:
        for line in fh:
            if "S:" in line: S = line.split("S:")[-1].split(",")[0].replace("%","").strip()
            if "D:" in line: D = line.split("D:")[-1].split(",")[0].replace("%","").strip()
            if "F:" in line: F = line.split("F:")[-1].split(",")[0].replace("%","").strip()
            if "M:" in line: M = line.split("M:")[-1].split(",")[0].replace("%","").strip()
    if S is None or D is None:
        return {}
    try:
        C = str(round(float(S) + float(D), 2))
    except Exception:
        C = ""
    return {
        "Complete":   C,
        "Single":     S or "",
        "Duplicated": D or "",
        "Fragmented": F or "",
        "Missing":    M or "",
        "Total":      ""
    }

def _find_busco_genes_table(busco_dir: str, busco_db: str) -> str:
    """
    Return the first existing BUSCO/Compleasm 'full_table' path.
    Supports both BUSCO and Compleasm layouts.
    """
    candidates = [
        os.path.join(busco_dir, f"run_{busco_db}_odb12", "full_table.tsv"),
        os.path.join(busco_dir, f"{busco_db}_odb12", "full_table_busco_format.tsv"),
        os.path.join(busco_dir, f"{busco_db}_odb12", "full_table.tsv"),
    ]
    for p in candidates:
        if os.path.exists(p):
            return p
    return ""

# ------------------------------ Main ------------------------------
def html_reporter(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """
    Render the EGAP HTML report for a sample.
    Returns: path to <sample_dir>/<sample_id>_EGAP_summary.html (or a stub if templates are missing)
    """
    # Resolve critical paths first (CWD-safe)
    input_csv_abs = _abs(input_csv)
    output_dir_abs = _abs(output_dir)

    # Locate or fetch templates
    templates_dir = _resolve_templates_dir(output_dir_abs)

    now = datetime.now().strftime("%Y%m%d-%H:%M:%S")

    # Metadata row
    input_df = pd.read_csv(input_csv_abs)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]

    # CSV columns (tolerate missing)
    illumina_sra         = current_series.get("ILLUMINA_SRA")
    illumina_f_raw_reads = current_series.get("ILLUMINA_RAW_F_READS")
    illumina_r_raw_reads = current_series.get("ILLUMINA_RAW_R_READS")
    ont_sra              = current_series.get("ONT_SRA")
    ont_raw_reads        = current_series.get("ONT_RAW_READS")
    pacbio_sra           = current_series.get("PACBIO_SRA")
    pacbio_raw_reads     = current_series.get("PACBIO_RAW_READS")
    first_busco_db       = current_series.get("BUSCO_1")
    second_busco_db      = current_series.get("BUSCO_2")
    ref_seq_gca          = current_series.get("REF_SEQ_GCA")
    ref_seq              = current_series.get("REF_SEQ")
    species_id           = current_series.get("SPECIES_ID")
    est_size             = current_series.get("EST_SIZE")

    # Handle the iNat column being typo'd in some sheets
    inat_id = None
    for key in ("INATRUALIST_ID", "INATURALIST_ID", "INAT_ID"):
        if key in current_series and pd.notna(current_series[key]):
            inat_id = current_series[key]
            break

    # Build project dirs
    species_dir = os.path.join(output_dir_abs, str(species_id))
    sample_dir  = os.path.join(species_dir, sample_id)
    ont_dir     = os.path.join(species_dir, "ONT")
    illumina_dir= os.path.join(species_dir, "Illumina")
    pacbio_dir  = os.path.join(species_dir, "PacBio")
    masurca_dir = os.path.join(sample_dir, "masurca_assembly")
    flye_dir    = os.path.join(sample_dir, "flye_assembly")
    spades_dir  = os.path.join(sample_dir, "spades_assembly")
    hifiasm_dir = os.path.join(sample_dir, "hifiasm_assembly")

    # Derive implied read/ref paths from SRA/GCA when raw path is empty
    if pd.notna(ont_sra) and pd.isna(ont_raw_reads):
        ont_raw_reads = os.path.join(ont_dir, f"{ont_sra}.fastq")
    if pd.notna(illumina_sra) and pd.isna(illumina_f_raw_reads) and pd.isna(illumina_r_raw_reads):
        illumina_f_raw_reads = os.path.join(illumina_dir, f"{illumina_sra}_1.fastq")
        illumina_r_raw_reads = os.path.join(illumina_dir, f"{illumina_sra}_2.fastq")
    if pd.notna(pacbio_sra) and pd.isna(pacbio_raw_reads):
        pacbio_raw_reads = os.path.join(pacbio_dir, f"{pacbio_sra}.fastq")
    if pd.notna(ref_seq_gca) and pd.isna(ref_seq):
        ref_seq = os.path.join(species_dir, "RefSeq", f"{species_id}_{ref_seq_gca}_RefSeq.fasta")

    # Metric headers
    inat_metrics = ["Collector","Observation Date","Latitude","Longitude","Country","State","County"]
    ont_metrics  = ["Mean read length","Mean read quality","Median read length","Median read quality",
                    "Number of reads","Read length N50","STDEV read length","Total bases",
                    ">Q10",">Q15",">Q20",">Q25",">Q30"]
    illumina_metrics = ["Filename","File type","Encoding","Total Sequences","Total Bases",
                        "Sequences flagged as poor quality","Sequence Length","%GC"]
    pacbio_metrics = list(ont_metrics)
    quast_metrics  = ["Number of contigs (>= 0 bp)","Number of contigs (>= 1000 bp)",
                      "Number of contigs (>= 5000 bp)","Number of contigs (>= 10000 bp)",
                      "Number of contigs (>= 25000 bp)","Number of contigs (>= 50000 bp)",
                      "Total length (>= 0 bp)","Total length (>= 1000 bp)","Total length (>= 5000 bp)",
                      "Total length (>= 10000 bp)","Total length (>= 25000 bp)","Total length (>= 50000 bp)",
                      "Number of contigs","Largest contig","Total length","GC (%)","N50","N90",
                      "auN","L50","L90","Number of N's per 100 kbp"]
    busco_metrics = ["Single","Duplicated","Fragmented","Missing"]

    # ------------------------------ iNaturalist (optional) ------------------------------
    if _HAS_INAT and pd.notna(inat_id):
        try:
            inat_html = f"https://www.inaturalist.org/observations/{int(inat_id)}"
        except Exception:
            inat_html = ""
        collector, observation_date, latitude, longitude, original_photo = get_inat_obs(inat_id)
        country, state, county = reverse_geocode(latitude, longitude)
        inat_df = pd.DataFrame({"Metric": inat_metrics,
                                "Value":  [collector, observation_date, latitude, longitude, country, state, county]})
        sample_stats_dict["INAT_DF"]    = inat_df
        sample_stats_dict["INAT_HTML"]  = inat_html
        sample_stats_dict["INAT_PHOTO"] = original_photo
    else:
        print(f"SKIP:\tNo iNaturalist ID provided or iNat disabled ({inat_id}).")

    # ------------------------------ ONT NanoPlot ------------------------------
    if os.path.isdir(ont_dir):
        for ndir in glob.glob(os.path.join(ont_dir, "*_ONT_nanoplot_analysis")):
            run_label = os.path.basename(ndir).split("_ONT")[0].upper()  # RAW / FILT / CORR
            html_matches  = glob.glob(os.path.join(ndir, "*NanoPlot-report.html"))
            stats_matches = glob.glob(os.path.join(ndir, "*NanoStats.txt"))
            if html_matches:
                sample_stats_dict[f"{run_label}_ONT_NANOPLOT_HTML"] = html_matches[0]
            if stats_matches:
                df_stats = parse_nanostats(stats_matches[0])
                sample_stats_dict[f"{run_label}_ONT_NANOSTATS_DF"] = df_stats
                print(f"\n=== {run_label} ONT NanoStats ===")
                try:
                    print(df_stats.to_markdown(index=False))
                except Exception:
                    print(df_stats.head())
    else:
        print(f"SKIP:\tNo ONT directory found at {ont_dir}.")

    # ------------------------------ Illumina FastQC ------------------------------
    if os.path.isdir(illumina_dir):
        for ndir in glob.glob(os.path.join(illumina_dir, "*fastqc_results")):
            base = os.path.basename(ndir)
            if base == "fastqc_results":
                run_label, f_pattern, r_pattern = "RAW", "*_1_fastqc.html", "*_2_fastqc.html"
            elif base == "dedup_fastqc_results":
                run_label, f_pattern, r_pattern = "DEDUP", "*_forward_dedup_fastqc.html", "*_reverse_dedup_fastqc.html"
            else:
                continue

            f_html_matches = glob.glob(os.path.join(ndir, f_pattern))
            r_html_matches = glob.glob(os.path.join(ndir, r_pattern))
            if f_html_matches:
                f_html = f_html_matches[0]
                sample_stats_dict[f"{run_label}_ILLU_F_FASTQC_HTML"] = f_html
                try:
                    sample_stats_dict[f"{run_label}_ILLU_F_BASIC_STATS_DF"] = parse_fastqc_basic_stats(f_html)
                except Exception as e:
                    print(f"WARN:\tFailed to parse FastQC (F) {f_html}: {e}")
            if r_html_matches:
                r_html = r_html_matches[0]
                sample_stats_dict[f"{run_label}_ILLU_R_FASTQC_HTML"] = r_html
                try:
                    sample_stats_dict[f"{run_label}_ILLU_R_BASIC_STATS_DF"] = parse_fastqc_basic_stats(r_html)
                except Exception as e:
                    print(f"WARN:\tFailed to parse FastQC (R) {r_html}: {e}")
    else:
        print(f"SKIP:\tNo Illumina directory found at {illumina_dir}.")

    # ------------------------------ PacBio NanoPlot ------------------------------
    if os.path.isdir(pacbio_dir):
        for ndir in glob.glob(os.path.join(pacbio_dir, "*_PacBio_nanoplot_analysis")):
            run_label = os.path.basename(ndir).split("_PacBio")[0].upper()  # RAW / FILT
            html_matches  = glob.glob(os.path.join(ndir, "*NanoPlot-report.html"))
            stats_matches = glob.glob(os.path.join(ndir, "*NanoStats.txt"))
            if html_matches:
                sample_stats_dict[f"{run_label}_PACBIO_NANOPLOT_HTML"] = html_matches[0]
            if stats_matches:
                df_stats = parse_nanostats(stats_matches[0])
                sample_stats_dict[f"{run_label}_PACBIO_NANOSTATS_DF"] = df_stats
                print(f"\n=== {run_label} PacBio NanoStats ===")
                try:
                    print(df_stats.to_markdown(index=False))
                except Exception:
                    print(df_stats.head())
    else:
        print(f"SKIP:\tNo PacBio directory found at {pacbio_dir}.")

    # ------------------------------ Per-assembler QUAST/BUSCO ------------------------------
    assembly_candidates = [masurca_dir, flye_dir, spades_dir, hifiasm_dir]
    for assembly_dir in assembly_candidates:
        if not os.path.isdir(assembly_dir):
            print(f"SKIP:\tNo Assembly directory found at {assembly_dir}.")
            continue

        asm_type = os.path.basename(assembly_dir).split("_")[0].upper()  # MASURCA/FLYE/SPADES/HIFIASM
        assembly_fasta = os.path.join(assembly_dir, f"{sample_id}_{asm_type.lower()}.fasta")

        # QUAST
        for quast_dir in glob.glob(os.path.join(assembly_dir, "*_quast")):
            try:
                quast_df, quast_html = process_quast_folder(quast_dir)
                sample_stats_dict[f"{asm_type}_QUAST_DF"]   = quast_df
                sample_stats_dict[f"{asm_type}_QUAST_HTML"] = quast_html
            except FileNotFoundError as e:
                print(f"QUAST processing failed for {quast_dir}: {e}")

        # BUSCO/Compleasm (skip if DB is NaN)
        for busco_index, busco_db in enumerate([first_busco_db, second_busco_db]):
            if pd.isna(busco_db):
                continue
            busco_type = "FIRST" if busco_index == 0 else "SECOND"
        
            # Accept both BUSCO and Compleasm output directories
            busco_dirs = glob.glob(os.path.join(assembly_dir, f"*_{busco_db}_busco")) + \
                         glob.glob(os.path.join(assembly_dir, f"*_{busco_db}_compleasm"))
            if not busco_dirs:
                print(f"SKIP:\tNo BUSCO/Compleasm outputs for {busco_db} in {assembly_dir}")
                continue
        
            for busco_dir in busco_dirs:
                # Prefer BUSCO short_summary, else Compleasm summary.txt
                metrics = {}
                summary_files = glob.glob(os.path.join(busco_dir, "short_summary*.txt"))
                if summary_files:
                    with open(summary_files[0], "r") as f:
                        for line in f:
                            if "%" in line and "C:" in line:
                                metrics = _parse_busco_summary_line(line)
                                break
                else:
                    metrics = _parse_compleasm_summary(os.path.join(busco_dir, "summary.txt"))
        
                if not metrics:
                    print(f"SKIP:\tNo parsable BUSCO/Compleasm summary in {busco_dir}")
                    continue
        
                metrics_df = pd.DataFrame(list(metrics.items()), columns=["Metric", "Value"])
        
                # Genes table (optional; supports BUSCO or Compleasm layouts)
                genes_file = _find_busco_genes_table(busco_dir, busco_db)
                busco_genes_df = pd.DataFrame()
                if genes_file:
                    try:
                        busco_genes_df = pd.read_csv(genes_file, sep="\t", comment="#")
                        for col in ["Busco id","Status","Sequence","Gene Start","Gene End","Strand","Score","Length"]:
                            if col not in busco_genes_df.columns:
                                busco_genes_df[col] = "-"
                    except Exception as e:
                        print(f"WARN:\tFailed to parse BUSCO/Compleasm genes table {genes_file}: {e}")
        
                outpath = os.path.join(busco_dir, f"{sample_id}_{busco_db}_EGAP_busco.html")
                sample_stats_dict[f"{asm_type}_{busco_type}_BUSCO_DF"]        = metrics_df
                sample_stats_dict[f"{asm_type}_{busco_type}_BUSCO_HTML"]      = outpath
                sample_stats_dict[f"{asm_type}_{busco_type}_BUSCO_GENES_DF"]  = busco_genes_df

    # ------------------------------ Final assembly QUAST/BUSCO (robust discovery) ------------------------------
    def _first_existing(paths):
        for p in paths:
            if p and os.path.exists(p):
                return p
        return None
    
    # Candidate basenames (without extension) we can accept
    candidate_bases = []
    if pd.notna(ref_seq) and pd.isna(est_size):
        candidate_bases.append(os.path.splitext(ref_seq)[0])
    candidate_bases += [
        os.path.join(sample_dir, f"{sample_id}_final_EGAP_assembly"),
        os.path.join(sample_dir, f"{sample_id}_final_curated"),
        os.path.join(sample_dir, f"{sample_id}_final_polish_assembly"),
        os.path.join(sample_dir, f"{sample_id}_EGAP_assembly"),
    ]
    
    # Pick the first base that has either a FASTA or a QUAST folder
    assembly_base = None
    for base in candidate_bases:
        if os.path.exists(base + ".fasta") or os.path.isdir(base + "_quast"):
            assembly_base = base
            break
    if assembly_base is None and candidate_bases:
        assembly_base = candidate_bases[0]
    
    assembly_fasta = (assembly_base + ".fasta") if assembly_base else ""
    print(f"DEBUG - final assembly base - {assembly_base}")
    
    # -------- QUAST: choose the first available report.tsv --------
    quast_dir_candidates = []
    if assembly_base:
        quast_dir_candidates = [assembly_base + "_quast"] + [b + "_quast" for b in candidate_bases if b != assembly_base]
    final_quast_dir = None
    for qd in quast_dir_candidates:
        if os.path.isdir(qd) and os.path.exists(os.path.join(qd, "report.tsv")):
            final_quast_dir = qd
            break
    
    if final_quast_dir:
        print(f"DEBUG - final_quast_dir - {final_quast_dir}")
        try:
            final_df, final_html = process_quast_folder(final_quast_dir)
            sample_stats_dict["FINAL_QUAST_DF"]   = final_df
            sample_stats_dict["FINAL_QUAST_HTML"] = final_html
        except FileNotFoundError as e:
            print(f"FINAL QUAST processing failed for {final_quast_dir}: {e}")
    else:
        if quast_dir_candidates:
            print(f"WARN:\tNo QUAST report found for final assembly; tried: {quast_dir_candidates}")
    
    # -------- BUSCO/Compleasm: accept either style and any of the basenames --------
    for busco_index, busco_db in enumerate([first_busco_db, second_busco_db]):
        if pd.isna(busco_db):
            continue
        busco_type = "FIRST" if busco_index == 0 else "SECOND"
    
        busco_dir_candidates = []
        for base in ([assembly_base] if assembly_base else []) + [b for b in candidate_bases if b and b != assembly_base]:
            busco_dir_candidates.append(f"{base}_{busco_db}_busco")
            busco_dir_candidates.append(f"{base}_{busco_db}_compleasm")
    
        chosen_busco_dir = _first_existing(busco_dir_candidates)
        if not chosen_busco_dir:
            print(f"WARN:\tNo BUSCO/Compleasm dir found for {busco_db}. Tried: {busco_dir_candidates}")
            continue
    
        metrics = {}
        short_sum = _first_existing(glob.glob(os.path.join(chosen_busco_dir, "short_summary*.txt")))
        if short_sum:
            with open(short_sum, "r") as f:
                for line in f:
                    if "%" in line and "C:" in line:
                        metrics = _parse_busco_summary_line(line)
                        break
        else:
            metrics = _parse_compleasm_summary(os.path.join(chosen_busco_dir, "summary.txt"))
    
        if not metrics:
            print(f"WARN:\tNo parsable BUSCO/Compleasm summary in {chosen_busco_dir}")
            continue
    
        metrics_df = pd.DataFrame(list(metrics.items()), columns=["Metric", "Value"])
    
        genes_file = _find_busco_genes_table(chosen_busco_dir, busco_db)
        busco_genes_df = pd.DataFrame()
        if genes_file:
            try:
                busco_genes_df = pd.read_csv(genes_file, sep="\t", comment="#")
                for col in ["Busco id","Status","Sequence","Gene Start","Gene End","Strand","Score","Length"]:
                    if col not in busco_genes_df.columns:
                        busco_genes_df[col] = "-"
            except Exception as e:
                print(f"WARN:\tFailed to parse BUSCO/Compleasm genes table {genes_file}: {e}")
    
        sample_stats_dict[f"FINAL_{busco_type}_BUSCO_DF"]        = metrics_df
        sample_stats_dict[f"FINAL_{busco_type}_BUSCO_HTML"]      = os.path.join(chosen_busco_dir, f"{sample_id}_{busco_db}_EGAP_busco.html")
        sample_stats_dict[f"FINAL_{busco_type}_BUSCO_GENES_DF"]  = busco_genes_df

    # ------------------------------ Render or stub ------------------------------
    os.makedirs(sample_dir, exist_ok=True)
    html_report_outpath = os.path.join(sample_dir, f"{sample_id}_EGAP_summary.html")

    if templates_dir is None:
        # No templates available; write a stub and succeed
        stub = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>{sample_id} — EGAP summary (templates missing)</title></head>
<body>
<h1>EGAP Summary (Templates Missing)</h1>
<p>A full report could not be generated because both template files are required:</p>
<ul>
  <li><code>EGAP_summary.html</code></li>
  <li><code>EGAP_busco.html</code></li>
</ul>
<p>Place these files somewhere discoverable (see logs), or set <code>EGAP_TEMPLATES_DIR</code>, then re-run.</p>
<p>The pipeline continues without failing.</p>
</body></html>"""
        with open(html_report_outpath, "w") as fh:
            fh.write(stub)
        print(f"WARN:\tTemplates unavailable; wrote stub report: {html_report_outpath}")
        return html_report_outpath

    # Templates are available: render real report
    print(f"DEBUG - templates_dir - {templates_dir}")
    env = Environment(
        loader=FileSystemLoader(str(templates_dir)),
        autoescape=select_autoescape(["html", "xml"])
    )
    main_tmpl  = env.get_template("EGAP_summary.html")
    busco_tmpl = env.get_template("EGAP_busco.html")

    # Helper to coerce DataFrame -> {metric: value}
    def _df_to_map(df):
        if df is None or (isinstance(df, pd.DataFrame) and df.empty):
            return {}
        if "Metric" in df.columns and "Value" in df.columns:
            return dict(zip(df["Metric"], df["Value"]))
        try:
            return dict(zip(df.index, df["Value"]))
        except Exception:
            return {}

    # iNat
    inat_df        = sample_stats_dict.get("INAT_DF")
    inat_link      = sample_stats_dict.get("INAT_HTML")
    inat_photo     = sample_stats_dict.get("INAT_PHOTO")
    inat_thumbnail = sample_stats_dict.get("INAT_THUMBNAIL") or inat_photo
    inat_map       = _df_to_map(inat_df)

    # ONT
    ont_raw_df  = sample_stats_dict.get("RAW_ONT_NANOSTATS_DF")
    ont_filt_df = sample_stats_dict.get("FILT_ONT_NANOSTATS_DF")
    ont_corr_df = sample_stats_dict.get("CORR_ONT_NANOSTATS_DF")
    ont_raw_map  = _df_to_map(ont_raw_df)
    ont_filt_map = _df_to_map(ont_filt_df)
    ont_corr_map = _df_to_map(ont_corr_df)
    ont_raw_link = sample_stats_dict.get("RAW_ONT_NANOPLOT_HTML", "")
    ont_filt_link= sample_stats_dict.get("FILT_ONT_NANOPLOT_HTML", "")
    ont_corr_link= sample_stats_dict.get("CORR_ONT_NANOPLOT_HTML", "")

    # Illumina
    illumina_raw_f_df   = sample_stats_dict.get("RAW_ILLU_F_BASIC_STATS_DF")
    illumina_raw_r_df   = sample_stats_dict.get("RAW_ILLU_R_BASIC_STATS_DF")
    illumina_dedup_f_df = sample_stats_dict.get("DEDUP_ILLU_F_BASIC_STATS_DF")
    illumina_dedup_r_df = sample_stats_dict.get("DEDUP_ILLU_R_BASIC_STATS_DF")

    illumina_raw_f_map   = _df_to_map(illumina_raw_f_df)
    illumina_raw_r_map   = _df_to_map(illumina_raw_r_df)
    illumina_dedup_f_map = _df_to_map(illumina_dedup_f_df)
    illumina_dedup_r_map = _df_to_map(illumina_dedup_r_df)
    illumina_raw_f_link  = sample_stats_dict.get("RAW_ILLU_F_FASTQC_HTML", "")
    illumina_raw_r_link  = sample_stats_dict.get("RAW_ILLU_R_FASTQC_HTML", "")
    illumina_dedup_f_link= sample_stats_dict.get("DEDUP_ILLU_F_FASTQC_HTML", "")
    illumina_dedup_r_link= sample_stats_dict.get("DEDUP_ILLU_R_FASTQC_HTML", "")

    # PacBio
    pacbio_raw_df  = sample_stats_dict.get("RAW_PACBIO_NANOSTATS_DF")
    pacbio_filt_df = sample_stats_dict.get("FILT_PACBIO_NANOSTATS_DF")
    pacbio_raw_map  = _df_to_map(pacbio_raw_df)
    pacbio_filt_map = _df_to_map(pacbio_filt_df)
    pacbio_raw_link = sample_stats_dict.get("RAW_PACBIO_NANOPLOT_HTML", "")
    pacbio_filt_link= sample_stats_dict.get("FILT_PACBIO_NANOPLOT_HTML", "")

    # Per-assembler QUAST/BUSCO maps + links
    def _maps(prefix):
        return (_df_to_map(sample_stats_dict.get(f"{prefix}_QUAST_DF")),
                _df_to_map(sample_stats_dict.get(f"{prefix}_FIRST_BUSCO_DF")),
                _df_to_map(sample_stats_dict.get(f"{prefix}_SECOND_BUSCO_DF")),
                sample_stats_dict.get(f"{prefix}_QUAST_HTML",""),
                sample_stats_dict.get(f"{prefix}_FIRST_BUSCO_HTML",""),
                sample_stats_dict.get(f"{prefix}_SECOND_BUSCO_HTML",""))

    masurca_quast_map, masurca_first_busco_map, masurca_second_busco_map, \
        masurca_quast_link, masurca_first_busco_link, masurca_second_busco_link = _maps("MASURCA")
    flye_quast_map, flye_first_busco_map, flye_second_busco_map, \
        flye_quast_link, flye_first_busco_link, flye_second_busco_link = _maps("FLYE")
    spades_quast_map, spades_first_busco_map, spades_second_busco_map, \
        spades_quast_link, spades_first_busco_link, spades_second_busco_link = _maps("SPADES")
    hifiasm_quast_map, hifiasm_first_busco_map, hifiasm_second_busco_map, \
        hifiasm_quast_link, hifiasm_first_busco_link, hifiasm_second_busco_link = _maps("HIFIASM")

    final_quast_map        = _df_to_map(sample_stats_dict.get("FINAL_QUAST_DF"))
    final_first_busco_map  = _df_to_map(sample_stats_dict.get("FINAL_FIRST_BUSCO_DF"))
    final_second_busco_map = _df_to_map(sample_stats_dict.get("FINAL_SECOND_BUSCO_DF"))
    final_quast_link       = sample_stats_dict.get("FINAL_QUAST_HTML","")
    final_first_busco_link = sample_stats_dict.get("FINAL_FIRST_BUSCO_HTML","")
    final_second_busco_link= sample_stats_dict.get("FINAL_SECOND_BUSCO_HTML","")

    # -------- Render BUSCO sections (now that we have busco_tmpl) --------
    # For per-assembler BUSCO
    for key in list(sample_stats_dict.keys()):
        if key.endswith("_BUSCO_HTML") and "FINAL" not in key and os.path.dirname(sample_stats_dict[key]):
            outpath = sample_stats_dict[key]
            parts = key.split("_")
            asm_type, busco_type = parts[0], parts[1]  # e.g., MASURCA_FIRST
            df = sample_stats_dict.get(f"{asm_type}_{busco_type}_BUSCO_DF", pd.DataFrame())
            html = busco_tmpl.render(
                assembly_name=os.path.basename(outpath).replace("_EGAP_busco.html", ""),
                generated_time=now,
                busco_data=[{
                    "assembly_fasta": "",
                    "busco_db": busco_type.capitalize(),
                    "svg_path": "",
                    "metrics": dict(zip(df["Metric"], df["Value"])) if not df.empty else {},
                    "busco_genes_df": pd.DataFrame()
                }]
            )
            os.makedirs(os.path.dirname(outpath), exist_ok=True)
            with open(outpath, "w") as fh:
                fh.write(html)

    # For FINAL BUSCO
    for busco_type in ("FIRST", "SECOND"):
        df = sample_stats_dict.get(f"FINAL_{busco_type}_BUSCO_DF", pd.DataFrame())
        outpath = sample_stats_dict.get(f"FINAL_{busco_type}_BUSCO_HTML", "")
        if not outpath:
            continue
        html = busco_tmpl.render(
            assembly_name=os.path.basename(outpath).replace("_EGAP_busco.html", ""),
            generated_time=now,
            busco_data=[{
                "assembly_fasta": "",
                "busco_db": busco_type.capitalize(),
                "svg_path": "",
                "metrics": dict(zip(df["Metric"], df["Value"])) if not df.empty else {},
                "busco_genes_df": pd.DataFrame()
            }]
        )
        os.makedirs(os.path.dirname(outpath), exist_ok=True)
        with open(outpath, "w") as fh:
            fh.write(html)

    # -------- Render MAIN summary --------
    os.makedirs(sample_dir, exist_ok=True)
    html_report_outpath = os.path.join(sample_dir, f"{sample_id}_EGAP_summary.html")

    render_ctx = {
        "sample_id": sample_id,
        "generated_time": now,

        "inat_metrics": inat_metrics, "inat_map": inat_map,
        "inat_link": sample_stats_dict.get("INAT_HTML"), "inat_thumbnail": inat_thumbnail, "inat_photo": inat_photo,

        "ont_metrics": ont_metrics,
        "ont_raw_map": ont_raw_map, "ont_filt_map": ont_filt_map, "ont_corr_map": ont_corr_map,
        "ont_raw_link": ont_raw_link, "ont_filtered_link": ont_filt_link, "ont_corrected_link": ont_corr_link,

        "illumina_metrics": illumina_metrics,
        "illumina_raw_f_map": illumina_raw_f_map, "illumina_raw_r_map": illumina_raw_r_map,
        "illumina_dedup_f_map": illumina_dedup_f_map, "illumina_dedup_r_map": illumina_dedup_r_map,
        "illumina_raw_f_link": illumina_raw_f_link, "illumina_raw_r_link": illumina_raw_r_link,
        "illumina_dedup_f_link": illumina_dedup_f_link, "illumina_dedup_r_link": illumina_dedup_r_link,

        "pacbio_metrics": pacbio_metrics,
        "pacbio_raw_map": pacbio_raw_map, "pacbio_filt_map": pacbio_filt_map,
        "pacbio_raw_link": pacbio_raw_link, "pacbio_filtered_link": pacbio_filt_link,

        "quast_metrics": quast_metrics, "busco_metrics": busco_metrics,

        "masurca_quast_map": masurca_quast_map,
        "masurca_first_busco_map": masurca_first_busco_map, "masurca_second_busco_map": masurca_second_busco_map,
        "masurca_quast_link": masurca_quast_link,
        "masurca_first_busco_link": masurca_first_busco_link, "masurca_second_busco_link": masurca_second_busco_link,

        "flye_quast_map": flye_quast_map,
        "flye_first_busco_map": flye_first_busco_map, "flye_second_busco_map": flye_second_busco_map,
        "flye_quast_link": flye_quast_link,
        "flye_first_busco_link": flye_first_busco_link, "flye_second_busco_link": flye_second_busco_link,

        "spades_quast_map": spades_quast_map,
        "spades_first_busco_map": spades_first_busco_map, "spades_second_busco_map": spades_second_busco_map,
        "spades_quast_link": spades_quast_link,
        "spades_first_busco_link": spades_first_busco_link, "spades_second_busco_link": spades_second_busco_link,

        "hifiasm_quast_map": hifiasm_quast_map,
        "hifiasm_first_busco_map": hifiasm_first_busco_map, "hifiasm_second_busco_map": hifiasm_second_busco_map,
        "hifiasm_quast_link": hifiasm_quast_link,
        "hifiasm_first_busco_link": hifiasm_first_busco_link, "hifiasm_second_busco_link": hifiasm_second_busco_link,

        "final_quast_map": final_quast_map,
        "final_first_busco_map": final_first_busco_map, "final_second_busco_map": final_second_busco_map,
        "final_quast_link": final_quast_link,
        "final_first_busco_link": final_first_busco_link,
        "final_second_busco_link": final_second_busco_link,
    }

    html = main_tmpl.render(**render_ctx)
    with open(html_report_outpath, "w") as fh:
        fh.write(html)

    print(f"PASS:\tWrote report: {html_report_outpath}")
    return html_report_outpath


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python3 html_reporter.py <sample_id> <input_csv> <output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)

    _ = html_reporter(
        sys.argv[1],              # sample_id
        sys.argv[2],              # input_csv
        sys.argv[3],              # output_dir
        str(sys.argv[4]),         # cpu_threads
        str(sys.argv[5])          # ram_gb
    )
