#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
html_reporter.py

Generates an EGAP HTML summary using Jinja templates. Robust to missing
artifacts (FastQC / NanoPlot / QUAST / BUSCO), CWD-safe, and tolerant of
metadata quirks.

Stage:
    Final HTML Report

Created on Wed Aug 16 2023

Updated on 2026-04-16

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os
import glob
import re
from io import StringIO
from datetime import datetime

import pandas as pd
from bs4 import BeautifulSoup
from jinja2 import Environment, FileSystemLoader, select_autoescape

from process_metadata import get_inat_obs, reverse_geocode
from record_provenance import (load_provenance, pipeline_steps, stage_timings,
                               meaningful_tools, resolve_versions,
                               versions_from_entries)
from file_operations import resolve_templates_dir
from qc_assessment import busco_output_paths, completeness_tool_order
from utilities import load_sample_context, run_stage, to_abs


# Human-readable labels for the pipeline stages recorded in provenance.
STAGE_LABELS = {
    "preprocess_refseq": "Reference preparation",
    "preprocess_illumina": "Illumina read preprocessing",
    "preprocess_ont": "ONT read preprocessing",
    "preprocess_pacbio": "PacBio read preprocessing",
    "decontaminate_reads": "Read decontamination (Kraken2)",
    "assemble_masurca": "MaSuRCA assembly",
    "assemble_flye": "Flye assembly",
    "assemble_spades": "SPAdes assembly",
    "assemble_hifiasm": "hifiasm assembly",
    "compare_assemblies": "Assembly comparison / best selection",
    "polish_assembly": "Polishing",
    "curate_assembly": "Curation / scaffolding",
    "decontaminate_assembly": "Assembly decontamination (Tiara)",
    "qc_assessment": "Quality assessment (QUAST / BUSCO)",
    "html_reporter": "Report generation",
    "process_metadata": "NCBI submission metadata",
}


def _file_link(path):
    """Return ``{'name': basename, 'href': abspath}`` if *path* exists, else ``None``.

    Used to surface clickable links to the actual reads/assembly files in the
    report so a reader can trace each artefact back to disk.
    """
    if path and isinstance(path, str) and os.path.exists(path):
        return {"name": os.path.basename(path), "href": os.path.abspath(path)}
    return None


# ------------------------------ Helpers ------------------------------
def lineage_dirname(busco_odb: str, db_version: str = "odb12") -> str:
    """
    Return a lineage directory name guaranteed to end with _odbNN,
    without doubling the suffix if it's already present.
    """
    s = str(busco_odb).strip()
    return s if re.search(r"_odb\d+$", s) else f"{s}_{db_version}"


def parse_compleasm_full_table(tsv_path: str) -> pd.DataFrame:
    """
    Parse a Compleasm/BUSCO 'full_table*.tsv' and normalize columns to a common schema:
      Required:  Busco id | Status | Sequence | Gene Start | Gene End | Strand | Score | Length
      Optional:  Identity | Fraction | Frameshift events | Best gene | Codons
    Returns an ordered DataFrame with required columns first, then any of the optional ones if present.
    """
    if not os.path.exists(tsv_path):
        return pd.DataFrame()

    # Robustly handle commented header lines (e.g., "# Busco id\tStatus\t...")
    header_fields = None
    data_lines = []
    with open(tsv_path, "r", encoding="utf-8") as fh:
        for line in fh:
            s = line.rstrip("\n")
            if s.startswith("#"):
                # pick the first commented header that contains 'busco id'
                if header_fields is None and "busco id" in s.lower():
                    header_fields = [c.strip() for c in s.lstrip("#").strip().split("\t")]
                # skip all commented lines (including non-header comments)
                continue
            if s.strip():
                data_lines.append(s)

    if header_fields is not None:
        buf = StringIO("\n".join(["\t".join(header_fields)] + data_lines))
        df = pd.read_csv(buf, sep="\t", dtype=str, keep_default_na=False)
    else:
        # Fallback: no commented header found; read normally and strip any leading '#'
        df = pd.read_csv(tsv_path, sep="\t", dtype=str, keep_default_na=False)
        df.columns = [str(c).lstrip("#").strip() for c in df.columns]

    if df is None or df.empty:
        return pd.DataFrame()

    # Build case/space-insensitive map: lower -> original (strip leading '#' just in case)
    lower_cols = {c.lower().strip().lstrip("#"): c for c in df.columns}

    # Canonical names and their common aliases (lowercased)
    alias_map = {
        "busco id":          ["busco id", "busco_id", "busco"],
        "status":            ["status"],
        "sequence":          ["sequence", "seq", "contig", "scaffold"],
        "gene start":        ["gene start", "start", "gene_start"],
        "gene end":          ["gene end", "end", "gene_end"],
        "strand":            ["strand"],
        "score":             ["score", "bitscore", "bit score"],
        "length":            ["length", "len"],
        # Compleasm extras
        "identity":          ["identity", "ident"],
        "fraction":          ["fraction", "frac", "coverage", "cov"],
        "frameshift events": ["frameshift events", "frameshifts", "frameshift_events"],
        "best gene":         ["best gene", "best_gene", "best hit", "best_hit"],
        "codons":            ["codons"]
    }

    # Build rename map from whatever we find -> canonical (fixed-case)
    proper_case = {
        "busco id": "Busco id",
        "status": "Status",
        "sequence": "Sequence",
        "gene start": "Gene Start",
        "gene end": "Gene End",
        "strand": "Strand",
        "score": "Score",
        "length": "Length",
        "identity": "Identity",
        "fraction": "Fraction",
        "frameshift events": "Frameshift Events",
        "best gene": "Best Gene",
        "codons": "Codons",
    }
    rename = {}
    for canonical, aliases in alias_map.items():
        for a in aliases:
            if a in lower_cols:
                rename[lower_cols[a]] = proper_case.get(canonical, canonical.title())
                break

    # Apply rename
    df = df.rename(columns=rename).copy()

    # Canonical required columns (exact wanted casing)
    required = ["Busco id","Status","Sequence","Gene Start","Gene End","Strand","Score","Length"]
    for col in required:
        if col not in df.columns:
            df[col] = "-"

    # Optional extras (only include if present after rename)
    optional = ["Identity","Fraction","Frameshift Events","Best Gene","Codons"]
    ordered_cols = required + [c for c in optional if c in df.columns]

    # Coerce numeric-ish fields that should remain strings but trimmed. Missing
    # cells become "-" like missing columns (pandas 3's astype(str) keeps them
    # missing; pandas 2 wrote "nan").
    for c in ordered_cols:
        if c in df.columns:
            df[c] = df[c].fillna("-").astype(str).str.strip()

    return df[ordered_cols]

def load_busco_or_compleasm_genes_table(busco_dir: str, busco_db: str) -> pd.DataFrame:
    """Parse the gene table of a BUSCO or Compleasm folder; empty when absent or unreadable."""
    path = find_busco_genes_table(busco_dir, busco_db) if busco_dir and busco_db else ""
    if not path:
        return pd.DataFrame()
    try:
        return parse_compleasm_full_table(path)
    except Exception as e:
        print(f"WARN:\tFailed to parse BUSCO/Compleasm genes table {path}: {e}")
        return pd.DataFrame()


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

def parse_busco_summary_line(line: str) -> dict:
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

def parse_compleasm_summary(path: str) -> dict:
    """
    Parse Compleasm summary.txt where S:/D:/F:/M: appear on separate lines.
    - Reads percentages for S/D/F/M (keeps your original format: no '%' in strings)
    - Tries to read counts if present (e.g., 'S:99.7%, 346' or 'S: 346 (99.7%)')
    - Computes Complete = S + D (rounded 2dp)
    - Computes Total = S_cnt + D_cnt + F_cnt + M_cnt, or uses 'n:' if present
    """
    if not os.path.exists(path):
        return {}

    import re

    perc = {"S": None, "D": None, "F": None, "M": None}
    cnt  = {"S": None, "D": None, "F": None, "M": None}
    total_n = None

    with open(path, "r") as fh:
        for line in fh:
            # capture 'n:' if present (some tools include the lineage size)
            if total_n is None:
                m_n = re.search(r"\bn\s*[:=]\s*(\d+)\b", line)
                if m_n:
                    try:
                        total_n = int(m_n.group(1))
                    except Exception:
                        total_n = None

            for key in ("S", "D", "F", "M"):
                if f"{key}:" not in line:
                    continue

                # Pattern A: "S:99.71%, 346"
                m = re.search(rf"{key}\s*:\s*([\d\.]+)\s*%\s*,\s*(\d+)", line)
                if m:
                    perc[key] = m.group(1).strip()
                    try:
                        cnt[key] = int(m.group(2))
                    except Exception:
                        cnt[key] = None
                    continue

                # Pattern B: "S: 346 (99.71%)"
                m = re.search(rf"{key}\s*:\s*(\d+)\s*\(([\d\.]+)\s*%\)", line)
                if m:
                    try:
                        cnt[key] = int(m.group(1))
                    except Exception:
                        cnt[key] = None
                    perc[key] = m.group(2).strip()
                    continue

                # Pattern C: "S: 99.71%"
                m = re.search(rf"{key}\s*:\s*([\d\.]+)\s*%", line)
                if m:
                    perc[key] = m.group(1).strip()
                    continue

                # Pattern D: "S: 346"
                m = re.search(rf"{key}\s*:\s*(\d+)\b", line)
                if m:
                    try:
                        cnt[key] = int(m.group(1))
                    except Exception:
                        cnt[key] = None
                    continue

    # Require at least S and D percentages to compute "Complete"
    if perc["S"] is None or perc["D"] is None:
        return {}

    try:
        complete = str(round(float(perc["S"]) + float(perc["D"]), 2))
    except Exception:
        complete = ""

    # Prefer explicit 'n:' if found; else sum available counts
    if total_n is None:
        counts = [c for c in (cnt["S"], cnt["D"], cnt["F"], cnt["M"]) if isinstance(c, int)]
        total_n = sum(counts) if counts else None

    return {
        "Complete":   complete,
        "Single":     perc["S"] or "",
        "Duplicated": perc["D"] or "",
        "Fragmented": perc["F"] or "",
        "Missing":    perc["M"] or "",
        "Total":      (str(total_n) if isinstance(total_n, int) else "")
    }


def extract_busco_db_info(busco_dir: str, busco_db: str) -> dict:
    """Extract the lineage database version actually used by BUSCO/Compleasm.

    The lineage name alone (e.g. ``basidiomycota``) omits the OrthoDB version
    (``_odb10`` / ``_odb12``), which *is* the database version. This reads the
    authoritative value from the tool's own output: the BUSCO
    ``short_summary*.txt`` ``"The lineage dataset is:"`` line (with its creation
    date and BUSCO version), or the Compleasm ``summary.txt`` lineage header.
    Falls back to the normalized lineage directory name.

    Parameters
    ----------
    busco_dir : str
        BUSCO/Compleasm output directory.
    busco_db : str
        Lineage as given in the TSV (may lack the ``_odbNN`` suffix).

    Returns
    -------
    dict
        ``{lineage, creation_date, tool, tool_version}`` (strings; empty when
        a field cannot be determined).
    """
    info = {"lineage": lineage_dirname(busco_db), "creation_date": "",
            "tool": "", "tool_version": ""}
    if not busco_dir or not os.path.isdir(busco_dir):
        return info

    short_summaries = glob.glob(os.path.join(busco_dir, "short_summary*.txt"))
    if short_summaries:
        info["tool"] = "BUSCO"
        try:
            with open(short_summaries[0], "r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if "lineage dataset is:" in line:
                        after = line.split("lineage dataset is:", 1)[1].strip()
                        lineage = after.split("(")[0].strip()
                        if lineage:
                            info["lineage"] = lineage
                        mdate = re.search(r"Creation date:\s*([0-9][0-9\-/]+)", line)
                        if mdate:
                            info["creation_date"] = mdate.group(1)
                    elif "BUSCO version is:" in line:
                        info["tool_version"] = line.split("BUSCO version is:", 1)[1].strip()
        except Exception:
            pass
        return info

    # Compleasm: no short_summary; lineage may appear in summary.txt or the dir name.
    info["tool"] = "Compleasm" if busco_dir.endswith("_compleasm") else "BUSCO/Compleasm"
    summary_txt = os.path.join(busco_dir, "summary.txt")
    if os.path.exists(summary_txt):
        try:
            with open(summary_txt, "r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    m = re.search(r"lineage\s*[:=]\s*(\S+)", line, re.IGNORECASE)
                    if m:
                        info["lineage"] = m.group(1).strip()
                        break
        except Exception:
            pass
    return info



def find_busco_genes_table(busco_dir: str, busco_db: str) -> str:
    """
    Return the first existing BUSCO/Compleasm 'full_table' path.
    Supports both BUSCO and Compleasm layouts.
    """
    lineage = lineage_dirname(busco_db)
    candidates = [
        # BUSCO; busco_output_paths also finds 6.1's run_<lineage>.2 folder
        busco_output_paths(busco_dir, lineage, os.path.basename(busco_dir))[1],
        os.path.join(busco_dir, lineage, "full_table_busco_format.tsv"),
        os.path.join(busco_dir, lineage, "full_table.tsv"),
    ]
    for p in candidates:
        if os.path.exists(p):
            return p
    return ""


# Read sets in the order the coverage table lists them: stats prefix -> label.
COVERAGE_LABELS = {
    "RAW_ILLU": "Illumina, raw",
    "TRIMMED_ILLU": "Illumina, trimmed",
    "DEDUPED_ILLU": "Illumina, deduplicated",
    "RAW_ONT": "ONT, raw",
    "FILT_ONT": "ONT, filtered",
    "CORRECT_ONT": "ONT, corrected",
    "RAW_PACBIO": "PacBio, raw",
    "FILT_PACBIO": "PacBio, filtered",
}


def coverage_rows(stats_csv: str):
    """Return ``(rows, basis)`` for the report's coverage table from the final stats CSV.

    *rows* holds one ``{read_set, total_bases, coverage}`` dict per read set
    that has a coverage value; *basis* says what the totals were divided by.
    Both are empty when the CSV is missing (QC failed) or has no coverage.
    """
    if not os.path.exists(stats_csv):
        return [], ""
    stats = pd.read_csv(stats_csv).iloc[0]
    rows = [{"read_set": label,
             "total_bases": f"{int(stats[f'{prefix}_TOTAL_BASES']):,}",
             "coverage": f"{stats[f'{prefix}_COVERAGE']:g}x"}
            for prefix, label in COVERAGE_LABELS.items()
            if pd.notna(stats.get(f"{prefix}_COVERAGE")) and pd.notna(stats.get(f"{prefix}_TOTAL_BASES"))]
    basis = stats.get("COVERAGE_BASIS")
    return rows, basis if isinstance(basis, str) else ""


def completeness_metrics(busco_dir: str) -> dict:
    """Return the summary metrics of a BUSCO or Compleasm folder, ``{}`` if none parse."""
    short_summaries = glob.glob(os.path.join(busco_dir, "short_summary*.txt"))
    if not short_summaries:
        return parse_compleasm_summary(os.path.join(busco_dir, "summary.txt"))
    with open(short_summaries[0], "r") as f:
        for line in f:
            if "%" in line and "C:" in line:
                return parse_busco_summary_line(line)
    return {}


def find_busco_svg(busco_dir: str, busco_db: str) -> str:
    """
    Try common SVG names emitted by your pipeline:
      - <base>_{lineage}_busco.svg
      - <base>_{lineage}_compleasm.svg
      - <busco_dir_basename>.svg
      - anything *{lineage}*busco.svg or *{lineage}*compleasm.svg in the parent dir
    """
    try:
        if not busco_dir or not os.path.isdir(busco_dir):
            return ""
        lineage = lineage_dirname(busco_db)
        parent = os.path.dirname(busco_dir)
        tag = os.path.basename(busco_dir)          # e.g., sample_final_{lineage}_busco

        candidates = [
            os.path.join(parent, f"{tag}.svg"),
        ]
        # add globs
        candidates += glob.glob(os.path.join(parent, f"*{lineage}*busco.svg"))
        candidates += glob.glob(os.path.join(parent, f"*{lineage}*compleasm.svg"))

        for c in candidates:
            if c and os.path.exists(c):
                return c
    except Exception:
        pass
    return ""

# ------------------------------ Main ------------------------------
def html_reporter(sample_id, input_tsv, output_dir, cpu_threads, ram_gb):
    """
    Render the EGAP HTML report for a sample.
    Returns: path to <sample_dir>/<sample_id>_EGAP_summary.html
    """
    # Resolve critical paths first (CWD-safe)
    input_tsvto_abs = to_abs(input_tsv)
    output_dirto_abs = to_abs(output_dir)

    # Fail now, before any work, if the install has no templates
    templates_dir = resolve_templates_dir()

    now = datetime.now().strftime("%Y%m%d-%H:%M:%S")

    # Metadata row
    ctx = load_sample_context(sample_id, input_tsvto_abs, output_dirto_abs, cpu_threads, ram_gb)
    sample_stats_dict = ctx.sample_stats_dict
    current_series = ctx.current_series

    # TSV columns (tolerate missing)
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
    species_dir = os.path.join(output_dirto_abs, str(species_id))
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
    busco_metrics = ["Complete", "Single", "Duplicated", "Fragmented", "Missing", "Total"]

    # ------------------------------ iNaturalist (optional) ------------------------------
    if pd.notna(inat_id):
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
            elif base == "trimmed_fastqc_results":
                run_label, f_pattern, r_pattern = "TRIMMED", "*_forward_mapped_fastqc.html", "*_reverse_mapped_fastqc.html"
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
        
            # Both tools' folders may exist. The last one with a usable summary
            # wins below, so list this run's tool last.
            this_tool, other_tool = completeness_tool_order()
            busco_dirs = glob.glob(os.path.join(assembly_dir, f"*_{busco_db}_{other_tool}")) + \
                         glob.glob(os.path.join(assembly_dir, f"*_{busco_db}_{this_tool}"))
            if not busco_dirs:
                print(f"SKIP:\tNo BUSCO/Compleasm outputs for {busco_db} in {assembly_dir}")
                continue
        
            for busco_dir in busco_dirs:
                svg_path = find_busco_svg(busco_dir, busco_db)
                
                metrics = completeness_metrics(busco_dir)
        
                if not metrics:
                    print(f"SKIP:\tNo parsable BUSCO/Compleasm summary in {busco_dir}")
                    continue
        
                metrics_df = pd.DataFrame(list(metrics.items()), columns=["Metric", "Value"])
        
                # Genes table (optional; supports BUSCO or Compleasm layouts)
                busco_genes_df = load_busco_or_compleasm_genes_table(busco_dir, busco_db)
        
                outpath = os.path.join(busco_dir, f"{sample_id}_{busco_db}_EGAP_busco.html")
                sample_stats_dict[f"{asm_type}_{busco_type}_BUSCO_DF"]        = metrics_df
                sample_stats_dict[f"{asm_type}_{busco_type}_BUSCO_HTML"]      = outpath
                sample_stats_dict[f"{asm_type}_{busco_type}_BUSCO_GENES_DF"]  = busco_genes_df

    # ------------------------------ Final assembly QUAST/BUSCO (robust discovery) ------------------------------
    
    # Candidate basenames (without extension) we can accept
    candidate_bases = []
    if pd.notna(ref_seq) and pd.isna(est_size):
        # QC-only mode: qc_assessment runs QUAST/BUSCO/Compleasm on the fetched
        # reference and writes them next to it, named
        # <sample_id>_EGAP_assembly_* in the reference's directory. Prefer that
        # base so the report finds the QC that actually ran; keep the
        # reference's own basename as a fallback.
        ref_dir = os.path.dirname(str(ref_seq))
        candidate_bases.append(os.path.join(ref_dir, f"{sample_id}_EGAP_assembly"))
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
    busco_db_details = []  # resolved lineage + DB version actually used (for the report)
    for busco_index, busco_db in enumerate([first_busco_db, second_busco_db]):
        if pd.isna(busco_db):
            continue
        busco_type = "FIRST" if busco_index == 0 else "SECOND"
    
        # This run's tool first; take the first folder with a usable summary.
        bases = ([assembly_base] if assembly_base else []) + [b for b in candidate_bases if b and b != assembly_base]
        busco_dir_candidates = [f"{base}_{busco_db}_{tool}" for base in bases
                                for tool in completeness_tool_order()]
    
        chosen_busco_dir, metrics = None, {}
        for candidate in busco_dir_candidates:
            if os.path.isdir(candidate):
                metrics = completeness_metrics(candidate)
                if metrics:
                    chosen_busco_dir = candidate
                    break
        if not chosen_busco_dir:
            print(f"WARN:\tNo BUSCO/Compleasm results found for {busco_db}. Tried: {busco_dir_candidates}")
            continue
    
        metrics_df = pd.DataFrame(list(metrics.items()), columns=["Metric", "Value"])
        busco_genes_df = load_busco_or_compleasm_genes_table(chosen_busco_dir, busco_db)

        sample_stats_dict[f"FINAL_{busco_type}_BUSCO_DF"]        = metrics_df
        sample_stats_dict[f"FINAL_{busco_type}_BUSCO_HTML"]      = os.path.join(chosen_busco_dir, f"{sample_id}_{busco_db}_EGAP_busco.html")
        sample_stats_dict[f"FINAL_{busco_type}_BUSCO_GENES_DF"]  = busco_genes_df

        # Record the lineage database version actually used (odbNN + creation date)
        _db_info = extract_busco_db_info(chosen_busco_dir, busco_db)
        _db_info["slot"] = "First" if busco_type == "FIRST" else "Second"
        busco_db_details.append(_db_info)

    # ------------------------------ Render ------------------------------
    os.makedirs(sample_dir, exist_ok=True)
    html_report_outpath = os.path.join(sample_dir, f"{sample_id}_EGAP_summary.html")

    print(f"DEBUG - templates_dir - {templates_dir}")
    env = Environment(
        loader=FileSystemLoader(str(templates_dir)),
        autoescape=select_autoescape(["html", "xml"])
    )
    main_tmpl  = env.get_template("EGAP_summary.html")
    busco_tmpl = env.get_template("EGAP_busco.html")

    # Helper to coerce DataFrame -> {metric: value}
    def df_to_map(df):
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
    inat_map       = df_to_map(inat_df)

    # ONT
    ont_raw_df  = sample_stats_dict.get("RAW_ONT_NANOSTATS_DF")
    ont_filt_df = sample_stats_dict.get("FILT_ONT_NANOSTATS_DF")
    ont_corr_df = sample_stats_dict.get("CORR_ONT_NANOSTATS_DF")
    ont_raw_map  = df_to_map(ont_raw_df)
    ont_filt_map = df_to_map(ont_filt_df)
    ont_corr_map = df_to_map(ont_corr_df)
    ont_raw_link = sample_stats_dict.get("RAW_ONT_NANOPLOT_HTML", "")
    ont_filt_link= sample_stats_dict.get("FILT_ONT_NANOPLOT_HTML", "")
    ont_corr_link= sample_stats_dict.get("CORR_ONT_NANOPLOT_HTML", "")

    # Illumina
    illumina_raw_f_df   = sample_stats_dict.get("RAW_ILLU_F_BASIC_STATS_DF")
    illumina_raw_r_df   = sample_stats_dict.get("RAW_ILLU_R_BASIC_STATS_DF")
    illumina_dedup_f_df = sample_stats_dict.get("DEDUP_ILLU_F_BASIC_STATS_DF")
    illumina_dedup_r_df = sample_stats_dict.get("DEDUP_ILLU_R_BASIC_STATS_DF")

    illumina_raw_f_map   = df_to_map(illumina_raw_f_df)
    illumina_raw_r_map   = df_to_map(illumina_raw_r_df)
    illumina_dedup_f_map = df_to_map(illumina_dedup_f_df)
    illumina_dedup_r_map = df_to_map(illumina_dedup_r_df)
    illumina_raw_f_link  = sample_stats_dict.get("RAW_ILLU_F_FASTQC_HTML", "")
    illumina_raw_r_link  = sample_stats_dict.get("RAW_ILLU_R_FASTQC_HTML", "")
    illumina_dedup_f_link= sample_stats_dict.get("DEDUP_ILLU_F_FASTQC_HTML", "")
    illumina_dedup_r_link= sample_stats_dict.get("DEDUP_ILLU_R_FASTQC_HTML", "")
    illumina_trimmed_f_map = df_to_map(sample_stats_dict.get("TRIMMED_ILLU_F_BASIC_STATS_DF"))
    illumina_trimmed_r_map = df_to_map(sample_stats_dict.get("TRIMMED_ILLU_R_BASIC_STATS_DF"))
    illumina_trimmed_f_link = sample_stats_dict.get("TRIMMED_ILLU_F_FASTQC_HTML", "")
    illumina_trimmed_r_link = sample_stats_dict.get("TRIMMED_ILLU_R_FASTQC_HTML", "")

    # Sequencing coverage per read set, from the stats CSV the final QC wrote
    coverage_table, coverage_basis = coverage_rows((assembly_base or "") + "_final_stats.csv")

    # PacBio
    pacbio_raw_df  = sample_stats_dict.get("RAW_PACBIO_NANOSTATS_DF")
    pacbio_filt_df = sample_stats_dict.get("FILT_PACBIO_NANOSTATS_DF")
    pacbio_raw_map  = df_to_map(pacbio_raw_df)
    pacbio_filt_map = df_to_map(pacbio_filt_df)
    pacbio_raw_link = sample_stats_dict.get("RAW_PACBIO_NANOPLOT_HTML", "")
    pacbio_filt_link= sample_stats_dict.get("FILT_PACBIO_NANOPLOT_HTML", "")

    # Per-assembler QUAST/BUSCO maps + links
    def maps(prefix):
        return (df_to_map(sample_stats_dict.get(f"{prefix}_QUAST_DF")),
                df_to_map(sample_stats_dict.get(f"{prefix}_FIRST_BUSCO_DF")),
                df_to_map(sample_stats_dict.get(f"{prefix}_SECOND_BUSCO_DF")),
                sample_stats_dict.get(f"{prefix}_QUAST_HTML",""),
                sample_stats_dict.get(f"{prefix}_FIRST_BUSCO_HTML",""),
                sample_stats_dict.get(f"{prefix}_SECOND_BUSCO_HTML",""))

    masurca_quast_map, masurca_first_busco_map, masurca_second_busco_map, \
        masurca_quast_link, masurca_first_busco_link, masurca_second_busco_link = maps("MASURCA")
    flye_quast_map, flye_first_busco_map, flye_second_busco_map, \
        flye_quast_link, flye_first_busco_link, flye_second_busco_link = maps("FLYE")
    spades_quast_map, spades_first_busco_map, spades_second_busco_map, \
        spades_quast_link, spades_first_busco_link, spades_second_busco_link = maps("SPADES")
    hifiasm_quast_map, hifiasm_first_busco_map, hifiasm_second_busco_map, \
        hifiasm_quast_link, hifiasm_first_busco_link, hifiasm_second_busco_link = maps("HIFIASM")

    final_quast_map        = df_to_map(sample_stats_dict.get("FINAL_QUAST_DF"))
    final_first_busco_map  = df_to_map(sample_stats_dict.get("FINAL_FIRST_BUSCO_DF"))
    final_second_busco_map = df_to_map(sample_stats_dict.get("FINAL_SECOND_BUSCO_DF"))
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
            metrics_df = sample_stats_dict.get(f"{asm_type}_{busco_type}_BUSCO_DF", pd.DataFrame())
            genes_df = sample_stats_dict.get(f"{asm_type}_{busco_type}_BUSCO_GENES_DF", pd.DataFrame())

            # Extract lineage (= busco_db) from filename: <sample_id>_<busco_db>_EGAP_busco.html
            fname = os.path.basename(outpath)
            busco_db_from_name = fname.replace(f"{sample_id}_", "").replace("_EGAP_busco.html", "")
            svg_path = find_busco_svg(os.path.dirname(outpath), busco_db_from_name)

            html = busco_tmpl.render(
                assembly_name=os.path.basename(outpath).replace("_EGAP_busco.html", ""),
                generated_time=now,
                busco_data=[{
                    "assembly_fasta": "",
                    "busco_db": str(busco_db_from_name),
                    "svg_path": svg_path,
                    "metrics": dict(zip(metrics_df["Metric"], metrics_df["Value"])) if not metrics_df.empty else {},
                    "busco_genes_df": genes_df
                }]
            )
            os.makedirs(os.path.dirname(outpath), exist_ok=True)
            with open(outpath, "w") as fh:
                fh.write(html)

    # For FINAL BUSCO
    for busco_type in ("FIRST", "SECOND"):
        metrics_df = sample_stats_dict.get(f"FINAL_{busco_type}_BUSCO_DF", pd.DataFrame())
        genes_df = sample_stats_dict.get(f"FINAL_{busco_type}_BUSCO_GENES_DF", pd.DataFrame())
        outpath = sample_stats_dict.get(f"FINAL_{busco_type}_BUSCO_HTML", "")
        if not outpath:
            continue

        # Extract lineage (= busco_db) from filename: <sample_id>_<busco_db>_EGAP_busco.html
        fname = os.path.basename(outpath)
        busco_db_from_name = fname.replace(f"{sample_id}_", "").replace("_EGAP_busco.html", "")
        svg_path = find_busco_svg(os.path.dirname(outpath), busco_db_from_name)

        html = busco_tmpl.render(
            assembly_name=os.path.basename(outpath).replace("_EGAP_busco.html", ""),
            generated_time=now,
            busco_data=[{
                "assembly_fasta": "",
                "busco_db": str(busco_db_from_name),
                "svg_path": svg_path,
                "metrics": dict(zip(metrics_df["Metric"], metrics_df["Value"])) if not metrics_df.empty else {},
                "busco_genes_df": genes_df
            }]
        )
        os.makedirs(os.path.dirname(outpath), exist_ok=True)
        with open(outpath, "w") as fh:
            fh.write(html)

    # -------- Render MAIN summary --------
    os.makedirs(sample_dir, exist_ok=True)
    html_report_outpath = os.path.join(sample_dir, f"{sample_id}_EGAP_summary.html")

    # -------- Reads / assembly file links (trace artefacts back to disk) --------
    dedup_f_path = os.path.join(illumina_dir, f"{species_id}_illu_forward_dedup.fastq")
    dedup_r_path = os.path.join(illumina_dir, f"{species_id}_illu_reverse_dedup.fastq")
    ont_highest_path = os.path.join(ont_dir, f"{species_id}_ONT_highest_mean_qual_long_reads.fastq")
    pacbio_highest_path = os.path.join(pacbio_dir, f"{species_id}_PacBio_highest_mean_qual_long_reads.fastq")
    reads_files = {
        "ont_raw": _file_link(ont_raw_reads if isinstance(ont_raw_reads, str) else None),
        "ont_highest": _file_link(ont_highest_path),
        "illumina_raw_f": _file_link(illumina_f_raw_reads if isinstance(illumina_f_raw_reads, str) else None),
        "illumina_raw_r": _file_link(illumina_r_raw_reads if isinstance(illumina_r_raw_reads, str) else None),
        "illumina_dedup_f": _file_link(dedup_f_path),
        "illumina_dedup_r": _file_link(dedup_r_path),
        "pacbio_raw": _file_link(pacbio_raw_reads if isinstance(pacbio_raw_reads, str) else None),
        "pacbio_highest": _file_link(pacbio_highest_path),
    }
    assembly_files = {
        "masurca": _file_link(os.path.join(masurca_dir, f"{sample_id}_masurca.fasta")),
        "flye": _file_link(os.path.join(flye_dir, f"{sample_id}_flye.fasta")),
        "spades": _file_link(os.path.join(spades_dir, f"{sample_id}_spades.fasta")),
        "hifiasm": _file_link(os.path.join(hifiasm_dir, f"{sample_id}_hifiasm.fasta")),
        "final": _file_link(assembly_fasta),
    }

    # -------- Pipeline provenance: steps, commands (settings), versions --------
    prov_source = "none"
    prov_steps = []
    prov_timings = []
    prov_programs = []
    prov_commands = []
    try:
        prov_entries, prov_source = load_provenance(output_dirto_abs, sample_id)
        prov_steps = [{"stage": s, "label": STAGE_LABELS.get(s, s)}
                      for s in pipeline_steps(prov_entries)]
        prov_timings = [{"label": STAGE_LABELS.get(t["stage"], t["stage"]),
                         "started": t.get("started", ""),
                         "duration": t["duration"],
                         "result": ("passed" if t.get("returncode") == 0
                                    else f"failed (exit code {t.get('returncode')})")}
                        for t in stage_timings(prov_entries)]
        _tools = meaningful_tools(prov_entries)
        # Prefer versions captured at run time (durable, match the run env);
        # probe/conda-resolve only the tools missing a captured version.
        _captured = versions_from_entries(prov_entries)
        _missing = [t for t in _tools if t not in _captured]
        _resolved = resolve_versions(_missing) if _missing else {}
        prov_programs = [{"tool": t, "version": _captured.get(t) or _resolved.get(t, "unknown")}
                         for t in _tools]
        prov_commands = [{"stage": e.get("stage", "unknown"),
                          "tool": e.get("tool", ""),
                          "command": e.get("command", "")}
                         for e in prov_entries if e.get("type") not in ("file", "stage")]
    except Exception as e:
        print(f"WARN:\tProvenance assembly failed: {e}")
    busco_dbs = [str(db) for db in (first_busco_db, second_busco_db) if pd.notna(db)]

    render_ctx = {
        "sample_id": sample_id,
        "generated_time": now,

        # Pipeline provenance + artefact links
        "prov_source": prov_source,
        "prov_steps": prov_steps,
        "prov_timings": prov_timings,
        "prov_programs": prov_programs,
        "prov_commands": prov_commands,
        "busco_dbs": busco_dbs,
        "busco_db_details": busco_db_details,
        "reads_files": reads_files,
        "assembly_files": assembly_files,

        "inat_metrics": inat_metrics, "inat_map": inat_map,
        "inat_link": sample_stats_dict.get("INAT_HTML"), "inat_thumbnail": inat_thumbnail, "inat_photo": inat_photo,

        "ont_metrics": ont_metrics,
        "ont_raw_map": ont_raw_map, "ont_filt_map": ont_filt_map, "ont_corr_map": ont_corr_map,
        "ont_raw_link": ont_raw_link, "ont_filt_link": ont_filt_link, "ont_corr_link": ont_corr_link,

        "illumina_metrics": illumina_metrics,
        "illumina_raw_f_map": illumina_raw_f_map, "illumina_raw_r_map": illumina_raw_r_map,
        "illumina_dedup_f_map": illumina_dedup_f_map, "illumina_dedup_r_map": illumina_dedup_r_map,
        "illumina_raw_f_link": illumina_raw_f_link, "illumina_raw_r_link": illumina_raw_r_link,
        "illumina_dedup_f_link": illumina_dedup_f_link, "illumina_dedup_r_link": illumina_dedup_r_link,
        "illumina_trimmed_f_map": illumina_trimmed_f_map, "illumina_trimmed_r_map": illumina_trimmed_r_map,
        "illumina_trimmed_f_link": illumina_trimmed_f_link, "illumina_trimmed_r_link": illumina_trimmed_r_link,
        "coverage_rows": coverage_table, "coverage_basis": coverage_basis,

        "pacbio_metrics": pacbio_metrics,
        "pacbio_raw_map": pacbio_raw_map, "pacbio_filt_map": pacbio_filt_map,
        "pacbio_raw_link": pacbio_raw_link, "pacbio_filt_link": pacbio_filt_link,

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
    run_stage(html_reporter)
