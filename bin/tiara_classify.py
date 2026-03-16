#!/usr/bin/env python3
"""
tiara_classify.py

Run Tiara on a genome assembly, optionally split output classes to FASTA,
then parse and visualize the Tiara TSV output in a single Plotly HTML report.

Used as an importable library module by decontaminate_assembly.py in the EGAP pipeline,
and can also be run standalone for ad-hoc classification and visualization.

Supports both Windows-style and native Linux/WSL paths:
    C:/Users/IanBollinger/Desktop/file.fasta
    /mnt/c/Users/IanBollinger/Desktop/file.fasta

Tiara command template:
    tiara --input INPUT_FASTA --output OUTPUT_TSV --min_len 3000 \\
          --prob_cutoff 0.65,0.65 --to_fasta all --threads 8 \\
          --probabilities --verbose --first_stage_kmer 6 \\
          --second_stage_kmer 7

First-stage classification labels:
    archaea, bacteria, prokarya, eukarya, organelle, unknown

Second-stage classification labels (for organelle-stage sequences):
    mitochondria, plastid, unknown, none

Created on Wed Feb 25 2026

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""

from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.io import to_html


VALID_TO_FASTA = {
    "mit",   # mitochondria
    "pla",   # plastid
    "bac",   # bacteria
    "arc",   # archaea
    "euk",   # eukarya
    "unk",   # unknown
    "pro",   # prokarya
    "all",   # all classes present in input fasta
}

VALID_CLASS_ORDER_FIRST = [
    "archaea",
    "bacteria",
    "prokarya",
    "eukarya",
    "organelle",
    "unknown",
]

VALID_CLASS_ORDER_SECOND = [
    "mitochondria",
    "plastid",
    "unknown",
    "none",
]


def windows_to_wsl_path(input_path: str) -> str:
    """
    Convert a Windows path to a WSL path if needed.

    Examples
    --------
    C:/Users/name/file.fasta -> /mnt/c/Users/name/file.fasta
    C:\\Users\\name\\file.fasta -> /mnt/c/Users/name/file.fasta

    If the path already appears to be a Linux path, it is returned unchanged.
    """
    if not input_path:
        raise ValueError("Input path is empty.")

    normalized = input_path.strip().replace("\\", "/")

    if normalized.startswith("/"):
        return normalized

    match_obj = re.match(r"^([A-Za-z]):/(.*)$", normalized)
    if match_obj:
        drive_letter = match_obj.group(1).lower()
        remainder = match_obj.group(2)
        return f"/mnt/{drive_letter}/{remainder}"

    return normalized


def ensure_parent_dir(file_path: str) -> None:
    """
    Create the parent directory for a file path if it does not already exist.
    """
    parent_dir = Path(file_path).parent
    parent_dir.mkdir(parents=True, exist_ok=True)


def ensure_dir(dir_path: str) -> None:
    """
    Create a directory if it does not already exist.
    """
    Path(dir_path).mkdir(parents=True, exist_ok=True)


def build_tiara_command(
    input_assembly: str,
    output_tsv: str,
    min_length_bp: int = 3000,
    probability_cutoff: str = "0.65,0.65",
    to_fasta: str = "all",
    cpu_threads: int = 8,
    first_stage_kmer: int = 6,
    second_stage_kmer: int = 7,
    tiara_executable: str = "tiara",
) -> List[str]:
    """
    Build the Tiara command as a subprocess-safe list.
    """
    if to_fasta not in VALID_TO_FASTA:
        raise ValueError(
            f"Invalid to_fasta value: {to_fasta}. "
            f"Valid options are: {', '.join(sorted(VALID_TO_FASTA))}"
        )

    if min_length_bp < 1:
        raise ValueError("min_length_bp must be >= 1.")

    if cpu_threads < 1:
        raise ValueError("cpu_threads must be >= 1.")

    if first_stage_kmer < 1 or second_stage_kmer < 1:
        raise ValueError("K-mer sizes must be >= 1.")

    cmd = [
        tiara_executable,
        "--input", input_assembly,
        "--output", output_tsv,
        "--min_len", str(min_length_bp),
        "--prob_cutoff", str(probability_cutoff),
        "--to_fasta", to_fasta,
        "--threads", str(cpu_threads),
        "--probabilities",
        "--verbose",
        "--first_stage_kmer", str(first_stage_kmer),
        "--second_stage_kmer", str(second_stage_kmer),
    ]
    return cmd


def run_subprocess_cmd(cmd: List[str], cwd: Optional[str] = None) -> int:
    """
    Run a subprocess command while streaming stdout/stderr in real time.

    Returns
    -------
    int
        Process return code.

    Raises
    ------
    RuntimeError
        If the subprocess exits with a non-zero return code.
    """
    printable_cmd = " ".join(shlex.quote(part) for part in cmd)
    print(f"\n[Tiara] Running command:\n{printable_cmd}\n")

    process = subprocess.Popen(
        cmd,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True,
    )

    if process.stdout is not None:
        for line in process.stdout:
            print(line, end="")

    process.wait()

    if process.returncode != 0:
        raise RuntimeError(
            f"Tiara failed with exit code {process.returncode}.\n"
            f"Command: {printable_cmd}"
        )

    return process.returncode


def collect_tiara_fasta_outputs(output_tsv: str) -> Dict[str, str]:
    """
    Collect Tiara-generated FASTA files in the same directory as the output TSV.

    Tiara creates FASTA files when --to_fasta is used. Naming can vary by class,
    so this function gathers FASTA-like files created alongside the TSV.

    Returns
    -------
    Dict[str, str]
        Dictionary mapping filename to full path.
    """
    output_dir = Path(output_tsv).parent
    output_stem = Path(output_tsv).stem
    fasta_dict: Dict[str, str] = {}

    for suffix in ("*.fa", "*.fasta", "*.fna"):
        for fasta_file in output_dir.glob(suffix):
            if fasta_file.name == Path(output_tsv).name:
                continue

            if output_stem in fasta_file.stem or fasta_file.stat().st_size > 0:
                fasta_dict[fasta_file.name] = str(fasta_file.resolve())

    return dict(sorted(fasta_dict.items()))


def simple_fasta_lengths(fasta_path: str) -> Dict[str, int]:
    """
    Parse a FASTA file and return {sequence_id: sequence_length}.
    """
    seq_len_dict: Dict[str, int] = {}
    current_id: Optional[str] = None
    current_len = 0

    with open(fasta_path, "r", encoding="utf-8", errors="replace") as in_handle:
        for raw_line in in_handle:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    seq_len_dict[current_id] = current_len

                current_id = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line)

    if current_id is not None:
        seq_len_dict[current_id] = current_len

    return seq_len_dict


def normalize_tiara_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize Tiara column names into predictable internal names.

    Required output columns:
        sequence_id
        first_stage_class
        second_stage_class
    """
    original_cols = list(df.columns)
    normalized_cols = []

    for col in original_cols:
        clean_col = (
            str(col)
            .strip()
            .lower()
            .replace(" ", "_")
            .replace("-", "_")
        )
        normalized_cols.append(clean_col)

    df.columns = normalized_cols
    rename_map = {}

    for col in df.columns:
        if col in {"sequence_id", "seq_id", "id"}:
            rename_map[col] = "sequence_id"
        elif "first" in col and "classification" in col:
            rename_map[col] = "first_stage_class"
        elif "second" in col and "classification" in col:
            rename_map[col] = "second_stage_class"

    df = df.rename(columns=rename_map)

    if "sequence_id" not in df.columns and len(df.columns) >= 1:
        df = df.rename(columns={df.columns[0]: "sequence_id"})

    if "first_stage_class" not in df.columns and len(df.columns) >= 2:
        df = df.rename(columns={df.columns[1]: "first_stage_class"})

    if "second_stage_class" not in df.columns and len(df.columns) >= 3:
        df = df.rename(columns={df.columns[2]: "second_stage_class"})

    required_cols = {"sequence_id", "first_stage_class", "second_stage_class"}
    missing_cols = required_cols - set(df.columns)

    if missing_cols:
        raise ValueError(
            "Could not normalize required Tiara columns. Missing: "
            f"{', '.join(sorted(missing_cols))}\n"
            f"Observed columns: {', '.join(map(str, original_cols))}"
        )

    df["sequence_id"] = df["sequence_id"].astype(str)
    df["first_stage_class"] = (
        df["first_stage_class"].astype(str).str.strip().str.lower()
    )
    df["second_stage_class"] = (
        df["second_stage_class"].astype(str).str.strip().str.lower()
    )
    df["second_stage_class"] = df["second_stage_class"].replace(
        {"nan": "none", "": "none", "na": "none", "none": "none"}
    )

    return df


def detect_probability_columns(df: pd.DataFrame) -> List[str]:
    """
    Detect likely probability columns among non-core Tiara columns.
    """
    core_cols = {"sequence_id", "first_stage_class", "second_stage_class"}
    prob_cols: List[str] = []

    for col in df.columns:
        if col in core_cols:
            continue

        numeric_series = pd.to_numeric(df[col], errors="coerce")
        valid_series = numeric_series.dropna()

        if len(valid_series) == 0:
            continue

        if ((valid_series >= 0) & (valid_series <= 1)).mean() >= 0.95:
            df[col] = numeric_series
            prob_cols.append(col)

    return prob_cols


def reorder_index(index_values: List[str], desired_order: List[str]) -> List[str]:
    """
    Put known classes first, preserving unexpected labels at the end.
    """
    observed = list(index_values)
    ordered = [item for item in desired_order if item in observed]
    extras = [item for item in observed if item not in ordered]
    return ordered + sorted(extras)


def save_dataframe(df: pd.DataFrame, out_path: str) -> None:
    """
    Save dataframe as TSV.
    """
    df.to_csv(out_path, sep="\t", index=False)


def build_plotly_report(
    df: pd.DataFrame,
    output_html: str,
    probability_columns: List[str],
) -> None:
    """
    Build a single self-contained HTML report with all applicable Plotly plots.
    """
    html_parts: List[str] = []

    html_header = """
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Tiara Classification Report</title>
<style>
body {
    font-family: Arial, sans-serif;
    margin: 24px;
    line-height: 1.5;
}
h1, h2, h3 {
    margin-top: 32px;
}
.section {
    margin-bottom: 36px;
}
.meta {
    background: #f5f5f5;
    padding: 12px 16px;
    border-radius: 8px;
    margin-bottom: 24px;
}
</style>
</head>
<body>
<h1>Tiara Classification Report</h1>
"""
    html_parts.append(html_header)

    html_parts.append(
        f"""
<div class="meta">
<p><strong>Total sequences:</strong> {len(df)}</p>
<p><strong>Columns:</strong> {", ".join(df.columns)}</p>
<p><strong>Detected probability columns:</strong> {", ".join(probability_columns) if probability_columns else "None detected"}</p>
</div>
"""
    )

    count_first = df["first_stage_class"].value_counts()
    count_first = count_first.reindex(
        reorder_index(count_first.index.tolist(), VALID_CLASS_ORDER_FIRST)
    ).dropna()

    if len(count_first) > 0:
        count_first_df = count_first.reset_index()
        count_first_df.columns = ["first_stage_class", "sequence_count"]
        fig = px.bar(
            count_first_df,
            x="first_stage_class",
            y="sequence_count",
            title="First-stage class counts",
        )
        html_parts.append('<div class="section"><h2>First-stage class counts</h2>')
        html_parts.append(to_html(fig, full_html=False, include_plotlyjs="cdn"))
        html_parts.append("</div>")

    count_second = df["second_stage_class"].value_counts()
    count_second = count_second.reindex(
        reorder_index(count_second.index.tolist(), VALID_CLASS_ORDER_SECOND)
    ).dropna()

    if len(count_second) > 0:
        count_second_df = count_second.reset_index()
        count_second_df.columns = ["second_stage_class", "sequence_count"]
        fig = px.bar(
            count_second_df,
            x="second_stage_class",
            y="sequence_count",
            title="Second-stage class counts",
        )
        html_parts.append('<div class="section"><h2>Second-stage class counts</h2>')
        html_parts.append(to_html(fig, full_html=False, include_plotlyjs=False))
        html_parts.append("</div>")

    matrix_df = pd.crosstab(
        df["first_stage_class"],
        df["second_stage_class"],
        dropna=False,
    )

    matrix_df = matrix_df.reindex(
        index=reorder_index(matrix_df.index.tolist(), VALID_CLASS_ORDER_FIRST),
        columns=reorder_index(matrix_df.columns.tolist(), VALID_CLASS_ORDER_SECOND),
        fill_value=0,
    )

    if not matrix_df.empty:
        fig = go.Figure(
            data=go.Heatmap(
                z=matrix_df.values,
                x=list(matrix_df.columns),
                y=list(matrix_df.index),
            )
        )
        fig.update_layout(title="First-stage vs second-stage classification matrix")
        html_parts.append(
            '<div class="section"><h2>First-stage vs second-stage matrix</h2>'
        )
        html_parts.append(to_html(fig, full_html=False, include_plotlyjs=False))
        html_parts.append("</div>")

    if "sequence_length_bp" in df.columns:
        length_df = df.copy()
        length_df["sequence_length_bp"] = pd.to_numeric(
            length_df["sequence_length_bp"],
            errors="coerce",
        )
        length_df = length_df[length_df["sequence_length_bp"].notna()].copy()

        if len(length_df) > 0:
            fig = px.histogram(
                length_df,
                x="sequence_length_bp",
                nbins=40,
                title="Sequence length distribution",
                log_x=True,
            )
            html_parts.append('<div class="section"><h2>Sequence length distribution</h2>')
            html_parts.append(to_html(fig, full_html=False, include_plotlyjs=False))
            html_parts.append("</div>")

            fig = px.box(
                length_df,
                x="first_stage_class",
                y="sequence_length_bp",
                title="Sequence length by first-stage class",
                points="outliers",
            )
            html_parts.append(
                '<div class="section"><h2>Sequence length by first-stage class</h2>'
            )
            html_parts.append(to_html(fig, full_html=False, include_plotlyjs=False))
            html_parts.append("</div>")

            bp_first = (
                length_df.groupby("first_stage_class")["sequence_length_bp"]
                .sum(min_count=1)
                .sort_values(ascending=False)
            )

            if len(bp_first) > 0:
                bp_first_df = bp_first.reset_index()
                bp_first_df.columns = ["first_stage_class", "total_bp"]
                fig = px.bar(
                    bp_first_df,
                    x="first_stage_class",
                    y="total_bp",
                    title="Total base pairs by first-stage class",
                )
                html_parts.append(
                    '<div class="section"><h2>Total base pairs by first-stage class</h2>'
                )
                html_parts.append(to_html(fig, full_html=False, include_plotlyjs=False))
                html_parts.append("</div>")

    for prob_col in probability_columns:
        plot_df = df.copy()
        plot_df[prob_col] = pd.to_numeric(plot_df[prob_col], errors="coerce")
        plot_df = plot_df[plot_df[prob_col].notna()].copy()

        if len(plot_df) == 0:
            continue

        fig = px.histogram(
            plot_df,
            x=prob_col,
            nbins=40,
            title=f"Probability distribution: {prob_col}",
        )
        html_parts.append(
            f'<div class="section"><h2>Probability distribution: {prob_col}</h2>'
        )
        html_parts.append(to_html(fig, full_html=False, include_plotlyjs=False))
        html_parts.append("</div>")

        fig = px.box(
            plot_df,
            x="first_stage_class",
            y=prob_col,
            title=f"{prob_col} by first-stage class",
            points="outliers",
        )
        html_parts.append(
            f'<div class="section"><h2>{prob_col} by first-stage class</h2>'
        )
        html_parts.append(to_html(fig, full_html=False, include_plotlyjs=False))
        html_parts.append("</div>")

        if "sequence_length_bp" in plot_df.columns:
            plot_df["sequence_length_bp"] = pd.to_numeric(
                plot_df["sequence_length_bp"],
                errors="coerce",
            )
            plot_df = plot_df[plot_df["sequence_length_bp"].notna()].copy()

            if len(plot_df) > 0:
                fig = px.scatter(
                    plot_df,
                    x="sequence_length_bp",
                    y=prob_col,
                    color="first_stage_class",
                    hover_data=["sequence_id"],
                    title=f"Sequence length vs probability: {prob_col}",
                    log_x=True,
                )
                html_parts.append(
                    f'<div class="section"><h2>Sequence length vs probability: {prob_col}</h2>'
                )
                html_parts.append(to_html(fig, full_html=False, include_plotlyjs=False))
                html_parts.append("</div>")

    html_parts.append("</body></html>")

    with open(output_html, "w", encoding="utf-8") as out_handle:
        out_handle.write("\n".join(html_parts))


def summarize_tiara(
    tiara_tsv: str,
    input_fasta: Optional[str] = None,
    output_dir: Optional[str] = None,
    output_html: Optional[str] = None,
) -> Dict[str, object]:
    """
    Parse Tiara TSV, optionally join FASTA lengths, generate summary TSVs,
    and build a single Plotly HTML report.
    """
    tiara_tsv = windows_to_wsl_path(tiara_tsv)

    if input_fasta is not None:
        input_fasta = windows_to_wsl_path(input_fasta)

    if output_dir is None:
        output_dir = str(Path(tiara_tsv).with_suffix("")) + "_viz"

    output_dir = windows_to_wsl_path(output_dir)
    ensure_dir(output_dir)

    if output_html is None:
        output_html = os.path.join(output_dir, "tiara_report.html")

    output_html = windows_to_wsl_path(output_html)

    if not os.path.exists(tiara_tsv):
        raise FileNotFoundError(f"Tiara TSV not found: {tiara_tsv}")

    df = pd.read_csv(tiara_tsv, sep="\t")
    df = normalize_tiara_columns(df)

    if input_fasta is not None:
        if not os.path.exists(input_fasta):
            raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")

        length_dict = simple_fasta_lengths(input_fasta)
        df["sequence_length_bp"] = df["sequence_id"].map(length_dict)

    prob_cols = detect_probability_columns(df)

    normalized_tsv = os.path.join(output_dir, "tiara_normalized.tsv")
    save_dataframe(df, normalized_tsv)

    generated_files = [normalized_tsv]

    count_first = df["first_stage_class"].value_counts()
    count_first = count_first.reindex(
        reorder_index(count_first.index.tolist(), VALID_CLASS_ORDER_FIRST)
    )
    count_first_df = count_first.reset_index()
    count_first_df.columns = ["first_stage_class", "sequence_count"]
    count_first_tsv = os.path.join(output_dir, "summary_first_stage_counts.tsv")
    save_dataframe(count_first_df, count_first_tsv)
    generated_files.append(count_first_tsv)

    count_second = df["second_stage_class"].value_counts()
    count_second = count_second.reindex(
        reorder_index(count_second.index.tolist(), VALID_CLASS_ORDER_SECOND)
    )
    count_second_df = count_second.reset_index()
    count_second_df.columns = ["second_stage_class", "sequence_count"]
    count_second_tsv = os.path.join(output_dir, "summary_second_stage_counts.tsv")
    save_dataframe(count_second_df, count_second_tsv)
    generated_files.append(count_second_tsv)

    matrix_df = pd.crosstab(
        df["first_stage_class"],
        df["second_stage_class"],
        dropna=False,
    )
    matrix_df = matrix_df.reindex(
        index=reorder_index(matrix_df.index.tolist(), VALID_CLASS_ORDER_FIRST),
        columns=reorder_index(matrix_df.columns.tolist(), VALID_CLASS_ORDER_SECOND),
        fill_value=0,
    )
    matrix_tsv = os.path.join(output_dir, "summary_first_vs_second_stage_matrix.tsv")
    matrix_df.to_csv(matrix_tsv, sep="\t")
    generated_files.append(matrix_tsv)

    if "sequence_length_bp" in df.columns:
        length_tsv = os.path.join(output_dir, "summary_lengths_with_classes.tsv")
        save_dataframe(df, length_tsv)
        generated_files.append(length_tsv)

        bp_first = (
            df.groupby("first_stage_class")["sequence_length_bp"]
            .sum(min_count=1)
            .sort_values(ascending=False)
        )
        bp_first_df = bp_first.reset_index()
        bp_first_df.columns = ["first_stage_class", "total_bp"]
        bp_first_tsv = os.path.join(output_dir, "summary_first_stage_total_bp.tsv")
        save_dataframe(bp_first_df, bp_first_tsv)
        generated_files.append(bp_first_tsv)

    if prob_cols:
        prob_summary_rows = []

        for prob_col in prob_cols:
            numeric_series = pd.to_numeric(df[prob_col], errors="coerce")
            valid_series = numeric_series.dropna()

            if len(valid_series) == 0:
                continue

            prob_summary_rows.append(
                {
                    "probability_column": prob_col,
                    "n": int(valid_series.notna().sum()),
                    "mean": float(valid_series.mean()),
                    "median": float(valid_series.median()),
                    "min": float(valid_series.min()),
                    "max": float(valid_series.max()),
                }
            )

        if prob_summary_rows:
            prob_summary_df = pd.DataFrame(prob_summary_rows)
            prob_summary_tsv = os.path.join(
                output_dir,
                "summary_probability_columns.tsv",
            )
            save_dataframe(prob_summary_df, prob_summary_tsv)
            generated_files.append(prob_summary_tsv)

    build_plotly_report(
        df=df,
        output_html=output_html,
        probability_columns=prob_cols,
    )
    generated_files.append(output_html)

    result_dict = {
        "tiara_tsv": tiara_tsv,
        "input_fasta": input_fasta,
        "output_dir": output_dir,
        "output_html": output_html,
        "probability_columns": prob_cols,
        "generated_files": generated_files,
    }

    print("\n[Tiara] Visualization complete.")
    print(f"[Tiara] Visualization directory: {output_dir}")
    print(f"[Tiara] HTML report: {output_html}")

    for file_path in generated_files:
        print(f"    {file_path}")

    return result_dict


def run_tiara(
    input_assembly: str,
    output_tsv: Optional[str] = None,
    min_length_bp: int = 3000,
    probability_cutoff: str = "0.65,0.65",
    to_fasta: str = "all",
    cpu_threads: int = 8,
    first_stage_kmer: int = 6,
    second_stage_kmer: int = 7,
    tiara_executable: str = "tiara",
    run_visualization: bool = True,
    visualization_output_dir: Optional[str] = None,
    visualization_output_html: Optional[str] = None,
) -> Dict[str, object]:
    """
    Run Tiara on an input assembly FASTA, then optionally parse and visualize
    the resulting TSV in a single Plotly HTML file.

    Parameters
    ----------
    input_assembly : str
        Input assembly FASTA file path.
    output_tsv : Optional[str]
        Output TSV path. If None, it will be derived from input_assembly.
    min_length_bp : int
        Minimum sequence length to classify. Use 1100 for Illumina-only
        assemblies (shorter contigs), 3000 for hybrid/long-read assemblies.
    probability_cutoff : str
        Tiara probability cutoff string, e.g. "0.65,0.65".
    to_fasta : str
        One of: mit, pla, bac, arc, euk, unk, pro, all.
    cpu_threads : int
        Number of CPU threads.
    first_stage_kmer : int
        Tiara first-stage k-mer length.
    second_stage_kmer : int
        Tiara second-stage k-mer length.
    tiara_executable : str
        Name or full path of the Tiara executable.
    run_visualization : bool
        Whether to parse the Tiara TSV and generate the HTML report.
    visualization_output_dir : Optional[str]
        Output directory for summary files and HTML report.
    visualization_output_html : Optional[str]
        Output HTML report path. If None, one is derived automatically.

    Returns
    -------
    Dict[str, object]
        Dictionary containing:
            input_assembly
            output_tsv
            fasta_outputs
            command
            visualization
    """
    input_assembly_wsl = windows_to_wsl_path(input_assembly)

    if output_tsv is None:
        if input_assembly.lower().endswith(".fasta"):
            output_tsv = input_assembly.replace(".fasta", "_tiara_out.tsv")
        elif input_assembly.lower().endswith(".fa"):
            output_tsv = input_assembly.replace(".fa", "_tiara_out.tsv")
        elif input_assembly.lower().endswith(".fna"):
            output_tsv = input_assembly.replace(".fna", "_tiara_out.tsv")
        else:
            output_tsv = f"{input_assembly}_tiara_out.tsv"

    output_tsv_wsl = windows_to_wsl_path(output_tsv)

    if not os.path.exists(input_assembly_wsl):
        raise FileNotFoundError(
            f"Input assembly was not found:\n{input_assembly_wsl}"
        )

    ensure_parent_dir(output_tsv_wsl)

    if shutil.which(tiara_executable) is None and not os.path.exists(tiara_executable):
        raise FileNotFoundError(
            f"Could not locate Tiara executable: {tiara_executable}\n"
            "Ensure Tiara is installed and available in your PATH."
        )

    cmd = build_tiara_command(
        input_assembly=input_assembly_wsl,
        output_tsv=output_tsv_wsl,
        min_length_bp=min_length_bp,
        probability_cutoff=probability_cutoff,
        to_fasta=to_fasta,
        cpu_threads=cpu_threads,
        first_stage_kmer=first_stage_kmer,
        second_stage_kmer=second_stage_kmer,
        tiara_executable=tiara_executable,
    )

    run_subprocess_cmd(cmd)

    if not os.path.exists(output_tsv_wsl):
        raise FileNotFoundError(
            "Tiara completed without raising a subprocess error, but the expected "
            f"TSV output was not found:\n{output_tsv_wsl}"
        )

    fasta_outputs = collect_tiara_fasta_outputs(output_tsv_wsl)

    visualization_result = None
    if run_visualization:
        visualization_result = summarize_tiara(
            tiara_tsv=output_tsv_wsl,
            input_fasta=input_assembly_wsl,
            output_dir=visualization_output_dir,
            output_html=visualization_output_html,
        )

    result_dict = {
        "input_assembly": input_assembly_wsl,
        "output_tsv": output_tsv_wsl,
        "fasta_outputs": fasta_outputs,
        "command": cmd,
        "visualization": visualization_result,
    }

    print("\n[Tiara] Completed successfully.")
    print(f"[Tiara] TSV output: {output_tsv_wsl}")

    if fasta_outputs:
        print("[Tiara] FASTA outputs:")
        for fasta_name, fasta_path in fasta_outputs.items():
            print(f"    {fasta_name}: {fasta_path}")

    return result_dict


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python tiara_classify.py <input_assembly.fasta> [output_tsv] [min_len] [threads]")
        print("\nExample:")
        print("  python tiara_classify.py /path/to/assembly.fasta")
        sys.exit(1)

    _input = sys.argv[1]
    _output_tsv = sys.argv[2] if len(sys.argv) > 2 else None
    _min_len = int(sys.argv[3]) if len(sys.argv) > 3 else 3000
    _threads = int(sys.argv[4]) if len(sys.argv) > 4 else 8

    run_tiara(
        input_assembly=_input,
        output_tsv=_output_tsv,
        min_length_bp=_min_len,
        probability_cutoff="0.65,0.65",
        to_fasta="all",
        cpu_threads=_threads,
        first_stage_kmer=6,
        second_stage_kmer=7,
        run_visualization=True,
        visualization_output_dir=None,
        visualization_output_html=None,
    )
