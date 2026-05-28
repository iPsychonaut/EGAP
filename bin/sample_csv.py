#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sample_csv.py

Per-sample CSV-row helpers, statistics dictionary, NanoPlot parsers,
long-read selection logic, and the typed ``SampleContext`` /
``AssemblerStage`` contracts.

Extracted from :mod:`utilities` in v3.4.1.  Hosts everything that
operates on a single row of the master sample-driver CSV.  Depends on
:mod:`file_operations` for :func:`to_abs`.

Stage:
    Cross-cutting (sample-row driver consumed by every pipeline stage)

Created on Wed Aug 16 2023

Updated on 2026-05-24

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os
import shutil
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Protocol, runtime_checkable

import pandas as pd

from file_operations import to_abs


# --------------------------------------------------------------
# Create a sample statistics dictionary from metadata
# --------------------------------------------------------------
def gen_sample_stats_dict(row):
    """Generate a sample statistics dictionary from a metadata row.

    Extracts key fields from a pandas Series and initializes placeholder
    values for all downstream quality metrics.

    Parameters
    ----------
    row : pandas.Series
        Single metadata row containing sample information columns.

    Returns
    -------
    dict
        Dictionary with all pipeline statistics fields initialized to
        ``None`` (except those extracted from *row*).
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


# --------------------------------------------------------------
# Extract and process sample metadata
# --------------------------------------------------------------
def get_current_row_data(input_df, sample_id):
    """Extract row data for a sample ID and generate a statistics dictionary.

    Filters *input_df* to the row matching *sample_id*, replaces any
    literal ``"None"`` strings with ``pd.NA`` so that downstream
    ``pd.isna()`` / ``pd.notna()`` guards work correctly, and builds
    the initial sample statistics dictionary.

    Parameters
    ----------
    input_df : pandas.DataFrame
        Full metadata DataFrame loaded from the input CSV.
    sample_id : str
        Sample identifier to filter on the ``SAMPLE_ID`` column.

    Returns
    -------
    tuple of (pandas.DataFrame, list, dict)
        ``(current_row, current_index, sample_stats_dict)`` - the filtered
        single-row DataFrame, its integer index list, and the initialized
        statistics dictionary.
    """
    # Filter the DataFrame for rows where the "SAMPLE_ID" column equals the provided sample_id
    current_row = input_df[input_df["SAMPLE_ID"] == sample_id].copy()

    # Replace literal string "None" with actual NaN so downstream pd.isna()
    # checks work correctly (some CSV editors write "None" instead of leaving
    # cells empty).
    current_row = current_row.replace(to_replace="None", value=pd.NA)

    sample_stats_dict = gen_sample_stats_dict(current_row)
    current_index = current_row.index.tolist()

    return current_row, current_index, sample_stats_dict


# --------------------------------------------------------------
# Per-sample context dataclass + loader
# --------------------------------------------------------------
@dataclass
class SampleContext:
    """Resolved per-sample inputs shared by every assembler / preprocessor.

    Bundles the artefacts that every ``bin/assemble_*.py`` and
    ``bin/preprocess_*.py`` module currently re-derives from the master
    sample-driver CSV: the single-row DataFrame, its first row as a
    ``pandas.Series``, the index list, and the initialized statistics
    dictionary.  Centralising the bundle behind a typed container makes
    the assembler/preprocessor contract checkable by tools like
    ``mypy --strict`` and removes the 6-line per-module re-derivation.

    Attributes
    ----------
    sample_id : str
        Sample identifier (matches ``SAMPLE_ID`` in *input_csv*).
    input_csv : str
        Absolute path to the master sample-driver CSV.
    output_dir : str
        Absolute path to the top-level results directory.
    cpu_threads : int
        CPU thread budget for this sample.
    ram_gb : int
        RAM budget (in GB) for this sample.
    current_row : pandas.DataFrame
        Single-row DataFrame filtered to the matching ``SAMPLE_ID``.
    current_series : pandas.Series
        ``current_row.iloc[0]`` - the row as a Series for ergonomic
        column access (``ctx.current_series["ILLUMINA_SRA"]``).
    current_index : list
        Integer index of *current_row* (``current_row.index.tolist()``).
    sample_stats_dict : dict
        Initialized per-sample statistics dictionary (populated as the
        pipeline advances).
    """

    sample_id: str
    input_csv: str
    output_dir: str
    cpu_threads: int
    ram_gb: int
    current_row: Any
    current_series: Any
    current_index: List[int] = field(default_factory=list)
    sample_stats_dict: Dict[str, Any] = field(default_factory=dict)


def load_sample_context(
    sample_id: str,
    input_csv: str,
    output_dir: str,
    cpu_threads: int,
    ram_gb: int,
) -> SampleContext:
    """Load and normalize the per-sample context expected by every stage.

    Performs the canonical first 6 lines of every assembler/preprocessor
    module: absolute-path normalization for *input_csv* / *output_dir*,
    ``pandas.read_csv`` of the driver, ``get_current_row_data`` lookup,
    and ``current_row.iloc[0]`` extraction.

    Parameters
    ----------
    sample_id : str
        Sample identifier matching the ``SAMPLE_ID`` column of *input_csv*.
    input_csv : str
        Path to the master sample-driver CSV (will be normalized to
        absolute).
    output_dir : str
        Path to the top-level results directory (will be normalized to
        absolute).
    cpu_threads : int
        CPU thread budget passed through into the returned context.
    ram_gb : int
        RAM budget (GB) passed through into the returned context.

    Returns
    -------
    SampleContext
        Fully populated context dataclass.
    """
    input_csv_abs = to_abs(input_csv)
    output_dir_abs = to_abs(output_dir)
    input_df = pd.read_csv(input_csv_abs)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]
    return SampleContext(
        sample_id=sample_id,
        input_csv=input_csv_abs,
        output_dir=output_dir_abs,
        cpu_threads=cpu_threads,
        ram_gb=ram_gb,
        current_row=current_row,
        current_series=current_series,
        current_index=current_index,
        sample_stats_dict=sample_stats_dict,
    )


@runtime_checkable
class AssemblerStage(Protocol):
    """Common entry-point signature for every EGAP assembler / preprocessor.

    Implementations live in ``bin/assemble_*.py`` and
    ``bin/preprocess_*.py``.  Every implementation accepts the same
    five-argument tuple the orchestrator hands out per sample, and
    returns a path to the produced FASTA/FASTQ (or ``None`` when no
    compatible inputs were available for that stage).

    Notes
    -----
    Declared with :func:`typing.runtime_checkable` so callers may use
    ``isinstance(fn, AssemblerStage)`` for diagnostics, but the contract
    is primarily intended for static type-checkers
    (``mypy --strict``).  Adding ``-> AssemblerStage`` annotations to
    the orchestrator's dispatch table is the cheapest way to lock in
    the convention identified by the v3.4.1 graphify pass.

    Parameters (of ``__call__``)
    ----------------------------
    sample_id : str
        Identifier matching ``SAMPLE_ID`` in *input_csv*.
    input_csv : str
        Path to the master sample-driver CSV.
    output_dir : str
        Top-level results directory; implementations anchor under it.
    cpu_threads : int
        Number of CPU threads the stage may use.
    ram_gb : int
        RAM budget in gigabytes.

    Returns
    -------
    Optional[str]
        Path to the produced artefact, or ``None`` when no compatible
        reads/inputs were available for this stage.
    """

    def __call__(
        self,
        sample_id: str,
        input_csv: str,
        output_dir: str,
        cpu_threads: int,
        ram_gb: int,
    ) -> Optional[str]: ...


# --------------------------------------------------------------
# Parse NanoPlot statistics for long reads
# --------------------------------------------------------------
def analyze_nanostats(READS_ORIGIN, nanoplot_out_file, sample_stats_dict):
    """Parse NanoPlot statistics and update the sample statistics dictionary.

    Reads NanoPlot output and extracts metrics based on read origin.

    Parameters
    ----------
    READS_ORIGIN : str
        Type and stage of reads (e.g., ``'Raw_ONT'``, ``'Filt_PacBio'``).
        Used to determine which keys in *sample_stats_dict* to populate.
    nanoplot_out_file : str
        Path to the NanoPlot ``NanoStats.txt`` output file.
    sample_stats_dict : dict
        Statistics dictionary to update in-place.

    Returns
    -------
    dict
        The updated *sample_stats_dict* with NanoPlot metrics filled in.
    """
    with open(nanoplot_out_file, "r") as nanostats:
        if "raw" in READS_ORIGIN.lower() and "ont" in READS_ORIGIN.lower():
            for line in nanostats:
                if "Number of reads:" in line:
                    sample_stats_dict["RAW_ONT_READS"] = float(line.split(":")[-1].replace(" ", "").replace(",", "").replace("\n", ""))
                elif "Mean read length:" in line:
                    sample_stats_dict["RAW_ONT_MEAN_LENGTH"] = float(line.split(":")[-1].replace(" ", "").replace(",", "").replace("\n", ""))
                elif "Mean read quality:" in line:
                    sample_stats_dict["RAW_ONT_MEAN_QUAL"] = float(line.split(":")[-1].replace(" ", "").replace(",", "").replace("\n", ""))
                elif "Total bases:" in line:
                    sample_stats_dict["RAW_ONT_TOTAL_BASES"] = float(line.split(":")[-1].replace(" ", "").replace(",", "").replace("\n", ""))
        elif "filt" in READS_ORIGIN.lower() and "ont" in READS_ORIGIN.lower():
            for line in nanostats:
                if "Number of reads:" in line:
                    sample_stats_dict["FILT_ONT_READS"] = float(line.split(":")[-1].replace(" ", "").replace(",", "").replace("\n", ""))
                elif "Mean read length:" in line:
                    sample_stats_dict["FILT_ONT_MEAN_LENGTH"] = float(line.split(":")[-1].replace(" ", "").replace(",", "").replace("\n", ""))
                elif "Mean read quality:" in line:
                    sample_stats_dict["FILT_ONT_MEAN_QUAL"] = float(line.split(":")[-1].replace(" ", "").replace(",", "").replace("\n", ""))
                elif "Total bases:" in line:
                    sample_stats_dict["FILT_ONT_TOTAL_BASES"] = float(line.split(":")[-1].replace(" ", "").replace(",", "").replace("\n", ""))
        elif "cor" in READS_ORIGIN.lower() and "ont" in READS_ORIGIN.lower():
            for line in nanostats:
                if "Number of reads:" in line:
                    sample_stats_dict["CORRECT_ONT_READS"] = float(line.split(":")[-1].replace(" ", "").replace(",", "").replace("\n", ""))
                elif "Mean read length:" in line:
                    sample_stats_dict["CORRECT_ONT_MEAN_LENGTH"] = float(line.split(":")[-1].replace(" ", "").replace(",", "").replace("\n", ""))
                elif "Mean read quality:" in line:
                    sample_stats_dict["CORRECT_ONT_MEAN_QUAL"] = float(line.split(":")[-1].replace(" ", "").replace(",", "").replace("\n", ""))
                elif "Total bases:" in line:
                    sample_stats_dict["CORRECT_ONT_TOTAL_BASES"] = float(line.split(":")[-1].replace(" ", "").replace(",", "").replace("\n", ""))
        elif "raw" in READS_ORIGIN.lower() and "pacbio" in READS_ORIGIN.lower():
            for line in nanostats:
                if "Number of reads:" in line:
                    sample_stats_dict["RAW_PACBIO_READS"] = float(line.split(":")[-1].replace(" ", "").replace(",", "").replace("\n", ""))
                elif "Mean read length:" in line:
                    sample_stats_dict["RAW_PACBIO_MEAN_LENGTH"] = float(line.split(":")[-1].replace(" ", "").replace(",", "").replace("\n", ""))
                elif "Mean read quality:" in line:
                    sample_stats_dict["RAW_PACBIO_MEAN_QUAL"] = float(line.split(":")[-1].replace(" ", "").replace(",", "").replace("\n", ""))
                elif "Total bases:" in line:
                    sample_stats_dict["RAW_PACBIO_TOTAL_BASES"] = float(line.split(":")[-1].replace(" ", "").replace(",", "").replace("\n", ""))
        elif "filt" in READS_ORIGIN.lower() and "pacbio" in READS_ORIGIN.lower():
            for line in nanostats:
                if "Number of reads:" in line:
                    sample_stats_dict["FILT_PACBIO_READS"] = float(line.split(":")[-1].replace(" ", "").replace(",", "").replace("\n", ""))
                elif "Mean read length:" in line:
                    sample_stats_dict["FILT_PACBIO_MEAN_LENGTH"] = float(line.split(":")[-1].replace(" ", "").replace(",", "").replace("\n", ""))
                elif "Mean read quality:" in line:
                    sample_stats_dict["FILT_PACBIO_MEAN_QUAL"] = float(line.split(":")[-1].replace(" ", "").replace(",", "").replace("\n", ""))
                elif "Total bases:" in line:
                    sample_stats_dict["FILT_PACBIO_TOTAL_BASES"] = float(line.split(":")[-1].replace(" ", "").replace(",", "").replace("\n", ""))
    return sample_stats_dict


# --------------------------------------------------------------
# Select highest quality long reads
# --------------------------------------------------------------
def select_long_reads(output_dir, input_csv, sample_id, cpu_threads):
    """Select the highest-mean-quality long reads from ONT or PacBio data.

    Parses NanoPlot statistics to compare quality across raw, filtered, and
    corrected read sets, then copies the best set to a canonical
    ``*_highest_mean_qual_long_reads.fastq`` path.

    Parameters
    ----------
    output_dir : str
        Root output directory containing per-species subdirectories.
    input_csv : str
        Path to the metadata CSV file.
    sample_id : str
        Sample identifier used to look up the row in *input_csv*.
    cpu_threads : int
        Number of threads available for compression tasks.

    Returns
    -------
    str or None
        Absolute path to the selected highest-quality reads file, or
        ``None`` if no suitable file can be found.
    """
    input_df = pd.read_csv(input_csv)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]  # Convert to Series (single row)

    # Identify read paths, reference, and BUSCO lineage info from CSV
    ont_sra = current_series["ONT_SRA"]
    ont_raw_reads = current_series["ONT_RAW_READS"]
    pacbio_raw_reads = current_series["PACBIO_RAW_READS"]
    pacbio_sra = current_series["PACBIO_SRA"]
    species_id = current_series["SPECIES_ID"]

    print(f"DEBUG - species_id - {species_id}")
    print(f"DEBUG - sample_id - {sample_id}")

    if pd.notna(ont_sra) and pd.isna(ont_raw_reads):
        ont_raw_reads = os.path.join(output_dir, species_id, "ONT", f"{ont_sra}.fastq")
    if pd.notna(pacbio_sra) and pd.isna(pacbio_raw_reads):
        pacbio_raw_reads = os.path.join(output_dir, species_id, "PacBio", f"{pacbio_sra}.fastq")

    print(f"DEBUG - ont_raw_reads - {ont_raw_reads}")
    print(f"DEBUG - pacbio_raw_reads - {pacbio_raw_reads}")

    if pd.notna(ont_raw_reads):
        print("DEBUG - PROCESSING ONT HIGHEST MEAN QUAL")
        reads_type = "ONT"
        reads_dir = os.path.dirname(ont_raw_reads)
        filtered_reads = os.path.join(reads_dir, f"{species_id}_{reads_type}_filtered.fastq")
        corrected_reads = os.path.join(reads_dir, f"{species_id}_{reads_type}_corrected.fastq")
        reads_origin_list = ["Raw_ONT_", "Filt_ONT_", "Corr_ONT_"]
    elif pd.notna(pacbio_raw_reads):
        print("DEBUG - PROCESSING PACBIO HIGHEST MEAN QUAL")
        reads_type = "PacBio"
        reads_dir = os.path.dirname(pacbio_raw_reads)
        filtered_reads = os.path.join(reads_dir, f"{species_id}_{reads_type}_filtered.fastq")
        corrected_reads = os.path.join(reads_dir, f"{species_id}_{reads_type}_corrected.fastq")
        reads_origin_list = ["Raw_PacBio_", "Filt_PacBio_"]
    else:
        print(f"ERROR:\tUNABLE TO PARSE LONG READS AS BOTH ONT AND PACBIO RAW READS ARE NONE: {ont_raw_reads} & {pacbio_raw_reads}")

    for reads_origin in reads_origin_list:
        sample_stats_dict = analyze_nanostats(reads_origin, os.path.join(reads_dir, f"{reads_origin}nanoplot_analysis", f"{reads_origin}NanoStats.txt"), sample_stats_dict)

    if pd.notna(ont_raw_reads):
        print("Selecting Highest Mean Quality Long reads...")
        highest_mean_qual_long_reads = corrected_reads
        if pd.notna(ont_raw_reads) and sample_stats_dict["CORRECT_ONT_MEAN_QUAL"] < sample_stats_dict["FILT_ONT_MEAN_QUAL"]:
            highest_mean_qual_long_reads = filtered_reads
            highest_mean_qual = sample_stats_dict["FILT_ONT_MEAN_QUAL"]
        else:
            highest_mean_qual = sample_stats_dict["CORRECT_ONT_MEAN_QUAL"]
    if pd.notna(pacbio_raw_reads):
        print("Selecting Highest Mean Quality Long reads...")
        highest_mean_qual_long_reads = pacbio_raw_reads
        if pd.notna(pacbio_raw_reads) and sample_stats_dict["RAW_PACBIO_MEAN_QUAL"] < sample_stats_dict["FILT_PACBIO_MEAN_QUAL"]:
            highest_mean_qual_long_reads = filtered_reads
            highest_mean_qual = sample_stats_dict["FILT_PACBIO_MEAN_QUAL"]
        else:
            highest_mean_qual = sample_stats_dict["RAW_ONT_MEAN_QUAL"]
    print(f"Highest Mean Quality Long reads: {highest_mean_qual_long_reads}")
    print(f"Mean Quality: {highest_mean_qual}")

    renamed_highest_mean_qual_long_reads = f"{species_id}_{reads_type}_highest_mean_qual_long_reads.fastq"
    if not os.path.exists(highest_mean_qual_long_reads):
        # try fallback: see if it's named like "Escherichia_coli_filtered.fastq"
        fallback_file = os.path.join(reads_dir, f"{species_id}_filtered.fastq")
        if os.path.exists(fallback_file):
            print(f"FALLBACK:\tFound fallback filtered file: {fallback_file}")
            highest_mean_qual_long_reads = fallback_file
        else:
            print("ERROR:\tNo usable highest-mean-quality long read file found.")
            return None

    renamed_highest_mean_qual_long_reads = os.path.join(reads_dir, f"{species_id}_{reads_type}_highest_mean_qual_long_reads.fastq")
    shutil.copy(highest_mean_qual_long_reads, renamed_highest_mean_qual_long_reads)

    print(f"NOTE:\tSelected highest quality long reads: {renamed_highest_mean_qual_long_reads} with mean quality {highest_mean_qual}")

    return renamed_highest_mean_qual_long_reads
