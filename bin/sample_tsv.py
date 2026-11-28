#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sample_tsv.py

Per-sample TSV-row helpers, statistics dictionary, NanoPlot parsers,
long-read selection logic, and the typed ``SampleContext`` contract.

Extracted from :mod:`utilities` in v3.4.1.  Hosts everything that
operates on a single row of the master sample-driver TSV.  Depends on
:mod:`file_operations` for :func:`to_abs`.

Stage:
    Cross-cutting (sample-row driver consumed by every pipeline stage)

Created on Wed Aug 16 2023

Updated on 2026-05-24

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import glob
import json
import os
import shutil
import sys
import zipfile
from dataclasses import dataclass, field
from typing import Any, Dict, List, NamedTuple

import pandas as pd

from file_operations import to_abs
from log import initialize_logging_environment


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
                         "COVERAGE_BASIS": None,
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
        Full metadata DataFrame loaded from the input TSV.
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

    # Replace literal string "None" with pd.NA so downstream pd.isna()
    # checks work correctly (some TSV editors write "None" instead of leaving
    # cells empty). Use mask() rather than DataFrame.replace(): on pandas
    # 1.4.x, replace() raises AttributeError when an object column already
    # holds pd.NA (the elementwise comparison collapses to a scalar bool).
    # mask() behaves identically on pandas 1.4.4, 2.0.3 and 3.x.
    current_row = current_row.mask(current_row == "None", other=pd.NA)

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
    sample-driver TSV: the single-row DataFrame, its first row as a
    ``pandas.Series``, the index list, and the initialized statistics
    dictionary.  Centralising the bundle behind a typed container makes
    the assembler/preprocessor contract checkable by tools like
    ``mypy --strict`` and removes the 6-line per-module re-derivation.

    Attributes
    ----------
    sample_id : str
        Sample identifier (matches ``SAMPLE_ID`` in *input_tsv*).
    input_tsv : str
        Absolute path to the master sample-driver TSV.
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
    input_tsv: str
    output_dir: str
    cpu_threads: int
    ram_gb: int
    current_row: Any
    current_series: Any
    current_index: List[int] = field(default_factory=list)
    sample_stats_dict: Dict[str, Any] = field(default_factory=dict)


# --------------------------------------------------------------
# The sample table's columns: one specification for every use
# --------------------------------------------------------------
SAMPLE_SHEET = "Sample Information"
METADATA_SHEET = "Sample Metadata"
SETTINGS_SHEET = "Sample Settings"   # reserved; see bin/sample_workbook.py


class Column(NamedTuple):
    """One column of the sample input: what it is called, where it lives, how to fill it."""
    name: str
    sheet: str               # SAMPLE_SHEET or METADATA_SHEET (a TSV holds both side by side)
    required: str            # "yes", "no", or the condition in words
    format: str              # expected cell format, shown in the workbook tooltip
    help: str                # what the value is for, shown in the workbook tooltip
    choices: tuple = ()      # allowed values, when the column is an enumeration
    expands_to: tuple = ()   # per-stage column names a comma-joined cell is split into


_ONE_OF_READS = "one of the SRA, RAW_DIR or RAW_READS columns per read type you have"
_KEYED = "yes; must match a SAMPLE_ID on the Sample Information sheet"

SAMPLE_COLUMNS = (
    # ---- Sample Information: the ground truth for which samples run ----
    Column("SPECIES_ID", SAMPLE_SHEET, "yes", "Genus_species, underscores, no spaces",
           "Species name. Names the output folder that samples of one species share."),
    Column("SAMPLE_ID", SAMPLE_SHEET, "yes", "<SPECIES_ID>-<identifier>, unique in the table",
           "Sample name. Names the sample's output folder, log and report."),
    Column("ORGANISM_KINGDOM", SAMPLE_SHEET, "yes", "one of the listed values",
           "Kingdom of the organism. Sets what Kraken2 and Tiara decontamination keep.",
           choices=("Bacteria", "Archaea", "Flora", "Funga", "Fauna")),
    Column("ORGANISM_KARYOTE", SAMPLE_SHEET, "yes", "one of the listed values",
           "Cell type of the organism. Sets QUAST and decontamination behaviour.",
           choices=("prokaryote", "eukaryote")),
    Column("PLOIDY", SAMPLE_SHEET, "no", "integer (1, 2, ...) or haploid / diploid",
           "Ploidy of the individual. Haploid samples skip the purge_dups haplotig removal."),
    Column("EST_SIZE", SAMPLE_SHEET, "yes, unless the sample is QC-only",
           "number followed by k, m or g (5m, 55m, 1.2g)",
           "Estimated genome size. Leave blank only to run QC on an existing assembly."),
    Column("BUSCO", SAMPLE_SHEET, "yes", "lineage1,lineage2 (agaricales,basidiomycota)",
           "The two Compleasm/BUSCO lineages completeness is measured against.",
           expands_to=("BUSCO_1", "BUSCO_2")),
    Column("ONT_SRA", SAMPLE_SHEET, _ONE_OF_READS, "SRA run accession (SRR...)",
           "Oxford Nanopore reads to download from SRA."),
    Column("ONT_RAW_DIR", SAMPLE_SHEET, _ONE_OF_READS, "folder path",
           "Folder of raw ONT FASTQ files to merge. Leave ONT_RAW_READS blank when set."),
    Column("ONT_RAW_READS", SAMPLE_SHEET, _ONE_OF_READS, "path to one .fastq or .fastq.gz",
           "Combined raw ONT reads."),
    Column("ILLUMINA_SRA", SAMPLE_SHEET, _ONE_OF_READS, "SRA run accession (SRR...)",
           "Paired-end Illumina reads to download from SRA."),
    Column("ILLUMINA_RAW_DIR", SAMPLE_SHEET, _ONE_OF_READS, "folder path",
           "Folder of raw paired-end Illumina FASTQ files. Leave ILLUMINA_RAW_READS blank when set."),
    Column("ILLUMINA_RAW_READS", SAMPLE_SHEET, _ONE_OF_READS,
           "forward_1.fastq,reverse_2.fastq (file names end _1 and _2; folder names must not)",
           "Raw forward and reverse Illumina reads.",
           expands_to=("ILLUMINA_RAW_F_READS", "ILLUMINA_RAW_R_READS")),
    Column("PACBIO_SRA", SAMPLE_SHEET, _ONE_OF_READS, "SRA run accession (SRR...)",
           "PacBio reads to download from SRA."),
    Column("PACBIO_RAW_DIR", SAMPLE_SHEET, _ONE_OF_READS, "folder path",
           "Folder of raw PacBio FASTQ files to merge. Leave PACBIO_RAW_READS blank when set."),
    Column("PACBIO_RAW_READS", SAMPLE_SHEET, _ONE_OF_READS, "path to one .fastq or .fastq.gz",
           "Combined raw PacBio reads."),
    Column("REF_SEQ_GCA", SAMPLE_SHEET, "no", "assembly accession (GCA_... or GCF_...)",
           "Reference genome to download: used for scaffolding and QC, or as the QC-only target."),
    Column("REF_SEQ", SAMPLE_SHEET, "no", "path to a FASTA file",
           "Local reference genome, with the same role as REF_SEQ_GCA."),
    Column("NT_ASSEMBLY_GCA", SAMPLE_SHEET, "no", "assembly accession (GCA_... or GCF_...)",
           "An existing assembly to carry into annotation, or to run QC on when no reads are given."),
    Column("NT_ASSEMBLY_PATH", SAMPLE_SHEET, "no", "path to a FASTA file",
           "Local existing assembly, with the same role as NT_ASSEMBLY_GCA."),
    Column("TX_EVI_GCA", SAMPLE_SHEET, "no", "accession, or a comma-separated list",
           "Transcript evidence for the annotation pipeline. EGAP carries it through unread."),
    Column("TX_EVI_PATH", SAMPLE_SHEET, "no", "path, or a comma-separated list",
           "Local transcript evidence for the annotation pipeline. EGAP carries it through unread."),
    Column("PT_EVI_GCA", SAMPLE_SHEET, "no", "accession, or a comma-separated list",
           "Protein evidence for the annotation pipeline. EGAP carries it through unread."),
    Column("PT_EVI_PATH", SAMPLE_SHEET, "no", "path, or a comma-separated list",
           "Local protein evidence for the annotation pipeline. EGAP carries it through unread."),
    Column("KRAKEN2_DB", SAMPLE_SHEET, "no", "folder path",
           "Kraken2 database for read decontamination. The KRAKEN2_DB environment variable "
           "wins when set; with neither, decontamination is skipped with a warning."),
    # ---- Sample Metadata: what the NCBI submission sheets need ----
    Column("SAMPLE_ID", METADATA_SHEET, _KEYED, "copy of a SAMPLE_ID from the first sheet",
           "Links this row to its sample. A row with no matching sample stops the run."),
    Column("BIOSAMPLE_ACCESSION", METADATA_SHEET, "no", "BioSample accession (SAMN...)",
           "BioSample the reads and the assembly belong to. Left blank in the sheets when empty."),
    Column("SAMPLE_TISSUE", METADATA_SHEET, "yes, for the SRA sheet", "free text (mycelium, fruiting body)",
           "Tissue or material sequenced. Goes into the SRA library title."),
    Column("SEQUENCING_METHODOLOGY", METADATA_SHEET, "yes, for the SRA sheet",
           "an SRA library strategy (WGS, WGA, RNA-Seq, ...)",
           "SRA library_strategy for this sample's libraries."),
    Column("LIBRARY_SOURCE", METADATA_SHEET, "no", "one of the listed values",
           "SRA library_source.",
           choices=("GENOMIC", "GENOMIC SINGLE CELL", "TRANSCRIPTOMIC", "TRANSCRIPTOMIC SINGLE CELL",
                    "METAGENOMIC", "METATRANSCRIPTOMIC", "SYNTHETIC", "VIRAL RNA", "OTHER")),
    Column("LIBRARY_SELECTION", METADATA_SHEET, "no", "an SRA library selection (RANDOM, PCR, ...)",
           "SRA library_selection: how the fragments were selected or enriched."),
    Column("INSTRUMENT_ILLUMINA", METADATA_SHEET, "no", "SRA instrument model (Illumina NovaSeq 6000)",
           "Sequencer that produced the Illumina reads."),
    Column("INSTRUMENT_ONT", METADATA_SHEET, "no", "SRA instrument model (MinION, PromethION)",
           "Sequencer that produced the ONT reads."),
    Column("INSTRUMENT_PACBIO", METADATA_SHEET, "no", "SRA instrument model (Sequel II, Revio)",
           "Sequencer that produced the PacBio reads."),
    Column("INATURALIST_ID", METADATA_SHEET, "no", "observation number from an iNaturalist URL",
           "iNaturalist observation of the specimen. Adds collector, date, place and photo to the report."),
)

# Columns every per-stage accessor expects to read off a row, derived from the
# specification: a collapsed cell contributes the names it expands to.
# ``read_sample_table`` adds any that are missing as NA, so stages can read
# ``current_series["BUSCO_1"]`` or ``["SAMPLE_TISSUE"]`` without ``KeyError``
# on a sparse or older table.
_LEGACY_COLUMNS = tuple(dict.fromkeys(
    name for column in SAMPLE_COLUMNS for name in (column.expands_to or (column.name,))))


def _split_collapsed_cell(value):
    """Split a comma-joined TSV cell into a list of stripped tokens.

    Returns ``[]`` for a blank/``NaN``/``None`` cell. A single value yields a
    one-element list.
    """
    if value is None:
        return []
    if isinstance(value, float) and pd.isna(value):
        return []
    text = str(value).strip()
    if text == "" or text.lower() in ("none", "nan"):
        return []
    return [tok.strip() for tok in text.split(",") if tok.strip() != ""]


def _expand_collapsed_pair(df, source_col, first_col, second_col):
    """Expand a collapsed ``source_col`` (``a,b``) into two legacy columns.

    Fills *first_col*/*second_col* from *source_col* only when the legacy
    columns are absent, so a table already shipping the expanded columns is
    left untouched. The collapsed column is dropped afterward.
    """
    if source_col not in df.columns:
        return df
    firsts, seconds = [], []
    for value in df[source_col]:
        tokens = _split_collapsed_cell(value)
        firsts.append(tokens[0] if len(tokens) >= 1 else pd.NA)
        seconds.append(tokens[1] if len(tokens) >= 2 else pd.NA)
    if first_col not in df.columns:
        df[first_col] = firsts
    if second_col not in df.columns:
        df[second_col] = seconds
    return df.drop(columns=[source_col])


def read_sample_table(input_path):
    """Read the tab-separated sample-driver table into a normalized DataFrame.

    EGAP's sample table moved from TSV to TSV: two fields are stored collapsed
    as comma-separated lists -- ``ILLUMINA_RAW_READS`` (forward,reverse) and
    ``BUSCO`` (lineage1,lineage2). This loader re-expands them into the legacy
    ``ILLUMINA_RAW_F_READS``/``_R_READS`` and ``BUSCO_1``/``BUSCO_2`` columns
    every downstream stage reads, and guarantees the full legacy/annotation
    column set is present, so the rest of the pipeline is untouched by the
    format change.

    Parameters
    ----------
    input_path : str
        Path to the master sample-driver TSV.

    Returns
    -------
    pandas.DataFrame
        Normalized table with the legacy per-stage column names present.
    """
    df = pd.read_csv(input_path, sep="\t")
    for column in SAMPLE_COLUMNS:
        if column.expands_to:
            df = _expand_collapsed_pair(df, column.name, *column.expands_to)
    for col in _LEGACY_COLUMNS:
        if col not in df.columns:
            df[col] = pd.NA
    return df


def load_sample_context(
    sample_id: str,
    input_tsv: str,
    output_dir: str,
    cpu_threads: int,
    ram_gb: int,
) -> SampleContext:
    """Load and normalize the per-sample context expected by every stage.

    Performs the canonical first 6 lines of every assembler/preprocessor
    module: absolute-path normalization for *input_tsv* / *output_dir*,
    ``pandas.read_csv`` of the driver, ``get_current_row_data`` lookup,
    and ``current_row.iloc[0]`` extraction.

    Parameters
    ----------
    sample_id : str
        Sample identifier matching the ``SAMPLE_ID`` column of *input_tsv*.
    input_tsv : str
        Path to the master sample-driver TSV (will be normalized to
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
    input_tsv_abs = to_abs(input_tsv)
    output_dir_abs = to_abs(output_dir)
    input_df = read_sample_table(input_tsv_abs)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]
    # preprocess_refseq copies a local REF_SEQ to <species>/RefSeq/..._RefSeq.fasta
    # (naming must match preprocess_refseq.place_assembly()). Hand every stage that
    # copy when the original extension is one tools reject (SPAdes exits 255 on
    # .fna). ``current_row`` keeps the raw cell, which preprocess_refseq reads.
    # ponytail: .gz is left alone because place_assembly copies bytes verbatim;
    # gunzip there if compressed local references need normalizing too.
    ref_seq = current_series["REF_SEQ"]
    if isinstance(ref_seq, str) and not ref_seq.lower().endswith((".fasta", ".fa", ".gz")):
        species_id, gca = current_series["SPECIES_ID"], current_series["REF_SEQ_GCA"]
        stem = f"{species_id}_{gca}_RefSeq" if pd.notna(gca) else f"{species_id}_RefSeq"
        canonical = os.path.join(output_dir_abs, str(species_id), "RefSeq", f"{stem}.fasta")
        if os.path.exists(canonical):
            current_series = current_series.copy()
            current_series["REF_SEQ"] = canonical
    return SampleContext(
        sample_id=sample_id,
        input_tsv=input_tsv_abs,
        output_dir=output_dir_abs,
        cpu_threads=cpu_threads,
        ram_gb=ram_gb,
        current_row=current_row,
        current_series=current_series,
        current_index=current_index,
        sample_stats_dict=sample_stats_dict,
    )


# --------------------------------------------------------------
# Shared __main__ for the stage scripts
# --------------------------------------------------------------
def run_stage(stage, cast=str, log=True):
    """Check ``sys.argv``, start the sample log, and call *stage*.

    Every stage script takes the same five arguments from the orchestrator:
    ``<sample_id> <input_tsv> <output_dir> <cpu_threads> <ram_gb>``.

    Parameters
    ----------
    stage : callable
        The stage function, called as
        ``stage(sample_id, input_tsv, output_dir, cpu_threads, ram_gb)``.
    cast : callable
        Applied to *cpu_threads* and *ram_gb* (``str`` or ``int``).
    log : bool
        ``False`` for the scripts that do not write to the sample log.

    Returns
    -------
    object
        Whatever *stage* returns, so the caller can turn it into an exit code.
    """
    if len(sys.argv) != 6:
        print(f"Usage: python3 {os.path.basename(sys.argv[0])} <sample_id> <input_tsv> "
              "<output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)
    sample_id, input_tsv, output_dir, cpu_threads, ram_gb = sys.argv[1:]
    if log:
        initialize_logging_environment(output_dir, sample_id)
    return stage(sample_id, input_tsv, output_dir, cast(cpu_threads), cast(ram_gb))


# --------------------------------------------------------------
# Parse NanoPlot statistics for long reads
# --------------------------------------------------------------
# NanoPlot folders written by preprocess_ont / preprocess_pacbio, per read set.
NANOSTATS_ORIGINS = {"ONT": ("Raw_ONT_", "Filt_ONT_", "Corr_ONT_"),
                     "PacBio": ("Raw_PacBio_", "Filt_PacBio_")}
# FastQC folders written by preprocess_illumina, keyed by their stats prefix.
FASTQC_FOLDERS = {"RAW_ILLU": "fastqc_results",
                  "TRIMMED_ILLU": "trimmed_fastqc_results",
                  "DEDUPED_ILLU": "dedup_fastqc_results"}
_BASE_UNITS = {"bp": 1, "kbp": 1e3, "mbp": 1e6, "gbp": 1e9}


def fastqc_total_bases(zip_path):
    """Return the "Total Bases" of one FastQC report zip as a number of bases, or ``None``.

    ``None`` covers an unreadable zip and FastQC releases before 0.12, which do
    not report the figure.
    """
    # ponytail: FastQC prints this to four significant digits ("396.4 Mbp"), so
    # coverage is good to about 0.05%. For exact totals, sum length x count from
    # the ">>Sequence Length Distribution" module (exact only when it is unbinned).
    try:
        with zipfile.ZipFile(zip_path) as archive:
            name = next(n for n in archive.namelist() if n.endswith("/fastqc_data.txt"))
            for line in archive.read(name).decode("utf-8", "replace").splitlines():
                if line.startswith("Total Bases\t"):
                    value, unit = line.split("\t")[1].split()
                    return float(value) * _BASE_UNITS[unit.lower()]
    except (OSError, zipfile.BadZipFile, StopIteration, ValueError, KeyError):
        pass
    return None


def load_read_stats(species_dir, sample_stats_dict):
    """Fill read counts and total bases from the NanoPlot and FastQC outputs on disk.

    Each stage is its own process with a fresh stats record, so nothing that
    preprocessing measured reaches the final QC unless it is read back here.
    Read sets with no output on disk keep ``None``.

    Parameters
    ----------
    species_dir : str
        ``<output_dir>/<SPECIES_ID>``, holding ``ONT/``, ``PacBio/`` and ``Illumina/``.
    sample_stats_dict : dict
        Statistics dictionary to update in-place.

    Returns
    -------
    dict
        The updated *sample_stats_dict*.
    """
    for technology, origins in NANOSTATS_ORIGINS.items():
        for origin in origins:
            nanostats = os.path.join(species_dir, technology, f"{origin}nanoplot_analysis",
                                     f"{origin}NanoStats.txt")
            if os.path.exists(nanostats):
                analyze_nanostats(origin, nanostats, sample_stats_dict)
    for prefix, folder in FASTQC_FOLDERS.items():
        totals = [fastqc_total_bases(report) for report in
                  sorted(glob.glob(os.path.join(species_dir, "Illumina", folder, "*_fastqc.zip")))]
        if len(totals) == 2 and None not in totals:  # forward + reverse
            sample_stats_dict[f"{prefix}_TOTAL_BASES"] = sum(totals)
    return sample_stats_dict


def add_read_coverage(sample_stats_dict, genome_bp):
    """Set ``<SET>_COVERAGE`` for every read set that has a ``<SET>_TOTAL_BASES``.

    Returns the number of read sets given a coverage (0 when *genome_bp* is
    unknown or no read set has a total, as for QC-only samples).
    """
    if genome_bp is None or pd.isna(genome_bp) or genome_bp <= 0:
        return 0
    done = 0
    for key, total in list(sample_stats_dict.items()):
        if key.endswith("_TOTAL_BASES") and total is not None and not pd.isna(total):
            sample_stats_dict[key.replace("_TOTAL_BASES", "_COVERAGE")] = round(total / genome_bp, 2)
            done += 1
    return done


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
    if not os.path.exists(nanoplot_out_file):
        # NanoPlot is run as a soft-guarded best-effort step by the
        # preprocessors, so a missing NanoStats file is not fatal here.
        # Skip rather than crash; downstream selection falls back to defaults.
        print(f"WARN:\tNanoStats file not found; skipping {READS_ORIGIN} stats: {nanoplot_out_file}")
        return sample_stats_dict
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
def select_long_reads(output_dir, input_tsv, sample_id, cpu_threads):
    """Select the highest-mean-quality long reads from ONT or PacBio data.

    Parses NanoPlot statistics to compare quality across raw, filtered, and
    corrected read sets, then copies the best set to a canonical
    ``*_highest_mean_qual_long_reads.fastq`` path.

    Parameters
    ----------
    output_dir : str
        Root output directory containing per-species subdirectories.
    input_tsv : str
        Path to the metadata TSV file.
    sample_id : str
        Sample identifier used to look up the row in *input_tsv*.
    cpu_threads : int
        Number of threads available for compression tasks.

    Returns
    -------
    str or None
        Absolute path to the selected highest-quality reads file, or
        ``None`` if no suitable file can be found.
    """
    input_df = read_sample_table(input_tsv)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]  # Convert to Series (single row)

    # Identify read paths, reference, and BUSCO lineage info from TSV
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
        reads_token = reads_type.lower()
        # NanoPlot analyses and filtered/corrected reads are written by
        # preprocess_ont into output_dir/species_id/ONT (ont_dir_abs), NOT
        # next to the raw input file.  Anchor here so an external
        # ONT_RAW_READS path does not break the lookup (see issue #17).
        reads_dir = os.path.join(output_dir, species_id, reads_type)
        filtered_reads = os.path.join(reads_dir, f"{species_id}_{reads_token}_filtered.fastq")
        corrected_reads = os.path.join(reads_dir, f"{species_id}_{reads_token}_corrected.fastq")
        reads_origin_list = ["Raw_ONT_", "Filt_ONT_", "Corr_ONT_"]
    elif pd.notna(pacbio_raw_reads):
        print("DEBUG - PROCESSING PACBIO HIGHEST MEAN QUAL")
        reads_type = "PacBio"
        reads_token = reads_type.lower()
        # See ONT branch: anchor to output_dir/species_id/PacBio
        # (pacbio_dirto_abs), where preprocess_pacbio writes its artifacts.
        reads_dir = os.path.join(output_dir, species_id, reads_type)
        filtered_reads = os.path.join(reads_dir, f"{species_id}_{reads_token}_filtered.fastq")
        corrected_reads = os.path.join(reads_dir, f"{species_id}_{reads_token}_corrected.fastq")
        reads_origin_list = ["Raw_PacBio_", "Filt_PacBio_"]
    else:
        print(f"ERROR:\tUNABLE TO PARSE LONG READS AS BOTH ONT AND PACBIO RAW READS ARE NONE: {ont_raw_reads} & {pacbio_raw_reads}")

    for reads_origin in reads_origin_list:
        sample_stats_dict = analyze_nanostats(reads_origin, os.path.join(reads_dir, f"{reads_origin}nanoplot_analysis", f"{reads_origin}NanoStats.txt"), sample_stats_dict)

    if pd.notna(ont_raw_reads):
        print("Selecting Highest Mean Quality Long reads...")
        highest_mean_qual_long_reads = corrected_reads
        if pd.notna(ont_raw_reads) and sample_stats_dict.get("CORRECT_ONT_MEAN_QUAL", 0.0) < sample_stats_dict.get("FILT_ONT_MEAN_QUAL", 0.0):
            highest_mean_qual_long_reads = filtered_reads
            highest_mean_qual = sample_stats_dict.get("FILT_ONT_MEAN_QUAL", 0.0)
        else:
            highest_mean_qual = sample_stats_dict.get("CORRECT_ONT_MEAN_QUAL", 0.0)
    if pd.notna(pacbio_raw_reads):
        print("Selecting Highest Mean Quality Long reads...")
        highest_mean_qual_long_reads = pacbio_raw_reads
        if pd.notna(pacbio_raw_reads) and sample_stats_dict.get("RAW_PACBIO_MEAN_QUAL", 0.0) < sample_stats_dict.get("FILT_PACBIO_MEAN_QUAL", 0.0):
            highest_mean_qual_long_reads = filtered_reads
            highest_mean_qual = sample_stats_dict.get("FILT_PACBIO_MEAN_QUAL", 0.0)
        else:
            highest_mean_qual = sample_stats_dict.get("RAW_PACBIO_MEAN_QUAL", 0.0)
    print(f"Highest Mean Quality Long reads: {highest_mean_qual_long_reads}")
    print(f"Mean Quality: {highest_mean_qual}")

    # What the chosen set is, for stages that must match a tool mode to it.
    # *_corrected.fastq is only a symlink to the filtered reads when Ratatosk
    # was skipped, so a symlink does not count as corrected.
    if highest_mean_qual_long_reads == corrected_reads and not os.path.islink(corrected_reads):
        reads_source = "corrected"
    elif highest_mean_qual_long_reads == pacbio_raw_reads:
        reads_source = "raw"
    else:
        reads_source = "filtered"

    renamed_highest_mean_qual_long_reads = f"{species_id}_{reads_type}_highest_mean_qual_long_reads.fastq"
    if not os.path.exists(highest_mean_qual_long_reads):
        reads_source = "filtered"
        # The chosen set may never have been written: when Ratatosk yields
        # nothing, preprocess_ont carries on with the filtered reads and no
        # *_corrected.fastq exists. Use the filtered set rather than returning
        # None, which makes every later stage fall back to the raw reads.
        # The second name is the older "<species>_filtered.fastq" layout.
        for fallback_file in (filtered_reads, os.path.join(reads_dir, f"{species_id}_filtered.fastq")):
            if os.path.exists(fallback_file) and os.path.getsize(fallback_file) > 0:
                print(f"FALLBACK:\t{highest_mean_qual_long_reads} does not exist; using filtered reads: {fallback_file}")
                highest_mean_qual_long_reads = fallback_file
                break
        else:
            print("ERROR:\tNo usable highest-mean-quality long read file found.")
            return None

    renamed_highest_mean_qual_long_reads = os.path.join(reads_dir, f"{species_id}_{reads_type}_highest_mean_qual_long_reads.fastq")
    shutil.copy(highest_mean_qual_long_reads, renamed_highest_mean_qual_long_reads)
    with open(renamed_highest_mean_qual_long_reads + LONG_READ_INFO_SUFFIX, "w") as fh:
        json.dump({"source": reads_source, "mean_quality": highest_mean_qual}, fh)

    print(f"NOTE:\tSelected highest quality long reads: {renamed_highest_mean_qual_long_reads} "
          f"({reads_source}) with mean quality {highest_mean_qual}")

    return renamed_highest_mean_qual_long_reads


# --------------------------------------------------------------
# Match Flye's ONT read-type mode to the selected reads
# --------------------------------------------------------------
LONG_READ_INFO_SUFFIX = ".info.json"

# Flye documents --nano-hq for reads under 5% error, which is about Q13.
NANO_HQ_MIN_MEAN_Q = 13.0


def long_read_info(selected_reads):
    """Return what select_long_reads recorded about *selected_reads*, or ``{}``.

    Keys: ``source`` (``corrected``, ``filtered`` or ``raw``) and
    ``mean_quality`` (NanoPlot mean read quality, may be ``None``).
    """
    try:
        with open(str(selected_reads) + LONG_READ_INFO_SUFFIX) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def flye_ont_mode(reads_in_use, selected_reads):
    """Choose ``--nano-corr``, ``--nano-hq`` or ``--nano-raw`` for Flye.

    ``--nano-corr`` assumes error-corrected reads (under 3% error). Given
    uncorrected reads of mean quality 11 it broke an E. coli chromosome into
    10 contigs where the other two modes gave 3; at mean quality 21 all three
    modes agreed. So it is used only for Ratatosk-corrected reads.

    Parameters
    ----------
    reads_in_use : str
        The file Flye will be given.
    selected_reads : str
        The canonical ``*_highest_mean_qual_long_reads.fastq`` path.
    """
    if os.path.abspath(str(reads_in_use)) != os.path.abspath(str(selected_reads)):
        return "--nano-raw"      # fell back to the raw reads; quality unknown
    info = long_read_info(selected_reads)
    if not info:
        return "--nano-corr"     # selected before this record existed; keep the old behaviour
    if info.get("source") == "corrected":
        return "--nano-corr"
    mean_q = info.get("mean_quality") or 0.0
    return "--nano-hq" if mean_q >= NANO_HQ_MIN_MEAN_Q else "--nano-raw"
