#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
process_metadata.py

Writes the NCBI submission sheets for a run: one row per read library for SRA
and one row per assembly for GenBank. Runs as a final step for every sample,
after the QC and the report, and appends to two files in the output folder, so
a run of many samples ends with one sheet of each kind.

Inputs come from the sample table: the read and reference columns, and the
"Sample Metadata" columns (BioSample, tissue, library strategy, instruments).
A sample with no metadata still gets its assembly row; it only misses the SRA rows.

Also holds the two iNaturalist helpers the HTML report uses.

Stage:
    Metadata (NCBI SRA / GenBank / iNaturalist)

Created on Wed Aug 16 2023

Updated on 2026-04-16

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""

import os
import sys
import warnings
from datetime import datetime
from pathlib import Path

import pandas as pd

from file_operations import resolve_templates_dir
from sample_tsv import load_sample_context, long_read_info, run_stage

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # EGAP.py sits beside bin/
from EGAP import VERSION  # noqa: E402


# ignore only the Data Validation extension warning from openpyxl
warnings.filterwarnings("ignore",
                        message="Data Validation extension is not supported and will be removed",
                        category=UserWarning,
                        module="openpyxl")

# --------------------------------------------------------------
# Reverse geocode latitude and longitude
# --------------------------------------------------------------
def reverse_geocode(lat, lon, user_agent="reverse_geocoder"):
    """Reverse geocode latitude and longitude to obtain location details.

    Uses Nominatim to convert coordinates into country, state, and county information,
    respecting usage policies with rate limiting.

    Args:
        lat (float): Latitude coordinate.
        lon (float): Longitude coordinate.
        user_agent (str): User agent identifier for Nominatim (default: 'reverse_geocoder').

    Returns:
        tuple: (country, state, county) as strings, or (None, None, None) if geocoding fails.
    """
    # Imported here so the metadata sheets below do not need geopy.
    from geopy.extra.rate_limiter import RateLimiter
    from geopy.geocoders import Nominatim

    geolocator = Nominatim(user_agent=user_agent)

    # RateLimiter to respect usage policy (1 request per second)
    geocode = RateLimiter(geolocator.reverse, min_delay_seconds=1)
    location = geocode((lat, lon), addressdetails=True)
    if location is None:
        return None, None, None
    address = location.raw.get("address", {})
    country = address.get("country")
    state = address.get("state") or address.get("region")
    county = address.get("county") or address.get("city_district")

    print(f"Country:\t{country}")
    print(f"State:\t{state}")
    print(f"County:\t{county}")

    return country, state, county


# --------------------------------------------------------------
# Retrieve iNaturalist observation data
# --------------------------------------------------------------
def get_inat_obs(search_id):
    """Retrieve observation data from iNaturalist by ID.

    Fetches observation details including collector, date, coordinates, and photo URLs
    using the pyinaturalist API.

    Args:
        search_id (str or int): iNaturalist observation ID.

    Returns:
        tuple: (collector, observation_date, latitude, longitude, photo_url)
               or None for photo URLs if no photos are available.
    """
    import requests  # here, so the metadata sheets below do not need it

    url = f"https://api.inaturalist.org/v1/observations/{search_id}"
    r = requests.get(url)
    r.raise_for_status()
    data = r.json()["results"][0]
    collector = data["user"]["name"]
    obs_date  = data["observed_on"]
    latitude  = data["geojson"]["coordinates"][1]
    longitude = data["geojson"]["coordinates"][0]
    photo = data["photos"][0]
    # Try the explicit original_url first
    photo_url = photo.get("original_url")
    # If it’s a size-template, swap in “original”
    if not photo_url and "{size}" in photo["url"]:
        photo_url = photo["url"].replace("{size}", "original")
    # Fallback to the raw url
    if not photo_url:
        photo_url = photo["url"]

    return collector, obs_date, latitude, longitude, photo_url


# --------------------------------------------------------------
# NCBI submission sheets
# --------------------------------------------------------------
# One entry per read technology: the sample-table columns that say a sample has
# it (the SRA column first), the processed reads EGAP leaves in the species
# folder, and how its SRA library row is described.
READ_SETS = (
    {"technology": "ONT", "platform": "OXFORD_NANOPORE", "layout": "SINGLE",
     "inputs": ("ONT_SRA", "ONT_RAW_DIR", "ONT_RAW_READS"),
     "files": ("ONT/{species_id}_ONT_highest_mean_qual_long_reads.fastq",),
     "design": "Filtlong length filtering, Ratatosk Illumina-Reads corrected"},
    {"technology": "Illumina", "platform": "ILLUMINA", "layout": "PAIRED",
     "inputs": ("ILLUMINA_SRA", "ILLUMINA_RAW_DIR", "ILLUMINA_RAW_F_READS"),
     "files": ("Illumina/{species_id}_illu_forward_dedup.fastq",
               "Illumina/{species_id}_illu_reverse_dedup.fastq"),
     "design": ("Trimmomatic adapter removal and quality trimming, BBDuk Decontamination "
                "using Kmers, Clumpify reads deduplication")},
    {"technology": "PacBio", "platform": "PACBIO_SMRT", "layout": "SINGLE",
     "inputs": ("PACBIO_SRA", "PACBIO_RAW_DIR", "PACBIO_RAW_READS"),
     "files": ("PacBio/{species_id}_PacBio_highest_mean_qual_long_reads.fastq",),
     "design": "Filtlong quality and length filtering"},
)
# select_long_reads records which stage its pick came from; the final stats use these prefixes.
_LONG_READ_STATS = {"corrected": "CORRECT", "filtered": "FILT", "raw": "RAW"}


def cell(row, column):
    """Return the stripped text of *row[column]*, or ``None`` for a blank or missing cell."""
    value = row.get(column)
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    return text if text and text.lower() not in ("none", "nan") else None


def template_columns(file_name, sheet_name):
    """Return the header of one NCBI template sheet, so newer templates can be dropped in."""
    return list(pd.read_excel(resolve_templates_dir() / file_name, sheet_name=sheet_name).columns)


def upsert_rows(path, columns, key_column, key_values, new_rows):
    """Replace the rows keyed by *key_values* in the TSV at *path* with *new_rows*.

    Creates the file when absent. A rerun therefore replaces a sample's rows
    rather than appending duplicates, and other samples' rows are kept.
    """
    # ponytail: read-modify-write of a file every sample shares. Safe because
    # samples run one at a time (EGAP.py's loop; maxForks = 1 under Nextflow).
    # Give each sample its own file and merge at the end if that ever changes.
    kept = []
    if os.path.exists(path):
        existing = pd.read_csv(path, sep="\t", dtype=str)
        kept = existing[~existing[key_column].isin(key_values)].to_dict("records")
    pd.DataFrame(kept + new_rows, columns=columns).to_csv(path, sep="\t", index=False)
    return path


def sample_read_sets(row):
    """Return the READ_SETS entries this sample has input for."""
    return [read_set for read_set in READ_SETS if any(cell(row, c) for c in read_set["inputs"])]


def genome_coverage(stats, read_sets, species_dir, species_id):
    """Return coverage of the reads the assembly was built from, as ``"87.5x"``, or ``None``.

    Uses the totals the final QC stored (deduplicated Illumina reads, and for
    long reads whichever of raw, filtered or corrected was selected), divided by
    the assembly length. Not the report's per-read-set coverage, whose
    denominator is the reference when there is one.
    """
    total = 0.0
    for read_set in read_sets:
        if read_set["technology"] == "Illumina":
            prefix = "DEDUPED_ILLU"
        else:
            selected = os.path.join(species_dir, read_set["files"][0].format(species_id=species_id))
            stage = _LONG_READ_STATS.get(long_read_info(selected).get("source"), "RAW")
            prefix = f"{stage}_{read_set['technology'].upper()}"
        bases = stats.get(f"{prefix}_TOTAL_BASES")
        if bases is not None and pd.notna(bases):
            total += float(bases)
    size = stats.get("GENOME_SIZE")
    if not total or size is None or pd.isna(size) or float(size) <= 0:
        return None
    return f"{total / float(size):.1f}x"


def sra_rows(row, sample_id, species_id, read_sets):
    """Return one SRA library row per read set that is not already in SRA."""
    tissue, strategy = cell(row, "SAMPLE_TISSUE"), cell(row, "SEQUENCING_METHODOLOGY")
    rows = []
    for read_set in read_sets:
        if cell(row, read_set["inputs"][0]):
            continue  # downloaded from SRA: nothing to submit
        files = [os.path.basename(f.format(species_id=species_id)) for f in read_set["files"]]
        rows.append({
            "biosample_accession": cell(row, "BIOSAMPLE_ACCESSION"),
            "library_ID": f"{sample_id}_{read_set['technology']}",
            "title": f"{strategy} of {species_id}: {tissue}",
            "library_strategy": strategy,
            "library_source": cell(row, "LIBRARY_SOURCE"),
            "library_selection": cell(row, "LIBRARY_SELECTION"),
            "library_layout": read_set["layout"],
            "platform": read_set["platform"],
            "instrument_model": cell(row, f"INSTRUMENT_{read_set['technology'].upper()}"),
            "design_description": read_set["design"],
            "filetype": "fastq",
            "filename": files[0],
            "filename2": files[1] if len(files) > 1 else None,
        })
    return rows


def process_metadata(sample_id, input_tsv, output_dir, cpu_threads, ram_gb):
    """Add one sample's rows to the run's NCBI SRA and assembly metadata sheets.

    Writes ``<output_dir>/<table stem>_EGAP_SRA_metadata.tsv`` and
    ``<output_dir>/<table stem>_EGAP_Assembly_metadata.tsv``, laid out like the
    NCBI templates in ``resources/templates``. Skips, without failing:

    * the whole step, when the sample has no final assembly statistics (a
      QC-only sample, or an assembly that did not finish);
    * the SRA rows, when SAMPLE_TISSUE or SEQUENCING_METHODOLOGY is blank, or
      every read set came from SRA.

    Args:
        sample_id (str): The SAMPLE_ID of the row to process.
        input_tsv (str): Path to the sample table the stages read.
        output_dir (str): Root output directory of the run.
        cpu_threads (int | str): Unused; every stage takes the same five arguments.
        ram_gb (int | str): Unused, as above.

    Returns:
        tuple[str | None, str | None]: ``(sra_path, assembly_path)``; ``None``
        for a sheet that was skipped.
    """
    ctx = load_sample_context(sample_id, input_tsv, output_dir, cpu_threads, ram_gb)
    row = ctx.current_series
    species_id = row["SPECIES_ID"]
    species_dir = os.path.join(ctx.output_dir, species_id)
    assembly_stem = f"{sample_id}_final_EGAP_assembly"
    stats_csv = os.path.join(species_dir, sample_id, f"{assembly_stem}_final_stats.csv")
    if not os.path.exists(stats_csv):
        print(f"SKIP:\tNo NCBI metadata for {sample_id}: no final assembly statistics at "
              f"{stats_csv} (a QC-only sample, or the assembly did not finish).")
        return None, None
    stats = pd.read_csv(stats_csv).iloc[0]
    read_sets = sample_read_sets(row)
    table_stem = os.path.splitext(os.path.basename(str(input_tsv)))[0]

    # ---- SRA: one row per read library EGAP processed from local files ----
    sra_path = None
    missing = [c for c in ("SAMPLE_TISSUE", "SEQUENCING_METHODOLOGY") if not cell(row, c)]
    libraries = [] if missing else sra_rows(row, sample_id, species_id, read_sets)
    if missing:
        print(f"SKIP:\tNo SRA rows for {sample_id}: {' and '.join(missing)} blank in the sample metadata.")
    elif not libraries:
        print(f"SKIP:\tNo SRA rows for {sample_id}: its reads were downloaded from SRA.")
    else:
        sra_path = upsert_rows(
            os.path.join(ctx.output_dir, f"{table_stem}_EGAP_SRA_metadata.tsv"),
            template_columns("SRA_metadata_acc.xlsx", "SRA_data"), "library_ID",
            [f"{sample_id}_{read_set['technology']}" for read_set in READ_SETS], libraries)
        print(f"PASS:\tWrote {len(libraries)} SRA row(s) for {sample_id}: {sra_path}")

    # ---- GenBank: one row for the assembly ----
    reference = cell(row, "REF_SEQ_GCA") or os.path.basename(cell(row, "REF_SEQ") or "") or None
    assembly_path = upsert_rows(
        os.path.join(ctx.output_dir, f"{table_stem}_EGAP_Assembly_metadata.tsv"),
        template_columns("GCA_metadata_acc.xlsx", "Genome_data"), "sample_name", [sample_id],
        [{"biosample_accession": cell(row, "BIOSAMPLE_ACCESSION"),
          "sample_name": sample_id,
          "assembly_date": datetime.now().strftime("%Y-%m-%d"),
          "assembly_name": assembly_stem,
          "assembly_method": "EGAP",
          "assembly_method_version": VERSION,
          "genome_coverage": genome_coverage(stats, read_sets, species_dir, species_id),
          "sequencing_technology": ", ".join(r["technology"] for r in read_sets) or None,
          "reference_genome": reference,
          "filename": f"{assembly_stem}.fasta"}])
    print(f"PASS:\tWrote the assembly row for {sample_id}: {assembly_path}")
    return sra_path, assembly_path


if __name__ == "__main__":
    run_stage(process_metadata, log=False)
