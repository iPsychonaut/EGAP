#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
process_metadata.py

Generates metadata for genome assemblies and Sequence Read Archive (SRA) submissions.
Processes sample data from a CSV to create TSV files with assembly metadata and SRA metadata,
and retrieves iNaturalist observation data with reverse geocoding for location details.

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""

import os, sys, warnings, requests
import pandas as pd
from datetime import datetime
from pathlib import Path
from geopy.geocoders import Nominatim
from geopy.extra.rate_limiter import RateLimiter
from utilities import get_current_row_data, calculate_genome_coverage

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
    geolocator = Nominatim(user_agent=user_agent)

    # RateLimiter to respect usage policy (1 request per second)
    geocode = RateLimiter(geolocator.reverse, min_delay_seconds=1)
    location = geocode((lat, lon), addressdetails=True)
    if location is None:
        return {"country": None, "state": None, "county": None}
    address = location.raw.get("address", {})
    country = address.get("country")
    state = address.get("state") or address.get("region"),
    county = address.get("county") or address.get("city_district"),

    print(f"Country:\t{country}")
    print(f"State:\t{state}")
    print(f"County:\t{county}")

    return country, state[0], county[0]


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
# Generate assembly metadata TSV
# --------------------------------------------------------------
def gen_assembly_metadata_tsv(input_csv, sample_id, output_dir, templates_dir):
    """Generate a TSV file with metadata for a genome assembly.

    Creates a TSV file containing assembly metadata (e.g., date, method, coverage)
    based on sample data from a CSV and the final assembly FASTA.

    Args:
        input_csv (str): Path to the input CSV with sample metadata.
        sample_id (str): Sample identifier.
        output_dir (str): Directory for output files.

    Returns:
        str: Path to the generated assembly metadata TSV file.
    """
    print(f"Generating Assembly Metadata TSV entry for {sample_id}...")

    # Load the input_csv into dataframe and collect current sample_id information
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
    species_id = current_series["SPECIES_ID"]

    # Versioning
    assembly_method = "EGAP"
    assembly_method_version = "3.1"
    
    # Determine where the final assembly FASTA lives
    species_id    = current_series["SPECIES_ID"]
    species_dir = os.path.join(output_dir, species_id)
    sample_dir = os.path.join(species_dir, sample_id)
    assembly_name = f"{sample_id}_final_EGAP_assembly.fasta"
    assembly_path = os.path.join(sample_dir, assembly_name)

    all_reads = []
    if pd.notna(ont_sra) or pd.notna(ont_raw_reads):
        all_reads.append(os.path.join(species_dir, "ONT", f"{species_id}_ONT_highest_mean_qual_long_reads.fastq"))
    if pd.notna(illumina_sra) or (pd.notna(illumina_f_raw_reads) and pd.notna(illumina_r_raw_reads)):
        all_reads.append(os.path.join(species_dir, "Illumina", f"{species_id}_illu_forward_dedup.fastq"))
        all_reads.append(os.path.join(species_dir, "Illumina", f"{species_id}_illu_reverse_dedup.fastq"))
    if pd.notna(pacbio_sra) or pd.notna(pacbio_raw_reads):
        all_reads.append(os.path.join(species_dir, "PacBio", f"{species_id}_PacBio_highest_mean_qual_long_reads.fastq"))    
    genome_coverage = str(calculate_genome_coverage(all_reads, assembly_path)) + "x" 

    print(f"DEBUG - genome_coverage - {genome_coverage}")
    
    # Comma‐separated tech list
    techs = []
    if pd.notna(illumina_sra): techs.append("Illumina")
    if pd.notna(ont_sra):      techs.append("ONT")
    if pd.notna(pacbio_sra):   techs.append("PacBio")
    sequencing_technology = ", ".join(techs) or None
    
    # Pick a reference if present
    reference_genome = None
    if pd.notna(current_series.get("REF_SEQ_GCA")):
        reference_genome = current_series["REF_SEQ_GCA"]
    elif pd.notna(current_series.get("REF_SEQ")):
        reference_genome = current_series["REF_SEQ"]
    
    # Final assembly date
    assembly_date = datetime.now().strftime("%Y-%m-%d")
    
    # Locate template and pull in the DataFrame
    gca_file  = os.path.join(templates_dir, "GCA_metadata_acc.xlsx")
    genome_df = pd.read_excel(gca_file, sheet_name="Genome_data")
    
    # Ensure one row exists for this sample
    mask = genome_df["sample_name"] == sample_id
    if not mask.any():
        blank = {col: None for col in genome_df.columns}
        blank["biosample_accession"] = species_id
        blank["sample_name"] = sample_id
        genome_df.loc[len(genome_df)] = blank
        mask = genome_df["sample_name"] == sample_id
    
    # Write eight metadata fields into that row
    genome_df.loc[mask, ["assembly_date",
                         "assembly_name",
                         "assembly_method",
                         "assembly_method_version",
                         "genome_coverage",
                         "sequencing_technology",
                         "reference_genome",
                         "filename"
    ]] = [assembly_date,
          assembly_name,
          assembly_method,
          assembly_method_version,
          genome_coverage,
          sequencing_technology,
          reference_genome,
          assembly_name]

    output_filename = f"{os.path.basename(input_csv).replace('.csv','')}_EGAP_Assembly_metadata.tsv"
    assembly_metadata_path = os.path.join(output_dir, output_filename)
    genome_df.to_csv(assembly_metadata_path, sep="\t", index=False)
    
    print(f"PASS:\tSuccessfully generated the {sample_id} Assembly Metadata TSV entry into: {assembly_metadata_path}")

    return assembly_metadata_path


# --------------------------------------------------------------
# Generate SRA metadata TSV
# --------------------------------------------------------------
def gen_sra_metadata_xlsx(input_csv, sample_id, output_dir, templates_dir):
    """Generate a TSV file with metadata for SRA submission.

    Creates a TSV file with SRA metadata (e.g., library details, read files) based on
    sample data from a CSV, supporting Illumina, ONT, and PacBio reads.

    Args:
        input_csv (str): Path to the input CSV with sample metadata.
        sample_id (str): Sample identifier.
        output_dir (str): Directory for output files.

    Returns:
        str or None: Path to the generated SRA metadata TSV file, or None if required fields are missing.
    """
    print("Generating SRA Metadata TSV entry...")
    
    # Load the input_csv into dataframe and collect current sample_id information
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
    species_id = current_series["SPECIES_ID"]
    search_id = current_series["INATRUALIST_ID"]
    sample_tissue = current_series["SAMPLE_TISSUE"]
    library_strategy = current_series["SEQUENCING_METHODOLOGY"]
    library_source = current_series["LIBRARY_SOURCE"]
    species_dir = os.path.join(output_dir, species_id)
        
    if pd.isna(sample_tissue):
        print("SKIP:\tNo SRA Metadata processing performed, missing SAMPLE_TISSUE: {sample_tissue}.")
    elif pd.isna(library_strategy):
        print("SKIP:\tNo SRA Metadata processing performed, missing SEQUENCING_METHODOLOGY: {library_strategy}.")
    else:        
        # Generate SRA Metadata Dataframe from either template or existing files
        existing_sra_metadata = os.path.join(os.path.dirname(input_csv), f"{os.path.basename(input_csv).replace('.csv','')}_EGAP_SRA_metadata.xlsx")
        if os.path.exists(existing_sra_metadata):
            base_file_path = existing_sra_metadata
        else:
            base_file_path = os.path.join(templates_dir, "SRA_metadata_acc.xlsx")
        sample_data = {}
        print(f"DEBUG - base_file_path - {base_file_path}")
        main_sra_df = pd.read_excel(base_file_path, sheet_name="SRA_data")
        sample_data[search_id] = main_sra_df.columns.tolist()

        # Initialize new SRA metadata DataFrame
        sra_columns = main_sra_df.columns.tolist()
        new_sra_df = pd.DataFrame(columns=sra_columns)

        reads = []
        library_ids = []
        filetype = "fastq"
        title = f"{library_strategy} of {species_id}: {sample_tissue}"
        if pd.notna(ont_sra) or pd.notna(ont_raw_reads):
            reads.append(os.path.join(species_dir, "ONT", f"{species_id}_ONT_highest_mean_qual_long_reads.fastq"))
            design_description = "Filtlong length filtering, Ratatosk Illumina-Reads corrected"
            library_ids.append(f"{species_id}_EGAP_Corrected_ONT")
        if pd.notna(illumina_sra) or (pd.notna(illumina_f_raw_reads) and pd.notna(illumina_r_raw_reads)):
            reads.append([os.path.join(species_dir, "Illumina", f"{species_id}_illu_forward_dedup.fastq"), os.path.join(species_dir, "Illumina", f"{species_id}_illu_reverse_dedup.fastq")])
            design_description = "Trimmomatic adapter removal and quality trimming, BBDuk Decontamination using Kmers, Clumpify reads deduplication"
            library_ids.append(f"{species_id}_EGAP_Forward_Deduplicated_Illumina")
            library_ids.append(f"{species_id}_EGAP_Reverse_Deduplicated_Illumina")
        if pd.notna(pacbio_sra) or pd.notna(pacbio_raw_reads):
            reads.append(os.path.join(species_dir, "PacBio", f"{species_id}_PacBio_highest_mean_qual_long_reads.fastq"))
            design_description = "Filtlong quality and length filtering"
            library_ids.append(f"{species_id}_EGAP_Filtered_PacBio")

        print(f"DEBUG - reads - {reads}")    
        print(f"DEBUG - library_ids - {library_ids}")
        print(f"DEBUG - title - {title}")
        print(f"DEBUG - design_description - {design_description}")

        for index, filename in enumerate(reads):
            library_id = library_ids[index]
        
            # split out one or two filenames
            if isinstance(filename, list):
                fname1, fname2 = filename
            else:
                fname1, fname2 = filename, None
        
            row = {"biosample_accession":    species_id,
                   "library_ID":             library_id,
                   "title":                  title.format(library_strategy=library_strategy,
                                                          species_id=species_id,
                                                          sample_tissue=sample_tissue),
                   "library_strategy":       library_strategy,
                   "library_source":         library_source,
                   "library_selection":      None,
                   "library_layout":         "PAIRED" if isinstance(filename, list) else "SINGLE",
                   "platform":               None,
                   "instrument_model":       None,
                   "design_description":     design_description,
                   "filetype":               filetype,
                   "filename":               fname1,
                   "filename2":              fname2,
                   "filename3":              None,
                   "filename4":              None,
                   "assembly":               None,
                   "fasta_file":             None}
                    
            # Append the new row using loc
            new_sra_df.loc[len(new_sra_df)] = row

        
        # Save the new SRA metadata workbook
        output_filename = f"{os.path.basename(input_csv).replace('.csv','')}_EGAP_SRA_metadata.tsv"
        sra_metadata_path = os.path.join(output_dir, output_filename)
        
        # Write single‐sheet TSV
        new_sra_df.to_csv(sra_metadata_path, sep="\t", index=False)

        print(f"PASS:\tSuccessfully generated the SRA Metadata TSV entry into: {sra_metadata_path}")

        return sra_metadata_path


def process_metadata(sample_id, input_csv, output_dir, cpu_threads, ram_gb):
    """
    PEP 8 Documentation
    """
    # Load input CSV and extract sample data
    input_df = pd.read_csv(input_csv)
    current_row, current_index, sample_stats_dict = get_current_row_data(input_df, sample_id)
    current_series = current_row.iloc[0]

    # Extract read paths and metadata from CSV
    species_id = current_series["SPECIES_ID"]
    search_id = current_series["INATRUALIST_ID"]
    species_dir = os.path.join(output_dir, species_id)

    print(f"DEBUG - species_id - {species_id}")
    print(f"DEBUG - species_dir - {species_dir}")
    print(f"DEBUG - search_id - {search_id}")

    # Path to this script
    script_file = Path(__file__).resolve()
    bin_dir = script_file.parent
    resources_dir = bin_dir.parent / "resources"
    templates_dir = resources_dir / "templates"
    print(f"DEBUG - templates_dir - {templates_dir}")

    # Generate SRA metadata TSV
    sra_output_path = gen_sra_metadata_xlsx(input_csv, sample_id, output_dir, templates_dir)

    print("PASS:\tSuccessfully generated NCBI SRA Metadata TSV: {sra_output_path}.")
    
    # Generate assembly metadata TSV
    assembly_metadata_path = gen_assembly_metadata_tsv(input_csv, sample_id, output_dir, templates_dir)

    print("PASS:\tSuccessfully generated NCBI Assembly Metadata TSV: {assembly_metadata_path}.")
        
    return sra_output_path, assembly_metadata_path


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python3 process_metadata.py <sample_id> <input_csv> "
            "<output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)
        
    sra_output_path, assembly_metadata_path = process_metadata(sys.argv[1],       # sample_id
                                                               sys.argv[2],       # input_csv
                                                               sys.argv[3],       # output_dir
                                                               str(sys.argv[4]),  # cpu_threads
                                                               str(sys.argv[5]))  # ram_gb