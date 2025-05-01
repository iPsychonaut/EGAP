#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
select_long_reads.py

Updated on Sat Mar 29 2025

This script selects the highest quality long reads (ONT or PacBio) based on NanoPlot mean quality.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import os
import sys
import subprocess

# --------------------------------------------------------------
# Parse NanoPlot summary for mean quality
# --------------------------------------------------------------
def get_mean_quality(nanoplot_dir, prefix):
    """Parse NanoPlot summary file to extract mean quality.

    Reads the NanoPlot summary file and extracts the mean quality value.

    Args:
        nanoplot_dir (str): Directory containing the NanoPlot summary file.
        prefix (str): Prefix for the NanoPlot summary file (e.g., 'Corrected_ONT_').

    Returns:
        float or None: Mean quality value if found, else None.
    """
    summary_file = f"{nanoplot_dir}/{prefix}NanoStats.txt"
    if not os.path.exists(summary_file):
        return None
    with open(summary_file, "r") as f:
        for line in f:
            if "Mean quality:" in line:
                return float(line.split(":")[1].strip())
    return None

# --------------------------------------------------------------
# Select highest quality long reads
# --------------------------------------------------------------
def select_long_reads(ont_file, pacbio_file, species_id):
    """Select the highest quality long reads from ONT or PacBio data.

    Compares mean quality scores from NanoPlot summaries and selects the best reads,
    creating a symbolic link to the selected file.

    Args:
        ont_file (str): Path to ONT reads file or 'None' if unavailable.
        pacbio_file (str): Path to PacBio reads file or 'None' if unavailable.
        species_id (str): Identifier for the species.

    Returns:
        None: Outputs the selected file path and creates a symbolic link.
    """
    ont_corrected_qual = get_mean_quality("corr_ont", "Corrected_ONT_") if ont_file != "None" else None
    ont_filtered_qual = get_mean_quality("filt_ont", "Filt_ONT_") if ont_file != "None" else None
    pacbio_raw_qual = get_mean_quality("raw_pacbio", "RawPacBio_") if pacbio_file != "None" else None
    pacbio_filtered_qual = get_mean_quality("filt_pacbio", "FiltPacBio_") if pacbio_file != "None" else None

    candidates = {}
    if ont_file != "None":
        if ont_corrected_qual and (not ont_filtered_qual or ont_corrected_qual >= ont_filtered_qual):
            candidates[ont_corrected_qual] = ont_file
        elif ont_filtered_qual:
            candidates[ont_filtered_qual] = ont_file.replace("corrected", "filtered")
    if pacbio_file != "None":
        if pacbio_filtered_qual and (not pacbio_raw_qual or pacbio_filtered_qual >= pacbio_raw_qual):
            candidates[pacbio_filtered_qual] = pacbio_file
        elif pacbio_raw_qual:
            candidates[pacbio_raw_qual] = pacbio_file.replace("filtered", "")

    if not candidates:
        print("SKIP:\tNo long reads available for selection")
        return

    highest_qual = max(candidates.keys())
    best_file = candidates[highest_qual]
    print(f"NOTE:\tSelected highest quality long reads: {best_file} with mean quality {highest_qual}")
    subprocess.run(["ln", "-sf", best_file, "highest_qual_long_reads.fastq.gz"], check=True)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python select_long_reads.py <ont_file>" 
              "<pacbio_file> <species_id>", file=sys.stderr)
        sys.exit(1)
    select_long_reads(sys.argv[1], sys.argv[2], sys.argv[3])