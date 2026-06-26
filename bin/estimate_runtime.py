#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
estimate_runtime.py

Pre-flight runtime calculator for the EGAP assemblers.

Given the resources available (CPU threads, RAM) and the data provided
(estimated genome size and the per-technology read volume), this module
produces an order-of-magnitude wall-clock estimate for each assembler EGAP can
run: Flye, MaSuRCA, SPAdes, and hifiasm. It also estimates peak RAM and flags
when an assembler is not applicable to the available read types or when the
estimated memory exceeds the provided budget.

These are heuristic models, not measurements. Each assembler's cost is modelled
as core-hours per Gbp of its dominant input, divided by the usable thread count.
The coefficients below are deliberately conservative ranges meant for relative
comparison and for catching pathological cases (e.g. a MaSuRCA mega-reads run
that will take days, or a SPAdes run that will exceed RAM). Calibrate the
COST_MODEL constants against real EGAP runs to tighten them over time.

Stage:
    Pre-assembly planning (informational; never fatal)

Created on 2026-06-13

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os
import re
import sys

import pandas as pd

from utilities import log_print, initialize_logging_environment, load_sample_context


# --------------------------------------------------------------
# Heuristic cost model (calibrate against real runs)
# --------------------------------------------------------------
# Per assembler:
#   stages       : list of (driver, core_hr_per_gbp_low, core_hr_per_gbp_high, parallel_eff)
#                  where *driver* names which read pool feeds the stage
#                  ("illumina" | "long" | "hifi"). Stages run sequentially, so
#                  their wall-clock contributions add.
#   graph_hr_per_mbp : small fixed overhead that scales with genome size
#                      (graph build, consensus, finalization), low/high.
#   ram_base_gb / ram_per_mbp / ram_per_gbp : rough peak-RAM model.
#   needs        : read types that must be present for the assembler to run.
COST_MODEL = {
    "flye": {
        "stages": [("long", 2.0, 6.0, 0.70)],
        "graph_hr_per_mbp": (0.004, 0.012),
        "ram_base_gb": 8.0, "ram_per_mbp": 0.04, "ram_per_gbp": 0.6,
        "needs": ("long",),
    },
    "masurca": {
        # super-reads from Illumina, then mega-reads from long reads. The
        # mega-reads (create_mega_reads) stage parallelizes poorly and is the
        # high-variance, hang-prone one, hence the low efficiency and wide band.
        "stages": [("illumina", 3.0, 8.0, 0.60),
                   ("long", 15.0, 60.0, 0.25)],
        "graph_hr_per_mbp": (0.02, 0.06),
        "ram_base_gb": 16.0, "ram_per_mbp": 0.05, "ram_per_gbp": 1.0,
        "needs": ("illumina",),
    },
    "spades": {
        "stages": [("illumina", 4.0, 12.0, 0.60)],
        "graph_hr_per_mbp": (0.01, 0.03),
        "ram_base_gb": 16.0, "ram_per_mbp": 0.20, "ram_per_gbp": 1.5,
        "needs": ("illumina",),
    },
    "hifiasm": {
        "stages": [("hifi", 1.5, 4.0, 0.80)],
        "graph_hr_per_mbp": (0.004, 0.010),
        "ram_base_gb": 8.0, "ram_per_mbp": 0.10, "ram_per_gbp": 1.0,
        "needs": ("hifi",),
    },
}

# Fallback coverage assumed when a read file is not yet on disk (estimate runs
# before preprocessing) so the calculator can still give a ballpark.
DEFAULT_COVERAGE = {"illumina": 50.0, "long": 30.0, "hifi": 30.0}

# FASTQ sequence-byte fraction: a 4-line record is header + seq + '+' + qual,
# with seq and qual roughly equal length, so sequence bases are ~47% of bytes.
_FASTQ_SEQ_FRACTION = 0.47
_GZ_EXPANSION = 4.0  # rough gzip ratio for FASTQ


# --------------------------------------------------------------
# Input measurement
# --------------------------------------------------------------
def parse_genome_size(est_size):
    """Parse an ``EST_SIZE`` value like ``"43m"`` or ``"1.2g"`` into base pairs.

    Parameters
    ----------
    est_size : str or float
        Estimated genome size with an ``m``/``g`` suffix.

    Returns
    -------
    int or None
        Genome size in base pairs, or ``None`` if it cannot be parsed.
    """
    m = re.match(r"^(\d+(?:\.\d+)?)(\D+)$", str(est_size).strip())
    if not m:
        return None
    multipliers = {"m": 10 ** 6, "g": 10 ** 9, "k": 10 ** 3}
    return int(float(m.group(1)) * multipliers.get(m.group(2).strip().lower(), 10 ** 6))


def estimate_fastq_bases(path):
    """Estimate the number of sequenced bases in a FASTQ file from its size.

    Uses a byte-size proxy rather than reading the file, so it is fast on
    multi-GB inputs. ``.gz`` inputs are expanded by a rough constant factor.

    Parameters
    ----------
    path : str or None
        Path to a FASTQ (optionally ``.gz``) file.

    Returns
    -------
    tuple of (int or None, str)
        ``(bases, source)`` where *source* describes how the figure was
        obtained (``"size-proxy"``, ``"gz-proxy"``, or ``"missing"``).
    """
    if not path or not isinstance(path, str) or not os.path.exists(path):
        return None, "missing"
    size = os.path.getsize(path)
    if size == 0:
        return None, "empty"
    if path.endswith(".gz"):
        return int(size * _GZ_EXPANSION * _FASTQ_SEQ_FRACTION), "gz-proxy"
    return int(size * _FASTQ_SEQ_FRACTION), "size-proxy"


def _first_existing(*paths):
    """Return the first path that exists and is non-empty, else the first non-empty string."""
    fallback = None
    for p in paths:
        if isinstance(p, str) and p.strip():
            if fallback is None:
                fallback = p
            if os.path.exists(p) and os.path.getsize(p) > 0:
                return p
    return fallback


def collect_metrics(sample_id, input_tsv, output_dir, cpu_threads, ram_gb):
    """Resolve genome size, resources, and per-technology read volume for a sample.

    Mirrors the read-path conventions used by the ``assemble_*`` modules
    (deduplicated Illumina, ``*_highest_mean_qual_long_reads`` ONT/PacBio, with
    raw/SRA fallbacks) so the estimate reflects the files the assemblers will
    actually consume.

    Returns
    -------
    dict
        Keys: ``sample_id``, ``threads``, ``ram_gb``, ``genome_bp``,
        ``illumina_gbp``, ``long_gbp``, ``hifi_gbp``, ``has`` (set of present
        read types), ``coverage`` (per type), and ``sources`` (how each volume
        was derived).
    """
    ctx = load_sample_context(sample_id, input_tsv, output_dir, cpu_threads, ram_gb)
    s = ctx.current_series
    species_id = s["SPECIES_ID"]
    species_dir = os.path.join(ctx.output_dir, str(species_id))

    def sra_path(col, subdir, suffix):
        sra = s[col]
        if pd.notna(sra):
            return os.path.join(species_dir, subdir, f"{sra}{suffix}")
        return None

    # Illumina: prefer dedup, fall back to raw, then SRA-derived.
    illu_f = _first_existing(
        os.path.join(species_dir, "Illumina", f"{species_id}_illu_forward_dedup.fastq"),
        s["ILLUMINA_RAW_F_READS"] if pd.notna(s["ILLUMINA_RAW_F_READS"]) else None,
        sra_path("ILLUMINA_SRA", "Illumina", "_1.fastq"),
    )
    illu_r = _first_existing(
        os.path.join(species_dir, "Illumina", f"{species_id}_illu_reverse_dedup.fastq"),
        s["ILLUMINA_RAW_R_READS"] if pd.notna(s["ILLUMINA_RAW_R_READS"]) else None,
        sra_path("ILLUMINA_SRA", "Illumina", "_2.fastq"),
    )
    ont = _first_existing(
        os.path.join(species_dir, "ONT", f"{species_id}_ONT_highest_mean_qual_long_reads.fastq"),
        s["ONT_RAW_READS"] if pd.notna(s["ONT_RAW_READS"]) else None,
        sra_path("ONT_SRA", "ONT", ".fastq"),
    )
    pacbio = _first_existing(
        os.path.join(species_dir, "PacBio", f"{species_id}_PacBio_highest_mean_qual_long_reads.fastq"),
        s["PACBIO_RAW_READS"] if pd.notna(s["PACBIO_RAW_READS"]) else None,
        sra_path("PACBIO_SRA", "PacBio", ".fastq"),
    )

    genome_bp = parse_genome_size(s["EST_SIZE"])
    genome_mbp = (genome_bp / 1e6) if genome_bp else None

    def gbp(*paths):
        total = 0
        src = "missing"
        for p in paths:
            bases, source = estimate_fastq_bases(p)
            if bases:
                total += bases
                src = source
        return (total / 1e9) if total else 0.0, src

    illumina_gbp, illu_src = gbp(illu_f, illu_r)
    long_gbp, long_src = gbp(ont)          # ONT drives the "long" pool for Flye/MaSuRCA/SPAdes
    hifi_gbp, hifi_src = gbp(pacbio)        # PacBio HiFi pool for hifiasm (and Flye if no ONT)

    # PacBio can also feed Flye/MaSuRCA as "long" when no ONT is present.
    if long_gbp == 0.0 and hifi_gbp > 0.0:
        long_gbp, long_src = hifi_gbp, hifi_src

    has = set()
    if illumina_gbp > 0 or pd.notna(s["ILLUMINA_SRA"]):
        has.add("illumina")
    if long_gbp > 0 or pd.notna(s["ONT_SRA"]) or pd.notna(s["PACBIO_SRA"]):
        has.add("long")
    if hifi_gbp > 0 or pd.notna(s["PACBIO_SRA"]):
        has.add("hifi")

    # Fill assumed volumes when reads are not yet on disk, so the estimate still
    # returns a ballpark (flagged via source = "assumed").
    sources = {"illumina": illu_src, "long": long_src, "hifi": hifi_src}
    if genome_bp:
        if "illumina" in has and illumina_gbp == 0.0:
            illumina_gbp = DEFAULT_COVERAGE["illumina"] * genome_bp / 1e9
            sources["illumina"] = "assumed"
        if "long" in has and long_gbp == 0.0:
            long_gbp = DEFAULT_COVERAGE["long"] * genome_bp / 1e9
            sources["long"] = "assumed"
        if "hifi" in has and hifi_gbp == 0.0:
            hifi_gbp = DEFAULT_COVERAGE["hifi"] * genome_bp / 1e9
            sources["hifi"] = "assumed"

    coverage = {}
    if genome_bp:
        coverage = {
            "illumina": illumina_gbp * 1e9 / genome_bp,
            "long": long_gbp * 1e9 / genome_bp,
            "hifi": hifi_gbp * 1e9 / genome_bp,
        }

    return {
        "sample_id": sample_id,
        "threads": _safe_int(cpu_threads, 1),
        "ram_gb": _safe_int(ram_gb, 0),
        "genome_bp": genome_bp,
        "genome_mbp": genome_mbp,
        "illumina_gbp": illumina_gbp,
        "long_gbp": long_gbp,
        "hifi_gbp": hifi_gbp,
        "has": has,
        "coverage": coverage,
        "sources": sources,
    }


def _safe_int(v, default):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return default


# --------------------------------------------------------------
# Per-assembler estimate
# --------------------------------------------------------------
def estimate_one(assembler, metrics):
    """Estimate wall-clock and peak RAM for a single assembler.

    Parameters
    ----------
    assembler : str
        One of ``"flye"``, ``"masurca"``, ``"spades"``, ``"hifiasm"``.
    metrics : dict
        Output of :func:`collect_metrics`.

    Returns
    -------
    dict
        Keys: ``assembler``, ``applicable`` (bool), ``reason`` (str when not
        applicable), ``low_hr``, ``high_hr``, ``peak_ram_gb``, ``ram_risk``
        (bool), ``notes`` (list of str).
    """
    model = COST_MODEL[assembler]
    notes = []

    # Applicability: every required read pool must be present.
    driver_gbp = {
        "illumina": metrics["illumina_gbp"],
        "long": metrics["long_gbp"],
        "hifi": metrics["hifi_gbp"],
    }
    for need in model["needs"]:
        if need not in metrics["has"]:
            label = {"illumina": "Illumina paired-end",
                     "long": "long (ONT/PacBio)",
                     "hifi": "PacBio HiFi"}[need]
            return {"assembler": assembler, "applicable": False,
                    "reason": f"no {label} reads", "low_hr": None, "high_hr": None,
                    "peak_ram_gb": None, "ram_risk": False, "notes": []}

    threads = max(1, metrics["threads"])
    genome_mbp = metrics["genome_mbp"] or 0.0

    low_hr = high_hr = 0.0
    for driver, ch_low, ch_high, eff in model["stages"]:
        gbp = driver_gbp.get(driver, 0.0)
        usable = max(1.0, threads * eff)
        low_hr += ch_low * gbp / usable
        high_hr += ch_high * gbp / usable
    g_low, g_high = model["graph_hr_per_mbp"]
    low_hr += g_low * genome_mbp
    high_hr += g_high * genome_mbp

    # Peak RAM estimate and budget check.
    total_input_gbp = sum(driver_gbp.get(d, 0.0) for d, *_ in model["stages"])
    peak_ram = (model["ram_base_gb"]
                + model["ram_per_mbp"] * genome_mbp
                + model["ram_per_gbp"] * total_input_gbp)
    ram_risk = bool(metrics["ram_gb"]) and peak_ram > metrics["ram_gb"]
    if ram_risk:
        notes.append(f"est. peak RAM {peak_ram:.0f} GB exceeds the {metrics['ram_gb']} GB "
                     f"budget; swapping can multiply runtime")

    # Per-assembler flavour notes.
    if assembler == "masurca" and metrics["long_gbp"] > 0:
        cov = metrics["coverage"].get("long", 0.0)
        notes.append(f"mega-reads stage is high-variance and hang-prone; "
                     f"long-read coverage ~{cov:.0f}x")
        if cov > 40:
            notes.append("consider reducing long-read coverage (<=30x) to cut mega-reads time")
    if metrics["sources"].get("illumina") == "assumed" or metrics["sources"].get("long") == "assumed":
        notes.append("read volume assumed from genome size (reads not yet on disk)")

    return {"assembler": assembler, "applicable": True, "reason": "",
            "low_hr": low_hr, "high_hr": high_hr,
            "peak_ram_gb": peak_ram, "ram_risk": ram_risk, "notes": notes}


def estimate_all(metrics):
    """Return a list of :func:`estimate_one` results for all four assemblers."""
    return [estimate_one(a, metrics) for a in ("flye", "masurca", "spades", "hifiasm")]


# --------------------------------------------------------------
# Formatting / logging
# --------------------------------------------------------------
def _fmt_hours(h):
    if h is None:
        return "-"
    if h < 1:
        return f"{h * 60:.0f} min"
    if h < 48:
        return f"{h:.1f} h"
    return f"{h / 24:.1f} days"


def format_estimate_line(est):
    """One-line human summary for a single assembler estimate."""
    a = est["assembler"]
    if not est["applicable"]:
        return f"{a:8s}: n/a ({est['reason']})"
    window = f"{_fmt_hours(est['low_hr'])} - {_fmt_hours(est['high_hr'])}"
    ram = f"~{est['peak_ram_gb']:.0f} GB RAM" + ("  [RAM RISK]" if est["ram_risk"] else "")
    return f"{a:8s}: est. {window}  ({ram})"


def format_table(metrics, estimates):
    """Render the full per-assembler comparison table as a string block."""
    lines = []
    lines.append("=" * 74)
    lines.append(f"RUNTIME ESTIMATE (heuristic) - sample {metrics['sample_id']}")
    gsize = f"{metrics['genome_mbp']:.1f} Mbp" if metrics["genome_mbp"] else "unknown"
    lines.append(f"  resources : {metrics['threads']} threads, {metrics['ram_gb']} GB RAM")
    cov = metrics["coverage"]
    lines.append(f"  genome    : {gsize}")
    lines.append(f"  data      : "
                 f"Illumina {metrics['illumina_gbp']:.2f} Gbp (~{cov.get('illumina', 0):.0f}x), "
                 f"long {metrics['long_gbp']:.2f} Gbp (~{cov.get('long', 0):.0f}x), "
                 f"HiFi {metrics['hifi_gbp']:.2f} Gbp (~{cov.get('hifi', 0):.0f}x)")
    lines.append("-" * 74)
    lines.append(f"  {'assembler':9s} {'estimate':>20s}   {'peak RAM':>9s}   notes")
    for est in estimates:
        if not est["applicable"]:
            lines.append(f"  {est['assembler']:9s} {'n/a':>20s}   {'-':>9s}   {est['reason']}")
            continue
        window = f"{_fmt_hours(est['low_hr'])} - {_fmt_hours(est['high_hr'])}"
        ram = f"~{est['peak_ram_gb']:.0f} GB" + ("*" if est["ram_risk"] else "")
        note = est["notes"][0] if est["notes"] else ""
        lines.append(f"  {est['assembler']:9s} {window:>20s}   {ram:>9s}   {note}")
        for extra in est["notes"][1:]:
            lines.append(f"  {'':9s} {'':>20s}   {'':>9s}   {extra}")
    lines.append("-" * 74)
    lines.append("  Heuristic order-of-magnitude figures for planning, not guarantees.")
    lines.append("  '*' = estimated peak RAM exceeds the provided budget.")
    lines.append("=" * 74)
    return "\n".join(lines)


def log_estimate_for(assembler, sample_id, input_tsv, output_dir, cpu_threads, ram_gb):
    """Compute and log a single assembler's estimate (for use inside assemble_*.py).

    Never raises: any failure is logged as a NOTE and swallowed so the estimate
    can never block an assembly.
    """
    try:
        metrics = collect_metrics(sample_id, input_tsv, output_dir, cpu_threads, ram_gb)
        est = estimate_one(assembler, metrics)
        log_print(f"NOTE:\t[runtime estimate] {format_estimate_line(est)}")
        for n in est["notes"]:
            log_print(f"NOTE:\t[runtime estimate] {assembler}: {n}")
    except Exception as exc:  # informational only; must not break assembly
        log_print(f"NOTE:\t[runtime estimate] unavailable for {assembler}: {exc}")


# --------------------------------------------------------------
# Standalone entry point: print the full comparison table
# --------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python3 estimate_runtime.py <sample_id> <input_tsv> "
              "<output_dir> <cpu_threads> <ram_gb>", file=sys.stderr)
        sys.exit(1)

    initialize_logging_environment(sys.argv[3], sys.argv[1])
    try:
        _metrics = collect_metrics(sys.argv[1], sys.argv[2], sys.argv[3],
                                   sys.argv[4], sys.argv[5])
        _estimates = estimate_all(_metrics)
        log_print("\n" + format_table(_metrics, _estimates))
    except Exception as exc:
        log_print(f"NOTE:\tRuntime estimate unavailable: {exc}")
    # Always succeed: this step is purely informational.
    sys.exit(0)
