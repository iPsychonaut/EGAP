#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
utilities.py

Backwards-compatibility re-export shim.

Historical ``utilities.py`` (Aug 2023 - Sep 2025) contained subprocess
helpers, logging globals, filesystem / FASTA-FASTQ helpers, and the
per-sample CSV driver in a single ~750-line module.  In v3.4.0 these
were split into four focused modules so each concern can evolve
independently and so the dependency direction is explicit:

    sh.py         -- subprocess execution (leaf)
    log.py        -- ANSI-coloured logging + file logging (leaf)
    fs.py         -- to_abs, validate_fasta, pigz wrappers, md5, etc.
                     depends on sh
    sample_csv.py -- SampleContext, AssemblerStage, get_current_row_data,
                     gen_sample_stats_dict, analyze_nanostats,
                     select_long_reads, load_sample_context
                     depends on fs

This shim re-exports every public symbol the old module exposed so
existing ``from utilities import run_subprocess_cmd`` lines across the
codebase keep working without modification.  New code should prefer
importing directly from the focused module.

Created on Wed Aug 16 2023

Updated on 2026-05-24 (split into sh / log / fs / sample_csv)

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
# flake8: noqa: F401  (the whole point of this module is to re-export)

from sh import (
    get_resource_values,
    run_subprocess_cmd,
)
from log import (
    generate_log_file,
    initialize_logging_environment,
    log_print,
)
from fs import (
    calculate_genome_coverage,
    md5_check,
    move_file_up,
    pigz_compress,
    pigz_decompress,
    sum_fasta_bases_with_pigz_safe,
    sum_fastq_bases_with_pigz_safe,
    to_abs,
    validate_fasta,
)
from sample_csv import (
    AssemblerStage,
    SampleContext,
    analyze_nanostats,
    gen_sample_stats_dict,
    get_current_row_data,
    load_sample_context,
    select_long_reads,
)

__all__ = [
    # sh
    "get_resource_values",
    "run_subprocess_cmd",
    # log
    "generate_log_file",
    "initialize_logging_environment",
    "log_print",
    # fs
    "calculate_genome_coverage",
    "md5_check",
    "move_file_up",
    "pigz_compress",
    "pigz_decompress",
    "sum_fasta_bases_with_pigz_safe",
    "sum_fastq_bases_with_pigz_safe",
    "to_abs",
    "validate_fasta",
    # sample_csv
    "AssemblerStage",
    "SampleContext",
    "analyze_nanostats",
    "gen_sample_stats_dict",
    "get_current_row_data",
    "load_sample_context",
    "select_long_reads",
]
