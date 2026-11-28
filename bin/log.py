#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
log.py

Logging helpers (file + ANSI-coloured stdout) for EGAP modules.

Extracted from :mod:`utilities` in v3.4.1.  Module-level
``DEFAULT_LOG_FILE`` singleton is set by
:func:`initialize_logging_environment` and read by :func:`log_print`,
so callers can ``log_print(msg)`` without threading a log handle
through every call site.

Stage:
    Cross-cutting (logging surface used by every pipeline stage)

Created on Wed Aug 16 2023

Updated on 2026-05-24

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os
from datetime import datetime


# --------------------------------------------------------------
# Module-level singleton populated by initialize_logging_environment.
# Declared at import time (rather than only inside the initializer) so
# ``global DEFAULT_LOG_FILE`` in log_print has a bound name even if a
# caller forgets to initialize first - they'll get a clear TypeError
# from the file open() instead of a NameError.
# --------------------------------------------------------------
DEFAULT_LOG_FILE = None


# --------------------------------------------------------------
# Log and print messages with color
# --------------------------------------------------------------
def log_print(input_message, log_file=None):
    """Log a message to a file and print it with ANSI-colored output.

    Prepends a timestamp, writes the message to the log file, and prints
    it in a color determined by the message prefix (e.g. ``ERROR`` ->
    red, ``PASS`` -> green).

    Parameters
    ----------
    input_message : str
        Message text to log and print.
    log_file : str, optional
        Path to the log file.  Defaults to the module-level
        ``DEFAULT_LOG_FILE`` set by ``initialize_logging_environment``.
    """
    global DEFAULT_LOG_FILE
    COLORS = {"grey": "\033[90m",
              "red": "\033[91m",
              "green": "\033[92m",
              "orange": "\033[38;5;208m",
              "yellow": "\033[93m",
              "blue": "\033[94m",
              "magenta": "\033[95m",
              "cyan": "\033[96m",
              "white": "\033[97m",
              "reset": "\033[0m"}
    if log_file is None:
        log_file = DEFAULT_LOG_FILE
    now = datetime.now()
    message = f"[{now:%Y-%m-%d %H:%M:%S}]\t{input_message}"
    message_type_dict = {"NOTE": "blue",
                         "CMD": "cyan",
                         "ERROR": "red",
                         "WARN": "yellow",
                         "PASS": "green",
                         "SKIP": "magenta",
                         "FAIL": "red"}
    print_color = "white"
    for key, value in message_type_dict.items():
        if key.lower() in input_message.lower():
            print_color = value
            break
    try:
        with open(log_file, "a") as file:
            print(message, file=file)
    except TypeError:
        # log_file was None (initialize_logging_environment never ran).
        print(f"UNLOGGED ERROR:\tUnable to load the log file provided: {log_file}")
    except (FileNotFoundError, OSError) as exc:
        # Defensive: a stale relative path (chdir'd subprocess) or a
        # transient FS hiccup should never crash a calling pipeline
        # stage that's just trying to log progress.  Fall back to stderr
        # so the message is at least visible in pipeline.log.
        print(
            f"UNLOGGED ERROR:\tCould not write to {log_file}: {exc}.  "
            f"Original message follows:\n  {message}"
        )
    color_code = COLORS.get(print_color, COLORS["white"])
    print(f"{color_code}{message}{COLORS['reset']}")


# --------------------------------------------------------------
# Set up logging environment
# --------------------------------------------------------------
def initialize_logging_environment(INPUT_FOLDER, sample_id=None):
    """Initialize the logging environment based on the input folder.

    Sets the module-level ``DEFAULT_LOG_FILE`` global and makes sure the log
    file exists. The file is never emptied: every stage of a sample is its own
    process and calls this, so truncating here erased the earlier stages.

    Parameters
    ----------
    INPUT_FOLDER : str
        Output folder path used to determine the log file location.
        When *sample_id* is ``None``, the log file is written as
        ``<INPUT_FOLDER>/<basename>_log.txt``.
    sample_id : str, optional
        When provided, the log file is written per-sample as
        ``<INPUT_FOLDER>/<sample_id>_log.txt`` so that each sample in
        a multi-sample TSV run gets its own log file.
    """
    global DEFAULT_LOG_FILE
    print(INPUT_FOLDER)
    input_file_path = f"{INPUT_FOLDER}/{sample_id or INPUT_FOLDER.split('/')[-1]}_log.txt"
    # Anchor the log path to its absolute form so a later subprocess that
    # ``chdir``s (e.g. ``assemble_masurca`` running ``assemble.sh`` from
    # the ``masurca_assembly`` subdir) can still resolve DEFAULT_LOG_FILE.
    # Without this, ``log_print`` raises FileNotFoundError when invoked
    # after any chdir, which masks the real per-stage error message and
    # crashes graceful-recovery paths that depend on logging.
    DEFAULT_LOG_FILE = os.path.abspath(input_file_path)
    print(DEFAULT_LOG_FILE)
    open(DEFAULT_LOG_FILE, "a").close()
