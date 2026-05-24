#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
log.py

Logging helpers (file + ANSI-coloured stdout) for EGAP modules.

Extracted from :mod:`utilities` in v3.4.0.  Module-level
``DEFAULT_LOG_FILE`` and ``ENVIRONMENT_TYPE`` singletons are set by
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
import platform
from datetime import datetime


# --------------------------------------------------------------
# Module-level singletons populated by initialize_logging_environment.
# Declared at import time (rather than only inside the initializer) so
# ``global DEFAULT_LOG_FILE`` in log_print has a bound name even if a
# caller forgets to initialize first - they'll get a clear TypeError
# from the file open() instead of a NameError.
# --------------------------------------------------------------
DEFAULT_LOG_FILE = None
ENVIRONMENT_TYPE = None


# --------------------------------------------------------------
# Create or manage a log file
# --------------------------------------------------------------
def generate_log_file(log_file_path, use_numerical_suffix=False):
    """Generate a log file, optionally with a numerical suffix if it exists.

    Parameters
    ----------
    log_file_path : str
        Desired path for the log file.
    use_numerical_suffix : bool
        If ``True``, append an incrementing integer suffix rather than
        overwriting an existing file.

    Returns
    -------
    str
        Path to the created or selected log file.
    """
    if os.path.exists(log_file_path) and use_numerical_suffix:
        counter = 1
        base, ext = os.path.splitext(log_file_path)
        new_log_file_path = f"{base}_{counter}{ext}"
        while os.path.exists(new_log_file_path):
            counter += 1
            new_log_file_path = f"{base}_{counter}{ext}"
        log_file_path = new_log_file_path
    else:
        open(log_file_path, "w").close()
    return log_file_path


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
        print(f"UNLOGGED ERROR:\tUnable to load the log file provided: {log_file}")
    color_code = COLORS.get(print_color, COLORS["white"])
    print(f"{color_code}{message}{COLORS['reset']}")


# --------------------------------------------------------------
# Set up logging environment
# --------------------------------------------------------------
def initialize_logging_environment(INPUT_FOLDER, sample_id=None):
    """Initialize the logging environment based on the input folder.

    Sets the module-level ``DEFAULT_LOG_FILE`` and ``ENVIRONMENT_TYPE``
    globals and creates the log file.  The log file path is adjusted for
    WSL/Linux (``/mnt/<drive>/...``) when running on a non-Windows OS.

    Parameters
    ----------
    INPUT_FOLDER : str
        Output folder path used to determine the log file location.
        When *sample_id* is ``None``, the log file is written as
        ``<INPUT_FOLDER>/<basename>_log.txt``.
    sample_id : str, optional
        When provided, the log file is written per-sample as
        ``<INPUT_FOLDER>/<sample_id>_log.txt`` so that each sample in
        a multi-sample CSV run gets its own log file.
    """
    global DEFAULT_LOG_FILE, ENVIRONMENT_TYPE
    print(INPUT_FOLDER)
    if sample_id:
        input_file_path = f"{INPUT_FOLDER}/{sample_id}_log.txt"
    else:
        input_file_path = f"{INPUT_FOLDER}/{INPUT_FOLDER.split('/')[-1]}_log.txt"
    os_name = platform.system()
    if os_name == "Windows":
        print("UNLOGGED:\tWINDOWS ENVIRONMENT")
        ENVIRONMENT_TYPE = "WIN"
    elif os_name in ["Linux", "Darwin"]:
        drive, path_without_drive = os.path.splitdrive(input_file_path)
        if drive:
            drive_letter = drive.strip(":\\/")
            path_without_drive_mod = path_without_drive.replace("\\", "/")
            input_file_path = f"/mnt/{drive_letter.lower()}{path_without_drive_mod}"
        print("UNLOGGED:\tLINUX/WSL/MAC ENVIRONMENT")
        ENVIRONMENT_TYPE = "LINUX/WSL/MAC"
    else:
        print(f"UNLOGGED ERROR:\tUnsupported OS: {os_name}")
        return
    print(input_file_path)
    run_log = generate_log_file(input_file_path, use_numerical_suffix=False)
    DEFAULT_LOG_FILE = run_log
