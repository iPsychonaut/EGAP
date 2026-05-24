#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sh.py

Subprocess execution helpers shared across EGAP modules.

Extracted from :mod:`utilities` in v3.4.0 so the subprocess concern is
isolated from logging, filesystem, and DataFrame helpers.  Callers that
still ``from utilities import run_subprocess_cmd`` continue to work via
the re-export shim in :mod:`utilities`.

Stage:
    Cross-cutting (used by every assembler / preprocessor / QC step)

Created on Wed Aug 16 2023

Updated on 2026-05-24

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import math
import subprocess


# --------------------------------------------------------------
# Calculate resource allocation
# --------------------------------------------------------------
def get_resource_values(percent_resources, total_cpu, total_ram):
    """Calculate CPU threads and RAM based on a percentage of total resources.

    Parameters
    ----------
    percent_resources : float
        Fraction of total resources to allocate (0.0-1.0).
    total_cpu : int
        Total available CPU threads.
    total_ram : int
        Total available RAM in GB.

    Returns
    -------
    tuple of (int, int)
        ``(cpu_threads, ram_gb)`` - number of threads and RAM in GB.
    """
    cpu_threads = int(math.floor(total_cpu * percent_resources))
    ram_gb = int(total_ram * percent_resources)
    return cpu_threads, ram_gb


# --------------------------------------------------------------
# Execute a subprocess command and log its output
# --------------------------------------------------------------
def run_subprocess_cmd(cmd_list, shell_check):
    """Execute a subprocess command and stream its output to stdout.

    Runs the command using ``subprocess.Popen``, captures combined
    stdout/stderr, and streams each line in real time.  Returns the
    process exit code.

    Parameters
    ----------
    cmd_list : str or list of str
        Command to execute, either as a shell string or as an argument list.
    shell_check : bool
        If ``True``, execute the command through the shell (required when
        *cmd_list* is a string with shell operators).

    Returns
    -------
    int
        The subprocess return code (0 on success, 127 if the executable
        was not found).
    """
    cmd_display = cmd_list if isinstance(cmd_list, str) else ' '.join(cmd_list)
    print(f"CMD:\t{cmd_display}")
    try:
        process = subprocess.Popen(cmd_list, shell=shell_check, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True)
    except (FileNotFoundError, PermissionError) as exc:
        # The executable was not found on PATH or is not executable.
        # This happens when the pipeline is launched from a conda env that
        # does not have the required tool installed.  Return 127 (the
        # conventional "command not found" exit code) instead of crashing.
        tool = cmd_list.split()[0] if isinstance(cmd_list, str) else cmd_list[0]
        print(f"ERROR:\tCould not launch '{tool}': {exc}. "
              f"Make sure the correct conda environment is active "
              f"(e.g. conda activate EGAP_env).")
        return 127
    for line in process.stdout:
        print(line, end="")
    process.wait()
    if process.returncode != 0:
        print(f"NOTE:\tCommand failed with return code {process.returncode}")
    else:
        print(f"PASS:\tSuccessfully processed command: {cmd_display}")
    return process.returncode
