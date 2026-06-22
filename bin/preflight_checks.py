#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
preflight_checks.py

Fail-fast startup checks for EGAP. Run once before any sample is processed so
misconfiguration aborts loudly with a clear message instead of cascading into
dozens of confusing per-step failures.

Two checks:
    1. Required tools are on PATH. The common failure is launching EGAP with
       the base conda Python instead of the EGAP_env Python, which leaves every
       bioinformatics tool unavailable. We probe a small set of stable-named
       sentinel tools that live only in EGAP_env; if they are absent the
       environment is wrong and we stop.
    2. The requested --ram_gb does not exceed the machine's actual RAM. Asking
       for more memory than exists drove the Linux OOM killer to terminate
       Pilon mid-polish; we refuse an over-commit up front.

Stage:
    Pre-flight (orchestrator startup)

Created on 2026-06-21

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os
import sys
import shutil


# Stable-named binaries that exist in EGAP_env and nowhere in base conda.
# If these are missing, the wrong environment is active -- abort.
SENTINEL_TOOLS = [
    "samtools", "minimap2", "masurca", "flye", "racon", "pilon",
    "kraken2", "datasets", "prefetch", "fasterq-dump", "busco", "compleasm",
]


def _actual_ram_gb():
    """Return the machine's total RAM in GB, or ``None`` if it can't be read."""
    try:
        import psutil
        return int(psutil.virtual_memory().total / (1024 ** 3))
    except Exception:
        pass
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemTotal:"):
                    return int(int(line.split()[1]) / (1024 ** 2))
    except Exception:
        pass
    return None


def run_preflight(ram_gb, required_tools=None, abort=True):
    """Run startup checks; abort loudly on failure.

    Parameters
    ----------
    ram_gb : int or str
        The ``--ram_gb`` value requested for the run.
    required_tools : list of str, optional
        Override the sentinel tool list (mainly for testing).
    abort : bool
        When ``True`` (default), ``sys.exit(1)`` on any failure. When ``False``,
        return the list of problem strings instead (used by the TUI to render
        the message in its own log before exiting).

    Returns
    -------
    list of str
        Problem descriptions. Empty when all checks pass. (Only returned when
        *abort* is ``False``; otherwise the process exits on failure.)
    """
    problems = []

    # 1) Required tools on PATH (wrong-environment guard).
    tools = required_tools if required_tools is not None else SENTINEL_TOOLS
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        problems.append(
            "Required tools are not on PATH: " + ", ".join(missing) + "\n"
            f"\t  Python in use : {sys.executable}\n"
            "\t  This usually means EGAP_env is not active. Launch with the env, e.g.\n"
            "\t    conda run -n EGAP_env --no-capture-output python EGAP.py ...\n"
            "\t  or activate it INSIDE your shell/tmux session first."
        )

    # 2) Requested RAM vs. actual RAM (OOM guard).
    try:
        req_gb = int(float(ram_gb))
    except (TypeError, ValueError):
        req_gb = None
    actual_gb = _actual_ram_gb()
    if req_gb is not None and actual_gb:
        if req_gb > actual_gb:
            problems.append(
                f"--ram_gb {req_gb} exceeds the machine's RAM (~{actual_gb} GB). "
                "Lower --ram_gb, or raise the VM memory (WSL: set memory in "
                ".wslconfig, then `wsl --shutdown`)."
            )
        elif req_gb > 0.9 * actual_gb:
            # Not fatal, but warn: leaves little headroom for the OS + tools.
            print(f"WARN:\t--ram_gb {req_gb} is over 90% of the ~{actual_gb} GB "
                  f"available; consider leaving more headroom.")

    if problems:
        msg = "ERROR:\tPre-flight checks failed:\n" + "\n".join(f"  - {p}" for p in problems)
        if abort:
            print("\n" + msg + "\n", file=sys.stderr)
            print(msg)
            sys.exit(1)
        return problems

    note = f"PASS:\tPre-flight OK: {len(tools)} sentinel tools found"
    if actual_gb:
        note += f"; --ram_gb {req_gb} within ~{actual_gb} GB RAM"
    print(note + ".")
    return []
