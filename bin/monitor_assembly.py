#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
monitor_assembly.py

Non-lethal runtime monitor for long-running assembly subprocesses.

Watches an assembler's working directory and reports elapsed time, output
throughput (MB/min), and a stall warning when no new bytes are written for a
configurable window. It never terminates the watched process; it only reports,
so the operator decides whether to keep waiting.

Motivation:
    MaSuRCA's create_mega_reads stage can peg a single core at ~100% CPU while
    producing no output for many hours (a pathological hang). CPU usage alone
    therefore cannot tell you whether an assembly is progressing. Watching the
    output directory's byte growth can: a healthy run keeps writing, a hung one
    flatlines even while the CPU stays busy.

Stage:
    Assembly (any) -- intended for MaSuRCA, Flye, SPAdes, hifiasm

Created on 2026-06-13

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os
import time
import threading

from utilities import log_print


# --------------------------------------------------------------
# Directory size sampling
# --------------------------------------------------------------
def _dir_size_and_newest(path):
    """Return ``(total_bytes, newest_mtime, newest_file)`` for *path*.

    Walks *path* recursively and sums file sizes. Files that disappear or are
    unreadable mid-walk are skipped rather than raising.

    Parameters
    ----------
    path : str
        Directory to scan.

    Returns
    -------
    tuple of (int, float, str or None)
        Total bytes under *path*, the most recent modification time seen, and
        the path of the most-recently modified file (``None`` if empty).
    """
    total = 0
    newest_mtime = 0.0
    newest_file = None
    for root, _dirs, files in os.walk(path):
        for fn in files:
            fp = os.path.join(root, fn)
            try:
                st = os.stat(fp)
            except OSError:
                continue
            total += st.st_size
            if st.st_mtime > newest_mtime:
                newest_mtime = st.st_mtime
                newest_file = fp
    return total, newest_mtime, newest_file


# --------------------------------------------------------------
# Background assembly monitor
# --------------------------------------------------------------
class AssemblyMonitor:
    """Background, non-lethal progress monitor for an assembly working directory.

    Samples *watch_dir* every *interval_s* seconds in a daemon thread and logs
    throughput and elapsed time. When no new bytes appear for *stall_warn_s*
    seconds it logs a ``WARNING`` so the operator can decide whether to keep
    waiting. It never signals or kills the watched process.

    Use as a context manager around the assembler subprocess call::

        with AssemblyMonitor(work_dir, label="S11 masurca"):
            run_subprocess_cmd(["bash", "assemble.sh"], False)

    Parameters
    ----------
    watch_dir : str or pathlib.Path
        Directory whose recursive byte growth is tracked.
    label : str, optional
        Tag included in every log line (e.g. ``"<sample> masurca"``).
    interval_s : int, optional
        Seconds between samples. Default 60.
    stall_warn_s : int, optional
        Seconds of zero growth before a stall WARNING is emitted. Default 7200
        (2 hours).
    """
    def __init__(self, watch_dir, label="assembly", interval_s=60,
                 stall_warn_s=2 * 60 * 60):
        self.watch_dir = str(watch_dir)
        self.label = label
        self.interval_s = max(5, int(interval_s))
        self.stall_warn_s = max(self.interval_s, int(stall_warn_s))
        self._stop = threading.Event()
        self._thread = None
        self._warned_stall = False

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exc_type, exc, tb):
        self.stop()
        return False  # never suppress exceptions from the watched block

    def start(self):
        """Launch the monitor daemon thread."""
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def stop(self):
        """Signal the monitor to stop and wait briefly for the thread to exit."""
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=self.interval_s + 5)

    def _run(self):
        start_t = time.monotonic()
        last_bytes, _, _ = _dir_size_and_newest(self.watch_dir)
        last_growth_t = start_t
        # _stop.wait returns True when stop() is called, ending the loop cleanly.
        while not self._stop.wait(self.interval_s):
            now = time.monotonic()
            cur_bytes, _newest_mtime, newest_file = _dir_size_and_newest(self.watch_dir)
            elapsed_min = (now - start_t) / 60.0
            delta = cur_bytes - last_bytes
            rate_mb_min = (delta / (1024 * 1024)) / (self.interval_s / 60.0)

            if delta > 0:
                last_growth_t = now
                self._warned_stall = False
                log_print(
                    f"NOTE:\t[{self.label} monitor] elapsed {elapsed_min:.0f} min, "
                    f"+{delta / 1024 / 1024:.1f} MB this interval "
                    f"({rate_mb_min:.1f} MB/min), total {cur_bytes / 1024 / 1024:.0f} MB"
                )
            else:
                quiet_min = (now - last_growth_t) / 60.0
                if (now - last_growth_t) >= self.stall_warn_s:
                    # Emit the WARNING once per stall, then a steady reminder each
                    # interval so a tailing operator keeps seeing it.
                    log_print(
                        f"WARNING:\t[{self.label} monitor] no output growth for "
                        f"{quiet_min:.0f} min (elapsed {elapsed_min:.0f} min). "
                        f"The assembler may be stalled even if CPU is busy. "
                        f"Newest file: {newest_file}. Not killing the process; "
                        f"decide whether to keep waiting or stop it manually."
                    )
                    self._warned_stall = True
                else:
                    log_print(
                        f"NOTE:\t[{self.label} monitor] elapsed {elapsed_min:.0f} min, "
                        f"no new output for {quiet_min:.0f} min (warns at "
                        f"{self.stall_warn_s / 60:.0f} min of no growth)"
                    )
            last_bytes = cur_bytes
