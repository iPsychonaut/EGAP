#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
EGAP_TUI.py

Textual TUI wrapper for the EGAP pipeline controller.

Layout goals:
- Banner along the top
- Parameters panel on the right (top-right)
- CPU/MEM panel on the right (bottom-right)
- Pipeline Progress in the bottom middle
- Streaming pipeline logs on the left

This script mirrors EGAP.py orchestration:
- preprocess_tsv()
- locate_bin_dir()
- processes list
- per-sample/per-process execution, then qc_assessment, then html_reporter


python EGAP_TUI.py -tsv /mnt/c/Users/theda/OneDrive/Desktop/EGAP/EGAP_TUI_test.tsv -o /mnt/c/Users/theda/OneDrive/Desktop/EGAP/EGAP_TUI_test -t 10 -r 24

Created on Sun Feb 22 14:56:27 2026

Updated on 2026-04-16

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""

from __future__ import annotations

import traceback
import argparse
import asyncio
import os
import sys
import time
import psutil
from collections import deque
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from rich.text import Text
from rich.table import Table
from rich.console import Group
from textual.widget import Widget
from textual.app import App, ComposeResult
from textual.containers import Horizontal, Vertical, VerticalScroll
from textual.widgets import Static, Header, Footer, RichLog, DataTable


# CPU Usage History display
HISTORY_LEN = 10  # number of blocks shown per core (__________)

def pct_to_block(p: float) -> Tuple[str, str]:
    """
    Return (char, style) for the most recent utilization sample.
    """
    if p <= 34:
        return "_", "bold bright_green"
    if p <= 67:
        return "▄", "rgb(255,215,0)"
    return "█", "rgb(255,69,0)"


@dataclass
class CoreHistory:
    blocks: Deque[Tuple[str, str]]  # (char, style)


class CpuThreadsHistoryWidget(Widget):
    """
    Two-column per-core utilization history.
    Newest sample is rightmost (closest to the percent), older slides left.
    """

    def __init__(
        self,
        refresh_s: float = 0.5,
        history_len: int = HISTORY_LEN,
        *,
        id: str | None = None,
        classes: str | None = None,
        name: str | None = None,
    ) -> None:
        super().__init__(id=id, classes=classes, name=name)

        self.refresh_s = refresh_s
        self.history_len = history_len
        self.n_cores = psutil.cpu_count(logical=True) or 1

        self.hist: List[CoreHistory] = [
            CoreHistory(
                blocks=deque(
                    [("_", "bold bright_green")] * self.history_len,
                    maxlen=self.history_len,
                )
            )
            for _ in range(self.n_cores)
        ]

        # warm-up sample
        psutil.cpu_percent(percpu=True)

    def on_mount(self) -> None:
        self.set_interval(self.refresh_s, self.tick)

    def tick(self) -> None:
        percpu = psutil.cpu_percent(percpu=True)
        for i in range(min(len(percpu), self.n_cores)):
            ch, style = pct_to_block(percpu[i])
            self.hist[i].blocks.append((ch, style))
        self.refresh()

    def render_core_line(self, core_idx: int, pct: float) -> Text:
        t = Text()
        t.append(f"C{core_idx:<2}  ", style="bright_white")

        for ch, style in self.hist[core_idx].blocks:
            t.append(ch, style)

        t.append(f"  {int(round(pct)):>3}%", style="bright_white")
        return t

    def render(self):
        percpu_now = psutil.cpu_percent(percpu=True)
        if len(percpu_now) < self.n_cores:
            percpu_now = (percpu_now + [0.0] * self.n_cores)[: self.n_cores]

        cpu_total = sum(percpu_now) / max(len(percpu_now), 1)

        cpu_title = Text(f"CPU Total Usage: {cpu_total:5.1f}%", style="bold bright_white")
        mid = (self.n_cores + 1) // 2

        cpu_table = Table.grid(padding=(0, 1))
        cpu_table.add_column(justify="left")
        cpu_table.add_column(justify="left")

        for row_i in range(mid):
            left_core = row_i
            right_core = row_i + mid

            left_txt = self.render_core_line(left_core, percpu_now[left_core])

            if right_core < self.n_cores:
                right_txt = self.render_core_line(right_core, percpu_now[right_core])
                cpu_table.add_row(left_txt, Text(" | ") + right_txt)
            else:
                cpu_table.add_row(left_txt, Text(""))

        vm = psutil.virtual_memory()
        sm = psutil.swap_memory()

        mem_text = Text(f"""\nRAM:  {vm.percent:5.1f}%  {fmt_bytes(vm.used)} / {fmt_bytes(vm.total)}\nSwap: {sm.percent:5.1f}%  {fmt_bytes(sm.used)} / {fmt_bytes(sm.total)}""",
                        style="bold bright_white")

        return Group(cpu_title, cpu_table, mem_text)


# ----------------------------
# Import EGAP controller helpers
# ----------------------------
# EGAP_TUI.py lives in bin/; EGAP.py lives one level up in the project root.
# Insert the project root into sys.path so `import EGAP` works regardless of
# which directory the user runs this script from.
_TUI_DIR     = Path(__file__).resolve().parent          # .../EGAP/bin
_PROJECT_DIR = _TUI_DIR.parent                          # .../EGAP
for _p in (str(_PROJECT_DIR), str(_TUI_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

try:
    import EGAP as egap
except Exception as e:
    raise RuntimeError(
        f"Could not import EGAP.py from {_PROJECT_DIR}.\n"
        f"Import error: {e}"
    )


# ----------------------------
# ANSI banner (pulled from EGAP.py banner style)
# ----------------------------
# Note: this is intentionally “banner-only” (top-of-screen). It uses ANSI sequences; we decode
# to Rich Text for consistent rendering inside Textual.
BANNER_ANSI = f"""
\033[91m.---.\033[92m________\\\033[38;5;208m=\033[96m/\033[92m________\033[91m.---.\033[0m
\033[91m|\033[38;5;208m[\033[93m_\033[38;5;208m]\033[91m|\033[38;2;0;55;200m--------\033[96m/\033[91m=\033[92m\\\033[38;2;0;55;200m--------\033[91m|\033[38;5;208m[\033[93m_\033[38;5;208m]\033[91m|   \033[38;2;0;55;200m.\033[96m---------.  \033[38;2;0;55;200m.\033[96m------.    \033[38;2;0;55;200m.\033[96m------.    \033[38;2;0;55;200m.\033[96m-------.\033[0m
\033[91m`---'\033[96m~~~~~~~(\033[38;5;208m===\033[92m)\033[96m~~~~~~~\033[91m`---'  \033[38;2;0;55;200m/\033[96m|         |\033[38;2;0;55;200m/\033[96m/        \\ \033[38;2;0;55;200m/\033[96m/        \\ \033[38;2;0;55;200m/\033[96m/         \\\033[0m
 \033[92m|\033[38;2;0;55;200m|\033[96m| \033[92m.--     \033[96m\\\033[91m=\033[92m/    \033[92m,--. \033[96m|\033[38;2;0;55;200m|\033[92m|  \033[38;2;0;55;200m| \033[96m|  .------\033[38;2;0;55;200m'\033[96m|   .------\033[38;2;0;55;200m'\033[96m|   .--\033[38;2;0;55;200m.   \033[96m. |   .--\033[38;2;0;55;200m.\033[96m   .\033[0m
 \033[92m|\033[38;2;0;55;200m|\033[96m| \033[38;2;0;55;200m|-      \033[92m/\033[38;5;208m=\033[96m\\    \033[38;2;0;55;200m|  _ \033[96m|\033[38;2;0;55;200m|\033[92m|  \033[38;2;0;55;200m| \033[96m|  |\033[38;2;0;55;200m----"| \033[96m|   |\033[38;2;0;55;200m----"| \033[96m|   |\033[38;2;0;55;200m-'\033[96m|   | |   |\033[38;2;0;55;200m-'\033[96m|   |\033[0m
 \033[92m|\033[38;2;0;55;200m|\033[96m| \033[96m`--    \033[92m(\033[91m===\033[96m)   \033[96m`--' \033[96m|\033[38;2;0;55;200m|\033[92m|  \033[38;2;0;55;200m| \033[96m|  `----.\033[38;2;0;55;200m|\033[96m |   |     \033[38;2;0;55;200m| \033[96m|   +--+   | |   +--+   |\033[0m
 \033[92m|\033[38;2;0;55;200m|\033[96m|         \033[92m\\\033[38;5;208m=\033[96m/         \033[96m|\033[38;2;0;55;200m|\033[92m|  \033[38;2;0;55;200m| \033[96m|       |\033[38;2;0;55;200m| \033[96m|   | \033[38;2;0;55;200m.\033[96m----.|          | |          |\033[0m
 \033[92m|\033[38;2;0;55;200m|\033[96m|         \033[96m/\033[91m=\033[92m\\         \033[96m|\033[38;2;0;55;200m|\033[92m|  \033[38;2;0;55;200m| \033[96m|  .----'\033[38;2;0;55;200m| \033[96m|   |\033[38;2;0;55;200m"\033[96m|    ||   +--+   | |   +------'\033[0m
 \033[92m|\033[38;2;0;55;200m|\033[96m|        \033[96m(\033[38;5;208m===\033[92m)        \033[96m|\033[38;2;0;55;200m|\033[92m|  \033[38;2;0;55;200m| \033[96m|  |\033[38;2;0;55;200m---" | \033[96m|   |\033[38;2;0;55;200m"\033[96m`-.  ||   |\033[38;2;0;55;200m| \033[96m|   | |   |\033[38;2;0;55;200m-----"\033[0m
 \033[92m|\033[38;2;0;55;200m|\033[96m|  \033[96m__     \033[96m\\\033[91m=\033[92m/    \033[96m,__. \033[96m|\033[38;2;0;55;200m|\033[92m|  \033[38;2;0;55;200m| \033[96m|  `------.|   `---'  ||   |\033[38;2;0;55;200m| \033[96m|   | |   |\033[0m
 \033[92m|\033[38;2;0;55;200m|\033[96m| \033[38;2;0;55;200m/__\\    \033[92m/\033[38;5;208m=\033[96m\\    \033[38;2;0;55;200m|__| \033[96m|\033[38;2;0;55;200m|\033[92m|  \033[38;2;0;55;200m`.\033[96m|         |`.        \033[96m/\033[38;2;0;55;200m.\033[96m|   |\033[38;2;0;55;200m| \033[96m|   |\033[38;2;0;55;200m.\033[96m|   |\033[0m
 \033[92m|\033[38;2;0;55;200m|\033[96m| \033[92m|  |   \033[92m(\033[91m===\033[96m)   \033[92m|    \033[96m|\033[38;2;0;55;200m|\033[92m|    \033[96m`---------'  `------'  `---' \033[38;2;0;55;200m`\033[96m'---' `---'\033[0m
\033[91m.---.\033[96m_____\033[92m[:\033[38;2;0;55;200m:\033[96m::::\033[38;2;0;55;200m:\033[92m]\033[96m_____\033[91m.---.\033[92m    \033[92m╔═══════════════════════════════════════════╗\033[0m
\033[91m|\033[38;5;208m[\033[93m_\033[38;5;208m]\033[91m|\033[38;2;0;55;200m------\033[92m|:\033[38;2;0;55;200m:\033[96m::\033[38;2;0;55;200m:\033[92m|\033[38;2;0;55;200m------\033[91m|\033[38;5;208m[\033[93m_\033[38;5;208m]\033[91m|    \033[92m║     \033[38;2;0;55;200mEnthe\033[96mome Ge\033[97mnome Assemb\033[96mly Pip\033[38;2;0;55;200meline     \033[92m║\033[0m
\033[91m`---'\033[92m~~~~~~\033[92m|::\033[38;2;0;55;200m:\033[96m:\033[38;2;0;55;200m:\033[92m|~~~~~~\033[91m`---'    \033[92m╚═══════════════════════════════════════════╝\033[0m

                    Curated & Maintained by Ian M Bollinger
                         (\033[38;2;0;55;200mian.bollinger@entheome.org\033[0m)

                                   \033[92mEGAP.py\033[0m
                                Version {egap.VERSION}

   Preprocess \033[38;2;0;55;200m-\033[92m>\033[0m Assemble \033[38;2;0;55;200m-\033[92m>\033[0m Compare \033[38;2;0;55;200m-\033[92m>\033[0m Polish \033[38;2;0;55;200m-\033[92m>\033[0m Curate \033[38;2;0;55;200m-\033[92m>\033[0m Assess \033[38;2;0;55;200m-\033[92m>\033[0m Report
"""


def fmt_bytes(n: float) -> str:
    for unit in ["B", "KiB", "MiB", "GiB", "TiB"]:
        if n < 1024.0:
            return f"{n:,.1f} {unit}"
        n /= 1024.0
    return f"{n:,.1f} PiB"


@dataclass
class StepStatus:
    sample_id: str
    step: str
    status: str  # PENDING, RUNNING, PASS, FAIL

class Panel(Static):
    def __init__(self, title: str, **kwargs):
        super().__init__("", **kwargs)
        self._title = title
        self._body = ""

    def set_body(self, body: str) -> None:
        self._body = body
        self.update(self.render())

    def render(self) -> str:
        return f"{self._title}\n\n{self._body}" if self._body else f"{self._title}\n"


class ENTHEOME_GENOME_ASSEMBLY_PIPELINE(App):
    CSS = """
    Screen { background: black; color: white; }

    #top_row {
        height: auto;
        margin: 1 2 0 2;
    }

    #banner {
        border: round #0037C8;
        padding: 0 1;
        height: auto;
        width: 84;
    }

    #params {
        border: round #0037C8;
        padding: 1 2;
        height: 26;
        width: 121;
    }

    #content { height: 1fr; margin: 0 2 1 2; }

    #left_col  { width: 80; min-width: 50; }
    #mid_col   { width: 73; min-width: 55; }
    #right_col { width: 52; min-width: 52; }

    #log_widget {
        border: round #0037C8;
        padding: 0 1;
        height: 1fr;
        width: 84;
    }

    #progress_table {
        border: round #0037C8;
        padding: 0 1;
        height: 1fr;
    }

    #cpu_threads {
        border: round #0037C8;
        padding: 1 2;
        height: 1fr;
        width: 52;
    }
    """

    BINDINGS = [
        ("q", "request_quit", "Quit"),
        ("ctrl+q", "request_quit", "Quit"),
        ("f2", "copy_logs", "Copy logs"),
    ]

    def __init__(self):
        super().__init__()

        self.current_process: Optional[asyncio.subprocess.Process] = None
        self.shutdown_requested = False

        # keep a copyable buffer of what you write to RichLog
        self.log_buffer_max_lines = 5000  # adjust as you like
        self.log_buffer: deque[str] = deque(maxlen=self.log_buffer_max_lines)

        # Persistent append log file (opened once output_dir is known in run_pipeline)
        self.log_file = None

        self.params_auto_scroll = True
        self.params_scroll_dir = 1  # 1 = down, -1 = up

        self.pipeline_task: Optional[asyncio.Task] = None

        self.steps: List[StepStatus] = []
        self.step_index: Dict[Tuple[str, str], int] = {}
        self.row_key_map: Dict[Tuple[str, str], str] = {}

        self.pipeline_start_time: str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        self.net_last = None
        self.net_last_t = None

        self.args = None

        self.COL_SAMPLE = 0
        self.COL_STEP = 1

    def auto_scroll_params(self) -> None:
        if not getattr(self, "params_auto_scroll", False):
            return
    
        try:
            params_container = self.query_one("#params", VerticalScroll)
        except Exception:
            return
    
        # Step size: 1 line per tick (adjust interval below for speed)
        next_y = params_container.scroll_y + self.params_scroll_dir
    
        # Bounce at ends (so it doesn’t just stop at bottom)
        max_y = params_container.max_scroll_y
        if next_y >= max_y:
            next_y = max_y
            self.params_scroll_dir = -1
        elif next_y <= 0:
            next_y = 0
            self.params_scroll_dir = 1
    
        params_container.scroll_to(y=next_y, animate=False)

    async def shutdown_pipeline(self):
        self.shutdown_requested = True
        self.log_line("Shutdown requested — terminating pipeline...")

        if self.pipeline_task and not self.pipeline_task.done():
            self.pipeline_task.cancel()

        if self.current_process is not None:
            try:
                # Send SIGTERM to the entire process group (MaSuRCA spawns children)
                pgid = os.getpgid(self.current_process.pid)
                os.killpg(pgid, 15)
                self.log_line("Sent SIGTERM to process group.")

                await asyncio.sleep(2)

                # If still around, force kill
                try:
                    os.killpg(pgid, 0)  # check still exists
                    os.killpg(pgid, 9)
                    self.log_line("Sent SIGKILL to process group.")
                except Exception:
                    pass

            except Exception as e:
                self.log_line(f"Shutdown error: {e}")

    async def action_request_quit(self):
        await self.shutdown_pipeline()
        await self.action_quit()

    async def on_unmount(self) -> None:
        # Safety: if the UI is closed by any means, attempt to stop pipeline
        if not self.shutdown_requested:
            await self.shutdown_pipeline()
        # Close the persistent log file
        if self.log_file is not None:
            try:
                self.log_file.close()
            except OSError:
                pass
            self.log_file = None

    async def action_copy_logs(self) -> None:
        text = "\n".join(self.log_buffer)

        if not text.strip():
            # Optional: give a small visual cue in the log
            self.log_widget.write("Nothing to copy (log buffer is empty).")
            return

        # Textual clipboard helper (terminal support varies; macOS Terminal is a known limitation)
        self.copy_to_clipboard(text)
        self.log_widget.write(f"Copied {len(self.log_buffer)} line(s) to clipboard.")

    def compose(self) -> ComposeResult:
        yield Header(show_clock=True)

        banner_text = Text.from_ansi(BANNER_ANSI.replace("\r\n", "\n"))

        with Horizontal(id="top_row"):
            self.banner = Static(banner_text, id="banner")
            yield self.banner

            with VerticalScroll(id="params"):
                self.params_text = Static("", id="params_text")
                yield self.params_text

        with Horizontal(id="content"):
            with Vertical(id="left_col"):
                self.log_widget = RichLog(id="log_widget", highlight=False, markup=False, wrap=True)
                yield self.log_widget

            with Vertical(id="mid_col"):
                self.progress_table = DataTable(id="progress_table", zebra_stripes=True)
                yield self.progress_table

            with Vertical(id="right_col"):
                self.cpu_threads_widget = CpuThreadsHistoryWidget(refresh_s=0.5, history_len=10, id="cpu_threads")
                yield self.cpu_threads_widget

        yield Footer()

    def on_mount(self) -> None:
        psutil.cpu_percent(interval=None)
        self.pipeline_task = asyncio.create_task(self.run_pipeline())

    def log_line(self, line: str) -> None:
        clean = line.rstrip("\n")
        self.log_buffer.append(clean)
        self.log_widget.write(clean)
        if self.log_file is not None:
            try:
                self.log_file.write(clean + "\n")
                self.log_file.flush()
            except OSError:
                pass

    def format_settings_block(self, title: str, items: Dict) -> Text:
        t = Text()
        # EGAP-like: section headers green
        t.append(f"{title}\n", style="bold bright_green")
    
        for k, v in items.items():
            if isinstance(v, dict):
                # Sub-section label
                t.append(f"  {k}\n", style="rgb(0,55,200)")
                for kk, vv in v.items():
                    t.append(f"    {str(kk):<34}", style="rgb(0,55,200)")
                    t.append(": ", style="bright_white")
                    t.append(f"{vv}\n", style="bright_white")
            else:
                t.append(f"  {str(k):<34}", style="rgb(0,55,200)")
                t.append(": ", style="bright_white")
                t.append(f"{v}\n", style="bright_white")
    
        return t
    
    def format_all_pipeline_settings(self, settings: Dict) -> Text:
        t = Text()
        # EGAP-like: red header bars, cyan dividers
        t.append("=" * 80 + "\n\n", style="bold bright_red")
    
        for section_name, section_data in settings.items():
            if isinstance(section_data, dict):
                t += self.format_settings_block(section_name, section_data)
            else:
                t.append(f"{section_name}\n", style="bold bright_green")
                t.append(f"  {section_data}\n", style="bright_white")
    
            t.append("\n", style="bright_white")
            t.append("-" * 71 + "\n\n", style="bright_cyan")
    
        t.append("=" * 80 + "\n", style="bold bright_red")
        return t

    def status_cell(self, status: str) -> Text:
        if status == "PASS":
            return Text("PASS", style="green")
        if status == "FAIL":
            return Text("FAIL", style="red")
        if status == "RUNNING":
            return Text("RUNNING", style="yellow")
        return Text("PENDING", style="dim")

    def init_step_plan(self, samples: List[Tuple[str, str]], processes: Optional[List[str]] = None,
                       sample_processes_map: Optional[Dict[str, List[str]]] = None) -> None:
        # Caller passes the (possibly skip-filtered) process list so the display
        # table matches what actually runs; fall back to the full list.
        # *sample_processes_map* overrides per-sample (e.g. QC-only samples) so
        # the plan only lists the steps each sample will actually run.
        if processes is None:
            processes = [
                "preprocess_refseq", "preprocess_illumina", "preprocess_ont",
                "preprocess_pacbio", "decontaminate_reads",
                "assemble_masurca", "assemble_flye",
                "assemble_spades", "assemble_hifiasm", "compare_assemblies",
                "polish_assembly", "curate_assembly", "decontaminate_assembly",
            ]
        sample_processes_map = sample_processes_map or {}

        self.steps = []
        for sample_id in samples:
            for p in sample_processes_map.get(sample_id, processes):
                self.steps.append(StepStatus(sample_id=sample_id, step=p, status="PENDING"))
            self.steps.append(StepStatus(sample_id=sample_id, step="qc_assessment(final)", status="PENDING"))
            self.steps.append(StepStatus(sample_id=sample_id, step="html_reporter", status="PENDING"))
    
        self.step_index = {(s.sample_id, s.step): i for i, s in enumerate(self.steps)}
    
        self.progress_table.clear(columns=True)
    
        # Use explicit keys so update_cell can target columns reliably
        self.progress_table.add_column("Sample", key="sample")
        self.progress_table.add_column("Step", key="step")
        self.progress_table.add_column("Status", key="status")
    
        self.row_key_map.clear()
        for s in self.steps:
            row_key = f"{s.sample_id}::{s.step}"
            self.row_key_map[(s.sample_id, s.step)] = row_key
            self.progress_table.add_row(
                s.sample_id,
                s.step,
                self.status_cell(s.status),
                key=row_key,
            )
    
        self.progress_table.cursor_type = "row"

    def update_table_row(self, sample_id: str, step: str) -> None:
        idx = self.step_index.get((sample_id, step))
        if idx is None:
            return
        s = self.steps[idx]
    
        row_key = self.row_key_map.get((sample_id, step))
        if not row_key:
            return
    
        try:
            self.progress_table.update_cell(row_key, "status", self.status_cell(s.status))
        except Exception as e:
            self.log_line(f"EXCEPTION updating table row: {sample_id}::{step} -> {e}")

    def set_step_status(
        self,
        sample_id: str,
        step: str,
        status: str,
        return_code: Optional[int] = None,
        started_at: Optional[str] = None,
        ended_at: Optional[str] = None,
    ) -> None:
        idx = self.step_index.get((sample_id, step))
        if idx is None:
            return
        s = self.steps[idx]
        s.status = status
        if return_code is not None:
            s.return_code = return_code
        if started_at is not None:
            s.started_at = started_at
        if ended_at is not None:
            s.ended_at = ended_at
        self.update_table_row(sample_id, step)

    def set_params_body(self) -> None:
        a = self.args
        current_moment = self.pipeline_start_time
    
        # Pull settings FROM EGAP.py (single source of truth)
        settings = egap.get_pipeline_settings(
            current_moment=current_moment,
            ram_gb=a.ram_gb,
            cpu_threads=a.cpu_threads,
            input_tsv=a.input_tsv,
            output_dir=a.output_dir,
        )
    
        params_renderable = self.format_all_pipeline_settings(settings)
    
        # This assumes you replaced the params panel with a VerticalScroll + Static
        self.params_text.update(params_renderable)

    # Lines containing any of these substrings are suppressed from the TUI log.
    # They are non-fatal noise from numpy C-extension API version mismatches
    # that occur when packages in the same env were compiled against different
    # numpy minor versions (e.g. tiara requires numpy 1.23 while other packages
    # were compiled against 1.24).  The mismatch is harmless in practice.
    _NOISE_SUBSTRINGS = (
        "module compiled against API version",
        "RuntimeError: module compiled against",
    )

    async def run_subprocess_stream(self, cmd: List[str], cwd: Optional[Path] = None) -> int:
        self.log_line(f"\n→ Running: {' '.join(cmd)}\n")

        # Inherit the current environment but force unbuffered Python output so
        # every print() in the child script appears in the TUI immediately rather
        # than being held in an 8 KB pipe buffer until the process exits.
        child_env = os.environ.copy()
        child_env["PYTHONUNBUFFERED"] = "1"

        try:
            proc = await asyncio.create_subprocess_exec(
                *cmd,
                cwd=str(cwd) if cwd else None,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.STDOUT,
                env=child_env,
                preexec_fn=os.setsid,  # critical for killing the whole tree
            )
        except Exception as e:
            self.log_line(f"EXCEPTION launching subprocess: {e}")
            return 127

        self.current_process = proc

        assert proc.stdout is not None
        while True:
            line = await proc.stdout.readline()
            if not line:
                break
            decoded = line.decode(errors="replace").rstrip("\n")
            # Suppress known non-fatal numpy API mismatch noise
            if not any(ns in decoded for ns in self._NOISE_SUBSTRINGS):
                self.log_line(decoded)

        rc = await proc.wait()
        self.current_process = None
        return rc

    async def run_pipeline(self) -> None:
        try:
            # ----------------------------------------------------------------
            # Argument resolution: two paths
            #   1. Launched via `EGAP.py --tui` → EGAP.py attaches its already-
            #      parsed argparse.Namespace as self._preloaded_args so we never
            #      re-parse sys.argv (which still contains --tui, --dry_run, etc.
            #      that this inner parser doesn't know about).
            #   2. Standalone launch (`python EGAP_TUI.py -tsv ...`) → parse
            #      sys.argv the normal way.
            # ----------------------------------------------------------------
            preloaded = getattr(self, "_preloaded_args", None)
            if preloaded is not None:
                self.args = preloaded
            else:
                parser = argparse.ArgumentParser(
                    description=f"Run EGAP with TUI (v{getattr(egap, 'VERSION', 'unknown')})"
                )
                parser.add_argument("--input_tsv", "-tsv", type=str, required=True)
                parser.add_argument("--output_dir", "-o", type=str, required=True)
                parser.add_argument("--cpu_threads", "-t", type=int, default=1)
                parser.add_argument("--ram_gb", "-r", type=int, default=8)
                parser.add_argument("--dry_run", action="store_true", default=False,
                                    help="Log file-management actions without executing them. "
                                         "Equivalent to EGAP_DRY_RUN=1.")
                # Per-assembler skip flags (mirror EGAP.py so standalone TUI
                # launches accept them too).
                parser.add_argument("--no-masurca", "-no_m", dest="no_masurca",
                                    action="store_true", default=False, help="Skip MaSuRCA.")
                parser.add_argument("--no-flye", "-no_f", dest="no_flye",
                                    action="store_true", default=False, help="Skip Flye.")
                parser.add_argument("--no-spades", "-no_s", dest="no_spades",
                                    action="store_true", default=False, help="Skip SPAdes.")
                parser.add_argument("--no-hifiasm", "-no_h", dest="no_hifiasm",
                                    action="store_true", default=False, help="Skip hifiasm.")
                self.args = parser.parse_args()

            # Propagate dry_run env var regardless of how args were obtained.
            if getattr(self.args, "dry_run", False):
                os.environ["EGAP_DRY_RUN"] = "1"

            # Per-sample log files are opened inside the sample loop below.

            self.set_params_body()
            settings = egap.get_pipeline_settings(
                current_moment=self.pipeline_start_time,
                ram_gb=self.args.ram_gb,
                cpu_threads=self.args.cpu_threads,
                input_tsv=self.args.input_tsv,
                output_dir=self.args.output_dir,
            )
            self.log_line("\nFULL PIPELINE SETTINGS:")
            for section, kv in settings.items():
                self.log_line(f"\n[{section}]")
                for k, v in kv.items():
                    self.log_line(f"  {k}: {v}")

            input_tsv = self.args.input_tsv
            output_dir = self.args.output_dir
            cpu_threads = self.args.cpu_threads
            ram_gb = self.args.ram_gb

            # _PROJECT_DIR and _TUI_DIR are resolved at module load time above.
            # egap.__file__ points to EGAP.py in the project root.
            project_dir = Path(egap.__file__).resolve().parent

            processes = [
                "preprocess_refseq", "preprocess_illumina", "preprocess_ont",
                "preprocess_pacbio", "decontaminate_reads",
                "assemble_masurca", "assemble_flye",
                "assemble_spades", "assemble_hifiasm", "compare_assemblies",
                "polish_assembly", "curate_assembly", "decontaminate_assembly",
            ]

            # Honour per-assembler skip flags (--no-masurca/-no_m, etc.), whether
            # they arrived via EGAP.py --tui (preloaded args) or a standalone launch.
            _skip_assemblers = {
                "assemble_masurca": getattr(self.args, "no_masurca", False),
                "assemble_flye":    getattr(self.args, "no_flye", False),
                "assemble_spades":  getattr(self.args, "no_spades", False),
                "assemble_hifiasm": getattr(self.args, "no_hifiasm", False),
            }
            _disabled = [proc for proc, skip in _skip_assemblers.items() if skip]
            if _disabled:
                processes = [p for p in processes if p not in _disabled]
                self.log_line("NOTE:\tSkipping assembler(s) by request: "
                              + ", ".join(p.replace("assemble_", "") for p in _disabled))
            if all(_skip_assemblers.values()):
                self.log_line("WARN:\tAll assemblers were disabled; no assembly will be produced.")

            bin_dir_candidate = egap.locate_bin_dir(processes, project_dir)
            if bin_dir_candidate is not None:
                bin_dir = bin_dir_candidate
            elif (project_dir / "bin").is_dir():
                bin_dir = project_dir / "bin"
            elif _TUI_DIR.name == "bin" and all(((_TUI_DIR / f"{p}.py").exists()) for p in processes):
                # EGAP_TUI.py is itself inside bin/ -- use that directory directly
                bin_dir = _TUI_DIR
            else:
                raise FileNotFoundError("Could not locate bin directory containing EGAP step scripts.")

            # ---- Fail-fast environment preflight (tools on PATH, RAM sanity) ----
            try:
                from preflight_checks import run_preflight
                _pf_problems = run_preflight(ram_gb, abort=False)
            except Exception as _pf_exc:
                self.log_line(f"WARN:\tPre-flight checks could not run ({_pf_exc}); continuing.")
                _pf_problems = []
            if _pf_problems:
                for _p in _pf_problems:
                    self.log_line("ERROR:\t" + _p)
                self.log_line("ERROR:\tAborting run; fix the above and relaunch.")
                return

            input_df = egap.preprocess_tsv(input_tsv)

            samples: List[Tuple[str, str]] = []
            for _, r in input_df.iterrows():
                sample_id = str(r["SAMPLE_ID"]).strip()
                samples.append((sample_id))

            # QC-only samples (reference, no reads, no EST_SIZE) run only the
            # reference fetch; QC + report run afterward. Build a per-sample
            # process map so both the execution loop and the step-plan display
            # reflect what actually runs.
            sample_proc_map = {}
            for _sid in samples:
                try:
                    _row = input_df[input_df["SAMPLE_ID"].astype(str).str.strip() == _sid].iloc[0]
                    _qc_only = egap.sample_is_qc_only(_row)
                except Exception:
                    _qc_only = False
                sample_proc_map[_sid] = (
                    [p for p in processes if p == "preprocess_refseq"] if _qc_only else processes
                )

            self.log_line(f"Loaded {len(samples)} sample(s) from TSV.")
            self.init_step_plan(samples, processes, sample_proc_map)

            for sample_id in samples:
                if self.shutdown_requested:
                    return

                # Open a per-sample log file and keep it active for this sample.
                _sample_log_path = Path(output_dir) / f"{sample_id}_log.txt"
                try:
                    _sample_log_path.parent.mkdir(parents=True, exist_ok=True)
                    self.log_file = open(_sample_log_path, "a", buffering=1)
                except OSError as e:
                    self.log_line(f"WARN:\tCould not open sample log file: {e}")
                    self.log_file = None

                self.log_line(f"\nNOTE:\t{'='*60}")
                self.log_line(f"NOTE:\tBeginning pipeline for sample: {sample_id}")
                self.log_line(f"NOTE:\tSample log: {_sample_log_path}")
                self.log_line(f"NOTE:\t{'='*60}")

                sample_step_failed = False

                sample_processes = sample_proc_map.get(sample_id, processes)
                if sample_processes is not processes:
                    self.log_line("NOTE:\tQC-only sample (reference, no reads, no EST_SIZE): "
                                  "running reference fetch + QC + report only.")

                for proc_name in sample_processes:
                    if self.shutdown_requested:
                        return

                    step_label = proc_name
                    script = bin_dir / f"{proc_name}.py"

                    started = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    self.set_step_status(sample_id, step_label, "RUNNING", started_at=started)

                    if not script.exists():
                        ended = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                        self.set_step_status(sample_id, step_label, "FAIL", return_code=127, ended_at=ended)
                        self.log_line(f"\nERROR:\tMissing script: {script}")
                        sample_step_failed = True
                        break

                    cmd = [
                        sys.executable,
                        str(script),
                        sample_id,
                        input_tsv,
                        output_dir,
                        str(cpu_threads),
                        str(ram_gb),
                    ]

                    rc = await self.run_subprocess_stream(cmd)
                    ended = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

                    if self.shutdown_requested:
                        self.set_step_status(sample_id, step_label, "FAIL", return_code=130, ended_at=ended)
                        return

                    if rc != 0:
                        self.set_step_status(sample_id, step_label, "FAIL", return_code=rc, ended_at=ended)
                        self.log_line(f"\nERROR:\t{proc_name} failed for {sample_id} (rc={rc})")
                        sample_step_failed = True
                        break

                    self.set_step_status(sample_id, step_label, "PASS", return_code=rc, ended_at=ended)

                if sample_step_failed:
                    self.log_line(f"\nWARN:\tOne or more steps failed for {sample_id}; "
                                  f"still running final QC and HTML report.")

                # ---- Final QC assessment (always runs per sample) ----
                qc_step_label = "qc_assessment(final)"
                qc_script = bin_dir / "qc_assessment.py"
                if not qc_script.exists():
                    qc_script = project_dir / "qc_assessment.py"

                started = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                self.set_step_status(sample_id, qc_step_label, "RUNNING", started_at=started)

                if not qc_script.exists():
                    ended = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    self.set_step_status(sample_id, qc_step_label, "FAIL", return_code=127, ended_at=ended)
                    self.log_line(f"ERROR:\tMissing script: {qc_script}")
                else:
                    qc_cmd = [
                        sys.executable,
                        str(qc_script),
                        "final",
                        input_tsv,
                        sample_id,
                        output_dir,
                        str(cpu_threads),
                        str(ram_gb),
                    ]
                    rc = await self.run_subprocess_stream(qc_cmd)
                    ended = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    if rc != 0:
                        self.set_step_status(sample_id, qc_step_label, "FAIL", return_code=rc, ended_at=ended)
                        self.log_line(f"WARN:\tqc_assessment returned non-zero for {sample_id} (rc={rc})")
                    else:
                        self.set_step_status(sample_id, qc_step_label, "PASS", return_code=rc, ended_at=ended)

                # ---- HTML report (always runs per sample) ----
                html_step_label = "html_reporter"
                html_script = bin_dir / "html_reporter.py"
                if not html_script.exists():
                    html_script = project_dir / "html_reporter.py"

                started = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                self.set_step_status(sample_id, html_step_label, "RUNNING", started_at=started)

                if not html_script.exists():
                    ended = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    self.set_step_status(sample_id, html_step_label, "FAIL", return_code=127, ended_at=ended)
                    self.log_line(f"ERROR:\tMissing script: {html_script}")
                else:
                    html_cmd = [
                        sys.executable,
                        str(html_script),
                        sample_id,
                        input_tsv,
                        output_dir,
                        str(cpu_threads),
                        str(ram_gb),
                    ]
                    rc = await self.run_subprocess_stream(html_cmd)
                    ended = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    if rc != 0:
                        self.set_step_status(sample_id, html_step_label, "FAIL", return_code=rc, ended_at=ended)
                        self.log_line(f"WARN:\thtml_reporter returned non-zero for {sample_id} (rc={rc})")
                    else:
                        self.set_step_status(sample_id, html_step_label, "PASS", return_code=rc, ended_at=ended)

                if sample_step_failed:
                    self.log_line(f"\nFAIL:\tSample {sample_id} completed with errors. "
                                  f"Continuing to next sample.")
                else:
                    self.log_line(f"\nPASS:\tSample {sample_id} completed successfully.")

                # Close this sample's log before moving to the next sample.
                if self.log_file is not None:
                    try:
                        self.log_file.close()
                    except OSError:
                        pass
                    self.log_file = None

            self.log_line("\nPASS:\tAll samples processed.")
        except Exception:
            self.log_line("\nFATAL:\tpipeline task crashed with exception:")
            self.log_line(traceback.format_exc())
            return


if __name__ == "__main__":
    ENTHEOME_GENOME_ASSEMBLY_PIPELINE().run()