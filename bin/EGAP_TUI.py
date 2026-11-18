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

Orchestration comes from EGAP.py (parser, stage list, stage commands,
preprocess_tsv); this file only adds the display and the async runner.


python EGAP_TUI.py -tsv /mnt/c/Users/theda/OneDrive/Desktop/EGAP/EGAP_TUI_test.tsv -o /mnt/c/Users/theda/OneDrive/Desktop/EGAP/EGAP_TUI_test -t 10 -r 24

Created on Sun Feb 22 14:56:27 2026

Updated on 2026-04-16

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""

from __future__ import annotations

import traceback
import asyncio
import os
import sys
import psutil
from collections import deque
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
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
    blocks: deque[Tuple[str, str]]  # (char, style)


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

    def init_step_plan(self, sample_proc_map: Dict[str, List[str]]) -> None:
        """Lay out one table row per step each sample will actually run."""
        self.steps = [
            StepStatus(sample_id=sample_id, step=p, status="PENDING")
            for sample_id, procs in sample_proc_map.items()
            for p in procs + egap.FINAL_STEPS + ["final_compress"]
        ]

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
            args=a,
        )
    
        params_renderable = self.format_all_pipeline_settings(settings)
    
        # This assumes you replaced the params panel with a VerticalScroll + Static
        self.params_text.update(params_renderable)

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
            if not any(ns in decoded for ns in egap._NOISE_SUBSTRINGS):
                self.log_line(decoded)

        rc = await proc.wait()
        self.current_process = None
        return rc

    async def run_pipeline(self) -> None:
        try:
            # Launched via `EGAP.py --tui`: EGAP.py attaches its already-parsed
            # argparse.Namespace as self._preloaded_args. Standalone launch
            # (`python EGAP_TUI.py -tsv ...`): parse sys.argv with the same parser.
            self.args = getattr(self, "_preloaded_args", None) or egap.build_parser().parse_args()
            egap.apply_run_flags(self.args)
            input_tsv = self.args.input_tsv
            output_dir = self.args.output_dir
            cpu_threads = self.args.cpu_threads
            ram_gb = self.args.ram_gb

            # Per-sample log files are opened inside the sample loop below.
            self.set_params_body()

            settings = egap.get_pipeline_settings(self.pipeline_start_time, ram_gb, cpu_threads,
                                                  input_tsv, output_dir, self.args)
            self.log_line("\nFULL PIPELINE SETTINGS:")
            for section, kv in settings.items():
                self.log_line(f"\n[{section}]")
                for k, v in kv.items():
                    self.log_line(f"  {k}: {v}")

            processes = egap.resolve_processes(self.args, say=self.log_line)
            bin_dir = egap.resolve_bin_dir(processes, Path(egap.__file__).resolve().parent)

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

            # QC-only samples (reference, no reads, no EST_SIZE) run only the
            # reference fetch; QC + report run afterward. The per-sample map
            # drives both the execution loop and the step-plan display.
            sample_proc_map = {
                str(r["SAMPLE_ID"]).strip(): ["preprocess_refseq"] if egap.sample_is_qc_only(r) else processes
                for _, r in input_df.iterrows()
            }
            self.log_line(f"Loaded {len(sample_proc_map)} sample(s) from TSV.")
            self.init_step_plan(sample_proc_map)

            def now() -> str:
                return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

            for sample_id, sample_processes in sample_proc_map.items():
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
                if sample_processes is not processes:
                    self.log_line("NOTE:\tQC-only sample (reference, no reads, no EST_SIZE): "
                                  "running reference fetch + QC + report only.")

                # Pipeline steps, then final QC and HTML report (those two always
                # run), then final_compress (only when nothing failed).
                sample_step_failed = False
                for proc in sample_processes + egap.FINAL_STEPS + ["final_compress"]:
                    if self.shutdown_requested:
                        return
                    final = proc in egap.FINAL_STEPS
                    if sample_step_failed and not final:
                        if proc == "final_compress":
                            self.log_line(f"\nSKIP:\tfinal_compress for {sample_id}: an earlier step failed.")
                        continue
                    if proc == "final_compress" and self.args.dry_run:
                        self.log_line(f"\nDRY-RUN:\tWould run final_compress for {sample_id}.")
                        continue

                    self.set_step_status(sample_id, proc, "RUNNING", started_at=now())
                    rc = await self.run_subprocess_stream(egap.stage_cmd(
                        bin_dir, proc, sample_id, input_tsv, output_dir, cpu_threads, ram_gb))
                    if self.shutdown_requested:
                        self.set_step_status(sample_id, proc, "FAIL", return_code=130, ended_at=now())
                        return
                    self.set_step_status(sample_id, proc, "PASS" if rc == 0 else "FAIL",
                                         return_code=rc, ended_at=now())
                    if rc == 0:
                        continue
                    if proc == "final_compress":
                        self.log_line(f"\nWARN:\tfinal_compress returned non-zero for {sample_id}; "
                                      f"files left uncompressed.")
                    elif final:
                        sample_step_failed = True
                        self.log_line(f"WARN:\t{proc} returned non-zero for {sample_id} (rc={rc})")
                    else:
                        sample_step_failed = True
                        self.log_line(f"\nERROR:\t{proc} failed for {sample_id} (rc={rc})")
                        self.log_line(f"\nWARN:\tOne or more steps failed for {sample_id}; "
                                      f"still running final QC and HTML report.")

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