"""Stage timing records, the per-sample log, and ``EGAP --estimate_only``."""
import importlib.util
import json
import subprocess
import sys
from datetime import datetime, timedelta
from pathlib import Path

import pytest

import log
import record_provenance as provenance

REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def egap():
    spec = importlib.util.spec_from_file_location("EGAP", REPO_ROOT / "EGAP.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _entries(tmp_path, sample_id="S1"):
    return json.loads((tmp_path / f"{sample_id}_provenance.json").read_text())


def test_only_the_latest_run_is_reported(tmp_path):
    first, second = "2026-11-01 09:00:00", "2026-11-02 10:00:00"
    provenance.record_stage(tmp_path, "S1", "preprocess_ont", first, "2026-11-01 09:00:05", 3725.4, 0)
    provenance.record_stage(tmp_path, "S1", "assemble_flye", first, "2026-11-01 10:02:10", 60, 1)
    provenance.record_stage(tmp_path, "S1", "preprocess_ont", second, "2026-11-02 10:00:05", 2.2, 0)
    provenance.record_stage(tmp_path, "S1", "qc_assessment", second, "2026-11-02 10:00:08", 90, 0)

    timings = provenance.stage_timings(_entries(tmp_path))
    assert [(t["stage"], t["duration"], t["returncode"]) for t in timings] == [
        ("preprocess_ont", "0:00:02", 0), ("qc_assessment", "0:01:30", 0)]
    assert provenance.stage_timings([]) == []
    assert provenance.stage_timings([{"type": "command", "stage": "x", "tool": "flye"}]) == []


def test_stage_records_are_not_mistaken_for_commands(tmp_path):
    provenance.record_stage(tmp_path, "S1", "compare_assemblies", "run", "2026-11-01 09:00:00", 5, 0)
    entries = _entries(tmp_path) + [{"type": "command", "stage": "assemble_flye", "tool": "flye",
                                     "command": "flye --version", "version": "2.9.6"}]
    assert provenance.meaningful_tools(entries) == ["flye"]
    assert provenance.versions_from_entries(entries) == {"flye": "2.9.6"}
    # A stage that ran no external command still appears in the step list.
    assert provenance.pipeline_steps(entries) == ["compare_assemblies", "assemble_flye"]


def test_finish_stage_reports_and_records(egap, tmp_path):
    said = []
    started = datetime.now() - timedelta(seconds=3725)
    elapsed = egap.finish_stage(tmp_path, "S1", "polish_assembly", "2026-11-01 09:00:00",
                                started, 0, say=said.append)
    assert elapsed in ("1:02:05", "1:02:06")
    assert said == [f"NOTE:\tpolish_assembly took {elapsed} (exit code 0)"]
    (entry,) = _entries(tmp_path)
    assert entry["type"] == "stage" and entry["stage"] == "polish_assembly"
    assert entry["run"] == "2026-11-01 09:00:00" and entry["returncode"] == 0
    assert entry["seconds"] == pytest.approx(3725, abs=2)


def test_starting_a_stage_keeps_the_sample_log(tmp_path, monkeypatch, capsys):
    """Every stage calls this on start-up; it used to empty the log each time."""
    monkeypatch.setattr(log, "DEFAULT_LOG_FILE", None)
    sample_log = tmp_path / "S1_log.txt"
    sample_log.write_text("earlier stage output\n")
    log.initialize_logging_environment(str(tmp_path), "S1")
    assert sample_log.read_text() == "earlier stage output\n"
    assert Path(log.DEFAULT_LOG_FILE) == sample_log.resolve()


def test_estimate_only_prints_a_table_and_touches_nothing(sample_tsv_factory, tmp_path):
    table = sample_tsv_factory(SPECIES_ID="Escherichia_coli", SAMPLE_ID="Escherichia_coli-Est",
                               EST_SIZE="5m", ONT_RAW_DIR="/data/ont_runs")
    out_dir = tmp_path / "never_created"
    done = subprocess.run([sys.executable, str(REPO_ROOT / "EGAP.py"), "--estimate_only",
                           "-tsv", str(table), "-o", str(out_dir), "-t", "8", "-r", "32"],
                          capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert "RUNTIME ESTIMATE (heuristic) - sample Escherichia_coli-Est" in done.stdout
    assert "8 threads, 32 GB RAM" in done.stdout
    # A folder of raw reads counts as long reads, so Flye gets an estimate.
    assert "no long (ONT/PacBio) reads" not in done.stdout
    assert not out_dir.exists()
