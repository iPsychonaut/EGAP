"""EGAP.py top-level helpers that do not launch subprocesses."""
import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def egap():
    """Import EGAP.py as a module without triggering the __main__ block."""
    spec = importlib.util.spec_from_file_location("EGAP", REPO_ROOT / "EGAP.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["EGAP"] = mod
    spec.loader.exec_module(mod)
    return mod


def _row(**kw):
    return pd.Series(kw)


def test_qc_only_when_ref_and_no_reads_no_est_size(egap):
    assert egap.sample_is_qc_only(_row(REF_SEQ_GCA="GCA_000005845.2", EST_SIZE=None,
                                       ILLUMINA_SRA=None))


def test_qc_only_via_nt_assembly_path(egap):
    assert egap.sample_is_qc_only(_row(NT_ASSEMBLY_PATH="/x/asm.fasta"))


def test_not_qc_only_when_reads_present(egap):
    assert not egap.sample_is_qc_only(_row(REF_SEQ_GCA="GCA_1", ILLUMINA_SRA="SRR1"))


def test_not_qc_only_when_est_size_present(egap):
    assert not egap.sample_is_qc_only(_row(REF_SEQ_GCA="GCA_1", EST_SIZE="5m"))


def test_blank_and_none_strings_do_not_count_as_values(egap):
    assert not egap._cell_has_value(_row(X=""), "X")
    assert not egap._cell_has_value(_row(X="None"), "X")
    assert not egap._cell_has_value(_row(X="nan"), "X")
    assert egap._cell_has_value(_row(X="GCA_1"), "X")


def test_version_string_matches_recipe(egap):
    """Guard against the 3.4.1 / 3.4.2 drift found in the root-file audit."""
    meta = (REPO_ROOT / "meta.yaml").read_text()
    assert f'version = "{egap.VERSION}"' in meta


def test_locate_bin_dir_prefers_repo_bin_over_hidden_worktree(egap, tmp_path):
    scripts = ["a", "b"]
    for d in (tmp_path / ".claude" / "worktrees" / "x" / "bin", tmp_path / "bin"):
        d.mkdir(parents=True)
        for s in scripts:
            (d / f"{s}.py").touch()
    assert egap.locate_bin_dir(scripts, tmp_path) == tmp_path / "bin"


def test_locate_bin_dir_walk_skips_hidden_dirs(egap, tmp_path):
    scripts = ["a", "b"]
    for d in (tmp_path / ".claude" / "worktrees" / "x" / "bin", tmp_path / "src" / "scripts"):
        d.mkdir(parents=True)
        for s in scripts:
            (d / f"{s}.py").touch()
    assert egap.locate_bin_dir(scripts, tmp_path) == tmp_path / "src" / "scripts"
    assert egap.locate_bin_dir(scripts, tmp_path / ".claude" / "nope") is None


def test_flags_stage_list_and_stage_commands(egap, tmp_path):
    args = egap.build_parser().parse_args(["-tsv", "s.tsv", "-o", "out", "-no_f", "--no-hifiasm"])
    said = []
    procs = egap.resolve_processes(args, say=said.append)
    assert [p for p in egap.PROCESSES if p not in procs] == ["assemble_flye", "assemble_hifiasm"]
    assert "flye, hifiasm" in said[0] and len(said) == 1
    settings = egap.get_pipeline_settings("now", 8, 1, "s.tsv", "out", args)["input-setup settings"]
    assert settings["assemblers skipped"] == "flye, hifiasm" and settings["dry run"] in ("True", "False")

    tail = ["out", "4", "16"]
    assert egap.stage_cmd(tmp_path, "polish_assembly", "S1", "s.tsv", "out", 4, 16)[2:] == ["S1", "s.tsv"] + tail
    assert egap.stage_cmd(tmp_path, "qc_assessment", "S1", "s.tsv", "out", 4, 16)[2:] == ["final", "s.tsv", "S1"] + tail


def test_nextflow_runs_every_egap_stage_in_order(egap):
    """draft_assembly.nf must not lag EGAP.py: same stage scripts, same order, same opt-outs."""
    import re
    nf = (REPO_ROOT / "draft_assembly.nf").read_text(encoding="utf-8")
    stages = egap.PROCESSES + egap.FINAL_STEPS + ["final_compress"]
    assert re.findall(r"^process (\w+) \{", nf, flags=re.M) == ["prepare"] + stages
    assert re.findall(r"stage\('(\w+)', sample_id\)", nf) == stages
    workflow = nf[nf.index("\nworkflow {"):]
    called = re.findall(r"\b(" + "|".join(stages) + r")\b", workflow)
    assert list(dict.fromkeys(called)) == stages
    for name in egap.ASSEMBLERS:
        assert f"params.no_{name} " in workflow
