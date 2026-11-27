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


TEMPLATES = ("EGAP_summary.html", "EGAP_busco.html",
             "GCA_metadata_acc.xlsx", "SRA_metadata_acc.xlsx")


def test_every_stage_has_a_script_beside_egap(egap):
    assert egap.BIN_DIR == REPO_ROOT / "bin"
    for stage in egap.PROCESSES + egap.FINAL_STEPS + ["final_compress"]:
        assert (egap.BIN_DIR / f"{stage}.py").is_file(), stage


def test_templates_are_bundled_and_overridable(monkeypatch, tmp_path):
    from file_operations import resolve_templates_dir

    monkeypatch.delenv("EGAP_TEMPLATES_DIR", raising=False)
    bundled = resolve_templates_dir()
    assert bundled == REPO_ROOT / "resources" / "templates"
    for name in TEMPLATES:
        assert (bundled / name).is_file(), name

    monkeypatch.setenv("EGAP_TEMPLATES_DIR", str(tmp_path))
    assert resolve_templates_dir() == tmp_path
    monkeypatch.setenv("EGAP_TEMPLATES_DIR", str(tmp_path / "missing"))
    with pytest.raises(FileNotFoundError):
        resolve_templates_dir()


def _usable_bash():
    """Path to a bash that understands this platform's paths, or None.

    On Windows ``bash`` may resolve to the WSL launcher in System32, which
    cannot take ``C:/`` paths; Git Bash can.
    """
    import shutil
    bash = shutil.which("bash")
    if bash and "system32" in bash.lower():
        return None
    return bash


@pytest.mark.skipif(_usable_bash() is None, reason="needs a bash that takes native paths")
def test_conda_build_script_installs_one_runnable_tree(tmp_path):
    """build.sh into a scratch PREFIX: the tree, the templates and a working launcher."""
    import os
    import subprocess
    prefix = tmp_path / "prefix"
    env = dict(os.environ, PREFIX=prefix.as_posix())
    subprocess.run([_usable_bash(), "build.sh"], cwd=REPO_ROOT, env=env, check=True,
                   capture_output=True)

    share = prefix / "share" / "egap"
    assert (share / "EGAP.py").is_file()
    assert (share / "bin" / "html_reporter.py").is_file()
    assert (share / "resources" / "TruSeq3-PE.fa").is_file()
    for name in TEMPLATES:
        assert (share / "resources" / "templates" / name).is_file(), name

    for command in ("EGAP", "EGAP_TUI"):
        assert (prefix / "bin" / command).is_file()
    done = subprocess.run([sys.executable, str(prefix / "bin" / "EGAP"), "--help"],
                          capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert "--input_tsv" in done.stdout


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
