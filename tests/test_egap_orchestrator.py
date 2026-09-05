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


def test_find_files_is_an_infinite_loop_on_any_match(tmp_path):
    """EGAP.find_files rebinds its accumulator to the os.walk file list and then
    appends to that same list while iterating it, so the first match never
    terminates. Run in a subprocess with a hard timeout to prove it without
    hanging the suite. No caller exists in the repo: delete the function."""
    import subprocess
    (tmp_path / "a_dedup.fastq").write_text("x")
    code = (
        "import importlib.util, sys;"
        f"spec=importlib.util.spec_from_file_location('EGAP', r'{REPO_ROOT / 'EGAP.py'}');"
        "m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m);"
        f"print(len(m.find_files('_dedup', folder=r'{tmp_path}')))"
    )
    with pytest.raises(subprocess.TimeoutExpired):
        subprocess.run([sys.executable, "-c", code], timeout=5, capture_output=True)
