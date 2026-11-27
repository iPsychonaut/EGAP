"""bin/qc_assessment.py: BUSCO output paths for BUSCO 5.x (odb12) and 6.1 (odb12.2)."""
import pytest

import qc_assessment as qc

TAG = "asm_enterobacterales_busco"
LINEAGE = "enterobacterales_odb12"


@pytest.mark.parametrize("dataset", ["enterobacterales_odb12", "enterobacterales_odb12.2"])
def test_finds_summary_and_table_for_both_busco_naming_schemes(tmp_path, dataset):
    summary = tmp_path / f"short_summary.specific.{dataset}.{TAG}.txt"
    table = tmp_path / f"run_{dataset}" / "full_table.tsv"
    table.parent.mkdir()
    summary.write_text("C:98.5%[S:98.5%,D:0.0%],F:0.4%,M:1.1%,n:542\n")
    table.write_text("")
    assert qc.busco_output_paths(str(tmp_path), LINEAGE, TAG) == (str(summary), str(table))


def test_falls_back_to_unversioned_names_before_busco_runs(tmp_path):
    summary, table = qc.busco_output_paths(str(tmp_path), LINEAGE, TAG)
    assert summary.endswith(f"short_summary.specific.{LINEAGE}.{TAG}.txt")
    assert table.endswith(f"run_{LINEAGE}/full_table.tsv".replace("/", __import__("os").sep))


# ---- --busco_only / EGAP_BUSCO_ONLY -------------------------------------------------

@pytest.fixture
def tool_calls(monkeypatch):
    """Replace both completeness tools with recorders; BUSCO reports metrics."""
    calls = []

    def fake_compleasm(_asm, _sample, stats, count, *_rest):
        calls.append("compleasm")
        if fake_compleasm.works:
            stats[f"{count.upper()}_BUSCO_C"] = 99.0

    def fake_busco(_asm, _sample, stats, count, *_rest):
        calls.append("busco")
        stats[f"{count.upper()}_BUSCO_C"] = 98.0

    fake_compleasm.works = True
    monkeypatch.setattr(qc, "compleasm_assembly", fake_compleasm)
    monkeypatch.setattr(qc, "busco_assembly", fake_busco)
    monkeypatch.setattr(qc, "_COMPLETENESS_FAILURES", [])
    monkeypatch.delenv("EGAP_BUSCO_ONLY", raising=False)
    return calls, fake_compleasm


def _evaluate():
    stats = {}
    return qc.run_lineage_eval("asm.fasta", "S1", stats, "first", "fungi", "final", 1), stats


def test_default_runs_compleasm_and_leaves_busco_alone(tool_calls):
    calls, _ = tool_calls
    measured, stats = _evaluate()
    assert measured and calls == ["compleasm"] and stats["FIRST_BUSCO_C"] == 99.0
    assert qc.completeness_tool_order() == ("compleasm", "busco")


def test_default_falls_back_to_busco_when_compleasm_gives_nothing(tool_calls):
    calls, compleasm = tool_calls
    compleasm.works = False
    measured, stats = _evaluate()
    assert measured and calls == ["compleasm", "busco"] and stats["FIRST_BUSCO_C"] == 98.0


def test_busco_only_never_calls_compleasm(tool_calls, monkeypatch):
    calls, _ = tool_calls
    monkeypatch.setenv("EGAP_BUSCO_ONLY", "1")
    measured, stats = _evaluate()
    assert measured and calls == ["busco"] and stats["FIRST_BUSCO_C"] == 98.0
    assert qc.completeness_tool_order() == ("busco", "compleasm")


def test_preflight_does_not_require_compleasm_under_busco_only(monkeypatch):
    import preflight_checks
    monkeypatch.setattr(preflight_checks.shutil, "which",
                        lambda tool: None if tool == "compleasm" else f"/bin/{tool}")
    monkeypatch.setattr(preflight_checks, "_actual_ram_gb", lambda: None)

    monkeypatch.delenv("EGAP_BUSCO_ONLY", raising=False)
    problems = preflight_checks.run_preflight(8, abort=False)
    assert problems and "compleasm" in problems[0]

    monkeypatch.setenv("EGAP_BUSCO_ONLY", "1")
    assert not preflight_checks.run_preflight(8, abort=False)
