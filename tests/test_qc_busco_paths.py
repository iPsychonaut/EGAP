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
