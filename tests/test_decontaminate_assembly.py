"""bin/decontaminate_assembly.py: kingdom keep-sets, Tiara parsing, done markers."""
import pandas as pd
import pytest

import decontaminate_assembly as da
import mockdata


@pytest.mark.parametrize("kingdom,expect_keep", [
    ("Funga", {"eukarya", "organelle", "unknown"}),
    ("flora", {"eukarya", "organelle", "unknown"}),
    ("FAUNA", {"eukarya", "organelle", "unknown"}),
    ("Bacteria", {"bacteria", "prokarya", "unknown"}),
    ("archaea", {"archaea", "prokarya", "unknown"}),
])
def test_keep_sets_by_kingdom(kingdom, expect_keep):
    keep, remove = da.get_decontamination_classes(kingdom, "eukaryote")
    assert keep == expect_keep
    assert keep | remove == da.ALL_TIARA_CLASSES
    assert not (keep & remove)


@pytest.mark.parametrize("bad", [None, "", "   ", pd.NA, "Plantae"])
def test_unknown_kingdom_falls_back_to_eukaryote(bad):
    keep, _ = da.get_decontamination_classes(bad, None)
    assert keep == {"eukarya", "organelle", "unknown"}


def test_unknown_always_kept_for_every_kingdom():
    for k in ("Funga", "Bacteria", "Archaea", None):
        keep, _ = da.get_decontamination_classes(k, None)
        assert "unknown" in keep


def test_parse_tiara_reads_first_stage_label(tmp_path):
    out = mockdata.write_tiara(tmp_path / "tiara.txt", {
        "contig_1": "eukarya", "contig_2": "Bacteria", "contam_1": "organelle"})
    parsed = da.parse_tiara(str(out))
    assert parsed["contig_1"] == "eukarya"
    assert parsed["contig_2"] == "bacteria"      # lower-cased
    assert parsed["contam_1"] == "organelle"


def test_parse_tiara_ingests_header_row_as_a_sequence(tmp_path):
    """Tiara writes a header line (sequence_id / class_fst_stage / ...). The
    parser has no header guard, so a phantom entry keyed 'sequence_id' with
    class 'class_fst_stage' appears in the mapping. Harmless today because
    that key never matches a FASTA record, but it inflates any count derived
    from the dict. Pin the behaviour; fix by skipping the header."""
    out = mockdata.write_tiara(tmp_path / "tiara.txt", {"contig_1": "eukarya"})
    parsed = da.parse_tiara(str(out))
    assert parsed.get("sequence_id") == "class_fst_stage"


def test_parse_tiara_skips_blank_and_comment_lines(tmp_path):
    p = tmp_path / "t.txt"
    p.write_text("# comment\n\ncontig_1\teukarya\tn/a\nbroken_line_no_tab\n")
    assert da.parse_tiara(str(p)) == {"contig_1": "eukarya"}


def test_done_marker_records_counts(tmp_path):
    m = tmp_path / "decontamination_done.txt"
    da.write_done_marker(str(m), "/x/asm.fasta", kept_count=7, removed_count=2,
                         kingdom="Funga", keep_classes={"eukarya"}, remove_classes={"bacteria"})
    text = m.read_text()
    assert "Sequences kept   : 7" in text
    assert "Sequences removed: 2" in text


def test_done_marker_does_not_verify_outputs(tmp_path):
    """Documents the checkpoint-integrity gap Matthew raised: the marker can be
    written even when the output it vouches for does not exist. The prototype
    in audit/prototypes/checkpoint.py is the proposed replacement."""
    m = tmp_path / "decontamination_done.txt"
    da.write_done_marker(str(m), "/does/not/exist.fasta", 0, 0)
    assert m.exists()  # nothing stopped it
