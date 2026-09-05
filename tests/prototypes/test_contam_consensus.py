import pytest

import mockdata
from audit.prototypes import contam_consensus as cc


@pytest.mark.parametrize("div,dom", [
    ("fung:basidiomycetes", "eukarya"), ("prok:g-proteobacteria", "bacteria"),
    ("arch:euryarchaeota", "archaea"), ("plnt:monocots", "eukarya"),
    ("anml:insects", "eukarya"), ("unassigned", "unknown"), ("", "unknown"), (None, "unknown")])
def test_fcs_div_mapping(div, dom):
    assert cc.fcs_div_to_domain(div) == dom


def test_segment_id_parsing():
    assert cc.split_segment_id("contig_1~1-5000") == ("contig_1", (1, 5000))
    assert cc.split_segment_id("contig_1~5001..9000") == ("contig_1", (5001, 9000))
    assert cc.split_segment_id("contig_1") == ("contig_1", None)


def _fixture(tmp_path):
    fcs = mockdata.write_fcsgx_rpt(tmp_path / "s.taxonomy.rpt", [
        {"seq_id": "contig_1", "seq_len": 20000, "taxid": 5341, "tax_name": "Psilocybe cubensis", "div": "fung:basidiomycetes"},
        {"seq_id": "contig_2", "seq_len": 15000, "taxid": 5341, "tax_name": "Psilocybe cubensis", "div": "fung:basidiomycetes"},
        {"seq_id": "contam_1", "seq_len": 3000, "taxid": 562, "tax_name": "Escherichia coli", "div": "prok:g-proteobacteria"},
        # chimera: two windows, two domains
        {"seq_id": "contig_3", "seq_len": 6000, "taxid": 5341, "tax_name": "Psilocybe cubensis", "div": "fung:basidiomycetes", "segment": (1, 6000)},
        {"seq_id": "contig_3", "seq_len": 4000, "taxid": 1280, "tax_name": "Staphylococcus aureus", "div": "prok:firmicutes", "segment": (6001, 10000)},
        # same-domain chimera: two species, both fungal
        {"seq_id": "contig_4", "seq_len": 5000, "taxid": 5341, "tax_name": "Psilocybe cubensis", "div": "fung:basidiomycetes", "segment": (1, 5000)},
        {"seq_id": "contig_4", "seq_len": 5000, "taxid": 5334, "tax_name": "Agaricus bisporus", "div": "fung:basidiomycetes", "segment": (5001, 10000)},
        # short contig
        {"seq_id": "contig_5", "seq_len": 600, "taxid": 562, "tax_name": "Escherichia coli", "div": "prok:g-proteobacteria"},
        # disagreement: FCS says bacteria, Tiara says eukarya
        {"seq_id": "contig_6", "seq_len": 8000, "taxid": 562, "tax_name": "Escherichia coli", "div": "prok:g-proteobacteria"},
        # FCS unknown
        {"seq_id": "contig_7", "seq_len": 8000, "taxid": 0, "tax_name": "unknown", "div": "unassigned"},
    ])
    tiara = mockdata.write_tiara(tmp_path / "tiara.txt", {
        "contig_1": "eukarya", "contig_2": "organelle", "contam_1": "bacteria",
        "contig_3": "eukarya", "contig_4": "eukarya", "contig_5": "bacteria",
        "contig_6": "eukarya", "contig_7": "unknown",
        "contig_8": "prokarya",   # Tiara-only contig (no FCS row)
    })
    return fcs, tiara


def test_collapse_detects_chimeras(tmp_path):
    fcs, _ = _fixture(tmp_path)
    col = cc.collapse_fcsgx(cc.parse_fcsgx_rpt(fcs))
    assert col["contig_3"]["species"] == "possible cross-domain chimera"
    assert col["contig_3"]["domain"] == "mixed"
    assert col["contig_3"]["seq_len"] == 10000
    assert col["contig_3"]["n_segments"] == 2
    assert col["contig_4"]["species"] == "possible chimera"
    assert col["contig_4"]["domain"] == "eukarya"
    assert col["contig_1"]["chimera"] is False
    assert col["contig_7"]["taxid"] == cc.UNKNOWN_TAXID


def test_compare_mirrors_fsp_tags(tmp_path):
    fcs, tiara = _fixture(tmp_path)
    cmp = cc.compare(cc.collapse_fcsgx(cc.parse_fcsgx_rpt(fcs)), cc.parse_tiara(tiara))
    assert cmp["contig_1"]["match"] == "match" and cmp["contig_1"]["blob_tag"] == "Psilocybe_cubensis"
    assert cmp["contig_2"]["match"] == "match"          # organelle counts as eukarya
    assert cmp["contam_1"]["match"] == "match" and cmp["contam_1"]["blob_tag"] == "Escherichia_coli"
    assert cmp["contig_3"]["blob_tag"] == "possible_chimera"
    assert cmp["contig_5"]["blob_tag"] == "unknown"     # < 1 kb override
    assert cmp["contig_6"]["match"] == "mismatch" and cmp["contig_6"]["blob_tag"] == "eukarya_bacteria"
    assert "contig_8" in cmp and cmp["contig_8"]["fcs_domain"] == "unknown"


def test_decision_rule_funga(tmp_path):
    fcs, tiara = _fixture(tmp_path)
    cmp = cc.compare(cc.collapse_fcsgx(cc.parse_fcsgx_rpt(fcs)), cc.parse_tiara(tiara))
    d = cc.partition("Funga", cmp)
    assert set(d["remove"]) == {"contam_1", "contig_5"}     # both classifiers say bacteria
    assert set(d["flag"]) == {"contig_3", "contig_4", "contig_6", "contig_8"}
    assert set(d["keep"]) == {"contig_1", "contig_2", "contig_7"}


def test_short_contig_with_agreement_is_removed_not_flagged():
    """The 1 kb override only affects the blob tag; the decision layer still
    sees both raw domain calls. Two agreeing bacterial calls on a fungal
    target remove the contig whatever its length."""
    assert cc.decide("Funga", "bacteria", "bacteria", False) == "remove"


@pytest.mark.parametrize("kingdom,t,f,chim,expect", [
    ("Funga", "eukarya", "eukarya", False, "keep"),
    ("Funga", "eukarya", "unknown", False, "keep"),
    ("Funga", "unknown", "unknown", False, "keep"),
    ("Funga", "bacteria", "bacteria", False, "remove"),
    ("Funga", "prokarya", "archaea", False, "remove"),
    ("Funga", "bacteria", "archaea", False, "flag"),      # both non-target but disagree
    ("Funga", "eukarya", "bacteria", False, "flag"),
    ("Funga", "bacteria", "unknown", False, "flag"),      # single-classifier evidence
    ("Funga", "eukarya", "eukarya", True, "flag"),
    ("Bacteria", "bacteria", "bacteria", False, "keep"),
    ("Bacteria", "eukarya", "eukarya", False, "remove"),
    ("Bacteria", "prokarya", "bacteria", False, "flag"),  # prokarya is ambiguous for a bacterial target
    ("Archaea", "archaea", "archaea", False, "keep"),
    ("", "eukarya", "eukarya", False, "keep"),            # blank kingdom -> eukaryote default
])
def test_decision_matrix(kingdom, t, f, chim, expect):
    assert cc.decide(kingdom, t, f, chim) == expect


def test_tiara_only_evidence_never_removes():
    """This is the behavioural change versus EGAP today, where Tiara alone
    removes. With FCS-GX absent (unknown) the contig is flagged, not removed."""
    assert cc.decide("Funga", "bacteria", "unknown", False) == "flag"


def test_review_table_written(tmp_path):
    fcs, tiara = _fixture(tmp_path)
    cmp = cc.compare(cc.collapse_fcsgx(cc.parse_fcsgx_rpt(fcs)), cc.parse_tiara(tiara))
    d = cc.partition("Funga", cmp)
    p = cc.write_review_table(tmp_path / "review.tsv", cmp, d)
    lines = p.read_text().splitlines()
    assert lines[0].split("\t")[-1] == "decision"
    assert len(lines) == 1 + len(cmp)


def test_prototype_tiara_parser_skips_header(tmp_path):
    t = mockdata.write_tiara(tmp_path / "t.txt", {"c1": "eukarya"})
    assert cc.parse_tiara(t) == {"c1": "eukarya"}
