import math

import pytest

import mockdata
from audit.prototypes import merqury_qv as mq


def test_parse_qv_and_completeness(tmp_path):
    mockdata.write_merquryfk(tmp_path / "spades", qv=38.2, completeness=97.5)
    qv = mq.parse_qv(tmp_path / "spades.qv")
    comp = mq.parse_completeness(tmp_path / "spades.completeness.stats")
    assert qv == {"spades": 38.2}
    assert comp == {"spades": 97.5}


def test_qv_consistent_with_error_rate(tmp_path):
    mockdata.write_merquryfk(tmp_path / "a", qv=30.0, completeness=90.0)
    rows = mq._read_table(tmp_path / "a.qv")
    err = float(rows[0]["Error %"]) / 100
    assert math.isclose(mq.qv_from_error_rate(err), 30.0, abs_tol=0.05)


def test_parse_qv_rejects_wrong_layout(tmp_path):
    p = tmp_path / "bad.qv"
    p.write_text("foo\tbar\n1\t2\n")
    with pytest.raises(ValueError):
        mq.parse_qv(p)


def _S(name, b1, b2, contigs, n50, qv=None, kc=None):
    return mq.AssemblyStats(name, b1, b2, contigs, n50, qv, kc)


def test_new_rule_and_egap_vote_agree_on_clear_winner():
    stats = [_S("SPAdes", 80, 78, 900, 20_000, qv=30), _S("MaSuRCA", 95, 94, 300, 80_000, qv=35)]
    assert mq.select_best(stats).name == "MaSuRCA"
    assert mq.egap_vote(stats) == "MaSuRCA"


def test_qv_breaks_busco_tie():
    stats = [_S("SPAdes", 95, 94, 300, 80_000, qv=32.0),
             _S("MaSuRCA", 95, 94, 300, 80_000, qv=36.5)]
    assert mq.select_best(stats).name == "MaSuRCA"
    # EGAP's vote: every metric ties, max() keeps the first -> SPAdes. Arbitrary.
    assert mq.egap_vote(stats) == "SPAdes"


def test_two_two_split_no_longer_depends_on_metric_order():
    """Same 2-2 splits pinned in tests/test_compare_assemblies.py. Under the
    explicit rule mean BUSCO decides first, so all three cases resolve to the
    assembly with the higher mean completeness, regardless of which lineage
    it won."""
    a = [_S("SPAdes", 90, 88, 100, 90_000), _S("MaSuRCA", 96, 95, 400, 30_000)]
    b = [_S("SPAdes", 96, 88, 400, 90_000), _S("MaSuRCA", 90, 95, 100, 30_000)]
    c = [_S("SPAdes", 90, 95, 400, 90_000), _S("MaSuRCA", 96, 88, 100, 30_000)]
    assert mq.select_best(a).name == "MaSuRCA"     # mean 95.5 vs 89.0
    assert mq.select_best(b).name == "MaSuRCA"     # mean 92.0 vs 92.5
    assert mq.select_best(c).name == "SPAdes"      # mean 92.5 vs 92.0
    # EGAP's current vote gives SPAdes for case b (BUSCO_1 winner). The two
    # rules disagree there; that disagreement is the point of making it explicit.
    assert mq.egap_vote(b) == "SPAdes"


def test_missing_qv_never_wins_a_qv_tiebreak():
    stats = [_S("A", 90, 90, 100, 50_000, qv=None), _S("B", 90, 90, 100, 50_000, qv=20.0)]
    assert mq.select_best(stats).name == "B"


def test_deterministic_last_resort_is_alphabetical():
    stats = [_S("Zeta", 90, 90, 100, 50_000), _S("Alpha", 90, 90, 100, 50_000)]
    assert mq.select_best(stats).name == "Zeta"   # max() of key tuple -> 'Zeta' > 'Alpha'


def test_rank_key_treats_none_busco_as_zero():
    stats = [_S("A", None, None, 10, 1_000_000, qv=50), _S("B", 60, 60, 5000, 5_000)]
    assert mq.select_best(stats).name == "B"
