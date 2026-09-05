import pytest

import mockdata
from audit.prototypes import kmer_sizing as ks


@pytest.mark.parametrize("median,expected", [
    (150, 99), (100, 65), (250, 165), (75, 49), (301, 199), (36, 23)])
def test_heuristic_k_is_two_thirds_forced_odd(median, expected):
    assert ks.heuristic_k(median) == expected


def test_force_odd_matches_fsp_awk():
    assert ks.force_odd(100) == 99
    assert ks.force_odd(99) == 99


def test_parse_kmergenie_html(tmp_path):
    html = mockdata.write_kmergenie_html(tmp_path / "report.html", 71)
    assert ks.parse_kmergenie_html(html.read_text()) == 71


def test_parse_kmergenie_html_missing_returns_none():
    assert ks.parse_kmergenie_html("<html>no prediction</html>") is None


@pytest.mark.parametrize("k,ok", [
    (21, True), (22, False), (13, False), (127, True), (129, False), (None, False)])
def test_validate_k_spades_bounds(k, ok):
    assert ks.validate_k(k) is ok


def test_validate_k_rejects_k_at_or_above_read_length():
    assert ks.validate_k(99, read_len=150)
    assert not ks.validate_k(99, read_len=99)
    assert not ks.validate_k(99, read_len=90)


def test_build_k_list_appends_valid_prediction():
    assert ks.build_k_list(ks.SPADES_DEFAULT_K, 99) == [21, 33, 55, 77, 99]


def test_build_k_list_ignores_duplicate_and_even():
    assert ks.build_k_list(ks.SPADES_DEFAULT_K, 55) == ks.SPADES_DEFAULT_K
    assert ks.build_k_list(ks.SPADES_DEFAULT_K, 56) == ks.SPADES_DEFAULT_K


def test_egap_hardcoded_list_survives_only_with_long_reads():
    """EGAP's [21,33,55,77,99] is fine for 150 bp reads. Trimmomatic CROP is
    145 and MINLEN 125 in EGAP settings, so 99 is safe there. It is NOT safe
    if reads are ever shorter than 100 bp (BBDuk floor is 100 bp per commit
    9924fd4), which is exactly the case the data-driven rule guards."""
    assert ks.build_k_list(ks.EGAP_SPADES_K, None, read_len=125) == [21, 33, 55, 77, 99]
    assert ks.build_k_list(ks.EGAP_SPADES_K, None, read_len=100) == [21, 33, 55, 77, 99]
    assert ks.build_k_list(ks.EGAP_SPADES_K, None, read_len=99) == [21, 33, 55, 77]
    assert ks.build_k_list(ks.EGAP_SPADES_K, None, read_len=76) == [21, 33, 55]


def test_median_from_mock_fastq(rng, tmp_path):
    recs = mockdata.make_genome(rng, n_contigs=2, min_len=5000, max_len=6000)
    r1, r2 = mockdata.write_paired_fastq(rng, recs, tmp_path / "r1.fastq.gz",
                                         tmp_path / "r2.fastq", n_pairs=300, read_len=150)
    assert ks.median_read_length(r1) == 150
    assert ks.median_read_length(r2) == 150
    assert ks.heuristic_k(ks.median_read_length(r1)) == 99


def test_median_variable_lengths(rng, tmp_path):
    recs = mockdata.make_genome(rng, n_contigs=2, min_len=5000, max_len=6000)
    r1, _ = mockdata.write_paired_fastq(rng, recs, tmp_path / "v1.fastq", tmp_path / "v2.fastq",
                                        n_pairs=1000, read_len=150, variable=True)
    med = ks.median_read_length(r1)
    assert 100 <= med <= 130          # uniform on [75,150] -> median near 112
    k = ks.heuristic_k(med)
    assert k % 2 == 1 and k < med


def test_k_dirs_derive_from_same_list():
    ks_ = ks.build_k_list(ks.SPADES_DEFAULT_K, 99)
    assert ks.spades_k_arg(ks_) == "21,33,55,77,99"
    assert ks.spades_k_dirs(ks_) == ["K21", "K33", "K55", "K77", "K99"]
