"""bin/compare_assemblies.py: ruling out assemblies that drop a replicon others agree on."""
import pytest

import compare_assemblies as ca


@pytest.mark.parametrize("intervals,length,expect", [
    ([], 100, 0.0),
    ([(0, 50)], 100, 0.5),
    ([(0, 40), (30, 60)], 100, 0.6),        # overlap counted once
    ([(70, 90), (0, 10), (5, 20)], 100, 0.4),
    ([(0, 100), (10, 20)], 100, 1.0),
])
def test_aligned_fraction(intervals, length, expect):
    assert ca.aligned_fraction(intervals, length) == pytest.approx(expect)


# Contigs of each assembly, and how much of each aligns to each other assembly:
# coverage[(target, query)][query contig] = fraction aligned to target.
CONTIGS = {"s": {"chr": 5_000_000, "p13": 13_000, "p160": 160_000, "junk": 8_000},
           "m": {"chr": 5_000_000, "p70": 70_000},
           "f": {"chr": 5_000_000, "p13x2": 27_000, "p99": 99_000}}
COVERAGE = {
    ("m", "s"): {"chr": 1.0, "p13": 0.05, "p160": 0.91, "junk": 0.0},
    ("f", "s"): {"chr": 1.0, "p13": 1.0, "p160": 0.95, "junk": 0.0},
    ("s", "m"): {"chr": 1.0, "p70": 1.0},
    ("f", "m"): {"chr": 1.0, "p70": 1.0},
    ("s", "f"): {"chr": 1.0, "p13x2": 0.5, "p99": 1.0},
    ("m", "f"): {"chr": 1.0, "p13x2": 0.09, "p99": 0.97},
}


@pytest.fixture
def assemblies(tmp_path, monkeypatch):
    files = {}
    for method, short in (("SPAdes", "s"), ("MaSuRCA", "m"), ("Flye", "f")):
        p = tmp_path / f"{short}.fasta"
        p.write_text("")
        files[method] = str(p)
    short_of = {v: v.rsplit("/", 1)[-1].rsplit("\\", 1)[-1][0] for v in files.values()}
    monkeypatch.setattr(ca, "contig_coverage", lambda t, q, *_: COVERAGE[(short_of[t], short_of[q])])
    monkeypatch.setattr(ca.SeqIO, "parse", lambda path, fmt: [
        type("R", (), {"id": c, "seq": "A" * n})() for c, n in CONTIGS[short_of[path]].items()])
    return files


def flagged(missing):
    return {m: sorted(c for _, c, _, _ in items) for m, items in missing.items()}


def test_with_reference_half_aligned_counts_as_present(assemblies, tmp_path):
    # MaSuRCA lacks the 13 kb plasmid that SPAdes and Flye share; junk that only
    # SPAdes made, and p160 that MaSuRCA holds 91% of, do not count.
    missing = ca.missing_shared_replicons(assemblies, 1, str(tmp_path / "w"))
    assert flagged(missing) == {"MaSuRCA": ["p13", "p13x2"]}


def test_without_reference_a_candidate_must_hold_95_percent(assemblies, tmp_path):
    missing = ca.missing_shared_replicons(assemblies, 1, str(tmp_path / "w"),
                                          ca.MIN_KEPT_FRACTION_NO_REF)
    assert flagged(missing) == {"MaSuRCA": ["p13", "p13x2", "p160"]}
