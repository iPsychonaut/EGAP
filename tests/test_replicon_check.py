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


def test_missing_needs_two_assemblers_to_agree(tmp_path, monkeypatch):
    # SPAdes and Flye both carry a plasmid; MaSuRCA lacks it. Flye also has a
    # contig nobody else has, which must not count against anyone.
    names = {"SPAdes": "s", "MaSuRCA": "m", "Flye": "f"}
    contigs = {"s": {"chr": 5_000_000, "p13": 13_000},
               "m": {"chr": 5_000_000},
               "f": {"chr": 5_000_000, "p13x2": 27_000, "junk": 8_000}}
    # absent[(target, query)]: query contigs not found in target
    absent = {("m", "s"): {"p13"}, ("m", "f"): {"p13x2", "junk"},
              ("s", "m"): set(), ("s", "f"): {"junk"},
              ("f", "m"): set(), ("f", "s"): set()}
    files = {}
    for m, short in names.items():
        p = tmp_path / f"{short}.fasta"
        p.write_text("".join(f">{c}\n{'A' * 10}\n" for c in contigs[short]))
        files[m] = str(p)
    path_to_short = {str(tmp_path / f"{s}.fasta"): s for s in names.values()}
    monkeypatch.setattr(ca, "contigs_absent_from",
                        lambda t, q, *_: absent[(path_to_short[t], path_to_short[q])])
    monkeypatch.setattr(ca.SeqIO, "parse", lambda path, fmt: [
        type("R", (), {"id": c, "seq": "A" * n})() for c, n in contigs[path_to_short[path]].items()])

    missing = ca.missing_shared_replicons(files, 1, str(tmp_path / "work"))

    assert set(missing) == {"MaSuRCA"}
    assert sorted((src, c) for src, c, _ in missing["MaSuRCA"]) == [("Flye", "p13x2"), ("SPAdes", "p13")]
