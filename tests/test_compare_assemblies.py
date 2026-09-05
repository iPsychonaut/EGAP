"""bin/compare_assemblies.py: best-assembly vote, with BUSCO/QUAST stubbed.

BUSCO and QUAST are monkeypatched so the vote logic runs on GitHub Actions
with no tools installed. Assembly discovery is also patched so the test
controls the candidate set and its insertion order.
"""
import os
from pathlib import Path

import pytest

import compare_assemblies as ca
import mockdata


def _setup(monkeypatch, tmp_path, sample_tsv_factory, rng, stats):
    """stats: dict assembler -> (busco1, busco2, contigs, n50), in insertion order."""
    tsv = sample_tsv_factory(SPECIES_ID="Sp", SAMPLE_ID="Sp-1", ORGANISM_KINGDOM="Funga",
                             ORGANISM_KARYOTE="eukaryote", EST_SIZE="50m",
                             BUSCO="agaricales,basidiomycota", ILLUMINA_SRA="SRR1")
    out = tmp_path / "out"
    sample_dir = out / "Sp" / "Sp-1"
    sample_dir.mkdir(parents=True)
    fastas = {}
    for name in stats:
        recs = mockdata.make_genome(rng, n_contigs=3, min_len=500, max_len=800)
        fastas[name] = str(mockdata.write_fasta(sample_dir / f"{name.lower()}.fasta", recs))

    # discovery order in compare_assemblies(): SPAdes, MaSuRCA, Flye, Hifiasm
    monkeypatch.setattr(ca, "discover_spades", lambda d, s: fastas.get("SPAdes"))
    monkeypatch.setattr(ca, "discover_masurca", lambda d, s: fastas.get("MaSuRCA"))
    monkeypatch.setattr(ca, "discover_flye", lambda d, s: fastas.get("Flye"))
    monkeypatch.setattr(ca, "discover_hifiasm", lambda d, s: fastas.get("Hifiasm"))

    by_path = {v: k for k, v in fastas.items()}
    monkeypatch.setattr(ca, "get_busco_score",
                        lambda asm, db, t, d: stats[by_path[asm]][0 if db == "agaricales" else 1])
    monkeypatch.setattr(ca, "get_quast_stats",
                        lambda asm, t, d, m: tuple(stats[by_path[asm]][2:4]))
    return tsv, out


def _winner(best_path, out):
    return Path(best_path).name


def test_clear_winner_selected(monkeypatch, tmp_path, sample_tsv_factory, rng):
    cwd = os.getcwd()
    try:
        tsv, out = _setup(monkeypatch, tmp_path, sample_tsv_factory, rng, {
            "SPAdes":  (80.0, 78.0, 900, 20_000),
            "MaSuRCA": (95.0, 94.0, 300, 80_000),
        })
        best = ca.compare_assemblies("Sp-1", str(tsv), str(out), 1, 4)
        assert Path(best).read_text() == Path(out / "Sp" / "Sp-1" / "masurca.fasta").read_text()
    finally:
        os.chdir(cwd)


@pytest.mark.parametrize("spades,masurca,expected", [
    # MaSuRCA wins both BUSCOs, SPAdes wins contigs + N50 -> 2-2 -> MaSuRCA
    ((90.0, 88.0, 100, 90_000), (96.0, 95.0, 400, 30_000), "masurca"),
    # SPAdes wins BUSCO_1 + N50, MaSuRCA wins BUSCO_2 + contigs -> 2-2 -> SPAdes
    ((96.0, 88.0, 400, 90_000), (90.0, 95.0, 100, 30_000), "spades"),
    # MaSuRCA wins BUSCO_1 + contigs, SPAdes wins BUSCO_2 + N50 -> 2-2 -> MaSuRCA
    ((90.0, 95.0, 400, 90_000), (96.0, 88.0, 100, 30_000), "masurca"),
])
def test_two_two_split_resolves_to_first_metric_winner(monkeypatch, tmp_path,
                                                        sample_tsv_factory, rng,
                                                        spades, masurca, expected):
    """With two assemblers and four metrics a 2-2 split is common. The vote is
    ``Counter(custom_stats).most_common(1)``; Counter preserves insertion order
    on ties, and custom_stats is built in METRIC order (BUSCO_1, BUSCO_2,
    contigs, N50), not assembler-discovery order. So the effective, undocumented
    tie-break is: whoever won BUSCO_1, else BUSCO_2, else contig count.

    Matthew's email described this as tie-breaking by insertion order. That is
    the mechanism; the observable rule is "first metric wins". The rule is
    defensible (completeness over contiguity) but implicit. These cases pin it
    so any change to the scoring is a deliberate, visible diff."""
    cwd = os.getcwd()
    try:
        tsv, out = _setup(monkeypatch, tmp_path, sample_tsv_factory, rng, {
            "SPAdes": spades, "MaSuRCA": masurca})
        best = ca.compare_assemblies("Sp-1", str(tsv), str(out), 1, 4)
        chosen = Path(best).read_text()
        assert chosen == (out / "Sp" / "Sp-1" / f"{expected}.fasta").read_text()
    finally:
        os.chdir(cwd)


def test_single_assembly_is_returned_without_vote(monkeypatch, tmp_path,
                                                  sample_tsv_factory, rng):
    cwd = os.getcwd()
    try:
        tsv, out = _setup(monkeypatch, tmp_path, sample_tsv_factory, rng, {
            "Flye": (70.0, 65.0, 50, 400_000),
        })
        best = ca.compare_assemblies("Sp-1", str(tsv), str(out), 1, 4)
        assert Path(best).name == "Sp_best_assembly.fasta"
    finally:
        os.chdir(cwd)


def test_failed_busco_on_one_assembly_does_not_crash(monkeypatch, tmp_path,
                                                     sample_tsv_factory, rng):
    cwd = os.getcwd()
    try:
        tsv, out = _setup(monkeypatch, tmp_path, sample_tsv_factory, rng, {
            "SPAdes":  (None, None, 100, 90_000),
            "MaSuRCA": (96.0, 95.0, 400, 30_000),
        })
        best = ca.compare_assemblies("Sp-1", str(tsv), str(out), 1, 4)
        assert best is not None
    finally:
        os.chdir(cwd)
