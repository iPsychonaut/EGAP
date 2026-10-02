"""bin/sample_tsv.py: TSV reading, collapsed-column expansion, row lookup."""
import pandas as pd
import pytest

import sample_tsv


def test_read_sample_table_expands_busco_pair(sample_tsv_factory):
    tsv = sample_tsv_factory(SPECIES_ID="Sp", SAMPLE_ID="Sp-1", ORGANISM_KINGDOM="Funga",
                             ORGANISM_KARYOTE="eukaryote", EST_SIZE="50m",
                             BUSCO="agaricales,basidiomycota")
    df = sample_tsv.read_sample_table(tsv)
    row = df.iloc[0]
    assert row["BUSCO_1"] == "agaricales"
    assert row["BUSCO_2"] == "basidiomycota"
    assert "BUSCO" not in df.columns


def test_read_sample_table_expands_illumina_pair(sample_tsv_factory):
    tsv = sample_tsv_factory(SPECIES_ID="Sp", SAMPLE_ID="Sp-1", BUSCO="a,b",
                             ILLUMINA_RAW_READS="/x/r_1.fastq,/x/r_2.fastq")
    row = sample_tsv.read_sample_table(tsv).iloc[0]
    assert row["ILLUMINA_RAW_F_READS"] == "/x/r_1.fastq"
    assert row["ILLUMINA_RAW_R_READS"] == "/x/r_2.fastq"


def test_single_busco_token_leaves_second_na(sample_tsv_factory):
    tsv = sample_tsv_factory(SPECIES_ID="Sp", SAMPLE_ID="Sp-1", BUSCO="fungi")
    row = sample_tsv.read_sample_table(tsv).iloc[0]
    assert row["BUSCO_1"] == "fungi"
    assert pd.isna(row["BUSCO_2"])


def test_literal_none_placeholder_treated_as_blank(sample_tsv_factory):
    tsv = sample_tsv_factory(SPECIES_ID="Sp", SAMPLE_ID="Sp-1", BUSCO="None")
    row = sample_tsv.read_sample_table(tsv).iloc[0]
    assert pd.isna(row["BUSCO_1"])


def test_local_fna_ref_seq_is_normalized_to_fasta(sample_tsv_factory, tmp_path):
    """A local .fna REF_SEQ reaches stages as the canonical RefSeq .fasta copy."""
    import preprocess_refseq

    fna = tmp_path / "GCA_000005845.2.fna"
    fna.write_text(">chr\nACGT\n")
    out = tmp_path / "out"
    tsv = sample_tsv_factory(SPECIES_ID="Sp", SAMPLE_ID="Sp-1", BUSCO="a,b", REF_SEQ=str(fna))

    # Before preprocess_refseq runs there is no copy, so the raw path is kept.
    assert sample_tsv.load_sample_context("Sp-1", tsv, out, 1, 1).current_series["REF_SEQ"] == str(fna)

    placed = preprocess_refseq.preprocess_refseq("Sp-1", str(tsv), str(out), 1, 1)
    canonical = out / "Sp" / "RefSeq" / "Sp_RefSeq.fasta"
    assert canonical.read_text() == ">chr\nACGT\n"
    assert placed == str(canonical.resolve())

    ref_seq = sample_tsv.load_sample_context("Sp-1", tsv, out, 1, 1).current_series["REF_SEQ"]
    assert ref_seq.endswith("_RefSeq.fasta") and canonical.samefile(ref_seq)

    # Re-running preprocess_refseq still sees the raw source, not its own copy.
    assert preprocess_refseq.preprocess_refseq("Sp-1", str(tsv), str(out), 1, 1) == placed


def test_get_current_row_data_returns_stats_schema(sample_tsv_factory):
    tsv = sample_tsv_factory(SPECIES_ID="Sp", SAMPLE_ID="Sp-1", BUSCO="a,b",
                             ONT_SRA="SRR000001")
    df = sample_tsv.read_sample_table(tsv)
    row, idx, stats = sample_tsv.get_current_row_data(df, "Sp-1")
    assert len(row) == 1
    # Fields Matthew referenced as present-but-unpopulated.
    assert "KMER_COMPLETENESS" in stats and stats["KMER_COMPLETENESS"] is None
    assert "QUAL_VAL" in stats and stats["QUAL_VAL"] is None


def test_get_current_row_data_unknown_sample_is_empty(sample_tsv_factory):
    tsv = sample_tsv_factory(SPECIES_ID="Sp", SAMPLE_ID="Sp-1", BUSCO="a,b")
    df = sample_tsv.read_sample_table(tsv)
    # Current behaviour: no error here; the empty frame only fails later at
    # ``current_row.iloc[0]`` inside each stage. A fail-fast check belongs here.
    row, idx, stats = sample_tsv.get_current_row_data(df, "does-not-exist")
    assert len(row) == 0 and idx == []


# ---- select_long_reads: which ONT read set the assemblers get ----

def _ont_layout(tmp_path, sample_tsv_factory, corrected=False, filtered=True):
    """<out>/Sp/ONT with raw reads, NanoStats for each stage, and optional read sets."""
    out = tmp_path / "out"
    ont = out / "Sp" / "ONT"
    ont.mkdir(parents=True)
    (ont / "SRRX.fastq").write_text("@raw\nACGT\n+\nIIII\n" * 50)
    if filtered:
        (ont / "Sp_ont_filtered.fastq").write_text("@filtered\nACGT\n+\nIIII\n")
    if corrected:
        (ont / "Sp_ont_corrected.fastq").write_text("@corrected\nACGT\n+\nIIII\n")
    # Ratatosk fell back to the filtered reads, so "Corr" stats equal "Filt" stats.
    for origin, qual in (("Raw_ONT_", 12.0), ("Filt_ONT_", 21.1), ("Corr_ONT_", 21.1)):
        d = ont / f"{origin}nanoplot_analysis"
        d.mkdir()
        (d / f"{origin}NanoStats.txt").write_text(f"Mean read quality:    {qual}\n")
    tsv = sample_tsv_factory(SPECIES_ID="Sp", SAMPLE_ID="Sp-1", BUSCO="a,b", ONT_SRA="SRRX")
    return out, ont, tsv


def test_select_long_reads_uses_filtered_when_corrected_was_never_written(sample_tsv_factory, tmp_path):
    """Ratatosk produced nothing: the selection must not return None (raw-read fallback)."""
    out, ont, tsv = _ont_layout(tmp_path, sample_tsv_factory, corrected=False)
    chosen = sample_tsv.select_long_reads(str(out), str(tsv), "Sp-1", 1)
    assert chosen == str(ont / "Sp_ONT_highest_mean_qual_long_reads.fastq")
    assert open(chosen).read().startswith("@filtered")


def test_select_long_reads_keeps_corrected_when_present(sample_tsv_factory, tmp_path):
    out, ont, tsv = _ont_layout(tmp_path, sample_tsv_factory, corrected=True)
    chosen = sample_tsv.select_long_reads(str(out), str(tsv), "Sp-1", 1)
    assert open(chosen).read().startswith("@corrected")


def test_select_long_reads_returns_none_with_no_candidate(sample_tsv_factory, tmp_path):
    out, ont, tsv = _ont_layout(tmp_path, sample_tsv_factory, corrected=False, filtered=False)
    assert sample_tsv.select_long_reads(str(out), str(tsv), "Sp-1", 1) is None
    assert not (ont / "Sp_ONT_highest_mean_qual_long_reads.fastq").exists()


# ---- flye_ont_mode: Flye read-type mode follows what was selected ----

def test_selection_records_source_and_quality(sample_tsv_factory, tmp_path):
    out, ont, tsv = _ont_layout(tmp_path, sample_tsv_factory, corrected=False)
    chosen = sample_tsv.select_long_reads(str(out), str(tsv), "Sp-1", 1)
    assert sample_tsv.long_read_info(chosen) == {"source": "filtered", "mean_quality": 21.1}

    out2, ont2, tsv2 = _ont_layout(tmp_path / "b", sample_tsv_factory, corrected=True)
    chosen2 = sample_tsv.select_long_reads(str(out2), str(tsv2), "Sp-1", 1)
    assert sample_tsv.long_read_info(chosen2)["source"] == "corrected"


def test_corrected_alias_symlink_is_not_corrected(sample_tsv_factory, tmp_path):
    """Without Illumina reads preprocess_ont symlinks *_corrected to the filtered reads."""
    out, ont, tsv = _ont_layout(tmp_path, sample_tsv_factory, corrected=False)
    try:
        (ont / "Sp_ont_corrected.fastq").symlink_to(ont / "Sp_ont_filtered.fastq")
    except OSError:
        pytest.skip("symlinks not permitted on this platform")
    chosen = sample_tsv.select_long_reads(str(out), str(tsv), "Sp-1", 1)
    assert sample_tsv.long_read_info(chosen)["source"] == "filtered"
    assert sample_tsv.flye_ont_mode(chosen, chosen) == "--nano-hq"


@pytest.mark.parametrize("info,expected", [
    ({"source": "corrected", "mean_quality": 30.0}, "--nano-corr"),
    ({"source": "filtered", "mean_quality": 21.1}, "--nano-hq"),
    ({"source": "filtered", "mean_quality": 13.0}, "--nano-hq"),
    ({"source": "filtered", "mean_quality": 11.0}, "--nano-raw"),
    ({"source": "filtered", "mean_quality": None}, "--nano-raw"),
    (None, "--nano-corr"),                      # selected by an older EGAP: no record
])
def test_flye_ont_mode(tmp_path, info, expected):
    import json
    reads = tmp_path / "Sp_ONT_highest_mean_qual_long_reads.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    if info is not None:
        (tmp_path / (reads.name + sample_tsv.LONG_READ_INFO_SUFFIX)).write_text(json.dumps(info))
    assert sample_tsv.flye_ont_mode(str(reads), str(reads)) == expected


def test_flye_ont_mode_raw_fallback(tmp_path):
    """If Flye is handed anything other than the selected file, assume raw reads."""
    selected = tmp_path / "Sp_ONT_highest_mean_qual_long_reads.fastq"
    raw = tmp_path / "SRRX.fastq"
    raw.write_text("@r\nACGT\n+\nIIII\n")
    assert sample_tsv.flye_ont_mode(str(raw), str(selected)) == "--nano-raw"
