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
