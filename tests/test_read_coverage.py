"""Read-set totals and coverage, read back from FastQC and NanoPlot outputs at final QC."""
import zipfile

import pandas as pd
import pytest

import sample_tsv

ROW_COLUMNS = ["SAMPLE_ID", "SPECIES_ID", "ONT_SRA", "ONT_RAW_READS", "ILLUMINA_SRA",
               "ILLUMINA_RAW_F_READS", "ILLUMINA_RAW_R_READS", "PACBIO_SRA",
               "PACBIO_RAW_READS", "REF_SEQ_GCA", "REF_SEQ"]


def fastqc_zip(folder, name, total_bases):
    """Write a minimal FastQC report zip; *total_bases* None mimics FastQC before 0.12."""
    folder.mkdir(parents=True, exist_ok=True)
    data = "##FastQC\t0.13.0\n>>Basic Statistics\tpass\nTotal Sequences\t10\n"
    if total_bases:
        data += f"Total Bases\t{total_bases}\n"
    path = folder / f"{name}_fastqc.zip"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(f"{name}_fastqc/fastqc_data.txt", data + ">>END_MODULE\n")
    return path


@pytest.mark.parametrize("text, expected", [
    ("950 bp", 950),
    ("12.5 kbp", 12_500),
    ("396.4 Mbp", 396_400_000),
    ("1.2 Gbp", 1_200_000_000),
])
def test_fastqc_total_bases_reads_each_unit(tmp_path, text, expected):
    assert sample_tsv.fastqc_total_bases(fastqc_zip(tmp_path, "r_1", text)) == pytest.approx(expected)


def test_fastqc_total_bases_is_none_when_absent_or_unreadable(tmp_path):
    assert sample_tsv.fastqc_total_bases(fastqc_zip(tmp_path, "old", None)) is None
    not_a_zip = tmp_path / "broken_fastqc.zip"
    not_a_zip.write_text("not a zip")
    assert sample_tsv.fastqc_total_bases(not_a_zip) is None
    assert sample_tsv.fastqc_total_bases(tmp_path / "missing_fastqc.zip") is None


def test_totals_are_read_from_disk_and_turned_into_coverage(tmp_path):
    species = tmp_path / "Escherichia_coli"
    illumina = species / "Illumina"
    for folder, forward, reverse in (("fastqc_results", "3 Mbp", "3 Mbp"),
                                     ("trimmed_fastqc_results", "2.5 Mbp", "2.5 Mbp")):
        fastqc_zip(illumina / folder, "reads_1", forward)
        fastqc_zip(illumina / folder, "reads_2", reverse)
    fastqc_zip(illumina / "dedup_fastqc_results", "only_forward", "2 Mbp")  # reverse missing
    nanoplot = species / "ONT" / "Raw_ONT_nanoplot_analysis"
    nanoplot.mkdir(parents=True)
    (nanoplot / "Raw_ONT_NanoStats.txt").write_text(
        "General summary:\nMean read length:   5,000.0\nMean read quality:     12.3\n"
        "Number of reads:    1,000.0\nTotal bases:    5,000,000.0\n")

    stats = sample_tsv.gen_sample_stats_dict(pd.Series({column: None for column in ROW_COLUMNS}))
    assert not any("HIFI" in key for key in stats)

    sample_tsv.load_read_stats(str(species), stats)
    assert stats["RAW_ILLU_TOTAL_BASES"] == pytest.approx(6_000_000)
    assert stats["TRIMMED_ILLU_TOTAL_BASES"] == pytest.approx(5_000_000)
    assert stats["DEDUPED_ILLU_TOTAL_BASES"] is None      # a lone report is not a pair
    assert stats["RAW_ONT_TOTAL_BASES"] == 5_000_000 and stats["RAW_ONT_READS"] == 1_000
    assert stats["FILT_ONT_TOTAL_BASES"] is None          # no NanoPlot output on disk

    assert sample_tsv.add_read_coverage(stats, 1_000_000) == 3
    assert stats["RAW_ILLU_COVERAGE"] == 6.0
    assert stats["TRIMMED_ILLU_COVERAGE"] == 5.0
    assert stats["RAW_ONT_COVERAGE"] == 5.0
    assert stats["DEDUPED_ILLU_COVERAGE"] is None


@pytest.mark.parametrize("genome_bp", [None, float("nan"), pd.NA, 0])
def test_coverage_is_left_empty_without_a_genome_size(genome_bp):
    stats = {"RAW_ONT_TOTAL_BASES": 5_000_000, "RAW_ONT_COVERAGE": None}
    assert sample_tsv.add_read_coverage(stats, genome_bp) == 0
    assert stats["RAW_ONT_COVERAGE"] is None


def test_report_coverage_table_and_fastqc_folders(tmp_path):
    """html_reporter needs jinja2 and bs4; CI installs them, a bare local Python may not."""
    for module in ("jinja2", "bs4", "geopy", "requests"):
        pytest.importorskip(module)
    import html_reporter

    stats_csv = tmp_path / "S1_final_EGAP_assembly_final_stats.csv"
    pd.DataFrame([{"RAW_ILLU_TOTAL_BASES": 6_000_000.0, "RAW_ILLU_COVERAGE": 6.0,
                   "TRIMMED_ILLU_TOTAL_BASES": None, "TRIMMED_ILLU_COVERAGE": None,
                   "RAW_ONT_TOTAL_BASES": 5_000_000.0, "RAW_ONT_COVERAGE": 5.25,
                   "COVERAGE_BASIS": "assembly length, 1,000,000 bp"}]).to_csv(stats_csv, index=False)
    rows, basis = html_reporter.coverage_rows(str(stats_csv))
    assert rows == [
        {"read_set": "Illumina, raw", "total_bases": "6,000,000", "coverage": "6x"},
        {"read_set": "ONT, raw", "total_bases": "5,000,000", "coverage": "5.25x"},
    ]
    assert basis == "assembly length, 1,000,000 bp"
    assert html_reporter.coverage_rows(str(tmp_path / "missing.csv")) == ([], "")
    # every stats prefix the table lists is one the final QC can fill
    known = set(sample_tsv.FASTQC_FOLDERS) | {"RAW_ONT", "FILT_ONT", "CORRECT_ONT",
                                              "RAW_PACBIO", "FILT_PACBIO"}
    assert set(html_reporter.COVERAGE_LABELS) == known
