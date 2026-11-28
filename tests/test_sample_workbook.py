"""The column specification, the sample workbook, and the NCBI metadata stage."""
import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

import sample_tsv
import sample_workbook

REPO_ROOT = Path(__file__).resolve().parents[1]
SAMPLE_NAMES = [c.name for c in sample_tsv.SAMPLE_COLUMNS if c.sheet == sample_tsv.SAMPLE_SHEET]
METADATA_NAMES = [c.name for c in sample_tsv.SAMPLE_COLUMNS if c.sheet == sample_tsv.METADATA_SHEET]


@pytest.fixture(scope="module")
def egap():
    spec = importlib.util.spec_from_file_location("EGAP", REPO_ROOT / "EGAP.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ---- column specification (no extra packages needed) --------------------------------

def test_every_column_says_what_is_required_and_in_what_format():
    for column in sample_tsv.SAMPLE_COLUMNS:
        assert column.required and column.format and column.help, column.name
        tooltip = sample_workbook.tooltip_text(column)
        assert tooltip.startswith(f"Required: {column.required}\nFormat: {column.format}\n")
        if column.choices:
            assert "Allowed: " + ", ".join(column.choices) in tooltip
    assert len(set(SAMPLE_NAMES)) == len(SAMPLE_NAMES)
    assert len(set(METADATA_NAMES)) == len(METADATA_NAMES)
    # Read-type prefixes mean "this sample has reads" to the Nextflow QC-only test.
    assert not [n for n in METADATA_NAMES if n.startswith(("ILLUMINA_", "ONT_", "PACBIO_"))]


def test_stage_columns_are_derived_from_the_specification():
    stage_columns = sample_tsv._LEGACY_COLUMNS
    for name in ("SAMPLE_ID", "BUSCO_1", "BUSCO_2", "ILLUMINA_RAW_F_READS", "ILLUMINA_RAW_R_READS",
                 "NT_ASSEMBLY_PATH", "KRAKEN2_DB", "SAMPLE_TISSUE", "INSTRUMENT_ONT"):
        assert name in stage_columns
    assert "BUSCO" not in stage_columns and "ILLUMINA_RAW_READS" not in stage_columns
    assert len(set(stage_columns)) == len(stage_columns)


def test_example_table_uses_only_known_columns():
    header = (REPO_ROOT / "resources" / "EGAP_test.tsv").read_text().splitlines()[0].split("\t")
    assert [name for name in header if name not in SAMPLE_NAMES] == []


def test_a_table_without_metadata_columns_still_reads(sample_tsv_factory):
    table = sample_tsv.read_sample_table(sample_tsv_factory(SAMPLE_ID="S1", BUSCO="a,b"))
    assert pd.isna(table.loc[0, "SAMPLE_TISSUE"]) and table.loc[0, "BUSCO_2"] == "b"


def test_command_line_requires_a_table_unless_writing_a_workbook(egap, capsys):
    with pytest.raises(SystemExit):
        egap.parse_args(["-o", "out"])
    assert egap.parse_args(["--new_workbook", "book.xlsx"]).new_workbook == "book.xlsx"
    args = egap.parse_args(["--input_workbook", "book.xlsx", "-o", "out"])
    assert args.input_tsv == "book.xlsx"


def test_a_tsv_input_is_left_alone(egap, sample_tsv_factory, tmp_path):
    table = str(sample_tsv_factory(SAMPLE_ID="S1"))
    args = egap.resolve_sample_table(egap.parse_args(["-tsv", table, "-o", str(tmp_path / "out")]))
    assert args.input_tsv == table and args.input_workbook is None
    assert not (tmp_path / "out").exists()


# ---- workbook (needs openpyxl) ------------------------------------------------------

def _samples():
    return pd.DataFrame([
        {"SPECIES_ID": "Escherichia_coli", "SAMPLE_ID": "Escherichia_coli-A",
         "ORGANISM_KINGDOM": "Bacteria", "EST_SIZE": "5m", "BUSCO": "enterobacterales,bacteria",
         "ILLUMINA_RAW_READS": "/d/a_1.fastq,/d/a_2.fastq", "PROJECT_CODE": "X1"},
        {"SPECIES_ID": "Escherichia_coli", "SAMPLE_ID": "Escherichia_coli-B",
         "ORGANISM_KINGDOM": "Bacteria", "EST_SIZE": "5m", "BUSCO": "enterobacterales,bacteria",
         "ONT_SRA": "SRR1", "PROJECT_CODE": None},
    ])


def test_generated_workbook_has_a_tooltip_on_every_column(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    path = sample_workbook.write_workbook(tmp_path / "samples.xlsx", _samples())
    book = openpyxl.load_workbook(path)
    assert book.sheetnames == [sample_tsv.SAMPLE_SHEET, sample_tsv.METADATA_SHEET]

    samples, metadata = book[sample_tsv.SAMPLE_SHEET], book[sample_tsv.METADATA_SHEET]
    assert [c.value for c in samples[1]] == SAMPLE_NAMES + ["PROJECT_CODE"]
    assert [c.value for c in metadata[1]] == METADATA_NAMES
    for sheet in (samples, metadata):
        assert sheet.freeze_panes == "A2"
        for header in sheet[1]:
            assert header.comment is not None and "Required:" in header.comment.text, header.value
    assert "Not an EGAP column" in samples.cell(row=1, column=len(SAMPLE_NAMES) + 1).comment.text

    dropdowns = {str(v.sqref).split(":")[0][:-1]: v.formula1 for v in samples.data_validations.dataValidation}
    kingdom = openpyxl.utils.get_column_letter(SAMPLE_NAMES.index("ORGANISM_KINGDOM") + 1)
    assert dropdowns[kingdom] == '"Bacteria,Archaea,Flora,Funga,Fauna"'
    # samples are copied in, and the metadata sheet starts with their IDs
    assert samples.cell(row=3, column=SAMPLE_NAMES.index("ONT_SRA") + 1).value == "SRR1"
    assert [metadata.cell(row=r, column=1).value for r in (2, 3)] == ["Escherichia_coli-A", "Escherichia_coli-B"]


def test_workbook_reads_back_as_the_same_table(tmp_path, egap):
    openpyxl = pytest.importorskip("openpyxl")
    book_path = sample_workbook.write_workbook(tmp_path / "run1.xlsx", _samples())
    book = openpyxl.load_workbook(book_path)
    metadata = book[sample_tsv.METADATA_SHEET]
    metadata.cell(row=3, column=METADATA_NAMES.index("SAMPLE_TISSUE") + 1, value="culture")
    book.create_sheet("Notes")["A1"] = "free text, no SAMPLE_ID column"
    book.save(book_path)

    args = egap.resolve_sample_table(egap.parse_args(["-tsv", book_path, "-o", str(tmp_path / "out")]))
    assert args.input_workbook == book_path
    assert Path(args.input_tsv) == tmp_path / "out" / "run1_samples.tsv"
    assert egap.resolve_sample_table(args).input_tsv == args.input_tsv   # safe to call twice

    table = sample_tsv.read_sample_table(args.input_tsv)
    assert list(table["SAMPLE_ID"]) == ["Escherichia_coli-A", "Escherichia_coli-B"]
    assert table.loc[0, "ILLUMINA_RAW_R_READS"] == "/d/a_2.fastq" and table.loc[1, "ONT_SRA"] == "SRR1"
    assert pd.isna(table.loc[0, "SAMPLE_TISSUE"]) and table.loc[1, "SAMPLE_TISSUE"] == "culture"
    assert table.loc[0, "PROJECT_CODE"] == "X1"


@pytest.mark.parametrize("sheet_name, sample_ids, message", [
    ("Sample Metadata", ["Escherichia_coli-A", "Escherichia_coli-Z"], "not on sheet"),
    ("Sample Settings", ["Escherichia_coli-Z"], "not on sheet"),
    ("Sample Metadata", ["Escherichia_coli-A", "Escherichia_coli-A"], "repeats SAMPLE_ID"),
])
def test_rows_for_unknown_or_repeated_samples_stop_the_run(tmp_path, sheet_name, sample_ids, message):
    pytest.importorskip("openpyxl")
    book_path = tmp_path / "bad.xlsx"
    with pd.ExcelWriter(book_path) as writer:
        _samples().drop(columns="PROJECT_CODE").to_excel(writer, sheet_name="Sample Information", index=False)
        pd.DataFrame({"SAMPLE_ID": sample_ids}).to_excel(writer, sheet_name=sheet_name, index=False)
    with pytest.raises(ValueError, match=message):
        sample_workbook.workbook_to_tsv(book_path, tmp_path / "out")


def test_repeated_sample_on_the_first_sheet_stops_the_run(tmp_path, egap):
    pytest.importorskip("openpyxl")
    book_path = tmp_path / "dup.xlsx"
    pd.concat([_samples(), _samples().head(1)]).to_excel(book_path, sheet_name="Sample Information", index=False)
    with pytest.raises(SystemExit, match="repeats SAMPLE_ID: Escherichia_coli-A"):
        egap.resolve_sample_table(egap.parse_args(["-tsv", str(book_path), "-o", str(tmp_path / "out")]))


# ---- NCBI metadata stage (reads the bundled .xlsx templates, so needs openpyxl) -----

def _finished_sample(out_dir, sample_id, stats, long_read_source=None):
    """Lay down what the final QC and read selection leave behind for one sample."""
    species_dir = out_dir / "Escherichia_coli"
    (species_dir / sample_id).mkdir(parents=True)
    pd.DataFrame([dict(GENOME_SIZE=5_000_000.0, **stats)]).to_csv(
        species_dir / sample_id / f"{sample_id}_final_EGAP_assembly_final_stats.csv", index=False)
    if long_read_source:
        (species_dir / "ONT").mkdir(exist_ok=True)
        selected = species_dir / "ONT" / "Escherichia_coli_ONT_highest_mean_qual_long_reads.fastq"
        Path(str(selected) + sample_tsv.LONG_READ_INFO_SUFFIX).write_text(
            json.dumps({"source": long_read_source, "mean_quality": 12.0}))


def _table(tmp_path, rows):
    path = tmp_path / "run1_samples.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return str(path)


def test_metadata_sheets_accumulate_and_a_rerun_does_not_duplicate(tmp_path):
    pytest.importorskip("openpyxl")
    import process_metadata as stage
    out = tmp_path / "out"
    table = _table(tmp_path, [
        {"SPECIES_ID": "Escherichia_coli", "SAMPLE_ID": "Escherichia_coli-A",
         "ILLUMINA_RAW_F_READS": "/d/a_1.fastq", "ILLUMINA_RAW_R_READS": "/d/a_2.fastq",
         "ONT_RAW_READS": "/d/a_ont.fastq", "REF_SEQ_GCA": "GCA_000005845.2",
         "BIOSAMPLE_ACCESSION": "SAMN01", "SAMPLE_TISSUE": "culture",
         "SEQUENCING_METHODOLOGY": "WGS", "LIBRARY_SOURCE": "GENOMIC", "INSTRUMENT_ONT": "MinION"},
        {"SPECIES_ID": "Escherichia_coli", "SAMPLE_ID": "Escherichia_coli-B",
         "ONT_SRA": "SRR1", "SAMPLE_TISSUE": "culture", "SEQUENCING_METHODOLOGY": "WGS"},
        {"SPECIES_ID": "Escherichia_coli", "SAMPLE_ID": "Escherichia_coli-C",
         "PACBIO_RAW_READS": "/d/c.fastq"},
        {"SPECIES_ID": "Escherichia_coli", "SAMPLE_ID": "Escherichia_coli-QcOnly",
         "REF_SEQ_GCA": "GCA_000005845.2"},
    ])
    _finished_sample(out, "Escherichia_coli-A", {"DEDUPED_ILLU_TOTAL_BASES": 300e6,
                                                   "FILT_ONT_TOTAL_BASES": 200e6,
                                                   "RAW_ONT_TOTAL_BASES": 900e6}, "filtered")
    # A and B are one species, so they share the ONT folder and its "filtered" selection.
    _finished_sample(out, "Escherichia_coli-B", {"FILT_ONT_TOTAL_BASES": 100e6})
    _finished_sample(out, "Escherichia_coli-C", {"RAW_PACBIO_TOTAL_BASES": 50e6,
                                                   "FILT_PACBIO_TOTAL_BASES": 40e6})

    run = lambda sample_id: stage.process_metadata(sample_id, table, str(out), 1, 1)
    sra_a, assembly = run("Escherichia_coli-A")
    assert run("Escherichia_coli-B")[0] is None            # reads came from SRA
    assert run("Escherichia_coli-C")[0] is None            # no tissue or strategy
    assert run("Escherichia_coli-QcOnly") == (None, None)  # no final assembly
    run("Escherichia_coli-A")                              # rerun replaces, never duplicates

    sra = pd.read_csv(sra_a, sep="\t", dtype=str)
    assert Path(sra_a).name == "run1_samples_EGAP_SRA_metadata.tsv"
    assert list(sra["library_ID"]) == ["Escherichia_coli-A_ONT", "Escherichia_coli-A_Illumina"]
    ont, illumina = sra.iloc[0], sra.iloc[1]
    assert (ont["platform"], ont["library_layout"], ont["instrument_model"]) == ("OXFORD_NANOPORE", "SINGLE", "MinION")
    assert ont["filename"] == "Escherichia_coli_ONT_highest_mean_qual_long_reads.fastq" and pd.isna(ont["filename2"])
    assert (illumina["platform"], illumina["library_layout"]) == ("ILLUMINA", "PAIRED")
    assert illumina["filename2"] == "Escherichia_coli_illu_reverse_dedup.fastq"
    assert "Trimmomatic" in illumina["design_description"] and "Filtlong" in ont["design_description"]
    assert set(sra["biosample_accession"]) == {"SAMN01"} and set(sra["title"]) == {"WGS of Escherichia_coli: culture"}

    genomes = pd.read_csv(assembly, sep="\t", dtype=str).set_index("sample_name")
    assert list(genomes.index) == ["Escherichia_coli-B", "Escherichia_coli-C", "Escherichia_coli-A"]
    a = genomes.loc["Escherichia_coli-A"]
    assert a["genome_coverage"] == "100.0x"                # (300 + 200 filtered) Mb over 5 Mb
    assert a["sequencing_technology"] == "ONT, Illumina" and a["reference_genome"] == "GCA_000005845.2"
    assert a["assembly_method"] == "EGAP" and a["assembly_method_version"] == stage.VERSION
    assert a["filename"] == "Escherichia_coli-A_final_EGAP_assembly.fasta"
    assert genomes.loc["Escherichia_coli-B", "genome_coverage"] == "20.0x"
    # No record of which PacBio reads were selected: raw is the fallback.
    assert genomes.loc["Escherichia_coli-C", "genome_coverage"] == "10.0x"
    assert pd.isna(genomes.loc["Escherichia_coli-C", "biosample_accession"])
