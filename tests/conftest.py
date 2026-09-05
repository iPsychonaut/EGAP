"""pytest configuration for EGAP.

Puts ``bin/`` and the repo root on sys.path so the pipeline modules import
exactly as they do when EGAP.py inserts bin/ at runtime, and so the audit
prototypes import as ``audit.prototypes.<name>``.

No external bioinformatics binaries are required. Tests that would need one
are marked ``requires_tools`` and skipped unless EGAP_RUN_TOOL_TESTS=1.
"""
import os
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"

for p in (str(REPO_ROOT), str(BIN_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

sys.path.insert(0, str(Path(__file__).resolve().parent))  # for `import mockdata`


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "requires_tools: needs conda-installed bioinformatics binaries")


def pytest_collection_modifyitems(config, items):
    if os.environ.get("EGAP_RUN_TOOL_TESTS") == "1":
        return
    skip = pytest.mark.skip(reason="set EGAP_RUN_TOOL_TESTS=1 to run tool-backed tests")
    for item in items:
        if "requires_tools" in item.keywords:
            item.add_marker(skip)


@pytest.fixture
def rng():
    import mockdata
    return mockdata.seeded(42)


@pytest.fixture
def mock_genome(rng, tmp_path):
    """A 5-contig genome with 2 planted high-GC contaminant contigs."""
    import mockdata
    recs = mockdata.make_genome(rng, n_contigs=5, contaminants=2)
    fasta = mockdata.write_fasta(tmp_path / "mock_genome.fasta", recs)
    return {"records": recs, "fasta": fasta, "ids": [r[0] for r in recs]}


@pytest.fixture
def sample_tsv_factory(tmp_path):
    """Build a one-row EGAP sample TSV using the collapsed v3.4.2 header."""
    header = ["SPECIES_ID", "SAMPLE_ID", "ORGANISM_KINGDOM", "ORGANISM_KARYOTE",
              "PLOIDY", "EST_SIZE", "BUSCO", "ONT_SRA", "ONT_RAW_DIR", "ONT_RAW_READS",
              "ILLUMINA_SRA", "ILLUMINA_RAW_DIR", "ILLUMINA_RAW_READS",
              "PACBIO_SRA", "PACBIO_RAW_DIR", "PACBIO_RAW_READS",
              "REF_SEQ_GCA", "REF_SEQ", "NT_ASSEMBLY_GCA", "NT_ASSEMBLY_PATH",
              "TX_EVI_GCA", "TX_EVI_PATH", "PT_EVI_GCA", "PT_EVI_PATH"]

    def _make(name="sample.tsv", **cells):
        row = {h: "" for h in header}
        row.update(cells)
        path = tmp_path / name
        with open(path, "w") as fh:
            fh.write("\t".join(header) + "\n")
            fh.write("\t".join(str(row[h]) for h in header) + "\n")
        return path

    return _make
