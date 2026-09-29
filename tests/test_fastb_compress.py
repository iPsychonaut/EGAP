"""bin/file_operations.py: FASTB for nucleotide FASTA, pigz for everything else."""
import os
import shutil

import pytest

import file_operations as fo

pytestmark = pytest.mark.skipif(
    not (shutil.which("fastb") and shutil.which("pigz")),
    reason="needs the fastb and pigz executables on PATH")

BARE = ">contig_1\nACGTNNNNacgtACGT\nACGT\n>contig_2\nGGGGCCCC\n"


def _write(path, text):
    path.write_text(text)
    return str(path)


def test_bare_headers_use_fastb_and_round_trip(tmp_path):
    fasta = _write(tmp_path / "asm.fasta", BARE)
    out = fo.pigz_compress(fasta, 1)
    assert out == fasta + ".fastb"
    assert os.path.exists(out) and not os.path.exists(fasta)

    assert fo.pigz_decompress(out, 1) == fasta
    assert not os.path.exists(out)
    with open(fasta, "rb") as fh:
        assert list(fo._fasta_records(fh)) == list(fo._fasta_records(BARE.encode().splitlines()))


@pytest.mark.parametrize("text", [
    ">contig_1 length=16 cov=30x\nACGTACGTACGTACGT\n",   # fastb drops the description
    ">prot_1\nMKFLILLFNILCLFPVLAADNHGVGPQGAS\n",         # fastb rejects protein
    "no header here\nACGT\n",                            # not FASTA at all
])
def test_lossy_or_rejected_input_falls_back_to_pigz(tmp_path, text):
    fasta = _write(tmp_path / "x.fasta", text)
    out = fo.pigz_compress(fasta, 1)
    assert out == fasta + ".gz"
    assert os.path.exists(out) and not os.path.exists(fasta + ".fastb")
    assert fo.pigz_decompress(out, 1) == fasta
    assert open(fasta).read() == text


def test_fastq_and_missing_fastb_use_pigz(tmp_path, monkeypatch):
    fastq = _write(tmp_path / "r.fastq", "@r1\nACGT\n+\nIIII\n")
    assert fo.pigz_compress(fastq, 1) == fastq + ".gz"

    monkeypatch.setattr(fo.shutil, "which", lambda name: None)
    fasta = _write(tmp_path / "asm.fasta", BARE)
    assert fo.pigz_compress(fasta, 1) == fasta + ".gz"


@pytest.mark.parametrize("ext", [".gz", ".fastb"])
def test_restore_if_compressed(tmp_path, ext, monkeypatch):
    fasta = _write(tmp_path / "asm.fasta", BARE)
    if ext == ".gz":
        monkeypatch.setattr(fo.shutil, "which", lambda name: None)
    assert fo.pigz_compress(fasta, 1) == fasta + ext
    monkeypatch.undo()
    assert fo.restore_if_compressed(fasta, 1) == fasta
    assert os.path.exists(fasta) and not os.path.exists(fasta + ext)
