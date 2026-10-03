"""bin/file_operations.py: INTERMEDIATE_FORMAT dispatch, FASTB round-trip, pigz fallback."""
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


@pytest.fixture()
def fmt(monkeypatch, request):
    monkeypatch.setenv("EGAP_INTERMEDIATE_FORMAT", request.param)
    return request.param


def test_default_format_is_pigz(monkeypatch):
    monkeypatch.delenv("EGAP_INTERMEDIATE_FORMAT", raising=False)
    assert fo.intermediate_format() == "pigz"
    monkeypatch.setenv("EGAP_INTERMEDIATE_FORMAT", "bogus")
    with pytest.raises(ValueError):
        fo.intermediate_format()


@pytest.mark.parametrize("fmt", ["pigz", "fastb"], indirect=True)
def test_compress_intermediate_follows_setting(tmp_path, fmt):
    fasta = _write(tmp_path / "asm.fasta", BARE)
    out = fo.compress_intermediate(fasta, 1)
    assert out == fasta + (".fastb" if fmt == "fastb" else ".gz")
    assert os.path.exists(out) and not os.path.exists(fasta)

    assert fo.decompress_intermediate(out, 2) == fasta
    assert os.path.exists(fasta) and not os.path.exists(out)
    # fastb decode wraps at 80 columns, so ignore line breaks
    with open(fasta) as fh:
        assert fh.read().replace("\n", "") == BARE.replace("\n", "")


@pytest.mark.parametrize("fmt", ["fastb"], indirect=True)
def test_fastq_always_uses_pigz(tmp_path, fmt):
    fastq = _write(tmp_path / "r.fastq", "@r1\nACGT\n+\nIIII\n")
    assert fo.compress_intermediate(fastq, 1) == fastq + ".gz"


@pytest.mark.parametrize("text", [
    ">contig_1 length=16 cov=30x\nACGTACGTACGTACGT\n",   # fastb drops the description
    ">prot_1\nMKFLILLFNILCLFPVLAADNHGVGPQGAS\n",         # fastb rejects protein
    "no header here\nACGT\n",                            # not FASTA at all
])
def test_lossy_or_rejected_input_falls_back_to_pigz(tmp_path, text):
    fasta = _write(tmp_path / "x.fasta", text)
    out = fo.fastb_compress(fasta, 1)
    assert out == fasta + ".gz"
    assert os.path.exists(out) and not os.path.exists(fasta + ".fastb")
    assert fo.pigz_decompress(out, 1) == fasta
    assert open(fasta).read() == text


def test_missing_fastb_executable_falls_back_to_pigz(tmp_path, monkeypatch):
    monkeypatch.setattr(fo.shutil, "which", lambda name: None)
    fasta = _write(tmp_path / "asm.fasta", BARE)
    assert fo.fastb_compress(fasta, 1) == fasta + ".gz"


@pytest.mark.parametrize("ext", [".gz", ".fastb"])
def test_restore_if_compressed(tmp_path, ext):
    fasta = _write(tmp_path / "asm.fasta", BARE)
    fn = fo.fastb_compress if ext == ".fastb" else fo.pigz_compress
    assert fn(fasta, 1) == fasta + ext
    assert fo.restore_if_compressed(fasta, 1) == fasta
    assert os.path.exists(fasta) and not os.path.exists(fasta + ext)
