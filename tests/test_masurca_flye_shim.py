"""bin/assemble_masurca.py: the <prefix>/Flye/bin shim MaSuRCA calls Flye through."""
import os
import shutil
import subprocess

import pytest

import assemble_masurca as am

pytestmark = pytest.mark.skipif(shutil.which("bash") is None or os.name == "nt",
                                reason="the shim is a bash script")


def _fake_env(tmp_path):
    bin_dir = tmp_path / "prefix" / "bin"
    bin_dir.mkdir(parents=True)
    for name, body in (("masurca", "exit 0"),
                       ("flye", 'printf "%s\\n" "$@"'),
                       ("flye-minimap2", "exit 0")):
        p = bin_dir / name
        p.write_text(f"#!/usr/bin/env bash\n{body}\n")
        p.chmod(0o755)
    return bin_dir


def test_shim_drops_kmer_size_and_links_minimap2(tmp_path, monkeypatch):
    bin_dir = _fake_env(tmp_path)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    vendor = tmp_path / "prefix" / "Flye" / "bin"
    vendor.mkdir(parents=True)
    os.symlink(bin_dir / "flye", vendor / "flye")  # an old-style link gets replaced

    am._masurca_vendor_flye_ensure_shim()

    out = subprocess.run([str(vendor / "flye"), "-t", "7", "--nano-corr", "r.fa",
                          "--kmer-size", "21", "-m", "2500", "-o", "out dir"],
                         capture_output=True, text=True, check=True).stdout.splitlines()
    assert out == ["-t", "7", "--nano-corr", "r.fa", "-m", "2500", "-o", "out dir"]
    assert os.path.realpath(vendor / "flye-minimap2") == os.path.realpath(bin_dir / "flye-minimap2")
