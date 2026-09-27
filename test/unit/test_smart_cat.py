"""Unit tests for bin/smart_cat.sh — decompress by content, not by file name.

The staged RVDB-prot release is named rvdb.fasta.gz but its URL ends .xz. Deciding the
decompressor from the extension made `gzip -dc` fail with "not in gzip format", which looks like
a corrupt download; the mobile stage went partial and n_mycovirus stayed NA on a full run.
"""
import gzip
import os
import shutil
import subprocess

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(os.path.dirname(HERE)), "bin", "smart_cat.sh")
PAYLOAD = b">seq1\nMKVLAA\n>seq2\nMTTQ\n"


def run(path):
    r = subprocess.run(["bash", SCRIPT, str(path)], capture_output=True)
    assert r.returncode == 0, r.stderr.decode()
    return r.stdout


def test_plain_file_passes_through(tmp_path):
    p = tmp_path / "a.fasta"; p.write_bytes(PAYLOAD)
    assert run(p) == PAYLOAD


def test_gzip_content(tmp_path):
    p = tmp_path / "a.fasta.gz"
    with gzip.open(p, "wb") as fh:
        fh.write(PAYLOAD)
    assert run(p) == PAYLOAD


def test_gzip_content_under_the_wrong_name(tmp_path):
    p = tmp_path / "a.fasta.xz"          # deliberately wrong extension
    with gzip.open(p, "wb") as fh:
        fh.write(PAYLOAD)
    assert run(p) == PAYLOAD


@pytest.mark.skipif(shutil.which("xz") is None, reason="xz not installed")
def test_xz_content_under_a_gz_name(tmp_path):
    """The exact shape that broke the mycovirus screen on the cluster."""
    plain = tmp_path / "a.fasta"; plain.write_bytes(PAYLOAD)
    p = tmp_path / "a.fasta.gz"
    with open(p, "wb") as out:
        subprocess.run(["xz", "-c", str(plain)], stdout=out, check=True)
    assert run(p) == PAYLOAD


@pytest.mark.skipif(shutil.which("bzip2") is None, reason="bzip2 not installed")
def test_bzip2_content(tmp_path):
    plain = tmp_path / "a.fasta"; plain.write_bytes(PAYLOAD)
    p = tmp_path / "a.fasta.bz2"
    with open(p, "wb") as out:
        subprocess.run(["bzip2", "-c", str(plain)], stdout=out, check=True)
    assert run(p) == PAYLOAD


def test_an_unreadable_file_fails_loudly(tmp_path):
    r = subprocess.run(["bash", SCRIPT, str(tmp_path / "missing.gz")], capture_output=True)
    assert r.returncode != 0 and b"cannot read" in r.stderr
