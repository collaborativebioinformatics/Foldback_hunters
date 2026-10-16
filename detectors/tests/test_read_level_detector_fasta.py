"""Tests for FASTA[.gz] input/output in detectors/read_level_detector.py."""

import gzip
import sys

import pytest
import read_level_detector as rld
from test_read_level_detector import _read_tsvs, _run, make_fastq


def fastq_to_fasta(fq, fa, wrap=None):
    lines = open(fq).read().split("\n")
    opener = gzip.open if str(fa).endswith(".gz") else open
    with opener(fa, "wt") as fh:
        for i in range(0, len(lines) - 3, 4):
            seq = lines[i + 1]
            if wrap:
                seq = "\n".join(seq[j:j + wrap] for j in range(0, len(seq), wrap))
            fh.write(f">{lines[i][1:]}\n{seq}\n\n")
    return fa


def read_fasta(path):
    return list(rld.iter_fasta(str(path)))


@pytest.mark.parametrize("name,fmt", [
    ("x.fasta", "fasta"), ("x.fa", "fasta"), ("x.fna", "fasta"),
    ("x.fasta.gz", "fasta.gz"), ("x.fa.gz", "fasta.gz"), ("x.fna.gz", "fasta.gz"),
])
def test_detect_format(name, fmt):
    assert rld.detect_format(name) == fmt


def test_detect_format_unknown():
    with pytest.raises(SystemExit):
        rld.detect_format("x.txt")


@pytest.mark.parametrize("suffix", [".fa", ".fa.gz"])
def test_iter_fasta(tmp_path, suffix):
    fa = tmp_path / f"r{suffix}"
    opener = gzip.open if suffix.endswith(".gz") else open
    with opener(fa, "wt") as fh:
        fh.write(">a comment here\nACGT\nacgt\n\n>b\nTTTT\n")
    assert read_fasta(fa) == [("a", "ACGTacgt", None, "comment here"), ("b", "TTTT", None, "")]


def test_fasta_matches_fastq_calls(tmp_path, monkeypatch):
    fq = make_fastq(tmp_path / "reads.fastq")
    fa = fastq_to_fasta(fq, tmp_path / "reads.fa", wrap=70)
    _run(fq, tmp_path / "o1", [], monkeypatch)
    _run(fa, tmp_path / "o2", [], monkeypatch)
    assert _read_tsvs(tmp_path / "o1", "reads") == _read_tsvs(tmp_path / "o2", "reads")


@pytest.mark.parametrize("action", ["tag", "clip", "split", "filter"])
@pytest.mark.parametrize("out", ["out.fa", "out.fa.gz"])
def test_fasta_output_actions(tmp_path, monkeypatch, action, out):
    fq = make_fastq(tmp_path / "reads.fastq")
    fa = fastq_to_fasta(fq, tmp_path / "reads.fa")
    _run(fq, tmp_path / "o1", ["--output", str(tmp_path / "ref.fastq"), "--action", action], monkeypatch)
    _run(fa, tmp_path / "o2", ["--output", str(tmp_path / out), "--action", action], monkeypatch)
    ref = [l.split("\n") for l in open(tmp_path / "ref.fastq").read().split("@")[1:]]
    got = read_fasta(tmp_path / out)
    assert [(r[0], len(r[1])) for r in got] == [(h.split()[0], len(s)) for h, s, *_ in ref]
    assert all(r[3].startswith("fb:i:") or "fb:i:" in r[3] for r in got)


def test_fasta_in_fastq_out_errors(tmp_path, monkeypatch):
    fa = fastq_to_fasta(make_fastq(tmp_path / "reads.fastq"), tmp_path / "reads.fa")
    with pytest.raises(SystemExit):
        _run(fa, tmp_path / "o", ["--output", str(tmp_path / "out.fastq")], monkeypatch)


def test_no_placeholder_warning_for_fasta_output(tmp_path, monkeypatch, capsys):
    fa = fastq_to_fasta(make_fastq(tmp_path / "reads.fastq"), tmp_path / "reads.fa")
    _run(fa, tmp_path / "o", ["--output", str(tmp_path / "out.fa")], monkeypatch)
    assert "without qualities" not in capsys.readouterr().err


def test_count_reads_fasta(tmp_path):
    fa = fastq_to_fasta(make_fastq(tmp_path / "reads.fastq"), tmp_path / "reads.fa")
    assert rld.count_reads(str(fa), "fasta") == 3
