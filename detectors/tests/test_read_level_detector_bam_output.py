"""Tests for unmapped-BAM output (Phase 03) in detectors/read_level_detector.py."""

import random

import pysam
import pytest
import read_level_detector as rld
from test_read_level_detector import _make_fq, _parse, _reads, _run
from test_read_level_detector_fasta import fastq_to_fasta


def _bam_records(path):
    with pysam.AlignmentFile(str(path), "rb", check_sq=False) as af:
        return [(r.query_name, r.query_sequence,
                 None if r.query_qualities is None else pysam.qualities_to_qualitystring(r.query_qualities),
                 r.flag, r.get_tag("fb")) for r in af]


@pytest.mark.parametrize("action", ["tag", "clip", "filter", "split"])
def test_bam_matches_fastq_output(tmp_path, monkeypatch, action):
    fq = _make_fq(tmp_path / "in.fastq", _reads())
    ref = tmp_path / "ref.fastq"
    out = tmp_path / "out.bam"
    _run(fq, tmp_path / "o1", ["--output", str(ref), "--action", action], monkeypatch)
    _run(fq, tmp_path / "o2", ["--output", str(out), "--action", action], monkeypatch)
    expected = [(i.split()[0], q, ql) for i, q, ql in _parse(ref)]
    got = _bam_records(out)
    assert [g[:3] for g in got] == expected
    assert all(g[3] == 4 for g in got)
    assert all(g[4] in (0, 1) for g in got)
    if action == "split":
        assert any(g[0].endswith("_1") for g in got) and any(g[0].endswith("_2") for g in got)


def test_gz_input(tmp_path, monkeypatch):
    import gzip
    fq = _make_fq(tmp_path / "in.fastq", _reads())
    gz = tmp_path / "in.fastq.gz"
    gz.write_bytes(gzip.compress(fq.read_bytes()))
    _run(gz, tmp_path / "o", ["--output", str(tmp_path / "out.bam")], monkeypatch)
    assert len(_bam_records(tmp_path / "out.bam")) == 3


def test_fasta_to_bam_no_quals_no_warning(tmp_path, monkeypatch, capsys):
    fq = _make_fq(tmp_path / "in.fastq", _reads())
    fa = fastq_to_fasta(fq, tmp_path / "in.fasta")
    _run(fa, tmp_path / "o", ["--output", str(tmp_path / "out.bam")], monkeypatch)
    recs = _bam_records(tmp_path / "out.bam")
    assert len(recs) == 3 and all(r[2] is None for r in recs)
    assert "WARNING" not in capsys.readouterr().err


def test_header_and_comment_dropped(tmp_path, monkeypatch):
    fq = _make_fq(tmp_path / "in.fastq", _reads())
    out = tmp_path / "out.bam"
    _run(fq, tmp_path / "o", ["--output", str(out)], monkeypatch)
    with pysam.AlignmentFile(str(out), "rb", check_sq=False) as af:
        h = af.header.to_dict()
        assert not h.get("SQ")
        assert h["PG"][0]["ID"] == "foldback_hunter"
        assert h["HD"]["SO"] == "unknown"
        names = [r.query_name for r in af]
    assert names == ["r1", "r2", "r3"]  # "runid=abc" comment is not in the name


def test_reingest(tmp_path, monkeypatch):
    fq = _make_fq(tmp_path / "in.fastq", _reads())
    out = tmp_path / "out.bam"
    _run(fq, tmp_path / "o", ["--output", str(out)], monkeypatch)
    again = [(r[0], r[1]) for r in rld.iter_bam(str(out))]
    assert again == [(r[0], r[1]) for r in rld.iter_fastq(str(fq))]


@pytest.mark.parametrize("par", [["--processes", "2"], ["--threads", "2"]])
def test_parallel_deterministic(tmp_path, monkeypatch, par):
    fq = _make_fq(tmp_path / "in.fastq", _reads() * 5)
    a, b = tmp_path / "a.bam", tmp_path / "b.bam"
    _run(fq, tmp_path / "o1", ["--output", str(a), "--action", "split"], monkeypatch)
    _run(fq, tmp_path / "o2", ["--output", str(b), "--action", "split"] + par, monkeypatch)
    assert _bam_records(a) == _bam_records(b)


def test_extension_inference_and_errors(tmp_path, monkeypatch):
    f = rld.resolve_output_format
    assert f("x.bam", None, "fastq") == "bam"
    assert f("x", "bam", "fasta") == "bam"
    for args in [("x.fastq", "bam", "fastq"), ("x.bam", "fastq", "fastq"),
                 ("x.bam.gz", None, "fastq"), ("x.bam", "bam", "bam"),
                 ("x.bam", None, "bam")]:
        with pytest.raises(ValueError):
            f(*args)


def test_near_end_fold_split(tmp_path, monkeypatch):
    # fold at ~0 / ~len: split must not crash (empty SEQ is legal in BAM)
    seq = "".join(random.Random(9).choice("ACGT") for _ in range(300))
    recs = [("e0", seq, 0), ("e1", seq, len(seq))]
    out = tmp_path / "out.bam"
    fq = _make_fq(tmp_path / "in.fastq", [(r, s, "") for r, s, _ in recs])
    # exercise the writer directly with edge positions
    with rld.open_bam_writer(str(out), ["t"]) as w:
        for rid, s, p in recs:
            for orid, oseq, oqual, tag in rld.apply_action("split", (rid, s, "I" * len(s), ""), True, p):
                w.write(rld.make_bam_record(w.header, orid, oseq, oqual, tag))
    got = _bam_records(out)
    assert len(got) == 4
    assert sum(len(g[1] or "") for g in got) == 2 * len(seq)
