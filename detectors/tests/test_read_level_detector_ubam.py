"""Tests for uBAM -> uBAM in detectors/read_level_detector.py."""

import array
import random

import pysam
import pytest
import read_level_detector as rld
from test_read_level_detector import _make_fq, _parse, _reads, _run


def _write_ubam(path, reads, with_sq=False, mapped=False):
    header = {"HD": {"VN": "1.6", "SO": "unknown"},
              "RG": [{"ID": "rg1", "SM": "s"}],
              "PG": [{"ID": "basecaller", "PN": "dorado"}],
              "CO": ["hello"]}
    if with_sq:
        header["SQ"] = [{"SN": "chr1", "LN": 100000}]
    with pysam.AlignmentFile(str(path), "wb", header=header) as w:
        for rid, seq, qual in reads:
            a = pysam.AlignedSegment(w.header)
            a.query_name, a.query_sequence = rid, seq
            a.query_qualities = pysam.qualitystring_to_array(qual)
            a.flag = 4
            if mapped:
                a.flag = 0  # mapped flag without @SQ
            a.set_tag("RG", "rg1")
            a.set_tag("mv", array.array("B", [5, 1, 0, 1]))
            a.set_tag("ts", 7, "i")
            a.set_tag("fb", 9, "i")
            c = [i for i, x in enumerate(seq) if x == "C"]
            if c:
                a.set_tag("MM", "C+m?," + ",".join(["0"] * len(c)) + ";")
                a.set_tag("ML", array.array("B", [i % 256 for i in range(1, len(c) + 1)]))
            w.write(a)


def _ubam(tmp_path):
    fq = _make_fq(tmp_path / "in.fastq", _reads())
    reads = [(i.split()[0], s, q) for i, s, q in _parse(fq)]
    bam = tmp_path / "in.bam"
    _write_ubam(bam, reads)
    return bam, reads, fq


def _recs(path):
    with pysam.AlignmentFile(str(path), "rb", check_sq=False) as af:
        return list(af)


@pytest.mark.parametrize("action", ["tag", "clip", "filter", "split"])
def test_ubam_matches_fastq_output(tmp_path, monkeypatch, action):
    bam, _reads_, fq = _ubam(tmp_path)
    ref, out = tmp_path / "ref.fastq", tmp_path / "out.bam"
    _run(fq, tmp_path / "o1", ["--output", str(ref), "--action", action], monkeypatch)
    _run(bam, tmp_path / "o2", ["--output", str(out), "--action", action], monkeypatch)
    got = _recs(out)
    assert [(r.query_name, r.query_sequence, pysam.qualities_to_qualitystring(r.query_qualities))
            for r in got] == [(i.split()[0], q, ql) for i, q, ql in _parse(ref)]
    for r in got:
        assert r.flag == 4 and r.get_tag("RG") == "rg1" and r.get_tag("fb") in (0, 1)
        cut = r.query_name.endswith(("_1", "_2")) or action == "clip" and r.get_tag("fb") == 1
        assert r.has_tag("mv") != cut and r.has_tag("ts") != cut


def test_header_preserved(tmp_path, monkeypatch):
    bam, _, _ = _ubam(tmp_path)
    out = tmp_path / "out.bam"
    _run(bam, tmp_path / "o", ["--output", str(out)], monkeypatch)
    with pysam.AlignmentFile(str(out), "rb", check_sq=False) as af:
        h = af.header.to_dict()
    assert h["RG"][0]["ID"] == "rg1" and h["CO"] == ["hello"] and not h.get("SQ")
    assert [p["ID"] for p in h["PG"]] == ["basecaller", "foldback_hunter"]
    assert h["PG"][1]["PP"] == "basecaller"


def test_uncut_tags_preserved(tmp_path, monkeypatch):
    bam, _, _ = _ubam(tmp_path)
    out = tmp_path / "out.bam"
    _run(bam, tmp_path / "o", ["--output", str(out)], monkeypatch)
    for a, b in zip(_recs(bam), _recs(out)):
        ta = {t: v for t, v in a.get_tags() if t != "fb"}
        assert {t: v for t, v in b.get_tags() if t != "fb"} == ta


@pytest.mark.parametrize("kw", [{"with_sq": True}, {"mapped": True}])
@pytest.mark.parametrize("action", ["tag", "split"])
def test_aligned_errors_before_output(tmp_path, monkeypatch, kw, action):
    bam = tmp_path / "al.bam"
    _write_ubam(bam, [("a", "ACGT" * 100, "I" * 400)], **kw)
    out = tmp_path / "out.bam"
    with pytest.raises(SystemExit, match="aligned BAM -> BAM is not supported"):
        _run(bam, tmp_path / "o", ["--output", str(out), "--action", action], monkeypatch)
    if "with_sq" in kw:
        assert not out.exists()


def test_reingest(tmp_path, monkeypatch):
    bam, reads, _ = _ubam(tmp_path)
    out = tmp_path / "out.bam"
    _run(bam, tmp_path / "o", ["--output", str(out), "--action", "tag"], monkeypatch)
    assert [(r[0], r[1]) for r in rld.iter_bam(str(out))] == [(r[0], r[1]) for r in reads]


@pytest.mark.parametrize("par", [["--processes", "2"], ["--threads", "2"]])
def test_parallel_deterministic(tmp_path, monkeypatch, par):
    bam, _, _ = _ubam(tmp_path)
    a, b = tmp_path / "a.bam", tmp_path / "b.bam"
    _run(bam, tmp_path / "o1", ["--output", str(a), "--action", "split"], monkeypatch)
    _run(bam, tmp_path / "o2", ["--output", str(b), "--action", "split"] + par, monkeypatch)
    assert [r.to_string() for r in _recs(a)] == [r.to_string() for r in _recs(b)]


# slice_modbase: seq = A C G C A C T C  (C at 1, 3, 5, 7)
SEQ = "ACGCACTC"


def test_slice_modbase_basic():
    mm, ml = rld.slice_modbase(SEQ, "C+m?,0,0,1;", [10, 20, 30], 3, 8)
    assert mm == "C+m?,0,1;" and list(ml) == [20, 30]  # calls at C pos 1,3,7 -> keep 3,7


def test_slice_modbase_first_call_skip():
    # calls at C ordinals 1,2,... skip 1 -> 2nd C (pos 3), then 0 -> pos 5
    mm, ml = rld.slice_modbase(SEQ, "C+m?,1,0;", [10, 20], 4, 8)
    assert mm == "C+m?,0;" and list(ml) == [20]


def test_slice_modbase_two_codes_and_all_dropped():
    mm, ml = rld.slice_modbase(SEQ, "C+hm?,0,0;A+a,0;", [1, 2, 3, 4, 9], 0, 4)
    assert mm == "C+hm?,0,0;A+a,0;" and list(ml) == [1, 2, 3, 4, 9]
    assert rld.slice_modbase(SEQ, "C+m?,0;", [1], 6, 8) == (None, None)


def test_slice_modbase_minus_strand_N_and_no_ml():
    assert rld.slice_modbase(SEQ, "C-m,0,0,0,0;", None, 2, 6) == ("C-m,0,0;", None)
    mm, ml = rld.slice_modbase(SEQ, "N+n?,2;", [5], 0, 4)
    assert mm == "N+n?,2;" and list(ml) == [5]


def test_slice_modbase_edge_cuts():
    assert rld.slice_modbase(SEQ, "C+m?,0,0,0,0;", [1, 2, 3, 4], 0, 0) == (None, None)
    mm, ml = rld.slice_modbase(SEQ, "C+m?,0,0,0,0;", [1, 2, 3, 4], 0, 8)
    assert mm == "C+m?,0,0,0,0;" and list(ml) == [1, 2, 3, 4]


def test_cut_record_mm_ml_and_mn(tmp_path):
    seq = "ACGCACTC"
    tags = [("MM", "C+m?,0,0,1;", "Z"), ("ML", array.array("B", [1, 2, 3]), "B"),
            ("MN", 8, "i"), ("mv", array.array("B", [1]), "B")]
    h = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}})
    a = rld.make_bam_record(h, "x_2", seq[4:], "I" * 4, "fb:i:1", (4, tags, seq), (4, 8))
    # calls at C pos 1,3,7 -> only 7 lies in [4,8): ordinal 1 in the slice "ACTC"
    assert a.get_tag("MM") == "C+m?,1;" and list(a.get_tag("ML")) == [3]
    assert a.get_tag("MN") == 4 and not a.has_tag("mv")
