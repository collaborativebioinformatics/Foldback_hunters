"""
Unit tests for detectors/read_level_detector.py's CLI: backend equivalence
(sequential / --processes / --threads), mutual exclusivity of --processes
and --threads, and the --threads + --mode seed warn-only behavior.
"""

import json
import random
import sys

import pytest
import read_level_detector as rld
from foldback_score import revcomp


def make_foldback_read(L, p, seed=0):
    """Build a noiseless foldback read: A + RC(B), true fold point p."""
    rng = random.Random(seed)
    seq = "".join(rng.choice("ACGT") for _ in range(L))
    t = min(p, L - p)
    first = seq[:p]
    tail_source = first[-t:]
    read = first + revcomp(tail_source)
    return read


def make_fastq(path):
    """A tiny FASTQ: two foldback reads, one unrelated (no-fold) read."""
    reads = [
        ("r1", make_foldback_read(2000, 1000, seed=1)),
        ("r2", make_foldback_read(2000, 1500, seed=2)),
        ("r3", "".join(random.Random(3).choice("ACGT") for _ in range(2000))),
    ]
    with open(path, "w") as fh:
        for read_id, seq in reads:
            fh.write(f"@{read_id}\n{seq}\n+\n{'I' * len(seq)}\n")
    return path


def _run(fastq_path, outdir, extra_args, monkeypatch, mode="probe"):
    argv = [
        "read_level_detector.py",
        "--input", str(fastq_path),
        "--outdir", str(outdir),
        "--mode", mode,
    ] + extra_args
    monkeypatch.setattr(sys, "argv", argv)
    rld.main()


def _read_tsvs(outdir, stem):
    calls_path = outdir / f"calls_foldback_hunter_{stem}.tsv"
    scores_path = outdir / f"scores_foldback_hunter_{stem}.tsv"
    return calls_path.read_text(), scores_path.read_text()


@pytest.mark.parametrize("extra_args", [[], ["--processes", "2"], ["--threads", "2"]])
def test_backends_produce_identical_output(tmp_path, monkeypatch, extra_args):
    fastq_path = make_fastq(tmp_path / "in.fastq")

    baseline_dir = tmp_path / "baseline"
    _run(fastq_path, baseline_dir, [], monkeypatch)
    baseline_calls, baseline_scores = _read_tsvs(baseline_dir, "in")

    variant_dir = tmp_path / "variant"
    _run(fastq_path, variant_dir, extra_args, monkeypatch)
    variant_calls, variant_scores = _read_tsvs(variant_dir, "in")

    assert variant_calls == baseline_calls
    assert variant_scores == baseline_scores


def test_threads_and_processes_are_mutually_exclusive(tmp_path, monkeypatch, capsys):
    fastq_path = make_fastq(tmp_path / "in.fastq")

    with pytest.raises(SystemExit):
        _run(
            fastq_path, tmp_path / "out",
            ["--processes", "2", "--threads", "2"],
            monkeypatch,
        )

    stderr = capsys.readouterr().err
    assert "not allowed with argument" in stderr


def test_threads_seed_mode_warns_but_completes(tmp_path, monkeypatch, capsys):
    fastq_path = make_fastq(tmp_path / "in.fastq")
    outdir = tmp_path / "out"

    _run(fastq_path, outdir, ["--threads", "2"], monkeypatch, mode="seed")

    stderr = capsys.readouterr().err
    assert "GIL" in stderr
    assert "--processes" in stderr
    assert (outdir / "calls_foldback_hunter_in.tsv").exists()
    assert (outdir / "scores_foldback_hunter_in.tsv").exists()


def test_summary_json_written(tmp_path, monkeypatch):
    fastq_path = make_fastq(tmp_path / "in.fastq")
    outdir = tmp_path / "out"

    _run(fastq_path, outdir, [], monkeypatch)

    calls_text, _ = _read_tsvs(outdir, "in")
    calls_rows = calls_text.strip().splitlines()[1:]
    flagged_count = sum(1 for row in calls_rows if row.split("\t")[2] == "True")

    summary = json.loads((outdir / "summary_foldback_hunter_in.json").read_text())

    assert summary["total_reads"] == 3
    assert summary["flagged_reads"] == flagged_count
    assert summary["clipped_reads"] == 0
    assert summary["filtered_reads"] == 0
    assert summary["foldback_rate"] == round(flagged_count / 3, 4)
    assert summary["mode"] == "probe"
    assert summary["threshold"] == 0.8
    assert summary["processes"] is None
    assert summary["threads"] is None
    assert isinstance(summary["runtime_sec"], (int, float))
    assert summary["runtime_sec"] >= 0


@pytest.mark.parametrize("path,expected", [
    ("reads.fastq", "fastq"),
    ("reads.fq", "fastq"),
    ("reads.fastq.gz", "fastq.gz"),
    ("reads.fq.gz", "fastq.gz"),
    ("reads.bam", "bam"),
])
def test_detect_format(path, expected):
    assert rld.detect_format(path) == expected


def test_detect_format_unrecognized_extension_raises():
    with pytest.raises(SystemExit):
        rld.detect_format("reads.sam")
