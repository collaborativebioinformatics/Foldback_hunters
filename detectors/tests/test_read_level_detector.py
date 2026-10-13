"""
Unit tests for detectors/read_level_detector.py's CLI: backend equivalence
(sequential / --processes / --threads), mutual exclusivity of --processes
and --threads, the --threads + --mode seed warn-only behavior, and the warn-only
oversubscription check.
"""

import json
import os
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


@pytest.mark.parametrize("flag", ["--threads", "--processes"])
def test_oversubscription_warns_but_completes(tmp_path, monkeypatch, capsys, flag):
    monkeypatch.setattr(rld, "available_cpus", lambda: 2)
    fastq_path = make_fastq(tmp_path / "in.fastq")
    outdir = tmp_path / "out"

    _run(fastq_path, outdir, [flag, "4"], monkeypatch)

    stderr = capsys.readouterr().err
    assert "exceeds" in stderr
    assert flag in stderr
    assert (outdir / "calls_foldback_hunter_in.tsv").exists()


def test_no_warning_at_cpu_count(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(rld, "available_cpus", lambda: 4)
    fastq_path = make_fastq(tmp_path / "in.fastq")

    _run(fastq_path, tmp_path / "out", ["--threads", "4"], monkeypatch)

    assert "exceeds" not in capsys.readouterr().err


def test_available_cpus_prefers_affinity(monkeypatch):
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: {0, 1, 2}, raising=False)
    monkeypatch.setattr(os, "cpu_count", lambda: 64)
    assert rld.available_cpus() == 3


def test_available_cpus_falls_back_to_cpu_count(monkeypatch):
    monkeypatch.delattr(os, "sched_getaffinity", raising=False)
    monkeypatch.setattr(os, "cpu_count", lambda: 8)
    assert rld.available_cpus() == 8


# --- --mode full confirmation ------------------------------------------------

def _stdin(monkeypatch, isatty, answer=""):
    class FakeStdin:
        def isatty(self):
            return isatty
    monkeypatch.setattr(sys, "stdin", FakeStdin())
    monkeypatch.setattr("builtins.input", lambda prompt="": answer)


@pytest.mark.parametrize("answer", ["y", "YES", " yes "])
def test_full_prompt_accept(tmp_path, monkeypatch, capsys, answer):
    fq = make_fastq(tmp_path / "in.fastq")
    _stdin(monkeypatch, True, answer)
    _run(fq, tmp_path / "out", [], monkeypatch, mode="full")
    assert "~3 reads (estimated), ETA ~0h 0m 4s" in capsys.readouterr().err
    assert (tmp_path / "out" / "calls_foldback_hunter_in.tsv").exists()


@pytest.mark.parametrize("answer", ["", "n", "no", "maybe"])
def test_full_prompt_decline(tmp_path, monkeypatch, answer):
    fq = make_fastq(tmp_path / "in.fastq")
    _stdin(monkeypatch, True, answer)
    with pytest.raises(SystemExit) as exc:
        _run(fq, tmp_path / "out", [], monkeypatch, mode="full")
    assert exc.value.code not in (0, None)
    assert not (tmp_path / "out" / "calls_foldback_hunter_in.tsv").exists()


@pytest.mark.parametrize("flag", ["--yes", "-y"])
def test_full_yes_skips_prompt(tmp_path, monkeypatch, flag):
    fq = make_fastq(tmp_path / "in.fastq")
    _stdin(monkeypatch, False)  # non-TTY; would exit without --yes
    monkeypatch.setattr("builtins.input", lambda prompt="": pytest.fail("prompted"))
    _run(fq, tmp_path / "out", [flag], monkeypatch, mode="full")
    assert (tmp_path / "out" / "calls_foldback_hunter_in.tsv").exists()


def test_full_non_tty_without_yes_exits(tmp_path, monkeypatch):
    fq = make_fastq(tmp_path / "in.fastq")
    _stdin(monkeypatch, False)
    with pytest.raises(SystemExit) as exc:
        _run(fq, tmp_path / "out", [], monkeypatch, mode="full")
    msg = str(exc.value.code)
    assert "~3 reads (estimated)" in msg and "ETA" in msg and "re-run with --yes" in msg
    assert not (tmp_path / "out").exists() or not list((tmp_path / "out").iterdir())


def test_probe_mode_never_prompts(tmp_path, monkeypatch):
    fq = make_fastq(tmp_path / "in.fastq")
    _stdin(monkeypatch, False)
    _run(fq, tmp_path / "out", [], monkeypatch, mode="probe")


def test_gz_count_is_estimated(tmp_path):
    import gzip
    fq = make_fastq(tmp_path / "in.fastq")
    gz = tmp_path / "in.fastq.gz"
    gz.write_bytes(gzip.compress(fq.read_bytes()))
    assert rld.count_reads(str(gz), "fastq.gz") == gz.stat().st_size // rld.FASTQ_GZ_BYTES_PER_READ


def test_format_eta():
    assert rld.format_eta(3723) == "1h 2m 3s"
