#!/usr/bin/env python3
"""
foldback_hunter: CLI for reference-free foldback read-level detection

    python detectors/read_level_detector.py --input <reads.fastq|.fastq.gz|.bam> \
        --outdir results/ \
        [--processes N | --threads N] [--mode probe|full|seed] [--yes]

The whole input is always scored. To subset it, use seqtk sample (FASTQ) or
samtools view -s (BAM). --mode full asks for confirmation (read count + ETA)
before scoring; pass --yes/-y to skip the prompt (required when stdin is not
a TTY).

Writes calls_foldback_hunter_<stem>.tsv (read_id, method, flagged)
and scores_foldback_hunter_<stem>.tsv (read_id, raw_score,
fold_position_bp, status). flagged is a boolean: raw_score >= --threshold
(default 0.8). <stem> is the input filename stem (e.g. foo.fastq.gz -> foo,
foo.bam -> foo).
"""

import argparse
import csv
import gzip
import itertools
import json
import math
import os
import sys
import time
from collections import deque
from contextlib import ExitStack
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Pool

from foldback_score import revcomp, score_read


def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def iter_fastq(path):
    """Yield (read_id, seq) from a FASTQ (plain or .gz)."""
    with open_maybe_gz(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip("\n")
            plus = fh.readline()
            fh.readline()  # qual, unused
            if not plus:
                break
            read_id = header.rstrip("\n").lstrip("@").split()[0]
            yield read_id, seq.upper()


def iter_bam(path):
    """Yield (read_id, seq) from a BAM (aligned or unaligned), primary records only."""
    import pysam  # deferred: only needed for BAM input, keeps FASTQ path pysam-free

    with pysam.AlignmentFile(path, "rb", check_sq=False) as af:
        for rec in af:
            if rec.is_secondary or rec.is_supplementary:
                continue
            seq = rec.query_sequence
            if seq is None:
                continue
            if rec.is_reverse:
                seq = revcomp(seq)
            yield rec.query_name, seq.upper()


def parse_stem(fastq_path):
    basename = os.path.basename(fastq_path)
    if basename.endswith(".gz"):
        basename = basename[: -len(".gz")]
    stem, _ext = os.path.splitext(basename)
    return stem


def detect_format(path):
    """Classify an input path by extension: fastq, fastq.gz, or bam."""
    lower = path.lower()
    if lower.endswith(".fastq.gz") or lower.endswith(".fq.gz"):
        return "fastq.gz"
    if lower.endswith(".fastq") or lower.endswith(".fq"):
        return "fastq"
    if lower.endswith(".bam"):
        return "bam"
    raise SystemExit(
        f"Unrecognized input extension: {path!r}. "
        "Expected .fastq, .fq, .fastq.gz, .fq.gz, or .bam."
    )


def _score_worker(task):
    read_id, seq, mode = task
    return score_read(read_id, seq, mode=mode)


CHUNKSIZE = 64  # probe/seed; --mode full uses 1 (~1.2 s/read, see FULL_SEC_PER_READ)


def _score_chunk(chunk):
    return [_score_worker(t) for t in chunk]


def _ordered_bounded(submit, wait, tasks, workers, chunksize):
    """Yield scored results in input order with a bounded number of chunks in flight.

    Pool.imap / executor.map drain the task generator eagerly, so memory still
    scales with read count; this keeps a bounded FIFO of submitted chunks. The
    window is sized in reads, max(4*workers, ceil(16*workers / chunksize))
    chunks, so small chunks still give enough lookahead that ordered yield does
    not stall on one slow read. At chunksize 64 this is 4*workers chunks.
    """
    window = max(4 * workers, math.ceil(16 * workers / chunksize))
    tasks = iter(tasks)
    inflight = deque()
    while True:
        chunk = list(itertools.islice(tasks, chunksize))
        if not chunk:
            break
        inflight.append(submit(chunk))
        if len(inflight) >= window:
            yield from wait(inflight.popleft())
    while inflight:
        yield from wait(inflight.popleft())


def threshold_type(value):
    threshold = float(value)
    if not 0 <= threshold <= 1:
        raise argparse.ArgumentTypeError(
            f"--threshold must be in [0, 1], got {value!r}"
        )
    return threshold


def available_cpus():
    """CPUs this process may use: the scheduler affinity mask where the OS
    exposes it (Linux, honors SLURM cpusets), else os.cpu_count()."""
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0))
    return os.cpu_count()


FULL_SEC_PER_READ = 1.2  # rough --mode full cost, used for the ETA
# Rough compressed bytes/read, measured on real HG002 ONT data (mean read length ~17-18 kb):
#   FASTQ.GZ: hg002_subsample.fastq.gz, 662 MB / 39,729 reads = 16.7 kB.
#   BAM: ONT open data giab_2025.01 R10 sup calls.sorted.bam (carries MM/ML modbase tags),
#        8 random 800-read windows = 23.1 kB aligned; same reads stripped to unaligned = 21.4 kB.
#        Unindexed BAMs are mostly unaligned, so the uBAM figure is used.
# Shorter reads or BAMs without modbase tags will be smaller per read.
FASTQ_GZ_BYTES_PER_READ = 16500
BAM_BYTES_PER_READ = 21000


def format_eta(seconds):
    seconds = int(round(seconds))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{h}h {m}m {s}s"


def count_reads(path, fmt):
    """Return an estimated read count. Cheap; see README for caveats."""
    if fmt == "fastq":
        with open(path, "rb") as fh:
            return sum(1 for _ in fh) // 4
    if fmt == "bam":
        import pysam

        with pysam.AlignmentFile(path, "rb", check_sq=False) as af:
            if af.has_index():
                # idxstats counts secondary/supplementary records too, so this
                # overestimates the primary reads actually scored.
                return sum(s.total for s in af.get_index_statistics()) + af.nocoordinate
    bytes_per_read = BAM_BYTES_PER_READ if fmt == "bam" else FASTQ_GZ_BYTES_PER_READ
    return os.path.getsize(path) // bytes_per_read


def confirm_full_run(path, fmt, assume_yes):
    """Print read count + ETA for --mode full; exit unless the user confirms."""
    n = count_reads(path, fmt)
    count = f"~{n} reads (estimated)"
    msg = f"[foldback_hunter] --mode full: {count}, ETA ~{format_eta(n * FULL_SEC_PER_READ)} at ~{FULL_SEC_PER_READ} s/read"
    if assume_yes:
        print(msg, file=sys.stderr)
        return
    if not sys.stdin.isatty():
        raise SystemExit(f"{msg}\nstdin is not a TTY; re-run with --yes to proceed.")
    print(msg, file=sys.stderr)
    if input("Proceed? [y/N] ").strip().lower() not in ("y", "yes"):
        raise SystemExit("Aborted.")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--input", required=True, help="Input reads: FASTQ (plain or .gz) or BAM.")
    ap.add_argument("--outdir", default="results")
    parallel_group = ap.add_mutually_exclusive_group()
    parallel_group.add_argument("--processes", type=int, default=None,
        help=("If set, score reads in parallel with this many worker processes. "
              "For --mode seed prefer --threads: similar speed, ~5x less memory "
              "(each worker process holds its own numba runtime, ~130 MiB). "
              "Warns (does not fail) if above the available CPU count."))
    parallel_group.add_argument("--threads", type=int, default=None,
        help=("If set, score reads in parallel with this many threads "
              "(concurrent.futures.ThreadPoolExecutor). Recommended for all "
              "modes, including --mode seed (its numba k-mer scan releases "
              "the GIL). Mutually exclusive with --processes. Warns (does not "
              "fail) if above the available CPU count."))
    ap.add_argument("--mode", choices=["probe", "full", "seed"], default="probe")
    ap.add_argument("-y", "--yes", action="store_true",
        help="Skip the --mode full confirmation prompt (required if stdin is not a TTY). "
             "The whole input is always scored; to subset it, use seqtk sample / "
             "samtools view -s.")
    ap.add_argument("--threshold", type=threshold_type, default=0.8,
        help="Minimum raw_score to flag a read as foldback, range [0,1]. "
             "Default 0.8 matches analysis/plot_Framing_A.py's published "
             "9-condition benchmark threshold — do not change the default "
             "without rerunning that analysis.")
    args = ap.parse_args()

    fmt = detect_format(args.input)
    if args.mode == "full":
        confirm_full_run(args.input, fmt, args.yes)

    if args.mode == "seed":
        # load the numba cache here before workers spawn.
        try:
            import seed_candidates
        except ImportError as e:
            raise SystemExit(
                f"--mode seed requires numba and numpy ({e}); "
                "pip install -r requirements.txt")

    requested = args.processes or args.threads  # mutually exclusive
    ncpu = available_cpus()
    if requested and ncpu and requested > ncpu:
        flag = "--processes" if args.processes else "--threads"
        print(
            f"[foldback_hunter] WARNING: {flag} {requested} exceeds the "
            f"{ncpu} CPUs available to this process; workers will "
            "oversubscribe cores. Continuing anyway.",
            file=sys.stderr,
        )

    stem = parse_stem(args.input)
    os.makedirs(args.outdir, exist_ok=True)

    if fmt in ("fastq", "fastq.gz"):
        reads = iter_fastq(args.input)
    elif fmt == "bam":
        reads = iter_bam(args.input)

    tasks = ((read_id, seq, args.mode) for read_id, seq in reads)

    calls_path = os.path.join(args.outdir, f"calls_foldback_hunter_{stem}.tsv")
    scores_path = os.path.join(args.outdir, f"scores_foldback_hunter_{stem}.tsv")

    chunksize = 1 if args.mode == "full" else CHUNKSIZE

    total_reads = 0
    flagged_reads = 0
    # runtime_sec covers scoring interleaved with TSV writing (streamed, so the
    # two cannot be timed separately).
    scoring_start = time.perf_counter()
    with open(calls_path, "w", newline="") as calls_fh, \
         open(scores_path, "w", newline="") as scores_fh, \
         ExitStack() as stack:
        if args.processes:
            pool = stack.enter_context(Pool(args.processes))
            results = _ordered_bounded(
                lambda c: pool.apply_async(_score_chunk, (c,)),
                lambda f: f.get(), tasks, args.processes, chunksize)
        elif args.threads:
            executor = stack.enter_context(ThreadPoolExecutor(max_workers=args.threads))
            results = _ordered_bounded(
                lambda c: executor.submit(_score_chunk, c),
                lambda f: f.result(), tasks, args.threads, chunksize)
        else:
            results = (_score_worker(t) for t in tasks)

        calls_writer = csv.writer(calls_fh, delimiter="\t")
        scores_writer = csv.writer(scores_fh, delimiter="\t")
        calls_writer.writerow(["read_id", "method", "flagged"])
        scores_writer.writerow(["read_id", "raw_score", "fold_position_bp", "status"])
        for r in results:
            raw_score = "" if r.raw_score is None else r.raw_score
            fold_pos = "" if r.fold_position_bp is None else r.fold_position_bp
            flagged = r.raw_score is not None and r.raw_score >= args.threshold
            calls_writer.writerow([r.read_id, "foldback_hunter", flagged])
            scores_writer.writerow([r.read_id, raw_score, fold_pos, r.status])
            total_reads += 1
            flagged_reads += flagged
    runtime_sec = time.perf_counter() - scoring_start

    foldback_rate = round(flagged_reads / total_reads, 4) if total_reads else 0.0

    summary = {
        "input": args.input,
        "mode": args.mode,
        "threshold": args.threshold,
        "processes": args.processes,
        "threads": args.threads,
        "total_reads": total_reads,
        "flagged_reads": flagged_reads,
        "clipped_reads": 0,
        "filtered_reads": 0,
        "foldback_rate": foldback_rate,
        "runtime_sec": round(runtime_sec, 4),
    }

    summary_path = os.path.join(args.outdir, f"summary_foldback_hunter_{stem}.json")
    with open(summary_path, "w") as summary_fh:
        json.dump(summary, summary_fh, indent=2)
        summary_fh.write("\n")

    print(f"[foldback_hunter] {total_reads} reads scored, mode={args.mode}")
    print(f"[foldback_hunter] wrote {calls_path}")
    print(f"[foldback_hunter] wrote {scores_path}")
    print(f"[foldback_hunter] wrote {summary_path}")


if __name__ == "__main__":
    main()
