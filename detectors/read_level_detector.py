#!/usr/bin/env python3
"""
foldback_hunter: CLI for reference-free foldback read-level detection

    python detectors/read_level_detector.py --fastq <input.fastq> \
        --outdir results/ \
        [--processes N | --threads N] [--mode probe|full|seed --max-reads N]

Writes calls_foldback_hunter_<stem>.tsv (read_id, method, flagged)
and scores_foldback_hunter_<stem>.tsv (read_id, raw_score,
fold_position_bp, status). flagged is the raw float score, not a threshold.
<stem> is the fastq filename stem (e.g. foo.fastq.gz -> foo).
"""

import argparse
import csv
import gzip
import itertools
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Pool

from foldback_score import score_read


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


def parse_stem(fastq_path):
    basename = os.path.basename(fastq_path)
    if basename.endswith(".gz"):
        basename = basename[: -len(".gz")]
    stem, _ext = os.path.splitext(basename)
    return stem


def _score_worker(task):
    read_id, seq, mode = task
    return score_read(read_id, seq, mode=mode)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--fastq", required=True, help="Input FASTQ (plain or .gz).")
    ap.add_argument("--outdir", default="results")
    parallel_group = ap.add_mutually_exclusive_group()
    parallel_group.add_argument("--processes", type=int, default=None,
        help="If set, score reads in parallel with this many worker processes.")
    parallel_group.add_argument("--threads", type=int, default=None,
        help=("If set, score reads in parallel with this many threads "
              "(concurrent.futures.ThreadPoolExecutor). Threads well for "
              "--mode full/probe (parasail/edlib release the GIL); does not "
              "thread well for --mode seed, whose k-mer matching loop is "
              "pure Python and holds the GIL — a warning is printed if you "
              "combine --threads with --mode seed, use --processes there "
              "instead for real speedup. Mutually exclusive with --processes."))
    ap.add_argument("--mode", choices=["probe", "full", "seed"], default="probe")
    ap.add_argument("--max-reads", type=int, default=None,
                     help="Cap the number of reads scored. Required when --mode full.")
    args = ap.parse_args()

    if args.mode == "full" and args.max_reads is None:
        raise SystemExit("--mode full requires --max-reads: no cap = no run.")

    if args.threads and args.mode == "seed":
        print(
            "[foldback_hunter] WARNING: --threads with --mode seed will not "
            "scale — seed mode's k-mer matching loop is pure Python and "
            "holds the GIL, so threads serialize on it (edlib's own GIL "
            "release inside that loop doesn't help). Use --processes for "
            "real speedup in seed mode. Continuing anyway.",
            file=sys.stderr,
        )

    stem = parse_stem(args.fastq)
    os.makedirs(args.outdir, exist_ok=True)

    reads = iter_fastq(args.fastq)
    if args.max_reads is not None:
        reads = itertools.islice(reads, args.max_reads)

    tasks = ((read_id, seq, args.mode) for read_id, seq in reads)

    if args.processes:
        with Pool(args.processes) as pool:
            results = pool.map(_score_worker, tasks)
    elif args.threads:
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            results = list(executor.map(_score_worker, tasks))
    else:
        results = [_score_worker(t) for t in tasks]

    calls_path = os.path.join(args.outdir, f"calls_foldback_hunter_{stem}.tsv")
    scores_path = os.path.join(args.outdir, f"scores_foldback_hunter_{stem}.tsv")

    with open(calls_path, "w", newline="") as calls_fh, \
         open(scores_path, "w", newline="") as scores_fh:
        calls_writer = csv.writer(calls_fh, delimiter="\t")
        scores_writer = csv.writer(scores_fh, delimiter="\t")
        calls_writer.writerow(["read_id", "method", "flagged"])
        scores_writer.writerow(["read_id", "raw_score", "fold_position_bp", "status"])
        for r in results:
            raw_score = "" if r.raw_score is None else r.raw_score
            fold_pos = "" if r.fold_position_bp is None else r.fold_position_bp
            calls_writer.writerow([r.read_id, "read_level_detector", raw_score])
            scores_writer.writerow([r.read_id, raw_score, fold_pos, r.status])

    print(f"[foldback_hunter] {len(results)} reads scored, mode={args.mode}")
    print(f"[foldback_hunter] wrote {calls_path}")
    print(f"[foldback_hunter] wrote {scores_path}")


if __name__ == "__main__":
    main()
