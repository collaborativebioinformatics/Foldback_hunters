#!/usr/bin/env python3
"""
foldback_hunter: CLI for reference-free foldback read-level detection

    python detectors/read_level_detector.py --input <reads.fastq|.fasta|.fa|.fna[.gz]|.bam> \
        --outdir results/ \
        [--processes N | --threads N] [--mode probe|full|seed] [--yes]

The whole input is always scored. To subset it, use seqtk sample (FASTQ) or
samtools view -s (BAM). --mode full asks for confirmation (read count + ETA)
before scoring; pass --yes/-y to skip the prompt (required when stdin is not
a TTY).

Writes calls_foldback_hunter_<stem>.tsv (read_id, method, flagged)
and scores_foldback_hunter_<stem>.tsv (read_id, raw_score,
fold_position_bp, status). flagged is a boolean: raw_score >= --threshold
(default 0.99). <stem> is the input filename stem (e.g. foo.fastq.gz -> foo).

Optional read output: --output <reads.fastq|.fasta[.gz]> with
--action {tag,clip,filter,split} (default tag). Every emitted record gets
" fb:i:1" (flagged) or " fb:i:0" appended to its header. clip keeps the longer
side of the fold; split emits <id>_1 / <id>_2; filter drops flagged reads;
clip/split cut at the detected fold position. FASTA output drops qualities. 
BAM input is written as FASTQ/FASTA with only fb:i:N
kept; records without qualities get '!' placeholders in FASTQ output.
"""

import argparse
import array
import bisect
import csv
import gzip
import itertools
import json
import math
import os
import re
import sys
import time
from collections import deque
from contextlib import ExitStack
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Pool

from foldback_score import revcomp, score_read


def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if path.lower().endswith(".gz") else open(path, mode)


def iter_fastq(path):
    """Yield (read_id, seq, qual, header_rest) from a FASTQ.
    seq keeps its original case; header_rest is the header text after the id"""
    with open_maybe_gz(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip("\n")
            plus = fh.readline()
            qual = fh.readline().rstrip("\n")
            if not plus:
                break
            fields = header.rstrip("\n").lstrip("@").split(None, 1)
            yield fields[0], seq, qual, (fields[1] if len(fields) > 1 else "")


def iter_fasta(path):
    """Yield (read_id, seq, None, header_rest) from a FASTA (multi-line records
    allowed). Holds one record at a time; seq keeps its original case."""
    def rec(header, parts):
        fields = header.split(None, 1)
        return (fields[0] if fields else ""), "".join(parts), None, (fields[1] if len(fields) > 1 else "")

    header = None
    parts = []
    with open_maybe_gz(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    yield rec(header, parts)
                header, parts = line[1:], []
            elif line and header is not None:
                parts.append(line)
    if header is not None:
        yield rec(header, parts)


ALIGNED_BAM_MSG = ("aligned BAM -> BAM is not supported; "
                   "use FASTQ/FASTA output or an unmapped BAM")


def read_ubam_header(path):
    """Return the header of an unmapped BAM; exit if it has @SQ (aligned)."""
    import pysam

    with pysam.AlignmentFile(path, "rb", check_sq=False) as af:
        if af.references:
            raise SystemExit(ALIGNED_BAM_MSG)
        return af.header


def iter_bam(path, keep_meta=False):
    """Yield (read_id, seq, qual, "") from a BAM (aligned or unaligned),
    primary records only, in original read orientation. 
    qual is None when the record has no qualities.
    keep_meta (uBAM -> uBAM): the 4th element is (flag, tags, seq) instead, and
    a mapped record is an error."""
    import pysam  # deferred: only needed for BAM input, keeps FASTQ path pysam-free

    with pysam.AlignmentFile(path, "rb", check_sq=False) as af:
        for rec in af:
            if keep_meta and not rec.is_unmapped:
                raise SystemExit(ALIGNED_BAM_MSG)
            if rec.is_secondary or rec.is_supplementary:
                continue
            seq = rec.query_sequence
            if seq is None:
                continue
            quals = rec.query_qualities
            qual = pysam.qualities_to_qualitystring(quals) if quals is not None else None
            if rec.is_reverse:
                seq = revcomp(seq)
                if qual is not None:
                    qual = qual[::-1]
            seq = seq.upper()
            yield rec.query_name, seq, qual, (
                (rec.flag, rec.get_tags(with_value_type=True), seq) if keep_meta else "")


def stored_fold_position(fold_position_bp, seq_len, is_reverse):
    return seq_len - fold_position_bp if is_reverse else fold_position_bp


def parse_stem(fastq_path):
    basename = os.path.basename(fastq_path)
    if basename.lower().endswith(".gz"):
        basename = basename[: -len(".gz")]
    stem, _ext = os.path.splitext(basename)
    return stem


FASTA_EXTS = (".fasta", ".fa", ".fna")


def detect_format(path):
    """Classify an input path by extension: fastq[.gz], fasta[.gz], or bam."""
    lower = path.lower()
    if lower.endswith(FASTA_EXTS):
        return "fasta"
    if lower.endswith(tuple(e + ".gz" for e in FASTA_EXTS)):
        return "fasta.gz"
    if lower.endswith(".fastq.gz") or lower.endswith(".fq.gz"):
        return "fastq.gz"
    if lower.endswith(".fastq") or lower.endswith(".fq"):
        return "fastq"
    if lower.endswith(".bam"):
        return "bam"
    raise SystemExit(
        f"Unrecognized input extension: {path!r}. "
        "Expected .fastq, .fq, .fasta, .fa, .fna (each optionally .gz), or .bam."
    )


def resolve_output_format(output, output_format, in_fmt):
    gz = output.lower().endswith(".gz")
    base = output.lower().removesuffix(".gz")
    ext_family = ("fasta" if base.endswith(FASTA_EXTS)
                  else "fastq" if base.endswith((".fastq", ".fq"))
                  else "bam" if base.endswith(".bam") else None)
    if output_format == "bam" or ext_family == "bam":
        if gz:
            raise ValueError("BAM output cannot be .gz")
    elif output_format and output_format.endswith(".gz") != gz:
        raise ValueError(f"--output-format {output_format} contradicts --output extension {output!r}")
    if output_format and ext_family and not output_format.startswith(ext_family):
        raise ValueError(f"--output-format {output_format} contradicts --output extension {output!r}")
    if output_format:
        fmt = output_format
    elif ext_family:
        fmt = ext_family + (".gz" if gz else "")
    elif gz:
        fmt = "fastq.gz"
    else:
        fmt = "fastq" if in_fmt == "bam" else in_fmt
    if in_fmt.startswith("fasta") and fmt.startswith("fastq"):
        raise ValueError("FASTA input has no qualities; use a FASTA or BAM output (--output-format fasta[.gz]|bam)")
    return fmt


def open_bam_writer(path, argv, template_header=None):
    """Unmapped-BAM writer: @HD SO:unknown + @PG with the command line, no @SQ.
    With template_header (uBAM input) its @HD/@RG/@PG/@CO are kept and our @PG
    is chained after the last input @PG."""
    import pysam  # deferred, like iter_bam

    header = template_header.to_dict() if template_header else {}
    header.setdefault("HD", {"VN": "1.6", "SO": "unknown"})
    pgs = header.setdefault("PG", [])
    ids = {pg["ID"] for pg in pgs}
    pg_id, n = "foldback_hunter", 0
    while pg_id in ids:
        n += 1
        pg_id = f"foldback_hunter.{n}"
    pg = {"ID": pg_id, "PN": "foldback_hunter", "CL": " ".join(argv)}
    if pgs:
        pg["PP"] = pgs[-1]["ID"]
    pgs.append(pg)
    return pysam.AlignmentFile(path, "wb", header=header)


_MM_ENTRY = re.compile(r"([ACGTUN])([+-])([a-z]+|\d+)([?.]?)((?:,\d+)*)$")


def slice_modbase(seq, mm, ml, a, b):
    """Rewrite MM/ML for the read slice seq[a:b]. Returns (mm, ml); both None
    when no call survives; ml is None if the input ml is."""
    mm_out, ml_out, off = [], [], 0
    for entry in filter(None, mm.split(";")):
        m = _MM_ENTRY.match(entry)
        if not m:
            raise ValueError(f"unparseable MM entry {entry!r}")
        base, strand, code, flag, skips = m.groups()
        n = len(code) if code.isalpha() else 1
        skips = [int(x) for x in skips.split(",")[1:]]
        pos = [i for i, c in enumerate(seq) if base == "N" or c == base]
        before = bisect.bisect_left(pos, a)
        new_skips, kept_ml, k, prev = [], [], -1, -1
        for j, s in enumerate(skips):
            k += s + 1
            if a <= pos[k] < b:
                new_k = k - before
                new_skips.append(new_k - prev - 1)
                prev = new_k
                kept_ml += ml[off + j * n: off + (j + 1) * n] if ml is not None else []
        off += len(skips) * n
        if new_skips:
            mm_out.append(f"{base}{strand}{code}{flag}" + "".join(f",{x}" for x in new_skips))
            ml_out += kept_ml
    if not mm_out:
        return None, None
    return ";".join(mm_out) + ";", (array.array("B", ml_out) if ml is not None else None)


def make_bam_record(header, rid, seq, qual, tag, meta=None, span=None):
    """meta = (flag, tags, full_seq) from iter_bam(keep_meta=True): copy flag and
    tags. span = (a, b) slice of the original read this record holds: MM/ML are
    rewritten to it, mv/ts dropped (not rewritten)."""
    import pysam

    a = pysam.AlignedSegment(header)
    a.query_name = rid
    a.flag = 4 if meta is None else meta[0]
    a.reference_id = -1
    a.reference_start = -1
    a.mapping_quality = 0
    a.query_sequence = seq  # must precede qualities: setting seq clears them
    if qual is not None:
        a.query_qualities = pysam.qualitystring_to_array(qual)
    if meta is not None:
        tags = {t: (v, vt) for t, v, vt in meta[1]}
        if span is not None:
            for t in ("mv", "ts"):
                tags.pop(t, None)
            if "MM" in tags:
                ml = tags["ML"][0] if "ML" in tags else None
                mm, ml = slice_modbase(meta[2], tags["MM"][0], ml, *span)
                tags.pop("MM"); tags.pop("ML", None)
                if mm is not None:
                    tags["MM"] = (mm, "Z")
                    if ml is not None:
                        tags["ML"] = (ml, "B")
            if "MN" in tags:
                tags["MN"] = (len(seq), "i")
        for t, (v, vt) in tags.items():
            if t != "fb":
                a.set_tag(t, v, value_type=None if vt == "B" else vt)  # arrays: type from typecode
    a.set_tag("fb", int(tag.rsplit(":", 1)[1]), value_type="i")
    return a


def format_fastq_record(rid, comment, seq, qual, tag):
    header = f"@{rid} {comment} {tag}" if comment else f"@{rid} {tag}"
    return f"{header}\n{seq}\n+\n{qual}\n"


def format_fasta_record(rid, comment, seq, tag):
    header = f">{rid} {comment} {tag}" if comment else f">{rid} {tag}"
    return f"{header}\n{seq}\n"


def _cut(qual, a, b):
    return None if qual is None else qual[a:b]


def apply_action(action, rec, flagged, fold_pos):
    """Return the output records [(id, seq, qual, tag, span)] (0, 1 or 2) for one read;
    span = (a, b) slice of the original read for cut records, else None.
    rec = (read_id, seq, qual, header_rest).
    Reads that are unflagged (or lack a fold position) pass through,
    tagged fb:i:0 (fb:i:1 if flagged without a position).
    """
    rid, seq, qual, _rest = rec
    tag = "fb:i:1" if flagged else "fb:i:0"
    if action == "filter":
        return [] if flagged else [(rid, seq, qual, tag, None)]
    if action == "tag" or not flagged or fold_pos is None:
        return [(rid, seq, qual, tag, None)]
    p = int(fold_pos)
    if action == "clip":
        if p >= len(seq) - p:
            return [(rid, seq[:p], _cut(qual, None, p), tag, (0, p))]
        return [(rid, seq[p:], _cut(qual, p, None), tag, (p, len(seq)))]
    # split
    return [(f"{rid}_1", seq[:p], _cut(qual, None, p), tag, (0, p)),
            (f"{rid}_2", seq[p:], _cut(qual, p, None), tag, (p, len(seq)))]


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
FASTA_GZ_BYTES_PER_READ = 4700


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
    if fmt == "fasta":
        with open(path, "rb") as fh:
            return sum(1 for line in fh if line.startswith(b">"))
    if fmt == "bam":
        import pysam

        with pysam.AlignmentFile(path, "rb", check_sq=False) as af:
            if af.has_index():
                # idxstats counts secondary/supplementary records too, so this
                # overestimates the primary reads actually scored.
                return sum(s.total for s in af.get_index_statistics()) + af.nocoordinate
    bytes_per_read = {"bam": BAM_BYTES_PER_READ, "fasta.gz": FASTA_GZ_BYTES_PER_READ}.get(
        fmt, FASTQ_GZ_BYTES_PER_READ)
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
    ap.add_argument("--input", required=True, help="Input reads: FASTQ or FASTA (plain or .gz) or BAM.")
    ap.add_argument("--outdir", default="results")
    parallel_group = ap.add_mutually_exclusive_group()
    parallel_group.add_argument("--processes", type=int, default=None,
        help=("If set, score reads in parallel with this many worker processes. "
              "Warns (does not fail) if above the available CPU count."))
    parallel_group.add_argument("--threads", type=int, default=None,
        help=("If set, score reads in parallel with this many threads "
              "(concurrent.futures.ThreadPoolExecutor). Threads well for "
              "--mode full/probe (parasail/edlib release the GIL); does not "
              "thread well for --mode seed, whose k-mer matching loop is "
              "pure Python and holds the GIL — a warning is printed if you "
              "combine --threads with --mode seed, use --processes there "
              "instead for real speedup. Mutually exclusive with --processes. Warns (does not fail) if "
              "above the available CPU count."))
    ap.add_argument("--mode", choices=["probe", "full", "seed"], default="probe")
    ap.add_argument("-y", "--yes", action="store_true",
        help="Skip the --mode full confirmation prompt (required if stdin is not a TTY). "
             "The whole input is always scored; to subset it, use seqtk sample / "
             "samtools view -s.")
    ap.add_argument("--threshold", type=threshold_type, default=0.99,
        help="Minimum raw_score to flag a read as foldback, range [0,1]. "
             "Default 0.99 from comprehensive benchmark analysis"
             " do not change without proper justification.")
    ap.add_argument("--output", default=None,
        help="Write reads to this FASTQ/FASTA (.gz ok) or unmapped BAM (.bam). FASTA output "
             "drops qualities (needed for FASTA input; FASTA->BAM stores none). BAM input to FASTQ/FASTA: "
             "tags/alignment dropped, only fb:i:N kept. FASTQ header comments are dropped in "
             "BAM output. Unmapped BAM -> BAM keeps the header, flag and tags (MM/ML "
             "rewritten, mv/ts dropped on clipped/split records); aligned BAM -> BAM errors. "
             "Without it only the TSVs/summary are written.")
    ap.add_argument("--output-format", choices=["fastq", "fastq.gz", "fasta", "fasta.gz", "bam"], default=None,
        help="Output read format. Default: from the --output extension, else the "
             "input format (.bam extension infers bam).")
    ap.add_argument("--action", choices=["tag", "clip", "filter", "split"], default=None,
        help="What to do with reads written to --output (default tag): tag appends "
             "fb:i:1|0 to the header; clip keeps the longer side of the fold; "
             "filter drops flagged reads; split emits <id>_1 and <id>_2")
    args = ap.parse_args()

    fmt = detect_format(args.input)
    if args.action and not args.output:
        ap.error("--action requires --output")
    action = args.action or "tag"
    if args.output:
        if os.path.abspath(args.output) == os.path.abspath(args.input):
            ap.error("--output must differ from --input")
        try:
            out_fmt = resolve_output_format(args.output, args.output_format, fmt)
        except ValueError as e:
            ap.error(str(e))
        # if action in ("clip", "split"):
        #     print(
        #         f"[foldback_hunter] WARNING: --action {action}"
        #         "cuts at a fold position that is unvalidated on "
        #         "noisy reads. Continuing anyway.",
        #         file=sys.stderr,
        #     )
    if args.mode == "full":
        confirm_full_run(args.input, fmt, args.yes)

    if args.threads and args.mode == "seed":
        print(
            "[foldback_hunter] WARNING: --threads with --mode seed will not "
            "scale — seed mode's k-mer matching loop is pure Python and "
            "holds the GIL, so threads serialize on it (edlib's own GIL "
            "release inside that loop doesn't help). Use --processes for "
            "real speedup in seed mode. Continuing anyway.",
            file=sys.stderr,
        )

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
    elif fmt in ("fasta", "fasta.gz"):
        reads = iter_fasta(args.input)
    elif fmt == "bam":
        reads = iter_bam(args.input, keep_meta=bool(args.output) and out_fmt == "bam")

    in_header = read_ubam_header(args.input) if (
        fmt == "bam" and args.output and out_fmt == "bam") else None

    pending = deque()  # records awaiting their result; bounded by _ordered_bounded

    warned_no_qual = False

    def _tasks():
        nonlocal warned_no_qual
        for rec in reads:
            if args.output:
                if rec[2] is None and out_fmt.startswith("fastq"):
                    if not warned_no_qual:
                        print("[foldback_hunter] WARNING: BAM record(s) without "
                              "qualities; writing '!' (Q0) placeholder",
                              file=sys.stderr)
                        warned_no_qual = True
                    rec = (rec[0], rec[1], "!" * len(rec[1]), rec[3])
                pending.append(rec)
            yield rec[0], rec[1].upper(), args.mode

    tasks = _tasks()

    calls_path = os.path.join(args.outdir, f"calls_foldback_hunter_{stem}.tsv")
    scores_path = os.path.join(args.outdir, f"scores_foldback_hunter_{stem}.tsv")

    chunksize = 1 if args.mode == "full" else CHUNKSIZE

    total_reads = 0
    flagged_reads = 0
    acted_reads = 0
    # runtime_sec covers scoring interleaved with TSV writing (streamed, so the
    # two cannot be timed separately).
    scoring_start = time.perf_counter()
    with open(calls_path, "w", newline="") as calls_fh, \
         open(scores_path, "w", newline="") as scores_fh, \
         ExitStack() as stack:
        out_fh = None
        if args.output and out_fmt == "bam":
            out_fh = stack.enter_context(open_bam_writer(args.output, sys.argv, in_header))
            bam_header = out_fh.header

            def emit(rid, comment, seq, qual, tag, span):
                meta = comment if isinstance(comment, tuple) else None
                out_fh.write(make_bam_record(bam_header, rid, seq, qual, tag, meta, span))
        elif args.output:
            out_fh = stack.enter_context(open_maybe_gz(args.output, "wt"))
            if out_fmt.startswith("fasta"):
                def emit(rid, comment, seq, qual, tag, span):
                    out_fh.write(format_fasta_record(rid, comment, seq, tag))
            else:
                def emit(rid, comment, seq, qual, tag, span):
                    out_fh.write(format_fastq_record(rid, comment, seq, qual, tag))
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
            if out_fh:
                rec = pending.popleft()
                for rid, seq, qual, tag, span in apply_action(action, rec, flagged, r.fold_position_bp):
                    emit(rid, rec[3], seq, qual, tag, span)
                if flagged and (action == "filter"
                                or (action in ("clip", "split") and r.fold_position_bp is not None)):
                    acted_reads += 1
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
        "clipped_reads": acted_reads if action == "clip" and args.output else 0,
        "filtered_reads": acted_reads if action == "filter" and args.output else 0,
        "foldback_rate": foldback_rate,
        "runtime_sec": round(runtime_sec, 4),
    }
    if args.output:
        summary["split_reads"] = acted_reads if action == "split" else 0
        summary["action"] = action

    summary_path = os.path.join(args.outdir, f"summary_foldback_hunter_{stem}.json")
    with open(summary_path, "w") as summary_fh:
        json.dump(summary, summary_fh, indent=2)
        summary_fh.write("\n")

    print(f"[foldback_hunter] {total_reads} reads scored, mode={args.mode}")
    print(f"[foldback_hunter] wrote {calls_path}")
    print(f"[foldback_hunter] wrote {scores_path}")
    print(f"[foldback_hunter] wrote {summary_path}")
    if args.output:
        print(f"[foldback_hunter] wrote {args.output}")


if __name__ == "__main__":
    main()
