"""
simulate_foldbacks.py

Takes a clean pbsim3-simulated ONT FASTQ (HG002 chr6) and converts a
controlled fraction of reads into foldback artifacts, emitting:
  - a modified FASTQ (clean reads untouched, foldback reads rewritten)
  - a truth table (TSV) with one row per read

WHAT A FOLDBACK IS
-------------------
A foldback read's second half is the reverse complement of its first half.
It happens naturally in ONT sequencing when a DNA molecule folds back on
itself inside the pore, or when the template and complement strands pass
through together. To fake one here: take a normal clean read, keep a
fraction of it from the start (the "fold point"), and paste on the reverse
complement of that same chunk. The result still looks like one continuous
read, exactly as it would if the sequencer produced it directly.

TWO FLAVORS OF FOLDBACK (why --adapter_present exists)
--------------------------------------------------------
1. Adapter-bridged: template and complement strands go through the pore
   back-to-back with the sequencing adapter still stuck between them.
   Dorado's duplex-splitting is *supposed* to catch these by searching for
   the adapter sequence inside the read.
2. Sequence-only: the molecule physically folds on itself, no adapter
   between the two halves -- direct sequence-to-reverse-complement
   junction. Adapter-search tools structurally cannot see these: there is
   no adapter sequence in the read for them to find.

Within a single run, foldback reads get adapter_present assigned per-read
(default: 20% get an adapter, 80% are sequence-only, via
--adapter_probability) so that a single dataset contains both flavors and
Task B's benchmark can facet recall by adapter_present. Use
--adapter_present true/false to force all foldback reads in a run to one
flavor instead, for targeted testing.

TRUTH TABLE SCHEMA (frozen -- do not change without telling the group)
------------------------------------------------------------------------
    read_id          matches the FASTQ header exactly, no leading '@'
    is_foldback       True / False
    fold_position     middle / quarter / near_end for foldbacks, empty for clean
    adapter_present   True / False for foldbacks, empty for clean
    source_locus      chr6 coordinates the original read came from, pulled
                       from the pbsim3 .maf file

OUTPUT NAMING
-------------
    truth_<fraction>_<position>_<adapterflag>.tsv
    sim_<fraction>_<position>_<adapterflag>.fastq

    e.g. truth_5pct_middle_mixed.tsv / sim_5pct_middle_mixed.fastq

USAGE
-----
    python simulate_foldbacks.py \\
        --source_fastq chr6_100x_0001.fq.gz \\
        --source_maf chr6_100x_0001.maf.gz \\
        --foldback_fraction 0.05 \\
        --fold_position middle \\
        --adapter_present mixed \\
        --outdir results/ \\
        --seed 1

Run this 9 times (3 fractions x 3 positions) -- see run_all_conditions.py
for a driver script that does this for you in one call.

*** PLACEHOLDER ADAPTER SEQUENCE -- READ THIS ***
The constant PLACEHOLDER_ADAPTER below is NOT verified against the actual
ONT LSK114 kit adapter sequence or against what duplex-tools searches for.
I could not confirm the exact current sequence. Before running the real
adapter_present-vs-duplex-tools comparison (the "half the paper" result),
replace PLACEHOLDER_ADAPTER with the verified sequence -- pull it from
duplex-tools' own source/config, or from ONT's community forum /
documentation for the SQK-LSK114 kit. Using an unverified sequence here
would silently invalidate that comparison, since duplex-tools would be
searching for a different string than the one actually embedded.
"""

import argparse
import csv
import gzip
import os
import random
import sys

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# *** VERIFY BEFORE REAL RUNS *** -- see module docstring above.
PLACEHOLDER_ADAPTER = "AATGTACTTCGTTCAGTTACGTATTGCT"  # UNVERIFIED placeholder

FOLD_POSITION_CHOICES = ("middle", "quarter", "near_end")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


# --------------------------------------------------------------------------
# I/O: pbsim3 FASTQ + MAF
# --------------------------------------------------------------------------

def iter_fastq(path):
    """Yield (read_id, seq, qual) from a FASTQ (plain or .gz)."""
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
            read_id = header.rstrip("\n").lstrip("@").split()[0]
            yield read_id, seq.upper(), qual


def parse_maf_source_loci(path):
    """
    Parse a pbsim3 MAF(.gz) into: read_id -> "ref_name:start-end"

    pbsim3 MAF blocks:
        a
        s ref_name  ref_start  ref_aln_size  ref_strand  ref_size  ref_aln_seq
        s read_id   q_start    q_aln_size    q_strand    q_size    q_aln_seq
    """
    loci = {}
    with open_maybe_gz(path) as fh:
        s_lines = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("a"):
                s_lines = []
            elif line.startswith("s"):
                s_lines.append(line.split())
                if len(s_lines) == 2:
                    ref_fields, q_fields = s_lines
                    ref_name = ref_fields[1]
                    ref_start = int(ref_fields[2])
                    ref_aln_size = int(ref_fields[3])
                    read_id = q_fields[1]
                    loci[read_id] = f"{ref_name}:{ref_start}-{ref_start + ref_aln_size}"
                    s_lines = []
    return loci


def write_fastq_record(fh, read_id, seq, qual):
    fh.write(f"@{read_id}\n{seq}\n+\n{qual}\n")


# --------------------------------------------------------------------------
# Foldback construction
# --------------------------------------------------------------------------

def choose_fold_point(seq_len, fold_position, rng):
    """
    Return the fold point (length of the untouched leading segment) for the
    given fold_position category.
      quarter   -> ~20-30% of read length
      middle    -> ~45-55% of read length
      near_end  -> fold lands in the last 500bp: RC tail length is a random
                   value between 50 and 500bp (so fold_point = L - tail_len)
    """
    if fold_position == "quarter":
        frac = rng.uniform(0.20, 0.30)
        return max(1, int(seq_len * frac))
    elif fold_position == "middle":
        frac = rng.uniform(0.45, 0.55)
        return max(1, int(seq_len * frac))
    elif fold_position == "near_end":
        tail_len = rng.randint(50, min(500, seq_len - 1))
        return max(1, seq_len - tail_len)
    else:
        raise ValueError(f"Unknown fold_position: {fold_position}")


def make_foldback(seq, qual, fold_position, adapter_present, adapter_seq, rng):
    """
    Build seq[:fold_point] + (adapter?) + revcomp(seq[:fold_point]).
    Returns (new_seq, new_qual).
    """
    L = len(seq)
    fold_point = choose_fold_point(L, fold_position, rng)

    first = seq[:fold_point]
    rc = revcomp(first)
    junction = adapter_seq if adapter_present else ""
    new_seq = first + junction + rc

    first_q = qual[:fold_point]
    rc_q = first_q[::-1]  # approximate: mirror quality with the folded segment
    junction_q = "I" * len(junction)  # placeholder high Phred for adapter bases
    new_qual = first_q + junction_q + rc_q

    return new_seq, new_qual


# --------------------------------------------------------------------------
# Main dataset construction
# --------------------------------------------------------------------------

def build_dataset(reads, loci, foldback_fraction, fold_position, adapter_mode,
                   adapter_seq, adapter_probability, rng):
    """
    reads: iterable of (read_id, seq, qual) from the clean pbsim3 FASTQ.
    Yields dicts with the truth-table fields plus final seq/qual.

    adapter_probability only applies when adapter_mode == "mixed": it's the
    per-foldback-read chance of getting an adapter at the fold junction
    (default 0.20, i.e. ~20% of foldbacks are adapter-bridged, ~80% are
    sequence-only).
    """
    for read_id, seq, qual in reads:
        source_locus = loci.get(read_id, "unknown")
        is_foldback = rng.random() < foldback_fraction

        if not is_foldback:
            yield {
                "read_id": read_id, "is_foldback": False, "fold_position": "",
                "adapter_present": "", "source_locus": source_locus,
                "seq": seq, "qual": qual,
            }
            continue

        if adapter_mode == "mixed":
            adapter_present = rng.random() < adapter_probability
        elif adapter_mode == "true":
            adapter_present = True
        else:  # "false"
            adapter_present = False

        new_seq, new_qual = make_foldback(
            seq, qual, fold_position, adapter_present, adapter_seq, rng
        )
        yield {
            "read_id": read_id, "is_foldback": True, "fold_position": fold_position,
            "adapter_present": adapter_present, "source_locus": source_locus,
            "seq": new_seq, "qual": new_qual,
        }


# --------------------------------------------------------------------------
# Sanity checks (run before anything gets committed)
# --------------------------------------------------------------------------

def run_sanity_checks(fastq_path, truth_path):
    """
    Checks:
      - row counts match between fastq and truth table
      - all read_ids unique
      - no empty fold_position when is_foldback == True
    Raises AssertionError with a clear message on the first failure.
    """
    fastq_ids = []
    with open(fastq_path) as fh:
        lines = fh.read().splitlines()
    for i in range(0, len(lines), 4):
        fastq_ids.append(lines[i][1:])  # strip '@'

    with open(truth_path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    assert len(fastq_ids) == len(rows), (
        f"Row count mismatch: FASTQ has {len(fastq_ids)} reads, "
        f"truth table has {len(rows)} rows."
    )

    truth_ids = [r["read_id"] for r in rows]
    assert len(truth_ids) == len(set(truth_ids)), "Duplicate read_id values found in truth table."
    assert set(fastq_ids) == set(truth_ids), "read_id sets differ between FASTQ and truth table."

    for r in rows:
        if r["is_foldback"] == "True":
            assert r["fold_position"] not in ("", None), (
                f"Row {r['read_id']} has is_foldback=True but empty fold_position."
            )
            assert r["adapter_present"] in ("True", "False"), (
                f"Row {r['read_id']} has is_foldback=True but adapter_present is not True/False."
            )
        else:
            assert r["fold_position"] == "", (
                f"Row {r['read_id']} has is_foldback=False but non-empty fold_position."
            )

    print(f"[sanity check] PASSED: {len(rows)} reads, "
          f"{sum(1 for r in rows if r['is_foldback']=='True')} foldback, "
          "row counts match, all read_ids unique, no missing fold_position.")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--source_fastq", required=True,
                     help="Clean pbsim3 output FASTQ (.fq.gz or plain).")
    ap.add_argument("--source_maf", required=True,
                     help="Clean pbsim3 output MAF (.maf.gz or plain); supplies source_locus.")
    ap.add_argument("--foldback_fraction", type=float, required=True,
                     help="Fraction of reads to convert to foldbacks, e.g. 0.05 for 5%%.")
    ap.add_argument("--fold_position", choices=FOLD_POSITION_CHOICES, required=True,
                     help="Where the fold happens: middle (~50%%), quarter (~25%%), "
                          "or near_end (in the last 500bp).")
    ap.add_argument("--adapter_present", choices=["mixed", "true", "false"], default="mixed",
                     help="mixed (default): each foldback read gets an adapter with "
                          "probability --adapter_probability, so both flavors appear in the "
                          "same dataset. true/false: force all foldback reads to one flavor.")
    ap.add_argument("--adapter_probability", type=float, default=0.20,
                     help="[mixed mode only] Per-foldback-read probability of being "
                          "adapter-bridged rather than sequence-only. Default 0.20 (20%%).")
    ap.add_argument("--adapter_seq", default=PLACEHOLDER_ADAPTER,
                     help="Adapter sequence inserted at the fold junction when adapter_present. "
                          "*** Default is an UNVERIFIED placeholder -- see module docstring. ***")
    ap.add_argument("--outdir", default="results", help="Output directory.")
    ap.add_argument("--seed", type=int, default=0, help="Random seed for reproducibility.")
    ap.add_argument("--skip_sanity_check", action="store_true",
                     help="Skip the post-write sanity checks (not recommended).")
    args = ap.parse_args()

    if args.adapter_seq == PLACEHOLDER_ADAPTER:
        print("WARNING: using the UNVERIFIED placeholder adapter sequence. "
              "Replace with the real LSK114 / duplex-tools search sequence "
              "before trusting the adapter_present comparison.", file=sys.stderr)

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    # Condition label for filenames, e.g. 5pct_middle_mixed
    pct = int(round(args.foldback_fraction * 100))
    condition = f"{pct}pct_{args.fold_position}_{args.adapter_present}"

    fastq_path = os.path.join(args.outdir, f"sim_{condition}.fastq")
    truth_path = os.path.join(args.outdir, f"truth_{condition}.tsv")

    loci = parse_maf_source_loci(args.source_maf)
    reads = iter_fastq(args.source_fastq)

    n_total = 0
    n_foldback = 0
    with open(fastq_path, "w") as fq_out, open(truth_path, "w") as truth_out:
        truth_out.write("read_id\tis_foldback\tfold_position\tadapter_present\tsource_locus\n")
        for row in build_dataset(
            reads, loci, args.foldback_fraction, args.fold_position,
            args.adapter_present, args.adapter_seq, args.adapter_probability, rng
        ):
            write_fastq_record(fq_out, row["read_id"], row["seq"], row["qual"])
            truth_out.write(
                f"{row['read_id']}\t{row['is_foldback']}\t{row['fold_position']}\t"
                f"{row['adapter_present']}\t{row['source_locus']}\n"
            )
            n_total += 1
            n_foldback += int(row["is_foldback"])

    rate = n_foldback / n_total if n_total else 0.0
    print(f"[{condition}] wrote {n_total} reads -> {fastq_path}")
    print(f"[{condition}] wrote truth table -> {truth_path}")
    print(f"[{condition}] foldback reads: {n_foldback} ({rate:.2%}, target {args.foldback_fraction:.2%})")

    if not args.skip_sanity_check:
        run_sanity_checks(fastq_path, truth_path)


if __name__ == "__main__":
    main()