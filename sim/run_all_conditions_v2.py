"""
run_all_conditions.py

Runs simulate_foldbacks.py across the full 3x3 sweep:
    foldback_fraction in {0.01, 0.05, 0.10}
    fold_position     in {middle, off_center, near_end}

adapter_present defaults to "mixed" for every condition (both flavors
present in each dataset), matching the design where Task B's benchmark
facets recall by adapter_present within a single condition rather than
needing separate adapter-forced datasets.

If --subsample_n is passed, it's forwarded to each simulate_foldbacks.py
run. The simulator uses a fixed subsampling seed internally, so all 9
conditions see the same clean reads (only the foldback assignments and
modifications differ across conditions).

Usage:
    python run_all_conditions.py \\
        --source_fastq chr6_100x_0001.fq.gz \\
        --source_maf chr6_100x_0001.maf.gz \\
        --subsample_n 230000 \\
        --outdir results/ \\
        --seed 1
"""

import argparse
import subprocess
import sys

FRACTIONS = [0.01, 0.05, 0.10]
POSITIONS = ["middle", "off_center", "near_end"]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source_fastq", required=True)
    ap.add_argument("--source_maf", required=True)
    ap.add_argument("--subsample_n", type=int, default=None,
                     help="If set, forwarded to simulate_foldbacks.py so every "
                          "condition uniformly subsamples the same N reads from the "
                          "input (via a fixed internal seed). Example: 230000 for "
                          "~20x on chr6.")
    ap.add_argument("--adapter_present", choices=["mixed", "true", "false"], default="mixed")
    ap.add_argument("--adapter_seq", default=None,
                     help="Pass through to simulate_foldbacks.py if you have the "
                          "verified adapter sequence ready.")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--dry_run", action="store_true",
                     help="Print the commands that would run without executing them.")
    args = ap.parse_args()

    failures = []
    for fraction in FRACTIONS:
        for position in POSITIONS:
            cmd = [
                sys.executable, "simulate_foldbacks.py",
                "--source_fastq", args.source_fastq,
                "--source_maf", args.source_maf,
                "--foldback_fraction", str(fraction),
                "--fold_position", position,
                "--adapter_present", args.adapter_present,
                "--outdir", args.outdir,
                "--seed", str(args.seed),
            ]
            if args.subsample_n:
                cmd += ["--subsample_n", str(args.subsample_n)]
            if args.adapter_seq:
                cmd += ["--adapter_seq", args.adapter_seq]

            print(f"\n=== fraction={fraction} position={position} ===")
            if args.dry_run:
                print("Would run:", " ".join(cmd))
                continue

            result = subprocess.run(cmd)
            if result.returncode != 0:
                failures.append((fraction, position))

    print("\n=== Sweep complete ===")
    if failures:
        print(f"{len(failures)} condition(s) FAILED: {failures}")
        sys.exit(1)
    else:
        print("All 9 conditions completed and passed sanity checks.")


if __name__ == "__main__":
    main()
