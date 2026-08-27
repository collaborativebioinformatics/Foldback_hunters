"""
run_all_conditions.py

Runs simulate_foldbacks_v2.py across the full 3x3 sweep:
    foldback_fraction in {0.01, 0.05, 0.10}
    fold_position     in {quarter, middle, near_end}

adapter_present defaults to "mixed" for every condition (both flavors
present in each dataset), matching the design where Task B's benchmark
facets recall by adapter_present within a single condition rather than
needing separate adapter-forced datasets.

Usage:
    python run_all_conditions.py \\
        --source_fastq chr6_100x_0001.fq.gz \\
        --source_maf chr6_100x_0001.maf.gz \\
        --outdir results/ \\
        --seed 1
"""

import argparse
import subprocess
import sys

FRACTIONS = [0.01, 0.05, 0.10]
POSITIONS = ["quarter", "middle", "near_end"]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source_fastq", required=True)
    ap.add_argument("--source_maf", required=True)
    ap.add_argument("--adapter_present", choices=["mixed", "true", "false"], default="mixed")
    ap.add_argument("--adapter_seq", default=None,
                     help="Pass through to simulate_foldbacks_v2.py if you have the "
                          "verified adapter sequence ready.")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    failures = []
    for fraction in FRACTIONS:
        for position in POSITIONS:
            cmd = [
                sys.executable, "simulate_foldbacks_v2.py",
                "--source_fastq", args.source_fastq,
                "--source_maf", args.source_maf,
                "--foldback_fraction", str(fraction),
                "--fold_position", position,
                "--adapter_present", args.adapter_present,
                "--outdir", args.outdir,
                "--seed", str(args.seed),
            ]
            if args.adapter_seq:
                cmd += ["--adapter_seq", args.adapter_seq]

            print(f"\n=== Running fraction={fraction} position={position} ===")
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