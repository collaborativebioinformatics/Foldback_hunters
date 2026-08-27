#!/usr/bin/env bash
#SBATCH --job-name=breakinator_clean
#SBATCH --output=logs/breakinator_clean_%j.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --partition=batch

set -eo pipefail

module load miniconda/3.13
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate foldback

FASTQ=~/foldback/PSIM-simulation/Test_sample_10k.fastq.gz
TRUTH=~/foldback/sim_data/truth_clean_control.tsv
REF=~/foldback/GRCh38.fasta
OUT_DIR=~/foldback/results/breakinator_results

SAM="${OUT_DIR}/clean_control.sam"
REPORT="${OUT_DIR}/clean_control_breakinator_raw.tsv"
CALLS="${OUT_DIR}/calls_breakinator_clean_control.tsv"

echo "=== breakinator: clean_control (aligning, name-sorted) ==="
minimap2 -ax map-ont -t "${SLURM_CPUS_PER_TASK}" "${REF}" "${FASTQ}" > "${SAM}"

echo "=== breakinator: running ==="
breakinator -i "${SAM}" -o "${REPORT}" --tabular --no-sym --threads "${SLURM_CPUS_PER_TASK}"

python3 - "${REPORT}" "${TRUTH}" "${CALLS}" <<'PYEOF'
import csv
import sys
report_path, truth_path, calls_path = sys.argv[1], sys.argv[2], sys.argv[3]
all_ids = set()
with open(truth_path) as f:
    next(f)
    for line in f:
        all_ids.add(line.split("\t", 1)[0])
flagged_ids = set()
with open(report_path) as f:
    reader = csv.DictReader(f, delimiter="\t")
    reader.fieldnames = [fn.lstrip("#") for fn in reader.fieldnames]
    for row in reader:
        if row.get("Classification", "").strip() == "Foldback":
            flagged_ids.add(row["Read_ID"])
with open(calls_path, "w") as out:
    out.write("read_id\tmethod\tflagged\n")
    for rid in sorted(all_ids):
        out.write(f"{rid}\tbreakinator\t{rid in flagged_ids}\n")
PYEOF

echo "Wrote ${CALLS}"
