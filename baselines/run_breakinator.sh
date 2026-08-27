#!/usr/bin/env bash
#SBATCH --job-name=breakinator
#SBATCH --output=logs/breakinator_%A_%a.log
#SBATCH --array=1-9
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=05:00:00
#SBATCH --partition=batch

set -eo pipefail
module load miniconda/3.13
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate foldback

CONDITIONS=(
    1pct_middle
    1pct_quarter
    1pct_near_end
    5pct_middle
    5pct_quarter
    5pct_near_end
    10pct_middle
    10pct_quarter
    10pct_near_end
)
CONDITION="${CONDITIONS[$((SLURM_ARRAY_TASK_ID - 1))]}"

DATA_DIR=~/foldback/sim_data
OUT_DIR=~/foldback/results/breakinator_results
REF=~/foldback/GRCh38.fasta
mkdir -p "${OUT_DIR}"

FASTQ="${DATA_DIR}/sim_${CONDITION}.fastq"
TRUTH="${DATA_DIR}/truth_${CONDITION}.tsv"
SAM="${OUT_DIR}/${CONDITION}.sam"
REPORT="${OUT_DIR}/${CONDITION}_breakinator_raw.tsv"
CALLS="${OUT_DIR}/calls_breakinator_${CONDITION}.tsv"


echo "=== ${CONDITION}: aligning (name-sorted, minimap2's natural order -- do NOT coordinate-sort) ==="
minimap2 -ax map-ont -t "${SLURM_CPUS_PER_TASK}" "${REF}" "${FASTQ}" > "${SAM}"

echo "=== ${CONDITION}: running breakinator ==="
breakinator -i "${SAM}" -o "${REPORT}" --tabular --no-sym --threads "${SLURM_CPUS_PER_TASK}"

echo "=== ${CONDITION}: converting to per-read call table ==="
# Confirmed real breakinator --tabular columns:
#   #Break1_chr  Break1_loc  Break_direction  Break2_chr  Break2_loc  MapQ  Read_ID  Classification
# Classification is "Foldback", "Chimeric", or presumably absent/other for
# clean reads. Only "Foldback" counts as flagged=True here (not "Chimeric").
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
    # normalize header: strip a leading '#' from the first column name if present
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
