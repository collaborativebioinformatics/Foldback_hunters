#!/usr/bin/env bash
#SBATCH --job-name=yacrd
#SBATCH --output=logs/yacrd_%A_%a.log
#SBATCH --array=1-9
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=01:00:00
#SBATCH --partition=batch

set -eo pipefail
module load miniconda/3.13
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate foldback

# --- one condition per array task ---
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
OUT_DIR=~/foldback/results/yacrd_results
mkdir -p "${OUT_DIR}"

FASTQ="${DATA_DIR}/sim_${CONDITION}.fastq"
TRUTH="${DATA_DIR}/truth_${CONDITION}.tsv"
PAF="${OUT_DIR}/${CONDITION}.paf"
REPORT="${OUT_DIR}/${CONDITION}.yacrd"
CALLS="${OUT_DIR}/calls_yacrd_${CONDITION}.tsv"

conda activate foldback

echo "=== ${CONDITION} ==="
minimap2 -x ava-ont -t "${SLURM_CPUS_PER_TASK}" "${FASTQ}" "${FASTQ}" > "${PAF}"
yacrd -i "${PAF}" -o "${REPORT}"

# yacrd's "Chimeric" label = coverage gap = does NOT fire on foldbacks
# (foldback's RC half is real sequence from the same locus, so other reads
# still overlap it). Near-zero recall here is the expected result.
python3 - "${REPORT}" "${TRUTH}" "${CALLS}" <<'PYEOF'
import sys
report_path, truth_path, calls_path = sys.argv[1], sys.argv[2], sys.argv[3]

all_ids = set()
with open(truth_path) as f:
    next(f)
    for line in f:
        all_ids.add(line.split("\t", 1)[0])

flagged = {}
with open(report_path) as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2:
            flagged[parts[1]] = (parts[0] == "Chimeric")

with open(calls_path, "w") as out:
    out.write("read_id\tmethod\tflagged\n")
    for rid in sorted(all_ids):
        out.write(f"{rid}\tyacrd\t{flagged.get(rid, False)}\n")
PYEOF

echo "Wrote ${CALLS}"
