#!/usr/bin/env bash
#SBATCH --job-name=yacrd_clean
#SBATCH --output=logs/yacrd_clean_%j.log
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
OUT_DIR=~/foldback/results/yacrd_results

PAF="${OUT_DIR}/clean_control.paf"
REPORT="${OUT_DIR}/clean_control.yacrd"
CALLS="${OUT_DIR}/calls_yacrd_clean_control.tsv"

echo "=== yacrd: clean_control ==="
minimap2 -x ava-ont -t "${SLURM_CPUS_PER_TASK}" "${FASTQ}" "${FASTQ}" > "${PAF}"
yacrd -i "${PAF}" -o "${REPORT}"

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
