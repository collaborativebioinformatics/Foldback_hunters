#!/usr/bin/env bash
#SBATCH --job-name=sniffles2
#SBATCH --output=sniffles2_%A_%a.log
#SBATCH --array=1-18
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=05:00:00
#SBATCH --partition=batch

set -euo pipefail
module load miniconda/3.13
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate foldback

# 18 tasks = 9 conditions x 2 filter variants (with_filter / no_filter).
# Each task does its OWN alignment rather than sharing a BAM between the two
# filter-variant tasks for the same condition -- two array tasks writing the
# same BAM file in parallel would race. This costs one extra alignment pass
# per condition versus the 9-job version, but avoids that race entirely and
# lets both variants run fully in parallel instead of one waiting on the
# other.
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
FILTER_TAGS=(with_filter no_filter)

IDX=$((SLURM_ARRAY_TASK_ID - 1))
CONDITION="${CONDITIONS[$((IDX / 2))]}"
FILTER_TAG="${FILTER_TAGS[$((IDX % 2))]}"

DATA_DIR=~/foldback/sim_data
OUT_DIR=~/foldback/results/sniffles2_${FILTER_TAG}
REF=~/foldback/GRCh38.fasta
mkdir -p "${OUT_DIR}"

FASTQ="${DATA_DIR}/sim_${CONDITION}.fastq"
TRUTH="${DATA_DIR}/truth_${CONDITION}.tsv"
BAM="${OUT_DIR}/${CONDITION}_${FILTER_TAG}.sorted.bam"
VCF="${OUT_DIR}/${CONDITION}_${FILTER_TAG}.vcf"
CALLS="${OUT_DIR}/calls_sniffles2_${FILTER_TAG}_${CONDITION}.tsv"


echo "=== ${CONDITION} / ${FILTER_TAG}: aligning ==="
minimap2 -ax map-ont -Y --MD -t "${SLURM_CPUS_PER_TASK}" "${REF}" "${FASTQ}" \
    | samtools sort -@ "${SLURM_CPUS_PER_TASK}" -o "${BAM}" -
samtools index "${BAM}"

echo "=== ${CONDITION} / ${FILTER_TAG}: sniffles2 ==="
if [ "${FILTER_TAG}" = "with_filter" ]; then
    sniffles --input "${BAM}" --vcf "${VCF}" --output-rnames --threads "${SLURM_CPUS_PER_TASK}"
else
    sniffles --input "${BAM}" --vcf "${VCF}" --output-rnames --minsvlen 0 --threads "${SLURM_CPUS_PER_TASK}"
fi

echo "=== ${CONDITION} / ${FILTER_TAG}: converting to per-read call table ==="
python3 - "${VCF}" "${TRUTH}" "${CALLS}" "sniffles2_${FILTER_TAG}" <<'PYEOF'
import re
import sys

vcf_path, truth_path, calls_path, method = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

all_ids = set()
with open(truth_path) as f:
    next(f)
    for line in f:
        all_ids.add(line.split("\t", 1)[0])

# any read named in any SV call's RNAMES is flagged True
supporting_reads = set()
with open(vcf_path) as f:
    for line in f:
        if line.startswith("#"):
            continue
        info = line.rstrip("\n").split("\t")[7]
        m = re.search(r"RNAMES=([^;]+)", info)
        if m:
            supporting_reads.update(m.group(1).split(","))

with open(calls_path, "w") as out:
    out.write("read_id\tmethod\tflagged\n")
    for rid in sorted(all_ids):
        out.write(f"{rid}\t{method}\t{rid in supporting_reads}\n")
PYEOF

echo "Wrote ${CALLS}"
