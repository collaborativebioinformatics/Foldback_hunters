#!/usr/bin/env bash
#SBATCH --job-name=sniffles2_with_filter_clean
#SBATCH --output=logs/sniffles2_with_filter_clean_%j.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=01:00:00
#SBATCH --partition=batch

set -eo pipefail

module load miniconda/3.13
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate foldback

FASTQ=~/foldback/PSIM-simulation/Test_sample_10k.fastq.gz
TRUTH=~/foldback/sim_data/truth_clean_control.tsv
REF=~/foldback/GRCh38.fasta
OUT_DIR=~/foldback/results/sniffles2_with_filter

BAM="${OUT_DIR}/clean_control_with_filter.sorted.bam"
VCF="${OUT_DIR}/clean_control_with_filter.vcf"
CALLS="${OUT_DIR}/calls_sniffles2_with_filter_clean_control.tsv"

echo "=== sniffles2_with_filter: clean_control (aligning) ==="
minimap2 -ax map-ont -Y --MD -t "${SLURM_CPUS_PER_TASK}" "${REF}" "${FASTQ}" \
    | samtools sort -@ "${SLURM_CPUS_PER_TASK}" -o "${BAM}" -
samtools index "${BAM}"

echo "=== sniffles2_with_filter: running (default settings) ==="
sniffles --input "${BAM}" --vcf "${VCF}" --output-rnames --threads "${SLURM_CPUS_PER_TASK}"

python3 - "${VCF}" "${TRUTH}" "${CALLS}" "sniffles2_with_filter" <<'PYEOF'
import re
import sys
vcf_path, truth_path, calls_path, method = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
all_ids = set()
with open(truth_path) as f:
    next(f)
    for line in f:
        all_ids.add(line.split("\t", 1)[0])
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
