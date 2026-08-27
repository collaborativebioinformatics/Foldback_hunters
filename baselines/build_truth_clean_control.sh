#!/usr/bin/env bash
# build_truth_clean_control.sh
# Builds ~/foldback/sim_data/truth_clean_control.tsv from the clean control
# FASTQ. Every read is_foldback=False. Idempotent -- skips if already built.
set -eo pipefail


FASTQ=~/foldback/PSIM-simulation/Test_sample_10k.fastq.gz
MAF=~/foldback/PSIM-simulation/Test_sample_10k.maf.gz   # optional
TRUTH=~/foldback/sim_data/truth_clean_control.tsv


if [ -f "${TRUTH}" ]; then
    echo "${TRUTH} already exists, skipping."
    exit 0
fi

python3 - "${FASTQ}" "${MAF}" "${TRUTH}" <<'PYEOF'
import gzip
import sys

fastq_path, maf_path, truth_path = sys.argv[1], sys.argv[2], sys.argv[3]

def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def read_ids_from_fastq(path):
    ids = []
    with open_maybe_gz(path) as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                ids.append(line.rstrip("\n").lstrip("@").split()[0])
    return ids

def parse_maf_loci(path):
    loci = {}
    try:
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
    except FileNotFoundError:
        pass
    return loci

read_ids = read_ids_from_fastq(fastq_path)
loci = parse_maf_loci(maf_path)

with open(truth_path, "w") as out:
    out.write("read_id\tis_foldback\tfold_position\tadapter_present\tsource_locus\n")
    for rid in read_ids:
        out.write(f"{rid}\tFalse\t\t\t{loci.get(rid, 'unknown')}\n")

print(f"Wrote {len(read_ids)} clean reads to {truth_path} "
      f"({len(loci)} had a MAF-derived source_locus)")
PYEOF
