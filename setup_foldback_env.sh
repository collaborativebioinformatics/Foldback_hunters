#!/usr/bin/env bash
#SBATCH --job-name=setup_foldback_env
#SBATCH --output=logs/setup_foldback_env_%j.log
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=batch      # CHANGE to your cluster's actual partition name

# add pandas, matplotlib

set -euo pipefail

ENV_NAME="foldback"

echo "=== Loading conda ==="
module load miniconda
source "$(conda info --base)/etc/profile.d/conda.sh"

echo "=== Removing existing '${ENV_NAME}' env (if present) ==="
conda env remove -n "${ENV_NAME}" -y || echo "No existing '${ENV_NAME}' env found, continuing."

echo "=== Creating '${ENV_NAME}' env ==="
conda create -n "${ENV_NAME}" -y -c bioconda -c conda-forge \
    python=3.10 \
    samtools \
    minimap2 \
    sniffles \
    yacrd \
    pbsim3 \
    breakinator=1.1.1 \
    seqtk \
    pysam \
    jupyterlab \
    ipykernel

echo "=== Activating '${ENV_NAME}' env ==="
conda activate "${ENV_NAME}"

echo "=== Installing edlib via pip ==="
pip install edlib

echo "=== Registering Jupyter kernel (display name: foldback) ==="
python -m ipykernel install --user --name "${ENV_NAME}" --display-name "foldback"

echo "=== Verifying installed tools ==="
python --version
samtools --version | head -1
minimap2 --version
sniffles --version
yacrd --version
pbsim --version || true
breakinator --version
seqtk 2>&1 | head -3 || true
python -c "import pysam; print('pysam', pysam.__version__)"
python -c "import edlib; print('edlib OK')"
jupyter kernelspec list

echo "=== Done. Env '${ENV_NAME}' created and kernel registered. ==="
