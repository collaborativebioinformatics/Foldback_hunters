#!/usr/bin/env bash
#SBATCH --job-name=framing_a
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=05:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --partition=batch

set -eo pipefail

module load miniconda/3.13
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate foldback

python compute_framing_a.py
