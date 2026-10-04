#!/bin/bash
#SBATCH --job-name=multiadd_data
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --account=one_sc_default

set -euo pipefail

module load r/4.3.2
export MULTIADDGP_R_LIBS="${MULTIADDGP_R_LIBS:-}"

cd "$RUN_DIR"

Rscript --vanilla \
  scripts/generate_dataset_task.R \
  "$MANIFEST" \
  "$SLURM_ARRAY_TASK_ID"
