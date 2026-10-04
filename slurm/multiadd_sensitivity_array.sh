#!/bin/bash
#SBATCH --job-name=multiadd_sens
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64GB
#SBATCH --account=one_sc_default

set -euo pipefail

module load r/4.3.2
export MULTIADDGP_R_LIBS="${MULTIADDGP_R_LIBS:-}"

cd "$RUN_DIR"

Rscript --vanilla \
  scripts/multiadd_sensitivity.R \
  fit \
  "$RUN_DIR" \
  "$MANIFEST" \
  "$SLURM_ARRAY_TASK_ID"
