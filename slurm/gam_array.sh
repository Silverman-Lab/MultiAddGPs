#!/bin/bash
#SBATCH --job-name=multiadd_gam
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=96GB
#SBATCH --account=one_sc_default

set -euo pipefail

module load r/4.3.2
export MULTIADDGP_R_LIBS="${MULTIADDGP_R_LIBS:-}"

cd "$RUN_DIR"

Rscript --vanilla \
  scripts/run_model_task.R \
  "$MANIFEST" \
  "$SLURM_ARRAY_TASK_ID"
