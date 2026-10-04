#!/bin/bash
#SBATCH --job-name=multiadd_lgpr
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=96GB
#SBATCH --account=one_sc_default

set -euo pipefail

if [[ "${SLURM_CPUS_PER_TASK:-0}" -ne 4 ]]; then
  echo "lgpr requires exactly 4 allocated CPUs" >&2
  exit 2
fi

module load r/4.3.2
export MULTIADDGP_R_LIBS="${MULTIADDGP_R_LIBS:-}"

cd "$RUN_DIR"

Rscript --vanilla \
  scripts/run_model_task.R \
  "$MANIFEST" \
  "$SLURM_ARRAY_TASK_ID"
