#!/bin/bash
# Create and submit the MultiAddGP hyperparameter sensitivity array.
set -euo pipefail

if [[ $# -gt 3 ]]; then
  echo "Usage: $0 [D] [N] [SEED_COUNT]" >&2
  exit 2
fi

RUN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
D="${1:-3}"
N="${2:-50}"
SEED_COUNT="${3:-10}"
MANIFEST="$RUN_DIR/manifests/multiadd_sensitivity_D${D}_N${N}_${SEED_COUNT}seed.csv"
LOG_DIR="$RUN_DIR/logs"

module load r/4.3.2
export MULTIADDGP_R_LIBS="${MULTIADDGP_R_LIBS:-}"
mkdir -p "$LOG_DIR"
cd "$RUN_DIR"

# Build tasks from completed baseline MultiAddGP results.
Rscript --vanilla scripts/multiadd_sensitivity.R prepare \
  "$RUN_DIR" "$D" "$N" "$SEED_COUNT"

count=$(($(wc -l < "$MANIFEST") - 1))
if [[ "$count" -le 0 ]]; then
  echo "No sensitivity tasks were created: $MANIFEST" >&2
  exit 1
fi

# Limit concurrent fits because each task allocates substantial memory.
job_id=$(sbatch --parsable \
  --array="1-${count}%10" \
  --output="$LOG_DIR/sensitivity_D${D}_N${N}_%A_%a.out" \
  --error="$LOG_DIR/sensitivity_D${D}_N${N}_%A_%a.err" \
  --export="ALL,RUN_DIR=$RUN_DIR,MANIFEST=$MANIFEST" \
  "$RUN_DIR/slurm/multiadd_sensitivity_array.sh")

echo "Sensitivity manifest: $MANIFEST"
echo "Sensitivity tasks: $count"
echo "Submitted sensitivity job: $job_id"
echo "After completion, run:"
echo "  Rscript --vanilla scripts/multiadd_sensitivity.R report $RUN_DIR $D $N $SEED_COUNT"
