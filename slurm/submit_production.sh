#!/bin/bash
set -euo pipefail

RUN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST_DIR="$RUN_DIR/manifests"
LOG_DIR="$RUN_DIR/logs"
APPROVAL_FILE="$RUN_DIR/config/PILOTS_APPROVED"

if [[ ! -f "$APPROVAL_FILE" ]]; then
  echo "Review pilots, then create $APPROVAL_FILE before production." >&2
  exit 1
fi

array_size() {
  local lines
  lines=$(wc -l < "$1")
  echo $((lines - 1))
}

submit_array() {
  local manifest="$1" concurrency="$2" script="$3" dependency="${4:-}" extra_export="${5:-}"
  local count
  count=$(array_size "$manifest")
  local export_vars="ALL,RUN_DIR=$RUN_DIR,MANIFEST=$manifest"
  if [[ -n "$extra_export" ]]; then
    export_vars="$export_vars,$extra_export"
  fi
  local command=(
    sbatch --parsable
    --array="1-${count}%${concurrency}"
    --output="$LOG_DIR/%x_%A_%a.out"
    --error="$LOG_DIR/%x_%A_%a.err"
    --export="$export_vars"
  )
  if [[ -n "$dependency" ]]; then
    command+=(--dependency="afterok:$dependency")
  fi
  command+=("$script")
  "${command[@]}"
}

mkdir -p "$LOG_DIR"
data_job=$(submit_array "$MANIFEST_DIR/datasets.csv" 20 "$RUN_DIR/slurm/dataset_array.sh")
echo "Datasets: $data_job"
echo "GAM: $(submit_array "$MANIFEST_DIR/gam.csv" 20 "$RUN_DIR/slurm/gam_array.sh" "$data_job")"
echo "MultiAddGP: $(submit_array "$MANIFEST_DIR/multiadd.csv" 10 "$RUN_DIR/slurm/multiadd_array.sh" "$data_job")"
echo "Gaussian lgpr: $(submit_array "$MANIFEST_DIR/lgpr_gaussian.csv" 5 "$RUN_DIR/slurm/lgpr_array.sh" "$data_job")"
echo "NB lgpr: $(submit_array "$MANIFEST_DIR/lgpr_nb.csv" 2 "$RUN_DIR/slurm/lgpr_array.sh" "$data_job" "LGPR_ADAPT_DELTA=0.9999,LGPR_MAX_TREEDEPTH=11")"
