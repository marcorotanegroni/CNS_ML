#!/bin/bash
# Submit the revised predictive benchmark (R1 major comment 2) to SLURM.
# Run from the repository root:   bash code/revision/submit_prediction_slurm.sh
# Optional: MAIL=address for end/failure e-mails; FRAMEWORKS="Tao" for a subset;
#           CPUS (default 32) and MEM (default 64gb) per job.
# Per framework: tuning on split 0, then evaluation of split 0 + 20 repeated
# splits (starts automatically when that tuning succeeds); the check tuning on
# split 1 runs in parallel. A final report job waits for all of them.
set -euo pipefail
cd "$(dirname "$0")/../.."
mkdir -p logs
FRAMEWORKS=${FRAMEWORKS:-"Drews Steele Tao"}
opts=(--cpus-per-task="${CPUS:-32}" --mem="${MEM:-64gb}")
if [ -n "${MAIL:-}" ]; then opts+=(--mail-type=END,FAIL --mail-user="$MAIL"); fi
job="code/revision/prediction_job.slurm"
waits=()
for f in $FRAMEWORKS; do
  tune0=$(sbatch --parsable "${opts[@]}" --job-name="tune0_$f" "$job" --framework "$f" --step tune)
  eval=$(sbatch --parsable "${opts[@]}" --job-name="eval_$f" --dependency="afterok:$tune0" \
    "$job" --framework "$f" --step evaluate)
  tune1=$(sbatch --parsable "${opts[@]}" --job-name="tune1_$f" "$job" --framework "$f" \
    --step tune --tune_split 1)
  echo "$f: tune split 0 = $tune0, evaluate = $eval (after $tune0), tune split 1 = $tune1"
  waits+=("$eval" "$tune1")
done
deps=$(IFS=:; echo "${waits[*]}")
report=$(sbatch --parsable --cpus-per-task=1 --mem=8gb --job-name=report \
  --dependency="afterany:$deps" "$job" --step report)
echo "report = $report (after all evaluation and split-1 tuning jobs)"
echo "Check with: squeue -u \$USER ; logs in logs/slurm_*.out"
