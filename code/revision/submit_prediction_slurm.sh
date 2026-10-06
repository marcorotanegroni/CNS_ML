#!/bin/bash
# Submit the revised predictive benchmark (R1 major comment 2) to SLURM.
# Run from the repository root:   bash code/revision/submit_prediction_slurm.sh
# Optional: MAIL=address for end/failure e-mails; FRAMEWORKS="Tao" for a subset;
#           CPUS (default 32) and MEM (default 64gb) per job; NODE=xen7 to run
#           every job on one node (then e.g. CPUS=20 fits the three main tunings
#           together on 64 cores).
# Per framework: tuning on split 0, then evaluation of split 0 + 20 repeated
# splits (starts automatically when that tuning succeeds); the check tuning on
# split 1 is independent. Main tunings are submitted first so that they start
# first when resources are limited. A final report job waits for all jobs.
set -euo pipefail
cd "$(dirname "$0")/../.."
mkdir -p logs
FRAMEWORKS=${FRAMEWORKS:-"Drews Steele Tao"}
opts=(--cpus-per-task="${CPUS:-32}" --mem="${MEM:-64gb}")
if [ -n "${MAIL:-}" ]; then opts+=(--mail-type=END,FAIL --mail-user="$MAIL"); fi
job="code/revision/prediction_job.slurm"
if [ -n "${NODE:-}" ]; then opts+=(--nodelist="$NODE"); fi
tune0=()
for f in $FRAMEWORKS; do
  id=$(sbatch --parsable "${opts[@]}" --job-name="tune0_$f" "$job" --framework "$f" --step tune)
  tune0+=("$id")
  echo "$f: tune split 0 = $id"
done
waits=()
i=0
for f in $FRAMEWORKS; do
  eval=$(sbatch --parsable "${opts[@]}" --job-name="eval_$f" --dependency="afterok:${tune0[$i]}" \
    "$job" --framework "$f" --step evaluate)
  echo "$f: evaluate = $eval (after ${tune0[$i]})"
  waits+=("$eval"); i=$((i + 1))
done
for f in $FRAMEWORKS; do
  tune1=$(sbatch --parsable "${opts[@]}" --job-name="tune1_$f" "$job" --framework "$f" \
    --step tune --tune_split 1)
  echo "$f: tune split 1 = $tune1"
  waits+=("$tune1")
done
deps=$(IFS=:; echo "${waits[*]}")
report_opts=(--cpus-per-task=1 --mem=8gb)
if [ -n "${NODE:-}" ]; then report_opts+=(--nodelist="$NODE"); fi
report=$(sbatch --parsable "${report_opts[@]}" --job-name=report \
  --dependency="afterany:$deps" "$job" --step report)
echo "report = $report (after all evaluation and split-1 tuning jobs)"
echo "Check with: squeue -u \$USER ; logs in logs/slurm_*.out"
