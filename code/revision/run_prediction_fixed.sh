#!/usr/bin/env bash
# Revised predictive benchmark (R1 major comment 2) on the workstation.
# Run from the repository root inside tmux:  bash code/revision/run_prediction_fixed.sh
# One pipeline per framework in parallel: one-off tuning on split 0,
# fixed-parameter evaluation of split 0 + 20 repeated splits, then a check
# tuning on split 1. Logs in logs/.
# Rerunning resumes from completed tuning and checkpointed splits.
set -u
THREADS=${THREADS:-4}
mkdir -p logs
for f in Drews Steele Tao; do
  (
    Rscript code/revision/07_prediction_run.R --framework "$f" --step tune --threads "$THREADS" \
      > "logs/tune_$f.log" 2>&1 &&
    Rscript code/revision/07_prediction_run.R --framework "$f" --step evaluate --threads "$THREADS" \
      > "logs/evaluate_$f.log" 2>&1
    echo "$(date '+%F %T') $f evaluation finished with exit code $?" >> logs/status.log
    # Robustness check, after the main results: repeat tuning on split 1.
    Rscript code/revision/07_prediction_run.R --framework "$f" --step tune --tune_split 1 \
      --threads "$THREADS" > "logs/tune_split1_$f.log" 2>&1
    echo "$(date '+%F %T') $f split-1 tuning finished with exit code $?" >> logs/status.log
  ) &
done
wait
Rscript code/revision/07_prediction_run.R --step report > logs/report.log 2>&1
echo "$(date '+%F %T') report written (exit code $?)" >> logs/status.log
