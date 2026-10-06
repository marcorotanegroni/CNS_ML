#!/usr/bin/env bash
# Revised predictive benchmark (R1 major comment 2) on the workstation.
# Run from the repository root inside tmux:  bash code/revision/run_prediction_fixed.sh
# One pipeline per framework in parallel: one-off tuning on split 0,
# fixed-parameter evaluation of split 0 + 20 repeated splits, then a check
# tuning on split 1. Logs in logs/.
# Rerunning resumes from completed tuning and checkpointed splits.
# FRAMEWORKS="Drews" restricts the run, e.g. to relaunch one framework.
set -u
THREADS=${THREADS:-4}
mkdir -p logs
FRAMEWORKS=${FRAMEWORKS:-"Drews Steele Tao"}
for f in $FRAMEWORKS; do
  (
    log() { echo "$(date '+%F %T') $f $1" >> logs/status.log; }
    if ! Rscript code/revision/07_prediction_run.R --framework "$f" --step tune \
        --threads "$THREADS" > "logs/tune_$f.log" 2>&1; then
      log "tuning FAILED (see logs/tune_$f.log)"; exit 1
    fi
    log "tuning finished"
    if ! Rscript code/revision/07_prediction_run.R --framework "$f" --step evaluate \
        --threads "$THREADS" > "logs/evaluate_$f.log" 2>&1; then
      log "evaluation FAILED (see logs/evaluate_$f.log)"; exit 1
    fi
    log "evaluation finished"
    # Robustness check, after the main results: repeat tuning on split 1.
    if Rscript code/revision/07_prediction_run.R --framework "$f" --step tune --tune_split 1 \
        --threads "$THREADS" > "logs/tune_split1_$f.log" 2>&1; then
      log "split-1 tuning finished"
    else
      log "split-1 tuning FAILED (see logs/tune_split1_$f.log)"
    fi
  ) &
done
wait
Rscript code/revision/07_prediction_run.R --step report > logs/report.log 2>&1
echo "$(date '+%F %T') report written (exit code $?)" >> logs/status.log
