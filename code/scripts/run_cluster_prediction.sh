#!/usr/bin/env bash
# Cluster 1 versus Cluster 2 classifiers of Figure 8 and Supplementary Figure 9,
# on a workstation. Run from the repository root (e.g. inside tmux):
#   THREADS=4 bash code/scripts/run_cluster_prediction.sh
# One pipeline per framework in parallel: tuning on split 0, fixed-parameter
# evaluation of split 0 + 20 repeated splits for the four predictor sets, and
# a check tuning on split 1. Then the FGA-only classifiers and the summary
# tables. Logs in logs/. Rerunning resumes from completed tuning and
# checkpointed splits. FRAMEWORKS="Drews" restricts the run to one framework.
set -u
THREADS=${THREADS:-4}
FRAMEWORKS=${FRAMEWORKS:-"Drews Steele Tao"}
VARIANTS=${VARIANTS:-"multiomic,expression_only,multiomic_no_purity,training_fold"}
run="code/scripts/cluster_prediction.R"
mkdir -p logs
for f in $FRAMEWORKS; do
  (
    log() { echo "$(date '+%F %T') $f $1" >> logs/status.log; }
    if ! Rscript "$run" --framework "$f" --step tune \
        --threads "$THREADS" > "logs/tune_$f.log" 2>&1; then
      log "tuning FAILED (see logs/tune_$f.log)"; exit 1
    fi
    log "tuning finished"
    if ! Rscript "$run" --framework "$f" --step evaluate --variants "$VARIANTS" \
        --threads "$THREADS" > "logs/evaluate_$f.log" 2>&1; then
      log "evaluation FAILED (see logs/evaluate_$f.log)"; exit 1
    fi
    log "evaluation finished"
    # Check of the tuning: repeat it on the training part of split 1.
    if Rscript "$run" --framework "$f" --step tune --tune_split 1 \
        --threads "$THREADS" > "logs/tune_split1_$f.log" 2>&1; then
      log "split-1 tuning finished"
    else
      log "split-1 tuning FAILED (see logs/tune_split1_$f.log)"
    fi
  ) &
done
wait
# FGA-only classifiers on the same patients and splits (notebook 06).
Rscript code/scripts/fga_prediction.R $FRAMEWORKS --threads "$THREADS" > logs/fga_prediction.log 2>&1
echo "$(date '+%F %T') FGA-only classifiers finished (exit code $?)" >> logs/status.log
Rscript "$run" --step report > logs/report.log 2>&1
echo "$(date '+%F %T') report written (exit code $?)" >> logs/status.log
