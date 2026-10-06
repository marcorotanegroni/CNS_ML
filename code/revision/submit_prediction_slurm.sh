#!/bin/bash
# Submit the revised predictive benchmark (R1 major comment 2) to SLURM.
# Run from the repository root:
#   NODES="xen7:64:40 xen5:38:40 xen3:32:28" bash code/revision/submit_prediction_slurm.sh
#
# Each inner-CV fold of each tuning is a separate job, checkpointed on
# completion; an assembly job then selects the configuration from the saved
# folds. Per framework:
#   tune0 folds 1-10 -> assemble0 -> evaluation (split 0 + 20 repeated splits)
#   tune1 folds 1-10 -> assemble1   (check tuning on split 1)
# Split-0 fold jobs are submitted first, longest framework first. A final
# report job waits for everything. Resubmitting skips completed folds/splits.
#
# Options (environment variables):
#   NODES       "name:cores[:memGB] ..." pins jobs to these nodes, in proportion
#               to their cores; the optional memGB is the memory per fold job
#               on that node (default: SLURM chooses nodes, MEM per job)
#   CPUS        CPUs per fold job (default 8); MEM memory per job (default 24gb)
#   CPUS_EVAL   CPUs for each evaluation job (default 16)
#   FRAMEWORKS  subset, e.g. "Tao" (default "Tao Steele Drews")
#   MAIL        address for end/failure e-mails
set -euo pipefail
cd "$(dirname "$0")/../.."
mkdir -p logs
FRAMEWORKS=${FRAMEWORKS:-"Tao Steele Drews"}
CPUS=${CPUS:-8}; MEM=${MEM:-24gb}; CPUS_EVAL=${CPUS_EVAL:-16}; FOLDS=10
job="code/revision/prediction_job.slurm"
mail=(); if [ -n "${MAIL:-}" ]; then mail=(--mail-type=FAIL --mail-user="$MAIL"); fi

# Node sequence: each node appears once per slot (cores / CPUS), interleaved,
# so consecutive jobs are spread across nodes in proportion to their size.
slots=(); slot_mem=()
if [ -n "${NODES:-}" ]; then
  while read -r name mem; do slots+=("$name"); slot_mem+=("$mem"); done < <(
    for spec in $NODES; do echo "$spec"; done | awk -F: -v cpus="$CPUS" -v mem="$MEM" '{
      n = int($2 / cpus); if (n < 1) n = 1
      m = ($3 != "") ? $3 "gb" : mem
      for (i = 0; i < n; i++) printf "%.6f %s %s\n", (i + 0.5) / n, $1, m }' |
      sort -n | awk '{print $2, $3}')
  first_node=${slots[0]}
fi
next=0
node_opt() {  # sets "node" (nodelist option) and "job_mem" for the next slot
  node=(); job_mem=$MEM
  if [ ${#slots[@]} -gt 0 ]; then
    node=(--nodelist="${slots[$next]}"); job_mem=${slot_mem[$next]}
    next=$(( (next + 1) % ${#slots[@]} ))
  fi
}
main_node=(); main_mem=$MEM
if [ -n "${first_node:-}" ]; then main_node=(--nodelist="$first_node"); main_mem=${slot_mem[0]}; fi

submit_folds() {  # $1 framework, $2 tuning split -> sets fold_ids (colon-joined)
  # Runs in the current shell (not $(...)) so the node counter advances.
  local f=$1 split=$2 k id ids=()
  for k in $(seq 1 $FOLDS); do
    node_opt
    id=$(sbatch --parsable --cpus-per-task="$CPUS" --mem="$job_mem" ${node[@]+"${node[@]}"} ${mail[@]+"${mail[@]}"} \
      --job-name="tune${split}_${f}_f$k" "$job" --framework "$f" --step tune \
      --tune_split "$split" --folds "$k")
    ids+=("$id")
  done
  fold_ids=$(IFS=:; echo "${ids[*]}")
}

waits=()
assemble0=()
for f in $FRAMEWORKS; do
  submit_folds "$f" 0; folds=$fold_ids
  a=$(sbatch --parsable --cpus-per-task=2 --mem="$main_mem" ${main_node[@]+"${main_node[@]}"} ${mail[@]+"${mail[@]}"} \
    --job-name="assemble0_$f" --dependency="afterok:$folds" "$job" --framework "$f" --step tune)
  assemble0+=("$a")
  echo "$f: split-0 fold jobs $folds; assembly $a"
done
i=0
for f in $FRAMEWORKS; do
  e=$(sbatch --parsable --cpus-per-task="$CPUS_EVAL" --mem="$main_mem" ${main_node[@]+"${main_node[@]}"} ${mail[@]+"${mail[@]}"} \
    --job-name="eval_$f" --dependency="afterok:${assemble0[$i]}" "$job" --framework "$f" --step evaluate)
  echo "$f: evaluation $e (after ${assemble0[$i]})"
  waits+=("$e"); i=$((i + 1))
done
for f in $FRAMEWORKS; do
  submit_folds "$f" 1; folds=$fold_ids
  a=$(sbatch --parsable --cpus-per-task=2 --mem="$main_mem" ${main_node[@]+"${main_node[@]}"} ${mail[@]+"${mail[@]}"} \
    --job-name="assemble1_$f" --dependency="afterok:$folds" "$job" --framework "$f" \
    --step tune --tune_split 1)
  echo "$f: split-1 fold jobs $folds; assembly $a"
  waits+=("$a")
done
deps=$(IFS=:; echo "${waits[*]}")
report=$(sbatch --parsable --cpus-per-task=1 --mem=8gb ${main_node[@]+"${main_node[@]}"} --job-name=report \
  --dependency="afterany:$deps" "$job" --step report)
echo "report $report (after all evaluations and split-1 assemblies)"
echo "Check with: squeue -u \$USER ; logs in logs/slurm_*.out"
