#!/bin/bash
# One step of the revised predictive benchmark as a SLURM job. Submitted by
# code/revision/submit_prediction_slurm.sh; arguments are passed to
# 07_prediction_run.R and the thread count follows --cpus-per-task.
#SBATCH --cpus-per-task=32
#SBATCH --mem=64gb
#SBATCH --time=95:00:00
#SBATCH --output=logs/slurm_%x_%j.out
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
echo "$(date '+%F %T') start $SLURM_JOB_NAME on $(hostname): $*"
Rscript code/revision/07_prediction_run.R "$@" --threads "${SLURM_CPUS_PER_TASK:-1}"
echo "$(date '+%F %T') end $SLURM_JOB_NAME"
