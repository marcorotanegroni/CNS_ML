#!/bin/bash
# One step of the cluster classifiers as a SLURM job. Submitted by
# code/scripts/slurm_submit_cluster_prediction.sh; arguments are passed to
# code/scripts/cluster_prediction.R and the thread count follows --cpus-per-task.
#SBATCH --cpus-per-task=32
#SBATCH --mem=64gb
#SBATCH --time=95:00:00
#SBATCH --output=logs/slurm_%x_%j.out
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
echo "$(date '+%F %T') start $SLURM_JOB_NAME on $(hostname): $*"
Rscript code/scripts/cluster_prediction.R "$@" --threads "${SLURM_CPUS_PER_TASK:-1}"
echo "$(date '+%F %T') end $SLURM_JOB_NAME"
