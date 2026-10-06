# Repeated held-out cluster classification

Status: workflow implemented and tested; results are written to `fixed/` by
the workstation run.

## Design (`code/revision/07_prediction_run.R`)

1. `tune`: once per framework, the original XGBoost grid with 10-fold
   class-stratified CV on the training part of split 0, with the archived
   final predictor matrices. The same tuning is repeated at the end on split 1
   as a check that the selected configuration does not depend on the
   partition (`fixed/tuning_comparison.csv`); evaluation always uses the
   split-0 configuration.
2. `evaluate`: the selected configuration is held fixed; one fit per outer
   split. Split 0 reconstructs the original partition with the original code
   and seed (complete archived matrix, `group_by(cancer_type) %>%
   sample_frac(.8)` after `set.seed(1234)`, restricted to the Cluster 1/2
   cohort). Splits 1-20 are new 80/20 splits stratified by cancer type (seeds
   20261001-20261020), so every cancer type with at least two patients appears
   in training and test. Partitions are saved in `fixed/manifests/` and reused
   by every run.
3. Metrics per split: Cluster 1 F1 with a 95% percentile bootstrap interval
   (2,000 resamples of test patients within cancer type, conditional on the
   fitted model), plus AUROC, average precision, MCC, balanced accuracy and
   Cluster 1 prevalence for reference. Across-split summaries describe the 20
   repeated splits; split 0 is reported separately.
4. `report`: `fixed/summary.csv`, `fixed/all_split_metrics.csv`,
   `fixed/tuning_comparison.csv`, `fixed/f1_repeated_splits.png/.pdf` (F1 only,
   for the manuscript) and `fixed/exploratory_metrics.png` (all metrics, not
   for the manuscript).

XGBoost uses `tree_method = "hist"`, set explicitly: the meaning of the
default `"auto"` changed in XGBoost 2.0 (exact greedy before, hist after), and
the version used for the manuscript was not recorded. On the reconstructed
split 0 of Drews, with the configuration selected when tuning was repeated
within two outer splits (eta 0.15, max_depth 2, 1,000 rounds), hist gave F1
0.932 (95% CI 0.915-0.947) and exact 0.933 (0.916-0.948), against 0.93 in the
manuscript; `--tree_method exact` is available. The reconstructed test set has
567 patients, including 100 BRCA as in Figure 8d.

Variants for the Reviewer 2 comments (filters refitted within the training
fold, no purity, scalar baseline) are implemented in the same script
(`--variants`) and use the same saved partitions; they are not part of this
run.

## Running

Input: `final_matrices.RData` (Zenodo record 20617051) in
`data/processed/zenodo/` or the repository root.

On the SLURM cluster (recommended), from the repository root:

```bash
bash code/revision/submit_prediction_slurm.sh
```

Per framework this submits the tuning on split 0, the evaluation (started
automatically when that tuning succeeds) and the check tuning on split 1 (in
parallel), plus a final report job. Each job uses 32 CPUs and 64 GB by default
(`CPUS`, `MEM`); `NODE=xen7` keeps every job on one node (with `CPUS=20` the
three main tunings run together on 64 cores, and are submitted first),
`NODE_CHECK=xen5 CPUS_CHECK=12` runs the split-1 check tunings on a second
node,
`FRAMEWORKS="Tao"` submits a subset and `MAIL=address` adds end/failure
e-mails. Monitor with `squeue -u $USER`; logs are in
`logs/slurm_*.out`, with a timestamp on every inner CV fold.

Without SLURM, `THREADS=4 bash code/revision/run_prediction_fixed.sh` runs the
same steps in parallel processes (logs in `logs/`, progress in
`logs/status.log`).

Each completed split is checkpointed and rerunning resumes. Checkpoints are
tied to the data, settings and code; the thread count can change between runs.
Saved tuning is reused only if training patients, folds, seed and tree method
match. XGBoost data are freed after every fit, so memory stays at a few GB per
framework.

## Input cohorts

| Framework | Cluster 1 | Cluster 2 | Total | Cancer types |
| --- | ---: | ---: | ---: | ---: |
| Drews | 2,127 | 614 | 2,741 | 32 |
| Steele | 1,299 | 3,012 | 4,311 | 32 |
| Tao | 1,630 | 4,406 | 6,036 | 32 |

These are distinct framework-specific targets with Cluster 1 positive; scores
do not rank the biological validity of the compendia. Archived cluster labels
are used without refitting. Predictors: 35,494 (Drews), 35,298 (Steele) and
35,155 (Tao), excluding identifiers, cancer type and the signature target.
The input checksum is verified against the Zenodo record.

## Output files per framework (`fixed/archived/<framework>/`)

`predictions.csv` (held-out probabilities and labels), `split_metrics.csv`,
`split_summary.csv`, `bootstrap.csv`, `manifests.csv`, `best_parameters.csv`,
`filters.csv`, `evaluation_notes.txt`. Fitted models are stored in ignored
checkpoints. Repeated predictions of the same patient across splits are never
pooled.
