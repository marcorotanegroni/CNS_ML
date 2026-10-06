# Predictive benchmark resampling

Status: revised fixed-parameter workflow implemented and tested; results are
written to `fixed/` when run. The earlier nested design (inner tuning inside
every outer split, `code/06_prediction_revision.Rmd`) was stopped after two
Drews splits: about 31 hours per split. Its two checkpoints are kept as
evidence that tuning is not a source of variation (both splits selected
eta 0.15, max_depth 2, 1000 rounds; mean CV accuracy spanned only about two
points across the 126 grid configurations).

## Design (`code/revision/07_prediction_run.R`)

1. `tune`: once per framework, the original XGBoost grid with 10-fold
   class-stratified CV on the training part of split 0, archived predictors.
   The same tuning is repeated afterwards on split 1 as a check that the
   selected configuration does not depend on the partition
   (`fixed/tuning_comparison.csv`); evaluation uses the split-0 configuration.
2. `evaluate`: the selected configuration is fixed; one fit per outer split.
   Split 0 reconstructs the historical partition (complete archived matrix,
   `group_by(cancer_type) %>% sample_frac(.8)` after `set.seed(1234)`,
   restricted to the Cluster 1/2 cohort). Splits 1-20 are new 80/20 splits
   stratified by cancer type (seeds 20261001-20261020), so every cancer type
   with at least two patients appears in training and test. All variants use
   the same saved partitions (`fixed/manifests/`). The R1 run uses the
   archived final predictor matrices (variant `archived`). The saved
   partitions allow the R2 variants (preprocessing within the training fold,
   no purity, scalar baseline) to be run later on identical splits; they are
   implemented but not part of this run.
3. Metrics per split: Cluster 1 F1 with a 95% bootstrap interval of test
   patients within cancer type, AUROC, average precision, MCC, balanced
   accuracy, Cluster 1 prevalence and the trivial F1 2p/(1+p) of predicting
   Cluster 1 for everyone.
4. `report`: `fixed/summary.csv`, `fixed/all_split_metrics.csv`,
   `fixed/summary.png`.

## Running on the workstation

Input in the repository root (ignored by Git): `final_matrices.RData`.
Run inside `tmux`; one pipeline per framework runs in parallel (tuning, then
evaluation), followed by the report. Logs go to `logs/`:

```bash
THREADS=4 bash code/revision/run_prediction_fixed.sh
```

XGBoost uses `tree_method = "auto"`, caret's historical default (exact greedy
algorithm for this data size), so the only change from the manuscript is the
outer partition. `--tree_method hist` is about three times faster and was used
only in a local check: on the reconstructed split 0 of Drews, with the
parameters selected in the two nested splits, exact gave F1 0.933 (95% CI
0.916-0.948) and hist 0.932 (0.915-0.947), against 0.93 in the manuscript.

Each completed split is checkpointed; rerunning a command resumes. The thread
count can change between runs without invalidating checkpoints.

## Verified input cohorts

| Framework | Cluster 1 | Cluster 2 | Total | Cancer types |
| --- | ---: | ---: | ---: | ---: |
| Drews | 2,127 | 614 | 2,741 | 32 |
| Steele | 1,299 | 3,012 | 4,311 | 32 |
| Tao | 1,630 | 4,406 | 6,036 | 32 |

There are no unmatched patient IDs in the normalized input. These are distinct
framework-specific targets, with Cluster 1 positive; scores do not rank the
biological validity of the compendia. Archived clusters are retained, not refitted.

The input checksums match Zenodo record 20617051. Analysis A has 35,494 Drews,
35,298 Steele and 35,155 Tao predictors, excluding identifiers, cancer type and
the original signature target. Analysis B's pre-filter omic input has
15,976 expression and 13,674 methylation columns, with duplicate gene names
across modalities preserved as distinct variables. Each framework adds 13,614
clinical/mutation/pathway columns. In B, feature selection is learned during training.
All retained omic features in the old final matrices map to this normalized
input; proportional values were checked on 80 patients per framework. In
particular all 21,880 Drews omic features matched, supporting that the larger
counts in the old preprocessing example are a documentation discrepancy.

## Files already present: pre-filter input audit for B

- `cohort_audit.csv`: cohort sizes and ID checks.
- `input_provenance.csv`: input filenames and MD5 checksums.
- `*_feature_manifest.csv`: complete candidate feature identities and sources.
  These are pre-filter manifests, not selected predictors or model results.

## Files produced by fitting

Within each analysis folder, each framework gets its own subdirectory containing held-out predictions,
split assignments, inner-CV tuning results, selected parameters, filter counts,
bootstrap draws and F1 summaries. Checkpoints store fitted models and learned
filters and are ignored by Git. Resume requires identical data, settings,
software and code; incompatible runs require a new output directory.

`F1_across_splits.csv` describes the distribution across repeated fits.
`F1_per_split_with_bootstrap.csv` reports separate intervals conditional on each
fitted model and its test mixture. Repeated predictions of the same patient are
never pooled as independent observations. The new intervals cannot be attached
to the historical single-split F1 values without the original predictions.

Synthetic end-to-end tests passed, including train/test and inner-fold isolation,
training-only filtering, bootstrap grouping, submodel predictions and checkpoint
resume. The Rmd renders successfully with training disabled. These checks do not
establish the empirical stability of the study F1 values; that requires the full run.
