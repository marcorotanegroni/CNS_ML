# Predictive benchmark resampling

Status: inputs and execution workflow validated; the repeated-cohort models
have not yet been fitted. No new study F1 estimates are present in this folder.

Open `code/06_prediction_revision.Rmd` in RStudio. Its analyses run independently:

- A: `RUN_SPLIT_STABILITY <- TRUE`, using `final_matrices.RData` unchanged.
  Results go to `split_stability/`; normalized omics are not loaded.
- B: `RUN_TRAINING_FOLD <- TRUE`, additionally using
  `exp_meth_post_normalization.RData`. Results go to `training_fold/`.

Leave the other flag `FALSE` to run one analysis at a time. The Zenodo inputs
belong in the repository root and are ignored by Git. Defaults for each analysis:
20 cancer-type-stratified 80/20 splits, 10-fold inner CV with the documented
XGBoost grid, and 2,000 stratified test-patient bootstrap replicates per split.
In B, feature filtering and RMS scaling are fitted separately inside each
training fold and again on the complete outer training set for the final refit.
A preserves the original preprocessing, including its global feature selection.

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
