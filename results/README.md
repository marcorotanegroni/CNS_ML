# Results

Tables produced by the notebooks and scripts in `code/`. Figure files are
written to `results/figures/` and are not tracked.

## `binary_concordance/` (notebook `01`, `code/R/binary_concordance.R`)

Binary agreement between the 58 signatures on the 5,881 matched TCGA tumors
under four activity rules: the zero/median rule of Figure 4 (`original`) and
the 25th, 50th, and 75th percentiles of the positive exposures of each
signature in its complete compendium (`positive_q25`, `positive_q50`,
`positive_q75`).

- `pair_metrics.csv`: one row per signature pair and rule. `n11`, `n00`,
  `n10`, `n01`: tumors active in both signatures, inactive in both, and active
  in only one; `state_jaccard`: patient-state Jaccard
  (n11 + n00) / (n11 + n00 + 2 (n10 + n01)), as in Figure 4;
  `active_share_of_agreements` and `inactive_share_of_agreements`:
  n11 / (n11 + n00) and n00 / (n11 + n00); `active_jaccard`:
  n11 / (n11 + n10 + n01); `cross_framework`: signatures from different
  compendia; `high_activity_pair`: both among the 12 signatures used for
  clustering; `active_prevalence1/2`: share of tumors classified as active.
  `expected_n11`, `observed_expected_n11`, and `cohen_kappa` are given for
  reference.
- `thresholds.csv`: threshold of each signature under each rule, the cohort on
  which it was computed, and the resulting number and share of active tumors.
- `principal_pairs.csv` and `principal_pair_ranks.csv`: the pairs discussed in
  the text, with their rank among the 48 cross-compendium pairs of the 12
  high-activity signatures (1 = most concordant).
- `provenance.txt`, `session_info.txt`: input checksum, conventions, and
  software versions.

## `clustering/` (notebook `02`)

- `cluster_summary.csv`: cluster sizes and mean signature activities
  (`code/scripts/figure5.R`).
- `supplementary_figure3_counts.csv`: counts shown in Supplementary Figure 3
  (`code/scripts/supplementary_figure3.R`).

## `survival/` (notebook `03`)

- `supplementary_figure5_numbers_at_risk.csv`: numbers at risk of the
  Kaplan-Meier curves of Figure 6 (`code/scripts/supplementary_figure5.R`).

## `cluster_prediction/` (notebook `04`)

Cluster 1 versus Cluster 2 classifiers, on the original train/test split
(split 0, seed 1234, notebook `00`) and on 20 repeated cancer-type-stratified
80/20 splits, with the hyperparameters selected on the split-0 training set.

| Folder | Predictors |
| --- | --- |
| `multiomic/` | All predictors of the final matrices (Figure 8) |
| `expression_only/` | Gene expression predictors only |
| `multiomic_no_purity/` | All predictors except tumor purity |
| `training_fold/` | Expression and methylation filters refitted on the training set of each split |
| `fga_only/` | Fraction genome altered only (notebook `06`) |

Each folder has one subfolder per compendium with `split_metrics.csv` (test-set
counts, F1 with its 95% bootstrap interval for split 0, AUROC, average
precision, MCC, balanced accuracy, Cluster 1 prevalence), `split_summary.csv`
(summary across the 20 repeated splits), `best_parameters.csv`, `filters.csv`,
and `evaluation_notes.txt`. The patient-level files `manifests.csv`,
`predictions.csv`, and `bootstrap.csv` are distributed through Zenodo.

- `summary.csv`, `all_split_metrics.csv`: all predictor sets and compendia.
- `tuning/`, `tuning_fga/`: cross-validation results of the hyperparameter
  search; `tuning_comparison.csv`: the search repeated on the training set of
  split 1.
- `inputs/`: input checksums and cohort sizes.
- `feature_importance/`: gain of every predictor of the split-0 multi-omic
  classifiers and the counts of Figure 8e (`code/scripts/figure8.R`).
- `reactome/`: gene lists and Reactome enrichment summary
  (`code/scripts/reactome_enrichment.R`).

## `fraction_genome_altered/` (notebook `06`)

- `clustering_fga_by_cluster.csv`: FGA median and quartiles by cluster.
- `clustering_fga_variance.csv`: share of the variation in FGA explained by
  cluster membership, overall and after accounting for cancer type.
- `survival_cox.csv`: Cox models of Figure 7 on the patients with an FGA value,
  without and with FGA (penalized spline; linear term for comparison), with
  Holm-adjusted P values.
