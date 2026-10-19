# Zenodo files

The following large processed files are distributed through the associated
[Zenodo record](https://doi.org/10.5281/zenodo.20617050) and are not tracked
by Git.

The Zenodo record is released under the Creative Commons Attribution 4.0
International license (CC BY 4.0). This applies to the authors' processed
outputs and does not replace the access, reuse, or citation requirements of
the underlying TCGA and third-party source data.

## `final_matrices.RData`

Final multi-omics matrices used for model development:

| Object | Dimensions |
| --- | --- |
| `mat_purity_mut_drews` | 5,024 x 35,497 |
| `mat_purity_mut_steele` | 7,893 x 35,301 |
| `mat_purity_mut_tao` | 7,823 x 35,158 |

## `exp_meth_post_normalization.RData`

Contains `ordinata`, the harmonized pan-cancer expression and methylation
matrix (8,988 x 29,650).

## `cluster_prediction_patient_level.zip`

Patient-level files of the Cluster 1 versus Cluster 2 classifiers (notebook
`04`). Extract the archive into `results/cluster_prediction/`:

- `manifests/<framework>.csv`: train/test assignment of every patient in the
  original split (split 0, seed 1234) and in the 20 repeated splits;
- `<predictor set>/<framework>/manifests.csv`: the same assignments with the
  Cluster 1/2 label of each patient;
- `<predictor set>/<framework>/predictions.csv`: test-set probability of
  Cluster 1 and predicted cluster;
- `<predictor set>/<framework>/bootstrap.csv`: the 2,000 bootstrap replicates
  of the test-set confusion matrix of each split.

Predictor sets: `multiomic`, `expression_only`, `multiomic_no_purity`,
`training_fold`, and `fga_only`.

DOI: https://doi.org/10.5281/zenodo.20617050

`code/R/prediction_inputs.R` verifies the checksums of the two `.RData` files:

- `final_matrices.RData`: MD5 `cd14cf9fd5eef77e3ea35d4022283eb5`.
- `exp_meth_post_normalization.RData`: MD5 `35574218fed4aad02c1557461be95c91`.
