# Zenodo files

The following large processed files are distributed through the associated
[Zenodo record](https://doi.org/10.5281/zenodo.20617051) and are not tracked
by Git.

The Zenodo record should use the Creative Commons Attribution 4.0
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

DOI: https://doi.org/10.5281/zenodo.20617051
