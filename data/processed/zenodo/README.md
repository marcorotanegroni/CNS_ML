# Zenodo files

The following large processed files are distributed through the associated
Zenodo record and are not tracked by Git.

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

SHA-256:

```text
68c8b08b7c227cf10c44adaed54618b2f71915bf1c5dda14bb7a9d11956b213e
```

## `exp_meth_post_normalization.RData`

Contains `ordinata`, the harmonized pan-cancer expression and methylation
matrix (8,988 x 29,650).

SHA-256:

```text
4dc56d9a91fe81e44ef6ca8f2d1621561d49e7f24994646436cf9efe5f1393d4
```

The Zenodo DOI will be added here after publication of the record.
