# Data

Place processed input files used by the analysis scripts in `data/processed/`.

Processed files currently used by the analysis notebooks:

```text
data/processed/
├── surv_data.RData
├── clustering_data.RData
├── signature_exploration_data.RData
├── list_varImp_class_80perc.RData
├── mut_finale.RData
├── methylation_batch_effect_data.RData
└── zenodo/
    ├── final_matrices.RData
    └── exp_meth_post_normalization.RData
```

Raw TCGA data and external source files should not be redistributed in this repository. If large processed files are stored outside Git, add the download link or accession information here.

Original processed data produced for this project are released under CC BY
4.0, subject to the third-party data conditions described in
`../LICENSE-DATA.md`.

The two files under `data/processed/zenodo/` are distributed through Zenodo
because they exceed the GitHub file-size limit. Their contents, dimensions,
and checksums are documented in `data/processed/zenodo/README.md`.

`methylation_batch_effect_data.RData` contains only the three LAML matrices
used to reproduce Supplementary Figure 6.

| File | Main contents | Used by |
| --- | --- | --- |
| `signature_exploration_data.RData` | Three signature matrices and comparison metadata | `01_signature_exploration.Rmd` |
| `clustering_data.RData` | Signature activities and final clustering objects | `02_clustering.Rmd` |
| `surv_data.RData` | Survival tables, purity, and final clustering objects | `03_survival_analysis.Rmd` |
| `list_varImp_class_80perc.RData` | Aggregated model variable-importance results | `04_machine_learning.Rmd` |
| `mut_finale.RData` | Tumor-specific mutation matrices | `00_data_preprocessing.Rmd` |
| `methylation_batch_effect_data.RData` | Paired LAML 27K, raw 450K, and normalized 450K matrices | `05_methylation_batch_effect.Rmd` |
