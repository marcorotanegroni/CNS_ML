# Graphical abstract (revision)

`Graphical_Abstract_revised.svg` is the submitted
`Graphical_Abstract_NAR_Cancer_final.svg` with four text changes; the PNG
(2000 x 800 px, as submitted) is rendered from it with the R package magick.
The SVG is the vector version.

- F1 values: Drews 0.93, Steele 0.64, Tao 0.64 (submitted: 0.93, 0.80, 0.24),
  from the reconstructed original split of the retrained classifiers
  (`results/revision/prediction/fixed/summary.csv`).
- "Variable test-set F1" became "F1 (distinct targets)": the three classifiers
  predict different framework-specific partitions.
- Bottom line: "Expression-only models retained F1 for Drews and Steele, not
  Tao" (expression-only F1 on the original split 0.93, 0.63 and 0.38;
  `variant = expression_only` in the same summary).
