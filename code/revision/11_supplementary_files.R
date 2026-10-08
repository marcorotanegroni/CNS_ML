# Supplementary items of the revision under their manuscript numbers.
# Copies the figures produced by 01_binary_concordance.R and
# 09_supplementary_figure5.R and writes Supplementary Table 2 (binary
# agreement per signature pair and the thresholds of each rule) as CSV.
# Exploratory quantities (expected n11, observed/expected n11, Cohen's kappa)
# stay in results/revision/binary_concordance/ and are not part of the table.
#
# Run from the repository root after those scripts:
#   Rscript code/revision/11_supplementary_files.R
# Output: results/revision/supplement/

out <- "results/revision/supplement"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
concordance <- "results/revision/binary_concordance"

copies <- c(
  Supplementary_Figure_5 = "results/revision/survival/Supplementary_Figure5",
  Supplementary_Figure_7 = file.path(concordance, "state_agreement_components"),
  Supplementary_Figure_8 = file.path(concordance, "threshold_sensitivity_principal_pairs"))
for (name in names(copies)) {
  for (ext in c("png", "tiff", "pdf")) {
    src <- paste0(copies[[name]], ".", ext)
    if (file.exists(src)) {
      stopifnot(file.copy(src, file.path(out, paste0(name, ".", ext)), overwrite = TRUE))
    }
  }
}

pairs <- read.csv(file.path(concordance, "pair_metrics.csv"))
pair_columns <- c("rule", "signature1", "signature2", "framework1", "framework2",
                  "cross_framework", "high_activity_pair", "n", "n11", "n00", "n10", "n01",
                  "state_jaccard", "active_share_of_agreements", "inactive_share_of_agreements",
                  "active_jaccard", "active_prevalence1", "active_prevalence2")
stopifnot(all(pair_columns %in% names(pairs)))
write.csv(pairs[, pair_columns], file.path(out, "Supplementary_Table_2_pair_metrics.csv"),
          row.names = FALSE, na = "NA")

thresholds <- read.csv(file.path(concordance, "thresholds.csv"))
threshold_columns <- c("rule", "rule_label", "signature", "framework", "reference_cohort",
                       "reference_n", "threshold", "comparison", "matched_n", "active_n",
                       "active_prevalence")
stopifnot(all(threshold_columns %in% names(thresholds)))
write.csv(thresholds[, threshold_columns], file.path(out, "Supplementary_Table_2_thresholds.csv"),
          row.names = FALSE, na = "NA")

writeLines(c(
  "Supplementary Table 2. Binary agreement between signatures under four activity thresholds.",
  "",
  "Supplementary_Table_2_pair_metrics.csv: one row per signature pair and rule (5,881 matched tumours).",
  "  rule: original (zero/median), positive_q25, positive_q50, positive_q75 (quantiles of positive exposure).",
  "  n11, n00, n10, n01: patients active in both, inactive in both, active only in signature1, active only in signature2.",
  "  state_jaccard: patient-state Jaccard (n11 + n00) / (n11 + n00 + 2 (n10 + n01)), as in Figure 4.",
  "  active_share_of_agreements, inactive_share_of_agreements: n11 / (n11 + n00) and n00 / (n11 + n00).",
  "  active_jaccard: n11 / (n11 + n10 + n01).",
  "  cross_framework: signatures from different compendia; high_activity_pair: both among the 12 signatures selected for clustering.",
  "  active_prevalence1/2: proportion of patients classified active for each signature.",
  "",
  "Supplementary_Table_2_thresholds.csv: threshold of each signature under each rule and the resulting number and proportion of active tumours.",
  "  comparison: how the threshold is applied (e.g. > 0 or >= median); reference_cohort and reference_n: cohort on which the threshold was computed."),
  file.path(out, "README.txt"))
message("Supplementary files written to ", normalizePath(out))
