Supplementary Table 2. Binary agreement between signatures under four activity thresholds.

Supplementary_Table_2_pair_metrics.csv: one row per signature pair and rule (5,881 matched tumours).
  rule: original (zero/median), positive_q25, positive_q50, positive_q75 (quantiles of positive exposure).
  n11, n00, n10, n01: patients active in both, inactive in both, active only in signature1, active only in signature2.
  state_jaccard: patient-state Jaccard (n11 + n00) / (n11 + n00 + 2 (n10 + n01)), as in Figure 4.
  active_share_of_agreements, inactive_share_of_agreements: n11 / (n11 + n00) and n00 / (n11 + n00).
  active_jaccard: n11 / (n11 + n10 + n01).
  cross_framework: signatures from different compendia; high_activity_pair: both among the 12 signatures selected for clustering.
  active_prevalence1/2: proportion of patients classified active for each signature.

Supplementary_Table_2_thresholds.csv: threshold of each signature under each rule and the resulting number and proportion of active tumours.
  comparison: how the threshold is applied (e.g. > 0 or >= median); reference_cohort and reference_n: cohort on which the threshold was computed.
