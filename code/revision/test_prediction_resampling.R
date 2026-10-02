# Run from the repository root: Rscript code/revision/test_prediction_resampling.R
source("code/revision/02_prediction_resampling.R")

known <- pr_f1(c("1", "1", "1", "2", "2", "2"),
               c("1", "1", "2", "1", "2", "2"))
stopifnot(identical(unlist(known[1:4], use.names = FALSE), c(2L, 1L, 1L, 2L)),
          abs(known$f1 - 2/3) < 1e-12,
          is.na(pr_f1("2", "2")$f1))
patients <- data.frame(patient_id = paste0("P", 1:14),
                       cancer_type = rep(c("A", "B", "Singleton"), c(10, 3, 1)))
splits <- suppressWarnings(pr_repeated_splits(patients, c(19L, 20L)))
stopifnot(identical(splits, suppressWarnings(
  pr_repeated_splits(patients[14:1, ], c(19L, 20L)))),
  all(splits$partition[splits$singleton_stratum] == "train"))
for (k in unique(splits$split_id)) {
  test <- subset(splits, split_id == k & partition == "test")
  stopifnot(identical(as.integer(table(test$cancer_type)), c(2L, 1L)))
}

predictions <- data.frame(framework = "Drews", split_id = "1",
  patient_id = paste0("P", 1:6), cancer_type = rep(c("A", "B"), each = 3),
  truth = c("1", "1", "1", "2", "2", "2"),
  predicted = c("1", "1", "2", "1", "2", "2"))
for (b in 1:20) {
  idx <- pr_bootstrap_indices(predictions$cancer_type)
  stopifnot(identical(table(predictions$cancer_type[idx]),
                      table(predictions$cancer_type)))
}
predictions <- rbind(predictions, transform(predictions, split_id = "2"))
result <- pr_bootstrap_test_f1(predictions, 100L, 19L)
stopifnot(identical(result, pr_bootstrap_test_f1(predictions[12:1, ], 100L, 19L)),
          nrow(result$summary) == 2L, nrow(result$bootstrap) == 200L,
          all(result$summary$n_test == 6L),
          all(result$summary$n_undefined == 0L),
          all(is.finite(result$summary$f1_lower_95)))
bad <- try(pr_bootstrap_test_f1(rbind(predictions, predictions[1, ]), 10L, 19L),
           silent = TRUE)
stopifnot(inherits(bad, "try-error"))
sparse <- predictions[1:2, ]
sparse$truth <- c("1", "2"); sparse$predicted <- c("2", "2")
sparse_result <- pr_bootstrap_test_f1(sparse, 100L, 19L)$summary
stopifnot(sparse_result$n_undefined > 0L,
          is.na(sparse_result$f1_lower_95), is.na(sparse_result$f1_upper_95))
cat("Verified F1 counts, stratified splits/bootstrap, reproducibility, separate\n",
    "split evaluations, duplicate rejection and undefined-bootstrap reporting.\n")
