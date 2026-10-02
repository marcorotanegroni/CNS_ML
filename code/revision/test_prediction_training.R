# Small synthetic end-to-end verification; never loads or trains the cohort.
# Run from the repository root: Rscript code/revision/test_prediction_training.R
source("code/revision/02_prediction_resampling.R")
source("code/revision/03_prediction_training.R")
pt_require_packages()
stopifnot(identical(pt_xgb_iteration_range(1450L, "2.0.3"), c(1L, 1451L)),
          identical(pt_xgb_iteration_range(1450L, "2.1.0"), c(1L, 1450L)),
          identical(pt_xgb_iteration_range(1L, "3.2.0"), c(1L, 1L)))
pt_check_xgboost()
# The blocked implementation must preserve caret's feature decisions, including
# negatively correlated features, duplicate columns/ties, and block boundaries.
set.seed(5721)
cor_x <- matrix(rnorm(200L * 120L), 200L, 120L)
cor_x[, 2L] <- cor_x[, 1L]
cor_x[, 3L] <- -cor_x[, 1L]
cor_x[, 55L] <- .99 * cor_x[, 54L] + .01 * cor_x[, 55L]
cor_x[, 119L] <- -.99 * cor_x[, 118L] + .01 * cor_x[, 119L]
expected_removed <- caret::findCorrelation(stats::cor(cor_x), cutoff = .90, exact = FALSE)
for (block in c(1L, 17L, 64L, 256L)) {
  stopifnot(identical(pt_correlated_columns(cor_x, .90, block), as.integer(expected_removed)))
}
stopifnot(identical(pt_correlated_columns(cor_x[, 1:90]),
                    caret::findCorrelation(stats::cor(cor_x[, 1:90]), cutoff = .90)),
          !length(pt_correlated_columns(cor_x[, 4:120], cutoff = 1)))
set.seed(9104)
n <- 80L
x <- matrix(rnorm(n * 15L), n, 15L,
            dimnames = list(sprintf("sample_%03d", seq_len(n)), paste0("feature_", seq_len(15L))))
x[, 1:10] <- x[, 1:10] + 5
x[, 2] <- x[, 1] * 2 # exercise correlation filtering
x[, 14] <- 1 # constant feature must disappear
x[, 15] <- rep(c(0, 1), n / 2L)
patients <- data.frame(patient_id = rownames(x), cancer_type = rep(c("A", "B"), each = n / 2L),
                       truth = ifelse(x[, 3] + x[, 7] > 10, "1", "2"))
types <- setNames(c(rep("expression", 5), rep("methylation", 5), rep("other", 5)), colnames(x))
prepared <- list(x = x, patients = patients, feature_type = types)
unchanged <- pt_fit_preprocessor(x, types, preprocessing = "archived_final")
stopifnot(identical(pt_apply_preprocessor(x, unchanged), x),
          identical(unchanged$retained, colnames(x)))
grid <- expand.grid(nrounds = c(3L, 6L), eta = .15, max_depth = c(1L, 2L), gamma = 0,
                    colsample_bytree = 1, min_child_weight = 1, subsample = .75)
out <- tempfile("prediction-training-test-")
on.exit_cleanup <- function() unlink(out, recursive = TRUE)
result <- pt_run_prediction(prepared, "synthetic", out, seeds = c(41L, 42L),
                             inner_folds = 2L, grid = grid, n_boot = 30L, threads = 1L)
stopifnot(nrow(result$predictions) == 32L, nrow(result$split_metrics) == 2L,
          nrow(result$bootstrap) == 60L, nrow(result$best_parameters) == 2L,
          nrow(result$cv_results) == 8L,
          identical(result$predictions$predicted,
                    ifelse(result$predictions$probability_cluster1 >= .5, "1", "2")))
checkpoints <- list.files(file.path(out, "checkpoints"), pattern = "^split_.*rds$", full.names = TRUE)
for (path in checkpoints) {
  fitted <- readRDS(path)
  stopifnot(!length(intersect(fitted$outer_training_ids, fitted$outer_test_ids)),
            !"feature_14" %in% fitted$outer_preprocessor$retained,
            identical(fitted$outer_preprocessor$training_ids, fitted$outer_training_ids))
  for (fold in fitted$inner_audit) {
    stopifnot(!length(intersect(fold$preprocessor$training_ids, fold$validation_ids)),
              !length(intersect(fold$preprocessor$training_ids, fitted$outer_test_ids)),
              setequal(c(fold$preprocessor$training_ids, fold$validation_ids), fitted$outer_training_ids))
    altered <- x
    excluded <- setdiff(rownames(x), fold$preprocessor$training_ids)
    altered[excluded, ] <- altered[excluded, ] * 1e4 + 1e6
    refitted <- pt_fit_preprocessor(altered[fold$preprocessor$training_ids, , drop = FALSE], types)
    stopifnot(identical(refitted, fold$preprocessor))
  }
}
# Model paths evaluated at a shorter nrounds must match a separately trained
# model with that nrounds and seed (the grid uses caret-style submodels).
filter <- pt_fit_preprocessor(x, types)
scaled <- pt_apply_preprocessor(x, filter)
assessment_missing <- x
assessment_missing[1L, filter$retained[1L]] <- Inf
assessment_scaled <- pt_apply_preprocessor(assessment_missing, filter)
stopifnot(is.na(assessment_scaled[1L, 1L]), identical(colnames(assessment_scaled), filter$retained))
training_missing <- x
training_missing[1L, 3L] <- NA_real_
missing_filter <- pt_fit_preprocessor(training_missing, types)
stopifnot("feature_3" %in% missing_filter$invalid_or_constant_removed)
long <- pt_xgb_fit(scaled, patients$truth, grid[2L, , drop = FALSE], 99L, 1L)
short <- pt_xgb_fit(scaled, patients$truth, grid[1L, , drop = FALSE], 99L, 1L)
stopifnot(isTRUE(all.equal(pt_xgb_predict(long, scaled, 3L), pt_xgb_predict(short, scaled), tolerance = 1e-7)))
# Verify that resume uses only exact input/settings/code matches.
resumed <- pt_run_prediction(prepared, "synthetic", out, seeds = c(41L, 42L),
                              inner_folds = 2L, grid = grid, n_boot = 30L, threads = 1L)
stopifnot(identical(result, resumed))
changed <- prepared
changed$x[1, 1] <- changed$x[1, 1] + .01
failure <- try(pt_run_prediction(changed, "synthetic", out, seeds = c(41L, 42L),
                                 inner_folds = 2L, grid = grid, n_boot = 30L, threads = 1L), silent = TRUE)
stopifnot(inherits(failure, "try-error"), grepl("different inputs/settings/code", as.character(failure)))
# The original-matrix sensitivity is a separate run, with no new filtering.
original_out <- tempfile("prediction-original-matrix-test-")
original_result <- pt_run_prediction(prepared, "synthetic", original_out, seeds = 41L,
  inner_folds = 2L, grid = grid, n_boot = 30L, threads = 1L,
  preprocessing = "archived_final")
original_checkpoint <- readRDS(list.files(file.path(original_out, "checkpoints"),
                                         full.names = TRUE)[1])
stopifnot(nrow(original_result$split_metrics) == 1L,
          identical(original_checkpoint$outer_preprocessor$retained, colnames(x)),
          all(original_result$filters$stage == "archived_unchanged"))
unlink(original_out, recursive = TRUE)
on.exit_cleanup()
cat("Synthetic nested training passed: split isolation, train-only filters, submodels, outputs, and strict resume.\n")
