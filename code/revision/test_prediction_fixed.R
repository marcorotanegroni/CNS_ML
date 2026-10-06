# Synthetic verification of the fixed-parameter workflow (split 0 + repeated
# splits supplied as manifests, one fit per split, missing-value baseline).
# Run from the repository root: Rscript code/revision/test_prediction_fixed.R
source("code/revision/02_prediction_resampling.R")
source("code/revision/03_prediction_training.R")
source("code/revision/04_prediction_inputs.R")
pt_require_packages()
set.seed(3301)
n <- 120L
ids <- sprintf("P%03d", seq_len(n))
x <- matrix(rnorm(n * 6L), n, 6L, dimnames = list(ids, paste0("f", 1:6)))
patients <- data.frame(patient_id = ids, cancer_type = rep(c("A", "B", "C"), each = n / 3L),
                       truth = ifelse(x[, 1] + x[, 2] > 0, "1", "2"))
archive <- data.frame(patient_id = c(ids, "X001", "X002"),
                      cancer_type = c(patients$cancer_type, "A", "D"))

# Split 0 is drawn on the whole archive, then restricted to eligible patients,
# and must be identical to dplyr's grouped sample_frac with seed 1234.
s0 <- pr_original_split(archive, ids)
set.seed(1234)
ref <- dplyr::sample_frac(dplyr::group_by(archive, cancer_type), .8)$patient_id
stopifnot(setequal(s0$patient_id[s0$partition == "train"], intersect(ref, ids)),
          setequal(s0$patient_id, ids), all(s0$split_id == "0"))
manifests <- rbind(s0, pr_repeated_splits(patients[, 1:2], c(11L, 12L))[, names(s0)])

fixed <- data.frame(nrounds = 5L, eta = .3, max_depth = 2L, gamma = 0,
                    colsample_bytree = 1, min_child_weight = 1, subsample = 1)
prepared <- list(x = x, patients = patients,
                 feature_type = setNames(rep("other", 6), colnames(x)))
out <- tempfile("prediction-fixed-test-")
res <- pt_run_prediction(prepared, "synthetic", out, n_boot = 20L, threads = 1L,
                         preprocessing = "archived_final", manifests = manifests,
                         fixed_parameters = fixed)
stopifnot(identical(sort(unique(res$split_metrics$split_id)), c("0", "1", "2")),
          res$split_summary$n_splits == 2L, is.null(res$cv_results),
          all(c("auroc", "mcc", "trivial_f1") %in% names(res$split_metrics)),
          all(res$best_parameters$nrounds == 5L))
# Test patients of split 0 are exactly its manifest test partition.
p0 <- res$predictions[res$predictions$split_id == "0", ]
stopifnot(setequal(p0$patient_id, s0$patient_id[s0$partition == "test"]))
# Resume with a different thread count reuses checkpoints.
res2 <- pt_run_prediction(prepared, "synthetic", out, n_boot = 20L, threads = 2L,
                          preprocessing = "archived_final", manifests = manifests,
                          fixed_parameters = fixed)
stopifnot(identical(res$predictions, res2$predictions))

# Missing values are allowed only in the unchanged (scalar) mode.
xm <- x; xm[1:10, 1] <- NA
prepared_na <- prepared; prepared_na$x <- xm
res_na <- pt_run_prediction(prepared_na, "synthetic", tempfile(), n_boot = 20L, threads = 1L,
                            preprocessing = "unchanged", manifests = manifests,
                            fixed_parameters = fixed)
stopifnot(nrow(res_na$split_metrics) == 3L,
          inherits(try(pt_fit_preprocessor(xm, prepared$feature_type, preprocessing = "archived_final"),
                       silent = TRUE), "try-error"))

# Dropping a predictor by source name.
prepared$feature_map <- data.frame(feature = colnames(x), source_name = c("age", "purity", paste0("g", 1:4)),
                                   feature_type = "other")
dropped <- pi_drop_features(prepared, "purity")
stopifnot(!"f2" %in% colnames(dropped$x), ncol(dropped$x) == 5L)
# Average precision does not depend on the order of tied probabilities and
# matches the step-wise definition on a hand-computed example.
y <- c("1", "2", "1", "2", "1"); pr <- c(.9, .5, .5, .5, .1)
m1 <- pr_classification_metrics(y, pr, ifelse(pr >= .5, "1", "2"))
o <- c(1, 4, 3, 2, 5)
m2 <- pr_classification_metrics(y[o], pr[o], ifelse(pr[o] >= .5, "1", "2"))
stopifnot(isTRUE(all.equal(m1$average_precision, m2$average_precision)),
          isTRUE(all.equal(m1$average_precision, 1 / 3 * 1 + 1 / 3 * (2 / 4) + 1 / 3 * (3 / 5))))
# Predictions all in one class: MCC is undefined but summaries still export.
one_class <- pr_classification_metrics(c("1", "2", "2"), c(.2, .3, .1), c("2", "2", "2"))
stopifnot(is.na(one_class$mcc))
fixed_low <- fixed; fixed_low$nrounds <- 1L; fixed_low$eta <- 1e-6
res_one <- pt_run_prediction(prepared, "synthetic", tempfile(), n_boot = 20L, threads = 1L,
                             preprocessing = "archived_final", manifests = manifests,
                             fixed_parameters = fixed_low)
stopifnot(nrow(res_one$split_summary) == 1L)

# One-off tuning on two different splits and the report's tuning comparison.
source("code/revision/07_prediction_run.R")
small_grid <- expand.grid(nrounds = c(3L, 6L), eta = .3, max_depth = c(1L, 2L), gamma = 0,
                          colsample_bytree = 1, min_child_weight = 1, subsample = 1)
root <- tempfile("prediction-root-")
prepared$feature_map <- NULL
for (split in c("0", "1")) {
  best <- pt_tune_once(prepared, manifests[manifests$split_id == split, ],
                       file.path(root, "tuning", "synthetic", paste0("split_", split)),
                       grid = small_grid, inner_folds = 3L, threads = 1L)
  stopifnot(nrow(best) == 1L)
}
again <- pt_tune_once(prepared, manifests[manifests$split_id == "0", ],
                      file.path(root, "tuning", "synthetic", "split_0"), grid = small_grid,
                      inner_folds = 3L)
# Saved tuning is not reused under different settings.
old_method <- getOption("cnsml.xgb_tree_method"); options(cnsml.xgb_tree_method = "exact")
stopifnot(inherits(try(pt_tune_once(prepared, manifests[manifests$split_id == "0", ],
                                    file.path(root, "tuning", "synthetic", "split_0"),
                                    grid = small_grid, inner_folds = 3L), silent = TRUE), "try-error"),
          inherits(try(pt_check_tuning(file.path(root, "tuning", "synthetic", "split_0"),
                                       manifests$patient_id[manifests$split_id == "1" &
                                                              manifests$partition == "train"],
                                       small_grid, 3L), silent = TRUE), "try-error"))
options(cnsml.xgb_tree_method = old_method)
# Notes written before hashes were recorded (current workstation run) pass the basic checks.
legacy <- tempfile("legacy-tuning-"); dir.create(legacy)
train0 <- manifests$patient_id[manifests$split_id == "0" & manifests$partition == "train"]
writeLines(c(paste("Training patients:", length(train0)),
             "Inner folds: 10 ; seed: 1234 ; threads: 4 ; tree_method: auto"),
           file.path(legacy, "tuning_notes.txt"))
stopifnot(isTRUE(pt_check_tuning(legacy, train0)))
file.copy(file.path(out), root, recursive = TRUE)
dir.create(file.path(root, "archived"))
file.rename(file.path(root, basename(out)), file.path(root, "archived", "synthetic"))
summary <- pr_report(root)
comparison <- read.csv(file.path(root, "tuning_comparison.csv"))
stopifnot(nrow(comparison) == 2L, setequal(comparison$tuned_on, c("split_0", "split_1")),
          nrow(summary) == 1L, summary$n_repeated == 2L, !is.na(summary$split0_f1))
unlink(c(out, root), recursive = TRUE)
cat("Fixed-parameter workflow passed: split 0 reconstruction, supplied manifests,",
    "single fit per split, extra metrics, thread-independent resume, missing-value mode,",
    "split-0/split-1 tuning and report, tie-invariant average precision,",
    "undefined MCC, tuning provenance checks.\n")
