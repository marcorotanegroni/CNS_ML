# Scalar copy number baseline for the cluster classifiers (Reviewer 2): Cluster
# 1 versus Cluster 2 predicted from fraction genome altered (FGA) alone, with
# the same patients, saved partitions (split 0 and 20 repeated splits),
# XGBoost workflow and bootstrap intervals as the multi-omic classifiers.
#
# Hyperparameters are tuned once per framework on the split-0 training set with
# the original grid and 10-fold CV, as for the multi-omic models (whose
# configuration was tuned for 35,000 predictors and is not reused for one).
# Patients without an FGA value (2-3%) are kept with a missing value, which
# XGBoost handles, so the cohorts and splits are identical.
#
# Run from the repository root after 13_fetch_fga.R and the multi-omic
# evaluation (it reads the saved partitions and labels):
#   Rscript code/revision/14_fga_prediction.R [Drews Steele Tao] [--threads 4]
# Outputs: results/revision/prediction/fixed/tuning_fga/<framework>/ and
#          results/revision/prediction/fixed/fga_baseline/<framework>/
# Then `Rscript code/revision/07_prediction_run.R --step report` adds the
# fga_baseline rows to fixed/summary.csv.

for (f in c("02_prediction_resampling.R", "03_prediction_training.R")) {
  source(file.path("code/revision", f))
}
args <- commandArgs(trailingOnly = TRUE)
threads <- 4L
if ("--threads" %in% args) {
  threads <- as.integer(args[which(args == "--threads") + 1L])
  args <- args[-c(which(args == "--threads"), which(args == "--threads") + 1L)]
}
frameworks <- if (length(args)) args else c("Drews", "Steele", "Tao")
stopifnot(all(frameworks %in% c("Drews", "Steele", "Tao")), !is.na(threads))
options(cnsml.xgb_tree_method = "hist")
root <- "results/revision/prediction/fixed"
fga <- utils::read.csv("data/processed/fga_pancanatlas.csv")
fga <- setNames(fga$fga, fga$patient)

for (framework in frameworks) {
  manifests <- utils::read.csv(file.path(root, "manifests", paste0(framework, ".csv")),
                               colClasses = c(split_id = "character"))
  labels <- utils::read.csv(file.path(root, "archived", framework, "manifests.csv"),
                            colClasses = c(split_id = "character", truth = "character"))
  patients <- labels[labels$split_id == "0", c("patient_id", "cancer_type", "truth")]
  stopifnot(setequal(patients$patient_id, manifests$patient_id[manifests$split_id == "0"]))
  x <- matrix(unname(fga[patients$patient_id]), ncol = 1L,
              dimnames = list(patients$patient_id, "fga"))
  message(framework, ": ", nrow(x), " patients, FGA missing for ", sum(is.na(x)))
  prepared <- list(x = x, patients = patients, feature_type = c(fga = "other"))

  best <- pt_tune_once(prepared, manifests[manifests$split_id == "0", ],
                       file.path(root, "tuning_fga", framework),
                       threads = threads, preprocessing = "unchanged")
  best <- best[, names(pt_default_grid()), drop = FALSE]
  pt_run_prediction(prepared, framework, file.path(root, "fga_baseline", framework),
                    threads = threads, preprocessing = "unchanged",
                    manifests = manifests, fixed_parameters = best)
}
