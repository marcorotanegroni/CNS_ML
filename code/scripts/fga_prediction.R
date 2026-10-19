# Cluster 1 versus Cluster 2 classifiers from fraction genome altered (FGA)
# alone (notebook 06, Supplementary Figure 9c), with the same patients, splits
# (split 0 and 20 repeated splits), XGBoost workflow and bootstrap intervals as
# the multi-omic classifiers of Figure 8.
#
# Hyperparameters are tuned once per framework on the split-0 training set with
# the grid and 10-fold CV of the multi-omic classifiers. Patients without an
# FGA value (about 2%) are kept with a missing value, which XGBoost handles, so
# the cohorts and splits are identical.
#
# Run from the repository root after code/scripts/fga_fetch.R and the
# multi-omic evaluation (it reads the saved splits and labels):
#   Rscript code/scripts/fga_prediction.R [Drews Steele Tao] [--threads 4]
# Outputs: results/cluster_prediction/tuning_fga/<framework>/ and
#          results/cluster_prediction/fga_only/<framework>/
# Then `Rscript code/scripts/cluster_prediction.R --step report` adds the
# fga_only rows to results/cluster_prediction/summary.csv.

for (f in c("prediction_resampling.R", "prediction_training.R")) {
  source(file.path("code/R", f))
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
root <- "results/cluster_prediction"
fga <- utils::read.csv("data/processed/fga_pancanatlas.csv")
fga <- setNames(fga$fga, fga$patient)

for (framework in frameworks) {
  manifests <- utils::read.csv(file.path(root, "manifests", paste0(framework, ".csv")),
                               colClasses = c(split_id = "character"))
  labels <- utils::read.csv(file.path(root, "multiomic", framework, "manifests.csv"),
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
  pt_run_prediction(prepared, framework, file.path(root, "fga_only", framework),
                    threads = threads, preprocessing = "unchanged",
                    manifests = manifests, fixed_parameters = best)
}
