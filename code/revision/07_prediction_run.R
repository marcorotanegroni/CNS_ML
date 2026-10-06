# Command-line entry point for the revised predictive benchmark (R1 major
# comment 2). Run from the repository root, one process per framework (they
# can run in parallel):
#
#   Rscript code/revision/07_prediction_run.R --framework Drews --step tune --threads 4
#   Rscript code/revision/07_prediction_run.R --framework Drews --step evaluate --threads 4
#   Rscript code/revision/07_prediction_run.R --step report
#
# tune      one-off 10-fold CV with the original grid on the training part of
#           the reconstructed historical split (split 0), archived predictors.
#           --tune_split 1 repeats the tuning on the training part of split 1
#           as a check that the selected configuration does not depend on the
#           partition; evaluation always uses the split-0 configuration.
# evaluate  fixed tuned parameters, split 0 + 20 repeated cancer-type-stratified
#           80/20 splits, archived predictors (variant "archived", default).
#           Optional variants for the later R2 comments run on the same saved
#           splits: training_fold, training_fold_no_purity, scalar_baseline
#           (not part of the R1 run; scalar source still to be decided).
# report    combine the completed variants into summary tables and a figure.
#
# Options: --variants a,b  --splits 0,1,2  --tune_split 0|1  --tree_method hist|exact (default hist,
# set explicitly because the meaning of "auto" depends on the XGBoost version)
# --params eta=0.15,max_depth=2,...
# (--params overrides tuning, for smoke tests only)  --out <dir>.
# Completed splits are checkpointed; rerunning the same command resumes.

pr_args <- function(args) {
  out <- list(framework = NULL, step = NULL, threads = 2L, variants = NULL,
              splits = NULL, params = NULL, out = NULL, tree_method = "hist",
              tune_split = "0")
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    if (!key %in% names(out) || i == length(args)) stop("Unknown or incomplete option: ", args[i])
    out[[key]] <- args[i + 1L]
    i <- i + 2L
  }
  out$threads <- as.integer(out$threads)
  out
}

run_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- pr_args(args)
  if (!opt$tree_method %in% c("auto", "exact", "hist", "approx")) stop("Invalid --tree_method.")
  options(cnsml.xgb_tree_method = opt$tree_method)
  project_dir <- normalizePath(".", mustWork = TRUE)
  if (!file.exists(file.path(project_dir, "code/revision/07_prediction_run.R"))) {
    stop("Run from the repository root.")
  }
  for (file in c("02_prediction_resampling.R", "03_prediction_training.R",
                 "04_prediction_inputs.R")) {
    source(file.path(project_dir, "code/revision", file))
  }
  root <- if (is.null(opt$out)) file.path(project_dir, "results/revision/prediction/fixed") else opt$out
  if (identical(opt$step, "report")) return(invisible(pr_report(root)))
  if (is.null(opt$framework) || !opt$framework %in% c("Drews", "Steele", "Tao")) {
    stop("--framework must be Drews, Steele or Tao.")
  }
  if (!opt$step %in% c("tune", "evaluate")) stop("--step must be tune, evaluate or report.")
  pt_require_packages()
  pt_check_xgboost()
  framework <- opt$framework
  # Zenodo inputs: data/processed/zenodo/ (as in the README) or repository root.
  zenodo_file <- function(name) {
    candidates <- file.path(project_dir, c("data/processed/zenodo", "."), name)
    found <- candidates[file.exists(candidates)]
    if (length(found)) found[1] else candidates[1]
  }
  final_path <- zenodo_file("final_matrices.RData")
  normalized_path <- zenodo_file("exp_meth_post_normalization.RData")
  cluster_path <- file.path(project_dir, "data/processed/clustering_data.RData")
  ascat_path <- file.path(project_dir, "Metadata_TCGA_ASCAT_penalty70.rds")
  variants <- "archived"
  if (!is.null(opt$variants)) {
    variants <- strsplit(opt$variants, ",")[[1]]
    if (any(!variants %in% c("archived", "training_fold", "training_fold_no_purity",
                             "scalar_baseline"))) stop("Unknown variant.")
  }

  message("[", format(Sys.time()), "] Loading archived inputs for ", framework)
  inputs <- pi_load_inputs(final_path, cluster_path = cluster_path,
                           output_dir = file.path(root, "inputs", framework),
                           input_mode = "archived_final", frameworks = framework)
  prepared <- pi_prepare_framework(inputs, framework)
  rm(inputs); invisible(gc())
  manifests <- pr_framework_manifests(prepared, file.path(root, "manifests", paste0(framework, ".csv")))
  # Training patients of each split, taken before any --splits filtering:
  # tuning and its provenance check always refer to the complete manifest.
  all_manifests <- manifests
  if (!is.null(opt$splits)) {
    keep <- strsplit(opt$splits, ",")[[1]]
    manifests <- manifests[manifests$split_id %in% keep, , drop = FALSE]
  }

  tuning_dir <- file.path(root, "tuning", framework)
  if (opt$step == "tune") {
    split <- opt$tune_split
    if (!split %in% all_manifests$split_id) stop("Unknown --tune_split: ", split)
    message("[", format(Sys.time()), "] One-off tuning on split ", split, " (", framework, ")")
    best <- pt_tune_once(prepared, all_manifests[all_manifests$split_id == split, ],
                         file.path(tuning_dir, paste0("split_", split)),
                         threads = opt$threads)
    print(best)
    return(invisible(best))
  }

  fixed <- if (!is.null(opt$params)) pr_parse_params(opt$params) else {
    best_path <- file.path(tuning_dir, "split_0", "best_parameters.csv")
    if (!file.exists(best_path)) stop("Run --step tune for ", framework, " first.")
    pt_check_tuning(dirname(best_path),
                    all_manifests$patient_id[all_manifests$split_id == "0" &
                                               all_manifests$partition == "train"])
    utils::read.csv(best_path)
  }
  fixed <- fixed[, names(pt_default_grid()), drop = FALSE]
  run_variant <- function(variant, prepared, preprocessing) {
    message("[", format(Sys.time()), "] ", framework, " / ", variant)
    pt_run_prediction(prepared, framework, file.path(root, variant, framework),
                      threads = opt$threads, preprocessing = preprocessing,
                      manifests = manifests, fixed_parameters = fixed)
    invisible(gc())
  }
  if ("archived" %in% variants) run_variant("archived", prepared, "archived_final")
  if ("scalar_baseline" %in% variants) {
    scalar <- pi_scalar_baseline(prepared, ascat_path)
    dir.create(file.path(root, "scalar_baseline"), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(cbind(framework = framework, scalar$coverage),
                     file.path(root, "scalar_baseline", paste0(framework, "_coverage.csv")),
                     row.names = FALSE)
    run_variant("scalar_baseline", scalar, "unchanged")
    rm(scalar)
  }
  rm(prepared); invisible(gc())
  if (any(c("training_fold", "training_fold_no_purity") %in% variants)) {
    message("[", format(Sys.time()), "] Loading pre-filter inputs for ", framework)
    inputs <- pi_load_inputs(final_path, normalized_path, cluster_path,
                             output_dir = file.path(root, "inputs_training_fold", framework),
                             input_mode = "training_fold", frameworks = framework)
    prepared <- pi_prepare_framework(inputs, framework)
    rm(inputs); invisible(gc())
    if ("training_fold" %in% variants) run_variant("training_fold", prepared, "training_fold")
    if ("training_fold_no_purity" %in% variants) {
      run_variant("training_fold_no_purity", pi_drop_features(prepared, "purity"), "training_fold")
    }
  }
  message("[", format(Sys.time()), "] Done: ", framework)
}

# Split 0 (historical seed-1234 partition) followed by 20 repeated splits.
# Saved on first use and reused, so every variant sees identical partitions.
pr_framework_manifests <- function(prepared, path, seeds = 20261001L + 0:19) {
  original <- pr_original_split(prepared$archive_frame, prepared$patients$patient_id)
  repeated <- pr_repeated_splits(prepared$patients, seeds)
  manifests <- rbind(original, repeated[, names(original)])
  manifests$split_id <- as.character(manifests$split_id)
  if (file.exists(path)) {
    saved <- utils::read.csv(path, colClasses = c(split_id = "character"))
    if (!identical(saved$patient_id, manifests$patient_id) ||
        !identical(saved$partition, manifests$partition)) {
      stop("Saved split manifest differs from the reconstructed one: ", path)
    }
  } else {
    # Atomic write: concurrent jobs for the same framework (e.g. tuning on
    # split 0 and split 1) may reach this point together.
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(pattern = "manifest-", tmpdir = dirname(path), fileext = ".csv")
    utils::write.csv(manifests, tmp, row.names = FALSE)
    if (!file.rename(tmp, path)) stop("Could not save split manifest: ", path)
  }
  manifests
}

pr_parse_params <- function(text) {
  values <- strsplit(strsplit(text, ",")[[1]], "=")
  out <- pt_default_grid()[1L, , drop = FALSE]
  for (v in values) {
    if (!v[1] %in% names(out)) stop("Unknown parameter: ", v[1])
    out[[v[1]]] <- as.numeric(v[2])
  }
  rownames(out) <- NULL
  out
}

pr_report <- function(root) {
  # Tuning check: configuration and CV surface on split 0 versus split 1.
  tuned <- Sys.glob(file.path(root, "tuning", "*", "split_*", "cv_results.csv"))
  if (length(tuned)) {
    tuning <- do.call(rbind, lapply(tuned, function(f) {
      cv <- utils::read.csv(f)
      best <- cv[which.max(cv$mean_accuracy), ]
      data.frame(framework = basename(dirname(dirname(f))), tuned_on = basename(dirname(f)),
                 best[, c("nrounds", "eta", "max_depth", "min_child_weight", "mean_accuracy", "sd_accuracy")],
                 min_mean_accuracy = min(cv$mean_accuracy),
                 n_within_1sd_of_best = sum(cv$mean_accuracy >= best$mean_accuracy - best$sd_accuracy))
    }))
    utils::write.csv(tuning, file.path(root, "tuning_comparison.csv"), row.names = FALSE)
    print(tuning, digits = 3)
  }
  files <- Sys.glob(file.path(root, "*", "*", "split_metrics.csv"))
  if (!length(files)) stop("No completed variants under ", root)
  metrics <- do.call(rbind, lapply(files, function(f) {
    x <- utils::read.csv(f, colClasses = c(split_id = "character"))
    cbind(variant = basename(dirname(dirname(f))), x)
  }))
  q <- function(v, p) unname(stats::quantile(v, p, na.rm = TRUE))
  med <- function(v) stats::median(v, na.rm = TRUE)
  rows <- split(metrics, list(metrics$variant, metrics$framework), drop = TRUE)
  summary <- do.call(rbind, lapply(rows, function(x) {
    s0 <- x[x$split_id == "0", ]; r <- x[x$split_id != "0", ]
    data.frame(variant = x$variant[1], framework = x$framework[1],
               prevalence_cluster1 = mean(x$prevalence_cluster1), trivial_f1 = mean(x$trivial_f1),
               split0_f1 = if (nrow(s0)) s0$f1 else NA, split0_f1_lower = if (nrow(s0)) s0$f1_lower_95 else NA,
               split0_f1_upper = if (nrow(s0)) s0$f1_upper_95 else NA,
               n_repeated = nrow(r), median_f1 = med(r$f1),
               f1_q025 = q(r$f1, .025), f1_q975 = q(r$f1, .975),
               median_auroc = med(r$auroc), auroc_q025 = q(r$auroc, .025),
               auroc_q975 = q(r$auroc, .975),
               median_average_precision = med(r$average_precision),
               median_mcc = med(r$mcc), n_undefined_mcc = sum(is.na(r$mcc)),
               median_balanced_accuracy = med(r$balanced_accuracy))
  }))
  rownames(summary) <- NULL
  utils::write.csv(metrics, file.path(root, "all_split_metrics.csv"), row.names = FALSE)
  utils::write.csv(summary, file.path(root, "summary.csv"), row.names = FALSE)
  print(summary, digits = 3)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    colours <- c(Drews = "#1f78b4", Steele = "#33a02c", Tao = "#e31a1c")
    metrics$framework <- factor(metrics$framework, c("Drews", "Steele", "Tao"))
    # Manuscript figure: F1 only, archived predictors. Repeated splits as
    # boxes and points; reconstructed original split as a cross with its 95%
    # bootstrap interval.
    a <- metrics[metrics$variant == "archived", ]
    if (nrow(a)) {
      p <- ggplot2::ggplot(a[a$split_id != "0", ], ggplot2::aes(framework, f1, colour = framework)) +
        ggplot2::geom_boxplot(outlier.shape = NA, width = .5) +
        ggplot2::geom_jitter(width = .1, height = 0, size = 1.2, alpha = .7) +
        ggplot2::geom_errorbar(data = a[a$split_id == "0", ],
                               ggplot2::aes(ymin = f1_lower_95, ymax = f1_upper_95),
                               width = .12, colour = "black", position = ggplot2::position_nudge(x = .35)) +
        ggplot2::geom_point(data = a[a$split_id == "0", ], shape = 4, size = 3, stroke = 1.1,
                            colour = "black", position = ggplot2::position_nudge(x = .35)) +
        ggplot2::scale_colour_manual(values = colours, guide = "none") +
        ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                       axis.text = ggplot2::element_text(colour = "black")) +
        ggplot2::labs(x = NULL, y = "Test-set F1 (Cluster 1)",
          caption = paste("Boxes and points: 20 repeated cancer-type-stratified 80/20 splits.",
                          "Cross: reconstructed original split with 95% bootstrap CI."))
      ggplot2::ggsave(file.path(root, "f1_repeated_splits.png"), p, width = 6, height = 4.5,
                      dpi = 300, bg = "white")
      ggplot2::ggsave(file.path(root, "f1_repeated_splits.pdf"), p, width = 6, height = 4.5)
    }
    # Exploratory overview (not for the manuscript): all variants and metrics.
    long <- do.call(rbind, lapply(c("f1", "auroc", "mcc"), function(m)
      data.frame(variant = metrics$variant, framework = metrics$framework,
                 split_id = metrics$split_id, metric = m, value = metrics[[m]])))
    long$metric <- factor(long$metric, c("f1", "auroc", "mcc"), c("F1 (Cluster 1)", "AUROC", "MCC"))
    trivial <- unique(metrics[, c("framework", "trivial_f1")])
    trivial$metric <- factor("F1 (Cluster 1)", levels(long$metric))
    p <- ggplot2::ggplot(long[long$split_id != "0", ], ggplot2::aes(variant, value, colour = framework)) +
      ggplot2::geom_boxplot(outlier.shape = NA, width = .6, position = ggplot2::position_dodge(.75)) +
      ggplot2::geom_point(data = long[long$split_id == "0", ], shape = 4, size = 2.5, stroke = 1,
                          position = ggplot2::position_dodge(.75)) +
      ggplot2::geom_hline(data = trivial, ggplot2::aes(yintercept = trivial_f1, colour = framework),
                          linetype = "dashed", linewidth = .4) +
      ggplot2::facet_wrap(~metric, scales = "free_y") +
      ggplot2::scale_colour_manual(values = colours) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
      ggplot2::labs(x = NULL, y = NULL, colour = NULL,
        caption = "Exploratory. Boxes: repeated splits. Crosses: reconstructed original split. Dashed: F1 of predicting Cluster 1 for every patient.")
    ggplot2::ggsave(file.path(root, "exploratory_metrics.png"), p, width = 12, height = 4.5,
                    dpi = 200, bg = "white")
  }
  invisible(summary)
}

if (sys.nframe() == 0L) run_cli()
