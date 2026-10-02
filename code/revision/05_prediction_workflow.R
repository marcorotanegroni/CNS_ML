# Shared entry point for RStudio and the revision notebook.

run_prediction_revision <- function(
    project_dir, final_path = file.path(project_dir, "final_matrices.RData"),
    normalized_path = file.path(project_dir, "exp_meth_post_normalization.RData"),
    output_dir = file.path(project_dir, "results/revision/prediction/split_stability"),
    frameworks = c("Drews", "Steele", "Tao"), seeds = 20261001L + 0:19,
    inner_folds = 10L, n_boot = 2000L, threads = 2L,
    input_mode = c("archived_final", "training_fold")) {
  input_mode <- match.arg(input_mode)
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  env <- new.env(parent = globalenv())
  for (file in c("02_prediction_resampling.R", "03_prediction_training.R",
                 "04_prediction_inputs.R")) {
    sys.source(file.path(project_dir, "code/revision", file), envir = env)
  }
  env$pt_require_packages()
  env$pt_check_xgboost()
  if (!length(frameworks) || anyDuplicated(frameworks) ||
      any(!frameworks %in% c("Drews", "Steele", "Tao"))) stop("Invalid frameworks.")
  inputs <- env$pi_load_inputs(final_path, normalized_path,
    file.path(project_dir, "data/processed/clustering_data.RData"), output_dir,
    input_mode = input_mode)
  print(inputs$audit)
  # Do not hold the three large prepared matrices in memory simultaneously.
  for (framework in frameworks) {
    message("Preparing ", framework, "...")
    prepared <- env$pi_prepare_framework(inputs, framework)
    result <- env$pt_run_prediction(
      prepared, framework, file.path(output_dir, framework), seeds = seeds,
      inner_folds = inner_folds, grid = env$pt_default_grid(), n_boot = n_boot,
      threads = threads, bootstrap_seed = 47001L, resume = TRUE,
      preprocessing = input_mode)
    utils::write.csv(prepared$feature_map,
      file.path(output_dir, framework, "feature_manifest.csv"), row.names = FALSE)
    rm(prepared, result)
    invisible(gc())
  }
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
  prediction_revision_report(output_dir, frameworks)
}

prediction_revision_report <- function(output_dir, frameworks = c("Drews", "Steele", "Tao")) {
  files <- file.path(output_dir, frameworks, "split_metrics.csv")
  available <- file.exists(files)
  if (!any(available)) return(NULL)
  frameworks <- frameworks[available]
  metrics <- do.call(rbind, lapply(files[available], utils::read.csv, stringsAsFactors = FALSE))
  summary <- do.call(rbind, lapply(file.path(output_dir, frameworks, "split_summary.csv"),
                                  utils::read.csv, stringsAsFactors = FALSE))
  utils::write.csv(summary, file.path(output_dir, "F1_across_splits.csv"), row.names = FALSE)
  utils::write.csv(metrics, file.path(output_dir, "F1_per_split_with_bootstrap.csv"), row.names = FALSE)
  plot <- NULL
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    metrics$framework <- factor(metrics$framework, levels = c("Drews", "Steele", "Tao"))
    metrics$split_id <- as.integer(metrics$split_id)
    plot <- ggplot2::ggplot(metrics, ggplot2::aes(split_id, f1, color = framework)) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = f1_lower_95, ymax = f1_upper_95),
                             linewidth = .6, na.rm = TRUE) +
      ggplot2::geom_point(size = 2) +
      ggplot2::facet_wrap(~framework, nrow = 1, drop = FALSE) +
      ggplot2::scale_color_manual(values = c(Drews = "#1f78b4", Steele = "#33a02c", Tao = "#e31a1c")) +
      ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
      ggplot2::labs(x = "Outer train/test split", y = "F1 (Cluster 1 positive)",
        title = "Prediction of framework-specific cluster assignments",
        subtitle = "Points: held-out F1. Bars: 95% test-patient bootstrap intervals.",
        caption = "Each framework defines a different target. Intervals condition on each fitted model; they are not intervals for the mean across splits.") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position = "none", strip.text = ggplot2::element_text(face = "bold"),
                     axis.text = ggplot2::element_text(color = "black"),
                     panel.grid.minor = ggplot2::element_blank(),
                     plot.caption = ggplot2::element_text(size = 8))
    ggplot2::ggsave(file.path(output_dir, "F1_repeated_splits.png"), plot,
                    width = 12, height = 4.5, dpi = 300)
    ggplot2::ggsave(file.path(output_dir, "F1_repeated_splits.pdf"), plot,
                    width = 12, height = 4.5)
  }
  list(summary = summary, metrics = metrics, plot = plot,
       available_frameworks = frameworks)
}
