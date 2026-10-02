# Nested, repeated held-out evaluation of the existing two-cluster targets.
# Source 02_prediction_resampling.R first. Sourcing this file never fits models.

pt_default_grid <- function() {
  expand.grid(nrounds = seq(1000L, 1500L, 75L), eta = c(.005, .01, .15),
              max_depth = c(2L, 4L, 6L), gamma = 0,
              colsample_bytree = 1, min_child_weight = c(3, 5),
              subsample = .75)
}

pt_require_packages <- function() {
  packages <- c("caret", "digest", "xgboost")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
  if (!exists("pr_repeated_splits", mode = "function")) {
    stop("Source code/revision/02_prediction_resampling.R before this file.")
  }
  invisible(TRUE)
}

pt_validate_prepared <- function(prepared) {
  if (!is.list(prepared) || !all(c("x", "patients", "feature_type") %in% names(prepared))) {
    stop("prepared must contain x, patients, and feature_type.")
  }
  x <- prepared$x; p <- prepared$patients
  if (!is.matrix(x) || !is.numeric(x) || !ncol(x) || nrow(x) < 4L ||
      is.null(rownames(x)) || is.null(colnames(x)) ||
      anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
    stop("x must be a numeric matrix with unique patient row names and feature column names.")
  }
  if (!is.data.frame(p) || !all(c("patient_id", "cancer_type", "truth") %in% names(p))) {
    stop("patients must contain patient_id, cancer_type, and truth.")
  }
  for (name in c("patient_id", "cancer_type", "truth")) {
    p[[name]] <- as.character(p[[name]])
    pr_check_text(p[[name]], name)
  }
  if (anyDuplicated(p$patient_id) || !setequal(p$patient_id, rownames(x)) ||
      !setequal(p$truth, c("1", "2"))) stop("Patient IDs or binary truth labels are invalid.")
  ft <- prepared$feature_type
  if (is.null(names(ft)) || anyDuplicated(names(ft)) ||
      !setequal(names(ft), colnames(x)) || anyNA(ft) ||
      any(!ft %in% c("expression", "methylation", "other"))) {
    stop("feature_type must name every x column as expression, methylation, or other.")
  }
  p <- p[match(rownames(x), p$patient_id), c("patient_id", "cancer_type", "truth"), drop = FALSE]
  rownames(p) <- NULL
  list(x = x, patients = p,
       feature_type = ft[colnames(x)])
}

pt_hash <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)

# Hash actual numeric input in small column blocks instead of serializing a
# multi-gigabyte matrix into one equally large intermediate raw vector.
pt_data_hash <- function(prepared, block_columns = 128L) {
  starts <- seq.int(1L, ncol(prepared$x), by = block_columns)
  blocks <- vapply(starts, function(start) {
    finish <- min(ncol(prepared$x), start + block_columns - 1L)
    pt_hash(prepared$x[, start:finish, drop = FALSE])
  }, character(1))
  pt_hash(list(blocks = blocks, patients = prepared$patients,
               feature_type = prepared$feature_type))
}

pt_seed <- function(seed, offset) {
  as.integer((as.double(seed) + as.double(offset)) %% (.Machine$integer.max - 1L))
}

# Match caret's default rule without retaining a full, very large correlation
# matrix. The fast rule (p >= 100) compares the mean absolute correlations of
# each above-cutoff pair. Each complete column mean is computed in one block,
# preserving the summation used by colMeans(abs(cor(x))). Only pair indices
# above the cutoff are retained. Small matrices use caret's exact algorithm.
pt_correlated_columns <- function(x, cutoff = .90, block_columns = 256L) {
  p <- ncol(x)
  if (p < 2L) return(integer())
  if (p < 100L) return(caret::findCorrelation(stats::cor(x), cutoff = cutoff))
  pr_check_integers(block_columns, "block_columns", minimum = 1L)
  if (length(block_columns) != 1L) stop("block_columns must be a scalar.")
  starts <- seq.int(1L, p, by = block_columns)
  average <- numeric(p)
  pair_rows <- pair_columns <- vector("list", length(starts))
  for (b in seq_along(starts)) {
    columns <- seq.int(starts[b], min(p, starts[b] + block_columns - 1L))
    absolute <- abs(stats::cor(x, x[, columns, drop = FALSE]))
    if (any(!is.finite(absolute))) stop("Nonfinite correlation in training features.")
    average[columns] <- colMeans(absolute)
    rows <- lapply(seq_along(columns), function(j) {
      if (columns[j] == 1L) return(integer())
      which(absolute[seq_len(columns[j] - 1L), j] > cutoff)
    })
    pair_rows[[b]] <- unlist(rows, use.names = FALSE)
    pair_columns[[b]] <- rep.int(columns, lengths(rows))
  }
  # Retain caret's ranking step and tie convention exactly.
  rank <- as.numeric(as.factor(average))
  rows <- unlist(pair_rows, use.names = FALSE)
  columns <- unlist(pair_columns, use.names = FALSE)
  discard_column <- rank[columns] > rank[rows]
  unique(c(columns[discard_column], rows[!discard_column]))
}

# All statistics are fitted exclusively on the rows supplied here.
# The CV is invariant to the historical positive RMS rescaling; computing it
# on unscaled values avoids doing that transformation twice.
pt_fit_preprocessor <- function(x, feature_type, drop_fraction = .2, cutoff = .90,
                                preprocessing = "training_fold") {
  preprocessing <- match.arg(preprocessing, c("archived_final", "training_fold"))
  if (nrow(x) < 2L) stop("At least two training rows are required.")
  feature_type <- feature_type[colnames(x)]
  if (preprocessing == "archived_final") {
    # R1 split sensitivity conditions on the unchanged historical predictor
    # matrix. No new selection or rescaling is fitted in this mode.
    if (any(!vapply(seq_len(ncol(x)), function(j) all(is.finite(x[, j])), logical(1)))) {
      stop("The archived final predictor matrix must contain finite numeric values.")
    }
    counts <- data.frame(stage = "archived_unchanged",
      feature_type = c("expression", "methylation", "other"),
      n_features = vapply(c("expression", "methylation", "other"),
        function(type) sum(feature_type == type), integer(1)))
    return(list(retained = colnames(x), rms = setNames(rep(1, ncol(x)), colnames(x)),
      training_ids = rownames(x), training_ids_hash = pt_hash(rownames(x)),
      counts = counts, preprocessing = preprocessing))
  }
  stages <- list(input = colnames(x))
  finite <- vapply(seq_len(ncol(x)), function(j) all(is.finite(x[, j])), logical(1))
  candidate_indices <- which(finite)
  sds <- vapply(candidate_indices, function(j) stats::sd(x[, j]), numeric(1))
  names(sds) <- colnames(x)[candidate_indices]
  valid <- names(sds)[is.finite(sds) & sds > 0]
  stages$finite_nonconstant <- valid
  if (!length(valid)) stop("No finite, nonconstant training features remain.")
  means <- vapply(valid, function(name) mean(x[, name]), numeric(1))
  cv <- sds[valid] / abs(means)
  # A variable with zero mean and positive SD has infinite CV and is retained
  # by the low-CV filter. Feature order resolves equal CV values deterministically.
  low_cv <- character()
  for (modality in c("expression", "methylation")) {
    features <- valid[feature_type[valid] == modality]
    n_drop <- floor(drop_fraction * length(features))
    if (n_drop > 0L) low_cv <- c(low_cv, features[order(cv[features], seq_along(features))][seq_len(n_drop)])
  }
  after_cv <- valid[!valid %in% low_cv]
  stages$after_low_cv <- after_cv
  omics <- after_cv[feature_type[after_cv] %in% c("expression", "methylation")]
  correlated <- character()
  if (length(omics) > 1L) {
    removed <- pt_correlated_columns(x[, omics, drop = FALSE], cutoff = cutoff)
    if (length(removed)) correlated <- omics[removed]
  }
  retained <- after_cv[!after_cv %in% correlated]
  stages$after_correlation <- retained
  if (!length(retained)) stop("No features remain after correlation filtering.")
  rms <- vapply(retained, function(name) sqrt(sum(x[, name]^2) / (nrow(x) - 1)), numeric(1))
  if (any(!is.finite(rms) | rms <= 0)) stop("Invalid training RMS scaling factors.")
  counts <- do.call(rbind, lapply(names(stages), function(stage) {
    data.frame(stage = stage, feature_type = c("expression", "methylation", "other"),
               n_features = vapply(c("expression", "methylation", "other"), function(type) {
                 sum(feature_type[stages[[stage]]] == type)
               }, integer(1)), row.names = NULL)
  }))
  list(retained = retained, rms = rms, low_cv_removed = low_cv,
       correlated_removed = correlated, invalid_or_constant_removed = setdiff(colnames(x), valid),
       training_ids = rownames(x), training_ids_hash = pt_hash(rownames(x)),
       counts = counts, cv = cv, drop_fraction = drop_fraction,
       correlation_cutoff = cutoff)
}

pt_apply_preprocessor <- function(x, fitted) {
  if (!all(fitted$retained %in% colnames(x))) stop("Prediction matrix is missing fitted features.")
  x <- x[, fitted$retained, drop = FALSE]
  # A feature can be complete in training but missing in assessment rows.
  # XGBoost uses the missing-value direction fitted on training; no assessment
  # statistic is estimated, and test missingness does not select features.
  for (j in seq_len(ncol(x))) {
    missing <- !is.finite(x[, j])
    if (any(missing)) x[missing, j] <- NA_real_
  }
  sweep(x, 2L, fitted$rms, "/")
}

pt_validate_grid <- function(grid) {
  columns <- names(pt_default_grid())
  if (!is.data.frame(grid) || !all(columns %in% names(grid)) || !nrow(grid)) {
    stop("grid must contain the xgbTree tuning columns.")
  }
  grid <- unique(grid[, columns, drop = FALSE])
  if (!all(vapply(grid, is.numeric, logical(1))) || any(!is.finite(as.matrix(grid)))) {
    stop("All grid values must be finite numbers.")
  }
  if (any(grid$nrounds < 1 | grid$nrounds != floor(grid$nrounds)) ||
      any(grid$max_depth < 1 | grid$max_depth != floor(grid$max_depth)) ||
      any(grid$eta <= 0 | grid$eta > 1) || any(grid$gamma < 0) ||
      any(grid$min_child_weight < 0) ||
      any(grid$colsample_bytree <= 0 | grid$colsample_bytree > 1) ||
      any(grid$subsample <= 0 | grid$subsample > 1)) stop("Invalid grid values.")
  # caret's xgbTree complexity ordering breaks exact Accuracy ties.
  grid <- grid[with(grid, order(nrounds, max_depth, eta, gamma,
                                colsample_bytree, min_child_weight, subsample)), , drop = FALSE]
  rownames(grid) <- NULL
  grid
}

pt_xgb_fit <- function(x, truth, params, seed, threads) {
  dtrain <- xgboost::xgb.DMatrix(x, label = as.integer(truth == "1"))
  args <- as.list(params[1L, setdiff(names(params), "nrounds"), drop = FALSE])
  args$objective <- "binary:logistic"
  args$nthread <- as.integer(threads)
  args$verbosity <- 0L
  pr_with_seed(seed, xgboost::xgb.train(params = args, data = dtrain,
                                      nrounds = as.integer(params$nrounds), verbose = 0))
}

pt_xgb_iteration_range <- function(nrounds, version = utils::packageVersion("xgboost")) {
  pr_check_integers(nrounds, "nrounds", minimum = 1L)
  if (length(nrounds) != 1L) stop("Supply one nrounds value.")
  # R's upper endpoint changed from exclusive to inclusive in XGBoost 2.1.0.
  # Both APIs use a one-based start; this is not Python's iteration_range.
  end <- as.integer(nrounds)
  if (numeric_version(as.character(version)) < numeric_version("2.1.0")) end <- end + 1L
  c(1L, end)
}

pt_xgb_predict <- function(model, x, nrounds = NULL) {
  dm <- xgboost::xgb.DMatrix(x)
  if (is.null(nrounds)) return(as.numeric(predict(model, dm)))
  method <- getS3method("predict", "xgb.Booster")
  if ("iterationrange" %in% names(formals(method))) {
    return(as.numeric(predict(model, dm, iterationrange = pt_xgb_iteration_range(nrounds))))
  }
  if ("ntreelimit" %in% names(formals(method))) {
    return(as.numeric(predict(model, dm, ntreelimit = as.integer(nrounds))))
  }
  stop("Unsupported xgboost prediction API: cannot evaluate nrounds submodels.")
}

pt_check_xgboost <- function() {
  pt_require_packages()
  message("Checking XGBoost ", utils::packageVersion("xgboost"), " prediction compatibility...")
  x <- cbind(a = seq(-2, 2, length.out = 32), b = rep(c(-1, 1), 16))
  truth <- ifelse(x[, 1] > 0, "1", "2")
  params <- data.frame(nrounds = 3L, eta = .3, max_depth = 2L, gamma = 0,
                       colsample_bytree = 1, min_child_weight = 1, subsample = 1)
  full <- pt_xgb_fit(x, truth, params, 1701L, 1L)
  for (n in 1:3) {
    shorter <- params; shorter$nrounds <- n
    fit <- if (n == 3L) full else pt_xgb_fit(x, truth, shorter, 1701L, 1L)
    if (!isTRUE(all.equal(pt_xgb_predict(full, x, n), pt_xgb_predict(fit, x), tolerance = 1e-7))) {
      stop("XGBoost round-selection self-check failed; do not start the benchmark.")
    }
  }
  invisible(TRUE)
}

pt_inner_tune <- function(x, truth, feature_type, grid, folds, seed, threads,
                          preprocessing = "training_fold") {
  if (any(table(factor(truth, levels = c("1", "2"))) < folds)) {
    stop("Each class needs at least inner_folds training patients; reduce folds only explicitly.")
  }
  validation <- pr_with_seed(seed, caret::createFolds(factor(truth, levels = c("1", "2")),
                                                    k = folds, list = TRUE, returnTrain = FALSE))
  accuracy <- matrix(NA_real_, nrow(grid), length(validation))
  colnames(accuracy) <- names(validation)
  group_columns <- setdiff(names(grid), "nrounds")
  group_key <- do.call(paste, c(grid[group_columns], sep = "|"))
  groups <- split(seq_len(nrow(grid)), factor(group_key, levels = unique(group_key)))
  audit <- vector("list", length(validation))
  for (f in seq_along(validation)) {
    message("  Inner fold ", f, "/", length(validation), ": ",
            if (preprocessing == "training_fold") "fit filters; " else "keep archived predictors; ",
            length(groups), " model paths")
    valid_idx <- validation[[f]]; train_idx <- setdiff(seq_len(nrow(x)), valid_idx)
    filter <- pt_fit_preprocessor(x[train_idx, , drop = FALSE], feature_type,
                                  preprocessing = preprocessing)
    train_x <- pt_apply_preprocessor(x[train_idx, , drop = FALSE], filter)
    valid_x <- pt_apply_preprocessor(x[valid_idx, , drop = FALSE], filter)
    audit[[f]] <- list(fold = names(validation)[f], preprocessor = filter,
                       validation_ids = rownames(x)[valid_idx],
                       validation_ids_hash = pt_hash(rownames(x)[valid_idx]))
    for (g in seq_along(groups)) {
      idx <- groups[[g]]
      longest <- idx[which.max(grid$nrounds[idx])]
      fit <- pt_xgb_fit(train_x, truth[train_idx], grid[longest, , drop = FALSE],
                        pt_seed(seed, f * 10000L + g), threads)
      for (candidate in idx) {
        prob <- pt_xgb_predict(fit, valid_x, grid$nrounds[candidate])
        accuracy[candidate, f] <- mean(ifelse(prob >= .5, "1", "2") == truth[valid_idx])
      }
      rm(fit)
    }
  }
  results <- cbind(grid, mean_accuracy = rowMeans(accuracy),
                   sd_accuracy = apply(accuracy, 1L, stats::sd), as.data.frame(accuracy))
  best <- which.max(results$mean_accuracy)
  list(best_parameters = grid[best, , drop = FALSE], cv_results = results,
       folds = audit, selection_metric = "Accuracy")
}

pt_atomic_save <- function(value, path) {
  tmp <- tempfile(pattern = "checkpoint-", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(value, tmp, compress = FALSE)
  if (!file.rename(tmp, path)) stop("Could not save checkpoint: ", path)
}

pt_run_prediction <- function(prepared, framework, output_dir, seeds = 1:20,
                              inner_folds = 10L, grid = pt_default_grid(),
                              n_boot = 2000L, bootstrap_seed = 47001L,
                              threads = 2L, resume = TRUE,
                              preprocessing = "training_fold") {
  pt_require_packages()
  preprocessing <- match.arg(preprocessing, c("archived_final", "training_fold"))
  prepared <- pt_validate_prepared(prepared)
  pr_check_text(framework, "framework")
  if (length(framework) != 1L) stop("Supply one framework per run.")
  for (name in c("inner_folds", "n_boot", "threads")) {
    value <- get(name)
    pr_check_integers(value, name, minimum = if (name == "threads") 1L else 2L)
    if (length(value) != 1L) stop(name, " must be a scalar.")
  }
  pr_check_integers(bootstrap_seed, "bootstrap_seed")
  if (length(bootstrap_seed) != 1L) stop("bootstrap_seed must be a scalar.")
  grid <- pt_validate_grid(grid)
  manifests <- pr_repeated_splits(prepared$patients, seeds)
  manifests$truth <- prepared$patients$truth[match(manifests$patient_id, prepared$patients$patient_id)]
  for (split in unique(manifests$split_id)) {
    m <- manifests[manifests$split_id == split, , drop = FALSE]
    if (any(table(factor(m$truth[m$partition == "train"], levels = c("1", "2"))) < inner_folds) ||
        !setequal(m$truth[m$partition == "test"], c("1", "2"))) {
      stop("Split ", split, " lacks sufficient training/test classes; inspect cohort and explicit seeds.")
    }
  }
  definitions <- sort(ls(envir = environment(pt_run_prediction), pattern = "^(pt|pr)_"))
  source_bodies <- lapply(definitions, function(name) {
    value <- get(name, envir = environment(pt_run_prediction))
    if (is.function(value)) paste(deparse(value), collapse = "\n") else NULL
  })
  packages <- c("caret", "digest", "xgboost")
  settings <- list(framework = framework, seeds = seeds, inner_folds = inner_folds,
                   grid = grid, n_boot = n_boot, bootstrap_seed = bootstrap_seed,
                   threads = threads, positive_class = "1", cutoff = .5,
                   selection_metric = "Accuracy", preprocessing = preprocessing,
                   low_cv_fraction = if (preprocessing == "training_fold") .2 else NA_real_,
                   correlation_cutoff = if (preprocessing == "training_fold") .90 else NA_real_,
                   scaling = if (preprocessing == "training_fold") "training RMS; no centering" else "archived predictors unchanged",
                   package_versions = setNames(vapply(packages, function(p) {
                     as.character(utils::packageVersion(p))
                   }, character(1)), packages), R_version = as.character(getRversion()))
  data_hash <- pt_data_hash(prepared)
  signature <- pt_hash(list(data_hash = data_hash, settings = settings, source = source_bodies))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- file.path(output_dir, "run_manifest.rds")
  if (file.exists(manifest_path)) {
    existing <- readRDS(manifest_path)
    if (!identical(existing$signature, signature)) {
      stop("Existing output belongs to different inputs/settings/code. Choose a new output_dir; stale checkpoints will not be reused.")
    }
    if (!resume) stop("Output directory already contains this run. Use resume=TRUE or a new output_dir.")
  } else {
    if (length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
      stop("Output directory is nonempty without a matching manifest. Choose a new output_dir.")
    }
    pt_atomic_save(list(signature = signature, settings = settings,
                        data_hash = data_hash, source_hash = pt_hash(source_bodies),
                        created = as.character(Sys.time()), session = capture.output(sessionInfo())), manifest_path)
  }
  utils::write.csv(manifests, file.path(output_dir, "split_manifest.csv"), row.names = FALSE)
  checkpoint_dir <- file.path(output_dir, "checkpoints")
  dir.create(checkpoint_dir, showWarnings = FALSE)
  completed <- vector("list", length(seeds))
  for (k in seq_along(seeds)) {
    path <- file.path(checkpoint_dir, sprintf("split_%03d_seed_%d.rds", k, seeds[k]))
    if (resume && file.exists(path)) {
      cached <- readRDS(path)
      if (!identical(cached$signature, signature)) stop("Checkpoint signature mismatch: ", path)
      message(framework, ": resume completed split ", k, "/", length(seeds))
      completed[[k]] <- cached
      next
    }
    message(framework, ": fitting split ", k, "/", length(seeds), " (seed ", seeds[k], ")")
    m <- manifests[manifests$split_id == as.character(k), , drop = FALSE]
    train <- m$patient_id[m$partition == "train"]
    test <- m$patient_id[m$partition == "test"]
    train_x <- prepared$x[train, , drop = FALSE]
    truth <- prepared$patients$truth[match(train, prepared$patients$patient_id)]
    tuning <- pt_inner_tune(train_x, truth, prepared$feature_type, grid,
                            inner_folds, pt_seed(seeds[k], 1000L), threads, preprocessing)
    filter <- pt_fit_preprocessor(train_x, prepared$feature_type, preprocessing = preprocessing)
    message("  Refit selected model on the complete outer training set")
    model <- pt_xgb_fit(pt_apply_preprocessor(train_x, filter), truth,
                        tuning$best_parameters, pt_seed(seeds[k], 1000000L), threads)
    probability <- pt_xgb_predict(model, pt_apply_preprocessor(prepared$x[test, , drop = FALSE], filter))
    patients <- prepared$patients[match(test, prepared$patients$patient_id), , drop = FALSE]
    predictions <- data.frame(framework = framework, split_id = as.character(k),
                              split_seed = seeds[k], patients,
                              probability_cluster1 = probability,
                              predicted = ifelse(probability >= .5, "1", "2"), row.names = NULL)
    boot <- pr_bootstrap_test_f1(predictions, n_boot, pt_seed(bootstrap_seed, k - 1L))
    filters <- do.call(rbind, c(lapply(tuning$folds, function(f) {
      cbind(framework = framework, split_id = as.character(k), fold = f$fold,
            training_ids_hash = f$preprocessor$training_ids_hash, f$preprocessor$counts)
    }), list(cbind(framework = framework, split_id = as.character(k), fold = "outer_refit",
                    training_ids_hash = filter$training_ids_hash, filter$counts))))
    result <- list(signature = signature, framework = framework, split_id = as.character(k),
                   split_seed = seeds[k], predictions = predictions,
                   best_parameters = cbind(framework = framework, split_id = as.character(k), tuning$best_parameters),
                   cv_results = cbind(framework = framework, split_id = as.character(k), tuning$cv_results),
                   inner_audit = tuning$folds, outer_preprocessor = filter,
                   outer_training_ids = train, outer_test_ids = test,
                   filters = filters, bootstrap = boot$bootstrap,
                   split_metrics = cbind(split_seed = seeds[k], boot$summary),
                   model_raw = xgboost::xgb.save.raw(model))
    pt_atomic_save(result, path)
    completed[[k]] <- result
  }
  bind <- function(name) do.call(rbind, lapply(completed, `[[`, name))
  metrics <- bind("split_metrics")
  f1 <- metrics$f1
  summary <- data.frame(framework = framework, n_splits = length(seeds),
                         mean_f1 = mean(f1), sd_f1 = if (length(f1) > 1L) stats::sd(f1) else NA_real_,
                         median_f1 = stats::median(f1), min_f1 = min(f1), max_f1 = max(f1))
  result <- list(predictions = bind("predictions"), split_metrics = metrics,
                 split_summary = summary, bootstrap = bind("bootstrap"),
                 manifests = manifests, filters = bind("filters"),
                 best_parameters = bind("best_parameters"), cv_results = bind("cv_results"))
  for (name in names(result)) {
    utils::write.csv(result[[name]], file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
  }
  preprocessing_notes <- if (preprocessing == "training_fold") c(
    "All feature filtering and RMS scaling are refitted inside each inner training fold.",
    "CV ranking removes floor(20% * eligible features) separately for expression/methylation.",
    "Correlation pruning is joint across expression/methylation, as in the original pipeline.",
    "Other features receive only finite/constant checks and training RMS scaling.",
    "Nonfinite assessment values are passed as missing to XGBoost; no test-based filtering/imputation is fitted."
  ) else c("Archived final predictor matrix used unchanged; no new feature filtering or scaling.",
           "Split stability is conditional on the original global preprocessing; this run does not resolve that limitation.")
  writeLines(c("Positive class: Cluster 1. Outer hold-out: 80/20 by cancer type.",
               "Inner CV: class-stratified; Accuracy selects parameters; probability cutoff = 0.5.",
               preprocessing_notes,
               "F1 intervals: 95% percentile bootstrap of test patients within cancer type, conditional on the fitted model.",
               "Across-split SD/range describe refitting variability; overlapping splits are not independent studies.",
               "Frameworks predict different cluster targets; their F1 values are not a biological quality ranking.",
               paste("Run signature:", signature)), file.path(output_dir, "evaluation_notes.txt"))
  result
}
