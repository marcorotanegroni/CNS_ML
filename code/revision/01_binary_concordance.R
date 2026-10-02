# Binary concordance sensitivity on a fixed, jointly complete matched cohort.
# Source this file to define functions; run Rscript from the repository root
# to write the analysis outputs. No input objects are modified.

bc_high_activity_signatures <- function() {
  c("Sig1", "Sig2", "Sig5", "Sig7", "CX1", "CX2", "CX3", "CX5",
    "CN1", "CN2", "CN9", "CN17")
}

bc_rule_labels <- function() {
  c(original = "Original full-cohort zero/median",
    positive_q25 = "Positive-exposure 25th percentile",
    positive_q50 = "Positive-exposure median",
    positive_q75 = "Positive-exposure 75th percentile")
}

bc_validate_matrix <- function(x, name) {
  if (!is.data.frame(x) || !identical(names(x)[1L], "cancer_type")) {
    stop(name, " must be a data frame with cancer_type first.")
  }
  ids <- rownames(x)
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop(name, " has absent, empty, or duplicated sample identifiers.")
  }
  if (anyDuplicated(names(x)) || ncol(x) < 2L ||
      !all(vapply(x[-1L], is.numeric, logical(1)))) {
    stop(name, " must have unique, numeric signature columns.")
  }
  values <- as.matrix(x[-1L])
  observed <- values[!is.na(values)]
  if (any(!is.finite(observed)) || any(observed < 0 | observed > 1)) {
    stop(name, " signature exposures must be finite fractions in [0, 1].")
  }
  invisible(TRUE)
}

bc_load_cohort <- function(input_path) {
  env <- new.env(parent = emptyenv())
  load(input_path, envir = env)
  full_names <- c(Tao = "tao", Drews = "drews", Steele = "steele")
  matched_names <- c(Tao = "tao_drews", Drews = "drews_tao",
                     Steele = "steele_drews")
  needed <- c(full_names, matched_names)
  if (!all(needed %in% ls(env))) stop("Missing required signature matrices.")
  full <- lapply(full_names, get, envir = env)
  matched <- lapply(matched_names, get, envir = env)
  for (i in seq_along(full)) {
    bc_validate_matrix(full[[i]], full_names[[i]])
    bc_validate_matrix(matched[[i]], matched_names[[i]])
  }
  ids <- rownames(matched$Drews)
  for (i in seq_along(matched)) {
    if (!setequal(ids, rownames(matched[[i]]))) {
      stop("Matched signature matrices have different sample sets.")
    }
    if (!all(ids %in% rownames(full[[i]]))) {
      stop("Matched samples are absent from a full reference cohort.")
    }
    matched[[i]] <- matched[[i]][ids, , drop = FALSE]
    if (!isTRUE(all.equal(matched[[i]], full[[i]][ids, , drop = FALSE],
                         tolerance = 1e-12))) {
      stop("A matched matrix differs from its full-cohort source.")
    }
  }
  cancer_type <- as.character(matched$Drews$cancer_type)
  for (x in matched) {
    if (!identical(cancer_type, as.character(x$cancer_type))) {
      stop("Cancer type is not aligned across matched signature matrices.")
    }
  }
  complete <- Reduce(`&`, lapply(matched, complete.cases)) &
    !is.na(cancer_type) & nzchar(cancer_type)
  if (!any(complete)) stop("No jointly complete matched samples remain.")
  matched <- lapply(matched, function(x) x[complete, , drop = FALSE])
  values <- do.call(cbind, lapply(matched, function(x) as.matrix(x[-1L])))
  if (anyDuplicated(colnames(values))) stop("Signature names overlap.")
  framework <- rep(names(matched), vapply(matched, ncol, integer(1)) - 1L)
  names(framework) <- colnames(values)
  if (!all(bc_high_activity_signatures() %in% colnames(values))) {
    stop("The manuscript's selected signatures are unavailable.")
  }
  list(full = full, values = values, framework = framework,
       samples = data.frame(sample_id = ids[complete],
                            cancer_type = cancer_type[complete]),
       initial_n = length(ids), excluded_n = sum(!complete))
}

bc_threshold_rule <- function(values, full, framework, rule) {
  binary <- matrix(0L, nrow(values), ncol(values), dimnames = dimnames(values))
  rows <- vector("list", ncol(values))
  for (i in seq_len(ncol(values))) {
    signature <- colnames(values)[i]
    reference <- full[[framework[[signature]]]][[signature]]
    observed <- reference[!is.na(reference)]
    positive <- observed[observed > 0]
    if (!length(observed)) stop("No observed reference exposures: ", signature)
    threshold <- switch(rule,
      original = median(observed),
      positive_q25 = if (length(positive)) unname(quantile(positive, .25)) else NA_real_,
      positive_q50 = if (length(positive)) unname(quantile(positive, .50)) else NA_real_,
      positive_q75 = if (length(positive)) unname(quantile(positive, .75)) else NA_real_)
    if (is.null(threshold)) stop("Unknown threshold rule: ", rule)
    positive_rule <- startsWith(rule, "positive_")
    strict <- !positive_rule && threshold == 0
    binary[, i] <- if (is.na(threshold)) 0L else if (strict)
      as.integer(values[, i] > threshold) else
      as.integer(values[, i] > 0 & values[, i] >= threshold)
    rows[[i]] <- data.frame(
      rule = rule, rule_label = unname(bc_rule_labels()[rule]),
      signature = signature, framework = framework[[signature]],
      estimand = if (strict) "detectable_positive_exposure" else
        if (positive_rule) "relative_exposure_among_positive_reference_samples" else
          "relative_exposure_in_reference_cohort",
      reference_cohort = "full_compendium",
      reference_n = length(reference), reference_observed_n = length(observed),
      reference_positive_n = length(positive), threshold = threshold,
      comparison = if (is.na(threshold)) "no positive reference exposures; all inactive" else
        if (strict) "x > threshold" else "x > 0 and x >= threshold",
      quantile_type = if (positive_rule) 7L else NA_integer_,
      reference_at_threshold_n = if (is.na(threshold)) NA_integer_ else sum(observed == threshold),
      matched_at_threshold_n = if (is.na(threshold)) NA_integer_ else sum(values[, i] == threshold),
      matched_n = nrow(values), active_n = sum(binary[, i]),
      active_prevalence = mean(binary[, i]), stringsAsFactors = FALSE)
  }
  list(binary = binary, thresholds = do.call(rbind, rows))
}

bc_safe_ratio <- function(numerator, denominator) {
  result <- rep(NA_real_, length(denominator))
  nonzero <- denominator > 0
  result[nonzero] <- numerator[nonzero] / denominator[nonzero]
  result
}

bc_pair_metrics <- function(binary, framework, rule) {
  if (!is.matrix(binary) || anyNA(binary) || any(!binary %in% c(0, 1)) ||
      nrow(binary) == 0L || ncol(binary) < 2L) stop("Invalid binary matrix.")
  pairs <- combn(seq_len(ncol(binary)), 2L)
  first <- pairs[1L, ]; second <- pairs[2L, ]
  active <- colSums(binary)
  n11 <- crossprod(binary)[cbind(first, second)]
  n10 <- active[first] - n11
  n01 <- active[second] - n11
  n00 <- nrow(binary) - n11 - n10 - n01
  discordant <- n10 + n01
  agreement <- n11 + n00
  names_first <- colnames(binary)[first]; names_second <- colnames(binary)[second]
  high <- bc_high_activity_signatures()
  data.frame(rule = rule, signature1 = names_first, signature2 = names_second,
    framework1 = unname(framework[names_first]),
    framework2 = unname(framework[names_second]),
    cross_framework = unname(framework[names_first] != framework[names_second]),
    high_activity_pair = names_first %in% high & names_second %in% high,
    n = nrow(binary), n_concordant = unname(agreement), n11 = unname(n11), n00 = unname(n00),
    n10 = unname(n10), n01 = unname(n01),
    state_jaccard = bc_safe_ratio(agreement, agreement + 2 * discordant),
    active_share_of_agreements = bc_safe_ratio(n11, agreement),
    inactive_share_of_agreements = bc_safe_ratio(n00, agreement),
    active_prevalence1 = unname(active[first] / nrow(binary)),
    active_prevalence2 = unname(active[second] / nrow(binary)),
    stringsAsFactors = FALSE)
}

bc_stability <- function(pairs) {
  baseline <- pairs[pairs$rule == "original", ]
  result <- list()
  for (rule in setdiff(unique(pairs$rule), "original")) {
    alternative <- pairs[pairs$rule == rule, ]
    stopifnot(identical(baseline$signature1, alternative$signature1),
              identical(baseline$signature2, alternative$signature2))
    include <- baseline$cross_framework & baseline$high_activity_pair
    for (metric in "state_jaccard") {
      a <- baseline[[metric]][include]; b <- alternative[[metric]][include]
      valid <- is.finite(a) & is.finite(b)
      rho <- if (sum(valid) >= 2L && length(unique(a[valid])) > 1L &&
                 length(unique(b[valid])) > 1L)
        cor(a[valid], b[valid], method = "spearman") else NA_real_
      result[[length(result) + 1L]] <- data.frame(
        reference_rule = "original", alternative_rule = rule,
        subset = "high_activity_cross_framework", metric = metric, n_pairs = sum(include),
        n_comparable = sum(valid), n_undefined_reference = sum(!is.finite(a)),
        n_undefined_alternative = sum(!is.finite(b)), spearman_rho = rho,
        median_absolute_change = if (any(valid)) median(abs(b[valid] - a[valid])) else NA_real_)
    }
  }
  do.call(rbind, result)
}

bc_verify_metrics <- function(original) {
  # Independent direct count reference, including pairs with no agreements.
  small <- cbind(a = c(1, 1, 0, 0), b = c(1, 0, 1, 0),
                 empty = c(0, 0, 0, 0), full = c(1, 1, 1, 1),
                 empty2 = c(0, 0, 0, 0), full2 = c(1, 1, 1, 1))
  for (x in list(small, original)) {
    framework <- setNames(rep("test", ncol(x)), colnames(x))
    metrics <- bc_pair_metrics(x, framework, "reference")
    for (i in seq_len(nrow(metrics))) {
      a <- x[, metrics$signature1[i]]; b <- x[, metrics$signature2[i]]
      counts <- c(sum(a == 1 & b == 1), sum(a == 0 & b == 0),
                  sum(a == 1 & b == 0), sum(a == 0 & b == 1))
      stopifnot(identical(as.numeric(metrics[i, c("n11", "n00", "n10", "n01")]),
                          as.numeric(counts)))
      same <- a == b
      direct <- c(sum(same) / (2 * length(a) - sum(same)),
                  if (any(same)) mean(a[same] == 1) else NA_real_,
                  if (any(same)) mean(a[same] == 0) else NA_real_)
      stopifnot(isTRUE(all.equal(
        as.numeric(metrics[i, c("state_jaccard", "active_share_of_agreements",
                                "inactive_share_of_agreements")]), direct)))
    }
  }
  invisible(TRUE)
}

bc_plot_heatmaps <- function(pairs, signatures, order, rules,
                             metrics = c(state_jaccard = "Patient-state Jaccard",
                                         inactive_share_of_agreements = "Inactive-inactive share of agreements"),
                             scale_name = "Value",
                             caption = "Fixed signature ordering. White diagonal: omitted self-comparisons. Grey: share undefined because no patients agree.") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for plots.")
  data <- pairs[pairs$rule %in% rules & pairs$signature1 %in% signatures &
                  pairs$signature2 %in% signatures, ]
  reverse <- data
  reverse$signature1 <- data$signature2; reverse$signature2 <- data$signature1
  data <- rbind(data, reverse)
  long <- do.call(rbind, lapply(names(metrics), function(metric) {
    data.frame(rule = data$rule, signature1 = data$signature1,
               signature2 = data$signature2, metric = metrics[[metric]], value = data[[metric]])
  }))
  long$signature1 <- factor(long$signature1, levels = order)
  long$signature2 <- factor(long$signature2, levels = rev(order))
  long$metric <- factor(long$metric, levels = unname(metrics))
  long$rule <- factor(long$rule, levels = rules, labels = unname(bc_rule_labels()[rules]))
  p <- ggplot2::ggplot(long, ggplot2::aes(signature1, signature2, fill = value)) +
    ggplot2::geom_tile(color = "grey25", linewidth = .1) +
    ggplot2::facet_grid(rule ~ metric, drop = FALSE) +
    ggplot2::scale_x_discrete(drop = FALSE) + ggplot2::scale_y_discrete(drop = FALSE) +
    # Use the blue family and tile outlines of Figure 4. Unlike its white
    # plateau below .25, this ramp also distinguishes low active agreement.
    # All panels keep the same mapping and fixed [0, 1] limits.
    ggplot2::scale_fill_gradientn(
      colors = c("white", "#deebf7", "lightblue", "#3182bd", "darkblue", "darkblue"),
      values = c(0, .1, .25, .5, .75, 1), limits = c(0, 1),
      breaks = seq(0, 1, .25), na.value = "grey50", name = scale_name) +
    ggplot2::coord_fixed() + ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = .5,
                                         size = if (length(signatures) > 12L) 5 else 8,
                                         color = "black"),
      axis.text.y = ggplot2::element_text(size = if (length(signatures) > 12L) 5 else 8,
                                         color = "black"),
      strip.text.y = ggplot2::element_text(angle = 0), legend.position = "bottom") +
    ggplot2::labs(x = NULL, y = NULL,
      title = paste0("Binary agreement: ", format(unique(data$n), big.mark = ","), " matched TCGA patients"),
      caption = caption)
  p
}

run_binary_concordance <- function(
    input_path = "data/processed/signature_exploration_data.RData",
    output_dir = "results/revision/binary_concordance", make_plots = TRUE) {
  cohort <- bc_load_cohort(input_path)
  rules <- names(bc_rule_labels())
  analyses <- lapply(rules, function(rule)
    bc_threshold_rule(cohort$values, cohort$full, cohort$framework, rule))
  names(analyses) <- rules
  # Reconstruct the notebook's dichotomization independently, before matching.
  reference <- do.call(cbind, lapply(cohort$full, function(x) {
    z <- sapply(x[-1L], function(v) {
      m <- median(v, na.rm = TRUE)
      if (m == 0) ifelse(v > m, 1L, 0L) else ifelse(v >= m, 1L, 0L)
    })
    rownames(z) <- rownames(x)
    z[rownames(cohort$values), , drop = FALSE]
  }))
  stopifnot(identical(analyses$original$binary, reference))
  bc_verify_metrics(analyses$original$binary)
  thresholds <- do.call(rbind, lapply(analyses, `[[`, "thresholds"))
  pairs <- do.call(rbind, lapply(rules, function(rule)
    bc_pair_metrics(analyses[[rule]]$binary, cohort$framework, rule)))
  rownames(thresholds) <- NULL; rownames(pairs) <- NULL
  stability <- bc_stability(pairs)
  principal <- list(CN1_Sig1 = c("CN1", "Sig1"), CN1_Sig2 = c("CN1", "Sig2"),
                    CX1_CN1 = c("CX1", "CN1"), CN1_CN2 = c("CN1", "CN2"))
  selected <- do.call(rbind, lapply(names(principal), function(label) {
    signatures <- principal[[label]]
    out <- pairs[pairs$signature1 %in% signatures & pairs$signature2 %in% signatures, ]
    cbind(selection = label, out)
  }))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(thresholds, file.path(output_dir, "thresholds.csv"), row.names = FALSE, na = "NA")
  write.csv(pairs, file.path(output_dir, "pair_metrics.csv"), row.names = FALSE, na = "NA")
  write.csv(stability, file.path(output_dir, "rank_stability.csv"), row.names = FALSE, na = "NA")
  write.csv(selected, file.path(output_dir, "principal_pairs.csv"), row.names = FALSE, na = "NA")
  all_order <- colnames(cohort$values)[hclust(dist(t(analyses$original$binary)), method = "complete")$order]
  high <- bc_high_activity_signatures()
  high_order <- high[hclust(dist(t(analyses$original$binary[, high, drop = FALSE])), method = "complete")$order]
  if (make_plots) {
    components <- c(active_share_of_agreements = "Active-active: n11 / (n11 + n00)",
                    inactive_share_of_agreements = "Inactive-inactive: n00 / (n11 + n00)")
    component_caption <- "Denominator: concordant patients for each pair. Shares sum to 1 when defined. Interpret alongside overall Jaccard."
    all_components <- bc_plot_heatmaps(pairs, all_order, all_order, "original",
      components, scale_name = "Share of agreements", caption = NULL) +
      ggplot2::labs(title = "All 58 signatures: original zero/median rule") +
      ggplot2::theme(strip.text.y = ggplot2::element_blank(), legend.position = "none")
    high_components <- bc_plot_heatmaps(pairs, high, high_order, "original",
      components, scale_name = "Share of agreements", caption = component_caption) +
      ggplot2::labs(title = "12 signatures selected for clustering: original zero/median rule") +
      ggplot2::theme(strip.text.y = ggplot2::element_blank())
    grDevices::png(file.path(output_dir, "state_agreement_components.png"),
                   width = 3600, height = 3600, res = 180, bg = "white")
    tryCatch({
      grid::grid.newpage()
      grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1)))
      print(all_components, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
      print(high_components, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    }, finally = grDevices::dev.off())
    sensitivity <- bc_plot_heatmaps(pairs, high, high_order, rules)
    ggplot2::ggsave(file.path(output_dir, "threshold_sensitivity_high_activity.png"),
                   sensitivity, width = 12, height = 20, units = "in", dpi = 180,
                   bg = "white")
  }
  writeLines(c(paste0("Input: ", normalizePath(input_path)),
    paste0("Input MD5: ", unname(tools::md5sum(input_path))),
    paste0("Matched samples before joint complete-case filtering: ", cohort$initial_n),
    paste0("Excluded samples: ", cohort$excluded_n),
    paste0("Analyzed samples: ", nrow(cohort$values)),
    paste0("Signatures: ", ncol(cohort$values)),
    "Verification: original dichotomization exactly reproduced; all pair metrics checked against direct counts.",
    "Positive quantiles use R quantile type 7; exposures tied at a positive threshold are active.",
    "Pairs with no concordant patients have undefined agreement shares (NA); an undefined positive quantile classifies all samples as inactive.",
    "State J = (n11+n00)/(n11+n00+2*(n10+n01)).",
    "Composition of observed agreements: active share = n11/(n11+n00); inactive share = n00/(n11+n00). Shares sum to 1 when defined.",
    "Rank stability summarizes the 48 cross-compendium pairs among the selected signatures, using finite pairs only.",
    "A constant vector or fewer than two comparable pairs yields NA for rank correlation.",
    "The selected 12 signatures reproduce Figure 4b; selection does not denote active-only agreement.",
    "Heatmap ordering uses complete-linkage clustering of Euclidean distances between original binary signature columns."),
    file.path(output_dir, "provenance.txt"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
  message("Binary concordance outputs: ", normalizePath(output_dir),
          " (", nrow(cohort$values), " patients; ", ncol(cohort$values),
          " signatures; ", length(rules), " rules).")
  invisible(list(thresholds = thresholds, pairs = pairs, stability = stability,
                 principal_pairs = selected, cohort = cohort, analyses = analyses))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (any(!args %in% "--no-plots")) stop("Usage: Rscript code/revision/01_binary_concordance.R [--no-plots]")
  run_binary_concordance(make_plots = !"--no-plots" %in% args)
}
