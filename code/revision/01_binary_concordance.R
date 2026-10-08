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
  # Chance reference: agreement expected if the two binary states were
  # independent given their observed active prevalences.
  p1 <- active[first] / nrow(binary); p2 <- active[second] / nrow(binary)
  expected_n11 <- nrow(binary) * p1 * p2
  observed <- agreement / nrow(binary)
  expected <- p1 * p2 + (1 - p1) * (1 - p2)
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
    active_jaccard = bc_safe_ratio(n11, n11 + discordant),
    expected_n11 = unname(expected_n11),
    observed_expected_n11 = bc_safe_ratio(n11, expected_n11),
    cohen_kappa = bc_safe_ratio(observed - expected, 1 - expected),
    active_prevalence1 = unname(p1), active_prevalence2 = unname(p2),
    stringsAsFactors = FALSE)
}

# Rank of each principal pair among the 48 cross-compendium pairs of the 12
# selected signatures (1 = strongest), per rule, for the patient-state and the
# active-state Jaccard.
bc_principal_ranks <- function(pairs, principal,
                               metrics = c("state_jaccard", "active_jaccard")) {
  subset <- pairs[pairs$cross_framework & pairs$high_activity_pair, ]
  key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "/")
  subset$pair <- key(subset$signature1, subset$signature2)
  wanted <- vapply(principal, function(x) key(x[1], x[2]), character(1))
  ranks <- list()
  for (metric in metrics) {
    for (rule in unique(subset$rule)) {
      x <- subset[subset$rule == rule, ]
      x$rank <- rank(-x[[metric]], ties.method = "min", na.last = "keep")
      hit <- x[x$pair %in% wanted, ]
      ranks[[length(ranks) + 1L]] <- data.frame(metric = metric, rule = rule,
        selection = names(wanted)[match(hit$pair, wanted)], pair = hit$pair,
        value = hit[[metric]], rank = hit$rank, n_pairs = nrow(x))
    }
  }
  do.call(rbind, ranks)
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
      union <- sum(a == 1 | b == 1)
      expected_n11 <- sum(a) * sum(b) / length(a)
      pe <- mean(a) * mean(b) + mean(1 - a) * mean(1 - b)
      direct <- c(sum(same) / (2 * length(a) - sum(same)),
                  if (any(same)) mean(a[same] == 1) else NA_real_,
                  if (any(same)) mean(a[same] == 0) else NA_real_,
                  if (union > 0) sum(a == 1 & b == 1) / union else NA_real_,
                  if (expected_n11 > 0) sum(a == 1 & b == 1) / expected_n11 else NA_real_,
                  if (pe < 1) (mean(same) - pe) / (1 - pe) else NA_real_)
      stopifnot(isTRUE(all.equal(
        as.numeric(metrics[i, c("state_jaccard", "active_share_of_agreements",
                                "inactive_share_of_agreements", "active_jaccard",
                                "observed_expected_n11", "cohen_kappa")]), direct)))
    }
  }
  invisible(TRUE)
}

# One agreement matrix in the layout of Figure 4: signatures in the given
# order on both axes (bottom to top), grey diagonal and grey undefined values.
# The colour ramp keeps the blue family and tile outlines of Figure 4 but,
# unlike its white plateau below .25, also distinguishes low values, which
# matters for active-state Jaccard. All panels share the [0, 1] mapping.
bc_matrix_panel <- function(pairs, order, metric, title, rule = "original", text_size = 6) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for plots.")
  data <- pairs[pairs$rule == rule & pairs$signature1 %in% order & pairs$signature2 %in% order, ]
  reverse <- data
  reverse$signature1 <- data$signature2; reverse$signature2 <- data$signature1
  long <- rbind(data.frame(x = data$signature1, y = data$signature2, value = data[[metric]]),
                data.frame(x = reverse$signature1, y = reverse$signature2, value = reverse[[metric]]),
                data.frame(x = order, y = order, value = NA_real_))
  long$x <- factor(long$x, levels = order); long$y <- factor(long$y, levels = order)
  ggplot2::ggplot(long, ggplot2::aes(x, y, fill = value)) +
    ggplot2::geom_tile(color = "grey25", linewidth = .1) +
    ggplot2::scale_fill_gradientn(
      colors = c("white", "#deebf7", "lightblue", "#3182bd", "darkblue", "darkblue"),
      values = c(0, .1, .25, .5, .75, 1), limits = c(0, 1),
      breaks = seq(0, 1, .25), na.value = "grey50", name = "Value") +
    ggplot2::coord_fixed(expand = FALSE) +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    ggplot2::theme_minimal(base_size = 7) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 7, hjust = .5),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = .5, size = text_size,
                                          colour = "black"),
      axis.text.y = ggplot2::element_text(size = text_size, colour = "black"))
}

# Compact threshold-sensitivity summary: distribution of each metric over the
# 48 cross-compendium pairs of the selected signatures, per rule, with the
# cross-compendium pairs discussed in the manuscript (CN1/Sig1, CN1/Sig2,
# CX1/CN1) highlighted and labelled by their rank among the 48.
bc_plot_threshold_sensitivity <- function(pairs, ranks, rules,
    metrics = c(state_jaccard = "Patient-state Jaccard",
                active_jaccard = "Active-state Jaccard")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for plots.")
  subset <- pairs[pairs$cross_framework & pairs$high_activity_pair & pairs$rule %in% rules, ]
  long <- do.call(rbind, lapply(names(metrics), function(m)
    data.frame(rule = subset$rule, metric = metrics[[m]], value = subset[[m]])))
  hits <- ranks[ranks$metric %in% names(metrics) & ranks$rule %in% rules, ]
  hits$metric <- metrics[hits$metric]
  rule_short <- c(original = "Original", positive_q25 = "Pos. Q25",
                  positive_q50 = "Pos. Q50", positive_q75 = "Pos. Q75")
  for (nm in c("long", "hits")) {
    d <- get(nm)
    d$rule <- factor(d$rule, levels = rules, labels = rule_short[rules])
    d$metric <- factor(d$metric, levels = unname(metrics))
    assign(nm, d)
  }
  ggplot2::ggplot(long, ggplot2::aes(rule, value)) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "grey92", colour = "grey55", width = .55) +
    ggplot2::geom_jitter(width = .12, height = 0, size = .7, colour = "grey60") +
    ggplot2::geom_line(data = hits, ggplot2::aes(group = pair, colour = pair), linewidth = .6) +
    ggplot2::geom_point(data = hits, ggplot2::aes(colour = pair), size = 2) +
    ggrepel::geom_text_repel(data = hits, ggplot2::aes(label = rank, colour = pair),
                             nudge_x = .32, direction = "y", size = 2.3, min.segment.length = Inf,
                             box.padding = .1, seed = 1, show.legend = FALSE) +
    ggplot2::facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    ggplot2::scale_colour_manual(values = c("#08519c", "#6baed6", "#e6550d"), name = NULL) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom",
                   axis.text = ggplot2::element_text(colour = "black"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0)) +
    # Grey: the 48 cross-compendium pairs of the 12 selected signatures;
    # numbers: rank of each principal pair among the 48 (1 = strongest).
    ggplot2::labs(x = NULL, y = NULL)
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
  principal <- list(CN1_Sig1 = c("CN1", "Sig1"), CN1_Sig2 = c("CN1", "Sig2"),
                    CX1_CN1 = c("CX1", "CN1"), CN1_CN2 = c("CN1", "CN2"),
                    CN4_CN10 = c("CN4", "CN10"))
  selected <- do.call(rbind, lapply(names(principal), function(label) {
    signatures <- principal[[label]]
    out <- pairs[pairs$signature1 %in% signatures & pairs$signature2 %in% signatures, ]
    cbind(selection = label, out)
  }))
  principal_ranks <- bc_principal_ranks(pairs, principal)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(thresholds, file.path(output_dir, "thresholds.csv"), row.names = FALSE, na = "NA")
  write.csv(pairs, file.path(output_dir, "pair_metrics.csv"), row.names = FALSE, na = "NA")
  write.csv(selected, file.path(output_dir, "principal_pairs.csv"), row.names = FALSE, na = "NA")
  write.csv(principal_ranks, file.path(output_dir, "principal_pair_ranks.csv"),
            row.names = FALSE, na = "NA")
  all_order <- colnames(cohort$values)[hclust(dist(t(analyses$original$binary)), method = "complete")$order]
  high <- bc_high_activity_signatures()
  high_order <- high[hclust(dist(t(analyses$original$binary[, high, drop = FALSE])), method = "complete")$order]
  if (make_plots) {
    # Supplementary figure, page width (6.8 in) at 600 dpi. The patient-state
    # Jaccard of all 58 signatures is Figure 4a; here its two components for
    # the 58 signatures (a, b) and all three measures for the 12 selected
    # signatures (c-e). Formulas are given in the caption.
    library_ok <- requireNamespace("patchwork", quietly = TRUE)
    if (!library_ok) stop("patchwork is required for the composite figure.")
    panels <- list(
      bc_matrix_panel(pairs, all_order, "inactive_share_of_agreements",
                      "Inactive-inactive share of agreements", text_size = 3.6),
      bc_matrix_panel(pairs, all_order, "active_jaccard", "Active-state Jaccard", text_size = 3.6),
      bc_matrix_panel(pairs, high_order, "state_jaccard", "Patient-state Jaccard"),
      bc_matrix_panel(pairs, high_order, "inactive_share_of_agreements", "Inactive-inactive share"),
      bc_matrix_panel(pairs, high_order, "active_jaccard", "Active-state Jaccard"))
    components <- (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]] | panels[[5]]) +
      patchwork::plot_layout(heights = c(1.45, 1), guides = "collect") +
      patchwork::plot_annotation(tag_levels = "a") &
      ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 9),
                     legend.position = "bottom", legend.key.height = grid::unit(5, "pt"),
                     legend.key.width = grid::unit(28, "pt"), legend.title.position = "top",
                     legend.title = ggplot2::element_text(hjust = .5))
    for (ext in c("png", "tiff", "pdf")) {
      args <- list(filename = file.path(output_dir, paste0("state_agreement_components.", ext)),
                   plot = components, width = 6.8, height = 6.7, units = "in", bg = "white")
      if (ext != "pdf") args$dpi <- 600
      if (ext == "tiff") args$compression <- "lzw"
      do.call(ggplot2::ggsave, args)
    }
    sensitivity <- bc_plot_threshold_sensitivity(pairs, principal_ranks, rules)
    old_heatmap <- file.path(output_dir, "threshold_sensitivity_high_activity.png")
    if (file.exists(old_heatmap)) file.remove(old_heatmap)
    for (ext in c("png", "tiff", "pdf")) {
      args <- list(filename = file.path(output_dir, paste0("threshold_sensitivity_principal_pairs.", ext)),
                   plot = sensitivity, width = 6.8, height = 3.6, units = "in", bg = "white")
      if (ext != "pdf") args$dpi <- 600
      if (ext == "tiff") args$compression <- "lzw"
      do.call(ggplot2::ggsave, args)
    }
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
    "Active-state J = n11/(n11+n10+n01). Expected n11 = N*p1*p2 under independence of the two binary states; Cohen kappa uses the same marginal prevalences.",
    "Principal-pair ranks are computed among the same 48 pairs (1 = strongest, ties share the minimum rank).",
    "The selected 12 signatures reproduce Figure 4b; selection does not denote active-only agreement.",
    "Heatmap ordering uses complete-linkage clustering of Euclidean distances between original binary signature columns."),
    file.path(output_dir, "provenance.txt"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
  message("Binary concordance outputs: ", normalizePath(output_dir),
          " (", nrow(cohort$values), " patients; ", ncol(cohort$values),
          " signatures; ", length(rules), " rules).")
  invisible(list(thresholds = thresholds, pairs = pairs,
                 principal_pairs = selected, principal_ranks = principal_ranks,
                 cohort = cohort, analyses = analyses))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (any(!args %in% "--no-plots")) stop("Usage: Rscript code/revision/01_binary_concordance.R [--no-plots]")
  run_binary_concordance(make_plots = !"--no-plots" %in% args)
}
