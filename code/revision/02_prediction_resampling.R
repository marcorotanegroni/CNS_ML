# Resampling utilities for existing binary cluster-classification targets.
# Source to define functions only. No data are loaded and no model is fitted.
# Labels must be explicitly encoded as "1" (Cluster 1, positive) and "2".

pr_check_text <- function(x, name) {
  if (!(is.character(x) || is.factor(x)) || !length(x) ||
      anyNA(x) || any(!nzchar(trimws(as.character(x))))) {
    stop(name, " must contain nonempty, nonmissing text values.")
  }
  invisible(TRUE)
}

pr_check_integers <- function(x, name, minimum = 0L) {
  if (!is.numeric(x) || !length(x) || anyNA(x) ||
      any(!is.finite(x) | x < minimum | x > .Machine$integer.max |
          x != floor(x))) stop(name, " must contain valid integer values.")
  invisible(TRUE)
}

pr_with_seed <- function(seed, code) {
  pr_check_integers(seed, "seed")
  if (length(seed) != 1L) stop("Provide exactly one seed.")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  })
  set.seed(seed)
  force(code)
}

# Call separately for each framework's confirmed eligible cohort.
# Singleton cancer types are retained in training and explicitly flagged.
pr_repeated_splits <- function(patients, seeds) {
  required <- c("patient_id", "cancer_type")
  if (!is.data.frame(patients) || !all(required %in% names(patients))) {
    stop("patients must contain patient_id and cancer_type.")
  }
  for (nm in required) pr_check_text(patients[[nm]], nm)
  if (anyDuplicated(patients$patient_id)) stop("Duplicated patient_id.")
  pr_check_integers(seeds, "seeds")
  if (anyDuplicated(seeds)) stop("Split seeds must be distinct.")
  patients <- as.data.frame(lapply(patients[required], as.character),
                            stringsAsFactors = FALSE)
  patients <- patients[order(patients$cancer_type, patients$patient_id), ]
  strata <- split(seq_len(nrow(patients)), patients$cancer_type)
  singleton <- lengths(strata) == 1L
  if (any(singleton)) {
    warning("Singleton cancer types retained in training: ",
            paste(names(strata)[singleton], collapse = ", "))
  }
  result <- lapply(seq_along(seeds), function(k) {
    train <- pr_with_seed(seeds[k], unlist(lapply(strata, function(idx) {
      n <- length(idx)
      n_train <- if (n == 1L) 1L else max(1L, min(n - 1L, round(.8 * n)))
      idx[sample.int(n, n_train)]
    }), use.names = FALSE))
    data.frame(split_id = as.character(k), split_seed = seeds[k], patients,
               partition = ifelse(seq_len(nrow(patients)) %in% train,
                                  "train", "test"),
               singleton_stratum = patients$cancer_type %in%
                 names(strata)[singleton], row.names = NULL)
  })
  do.call(rbind, result)
}

pr_f1 <- function(truth, predicted) {
  if (!length(truth) || length(truth) != length(predicted) ||
      anyNA(truth) || anyNA(predicted) ||
      any(!as.character(truth) %in% c("1", "2")) ||
      any(!as.character(predicted) %in% c("1", "2"))) {
    stop("truth and predicted must have equal lengths and only labels 1/2.")
  }
  truth <- as.character(truth); predicted <- as.character(predicted)
  tp <- sum(truth == "1" & predicted == "1")
  fn <- sum(truth == "1" & predicted == "2")
  fp <- sum(truth == "2" & predicted == "1")
  tn <- sum(truth == "2" & predicted == "2")
  denominator <- 2 * tp + fp + fn
  data.frame(tp = tp, fp = fp, fn = fn, tn = tn,
             f1 = if (denominator == 0) NA_real_ else 2 * tp / denominator)
}

# Sample complete patient rows with replacement within each cancer type.
pr_bootstrap_indices <- function(cancer_type) {
  pr_check_text(cancer_type, "cancer_type")
  strata <- split(seq_along(cancer_type), as.character(cancer_type))
  unlist(lapply(strata, function(idx) {
    idx[sample.int(length(idx), length(idx), replace = TRUE)]
  }), use.names = FALSE)
}

# Each input row is one held-out patient prediction. The same patient can
# appear in different splits/frameworks, but is never pooled across them.
# Confidence intervals condition on each fitted model and its test mixture;
# they do not quantify variation from refitting across train/test splits.
pr_bootstrap_test_f1 <- function(predictions, n_boot, seed) {
  required <- c("framework", "split_id", "patient_id", "cancer_type",
                "truth", "predicted")
  if (!is.data.frame(predictions) || !all(required %in% names(predictions))) {
    stop("predictions must contain: ", paste(required, collapse = ", "))
  }
  for (nm in required[1:4]) pr_check_text(predictions[[nm]], nm)
  if (anyDuplicated(predictions[c("framework", "split_id", "patient_id")])) {
    stop("Duplicated held-out patient within a framework/split.")
  }
  pr_f1(predictions$truth, predictions$predicted)
  pr_check_integers(n_boot, "n_boot", minimum = 2L)
  if (length(n_boot) != 1L) stop("Provide exactly one n_boot value.")
  groups <- unique(predictions[c("framework", "split_id")])
  groups <- groups[order(groups$framework, groups$split_id), , drop = FALSE]
  rownames(groups) <- NULL
  pr_with_seed(seed, {
    summaries <- draws <- vector("list", nrow(groups))
    for (i in seq_len(nrow(groups))) {
      key <- groups[i, , drop = FALSE]
      rows <- predictions$framework == key$framework &
        predictions$split_id == key$split_id
      test <- predictions[rows, , drop = FALSE]
      test <- test[order(test$cancer_type, test$patient_id), , drop = FALSE]
      if (!setequal(as.character(test$truth), c("1", "2"))) {
        stop("Both observed classes are required in every test framework/split.")
      }
      observed <- pr_f1(test$truth, test$predicted)
      boot <- lapply(seq_len(n_boot), function(b) {
        idx <- pr_bootstrap_indices(test$cancer_type)
        cbind(replicate = b, pr_f1(test$truth[idx], test$predicted[idx]))
      })
      boot <- do.call(rbind, boot)
      n_undefined <- sum(is.na(boot$f1))
      # Never silently discard a resample where F1 has no denominator.
      ci <- if (n_undefined > 0L) c(NA_real_, NA_real_) else
        unname(stats::quantile(boot$f1, c(.025, .975), type = 7))
      summaries[[i]] <- cbind(key, n_test = nrow(test), observed,
                              f1_lower_95 = ci[1], f1_upper_95 = ci[2],
                              n_boot = n_boot, n_undefined = n_undefined,
                              bootstrap_seed = seed)
      draws[[i]] <- cbind(key[rep(1L, n_boot), , drop = FALSE], boot)
    }
    list(summary = do.call(rbind, summaries),
         bootstrap = do.call(rbind, draws))
  })
}

# Reconstruct the historical 80/20 partition of notebook 00: the complete
# archived matrix, in stored row order, grouped by cancer type and sampled with
# dplyr::sample_frac(.80) after set.seed(1234). Only eligible patients are kept
# afterwards, so the held-out set is the historical one restricted to the
# Cluster 1/2 target cohort.
pr_original_split <- function(archive_frame, eligible_ids, seed = 1234L,
                              fraction = .80) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required.")
  required <- c("patient_id", "cancer_type")
  if (!is.data.frame(archive_frame) || !all(required %in% names(archive_frame))) {
    stop("archive_frame must contain patient_id and cancer_type in stored order.")
  }
  if (anyDuplicated(archive_frame$patient_id)) stop("Duplicated archived patient_id.")
  if (!all(eligible_ids %in% archive_frame$patient_id)) {
    stop("Eligible patients are missing from the archived matrix.")
  }
  train <- pr_with_seed(seed, {
    sampled <- dplyr::sample_frac(dplyr::group_by(archive_frame, cancer_type), fraction)
    as.character(sampled$patient_id)
  })
  eligible <- archive_frame[archive_frame$patient_id %in% eligible_ids, required]
  eligible <- eligible[order(eligible$cancer_type, eligible$patient_id), ]
  data.frame(split_id = "0", split_seed = as.integer(seed),
             patient_id = as.character(eligible$patient_id),
             cancer_type = as.character(eligible$cancer_type),
             partition = ifelse(eligible$patient_id %in% train, "train", "test"),
             singleton_stratum = FALSE, row.names = NULL)
}

# Threshold-free and prevalence-aware companions to the positive-class F1.
# AUROC uses the Mann-Whitney rank formulation; average precision is the
# step-wise area under the precision-recall curve. MCC is NA when all
# predictions fall in one class. The trivial F1 is the F1 of a
# classifier that labels every patient Cluster 1, i.e. 2p / (1 + p).
pr_classification_metrics <- function(truth, probability, predicted) {
  truth <- as.character(truth); predicted <- as.character(predicted)
  if (length(truth) != length(probability) || anyNA(probability)) {
    stop("probability must be complete and aligned with truth.")
  }
  positive <- truth == "1"
  n_pos <- sum(positive); n_neg <- sum(!positive)
  ranks <- rank(probability)
  auroc <- if (n_pos && n_neg) (sum(ranks[positive]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg) else NA_real_
  # Average precision over distinct score thresholds (tied probabilities form
  # one step), so the value does not depend on the order of tied patients.
  ord <- order(-probability)
  score <- probability[ord]; hits <- positive[ord]
  last <- which(c(diff(score) != 0, TRUE))
  tp_at <- cumsum(hits)[last]; fp_at <- cumsum(!hits)[last]
  average_precision <- if (n_pos) {
    sum(diff(c(0, tp_at / n_pos)) * tp_at / (tp_at + fp_at))
  } else NA_real_
  counts <- pr_f1(truth, predicted)
  tp <- as.numeric(counts$tp); fp <- as.numeric(counts$fp)
  fn <- as.numeric(counts$fn); tn <- as.numeric(counts$tn)
  mcc_den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  prevalence <- n_pos / length(truth)
  data.frame(prevalence_cluster1 = prevalence,
             trivial_f1 = 2 * prevalence / (1 + prevalence),
             auroc = auroc, average_precision = average_precision,
             balanced_accuracy = mean(c(tp / (tp + fn), tn / (tn + fp))),
             mcc = if (mcc_den > 0) (tp * tn - fp * fn) / mcc_den else NA_real_)
}
