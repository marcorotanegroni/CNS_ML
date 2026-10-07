# Adapter for the archived final matrices in Zenodo record 20617051.
# Source defines functions only; load large files explicitly from the notebook.
# training_fold additionally uses the normalized pre-filter archive.

pi_archive_spec <- function() {
  list(
    final_md5 = "cd14cf9fd5eef77e3ea35d4022283eb5",
    normalized_md5 = "35574218fed4aad02c1557461be95c91",
    n_expression = 15976L, n_methylation = 13674L,
    final_expression = c(Drews = 12653L, Steele = 12611L, Tao = 12599L),
    final_methylation = c(Drews = 9227L, Steele = 9073L, Tao = 8942L),
    objects = c(Drews = "mat_purity_mut_drews",
                Steele = "mat_purity_mut_steele", Tao = "mat_purity_mut_tao"),
    cluster_suffix = c(Drews = "drew", Steele = "steele", Tao = "tao"),
    targets = c(Drews = "CX1", Steele = "CN1", Tao = "Sig1"),
    dimensions = list(Drews = c(5024L, 35497L),
                      Steele = c(7893L, 35301L), Tao = c(7823L, 35158L)),
    mutation_start = c(Drews = 21886L, Steele = 21690L, Tao = 21547L)
  )
}

pi_check_ids <- function(ids, label) {
  if (is.null(ids) || anyNA(ids) || anyDuplicated(ids) ||
      any(!grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$", ids))) {
    stop(label, ": expected unique TCGA patient identifiers (12 characters).")
  }
}

pi_load_inputs <- function(final_path, normalized_path = NULL, cluster_path,
                           output_dir = NULL,
                           input_mode = c("archived_final", "training_fold"),
                           frameworks = c("Drews", "Steele", "Tao")) {
  input_mode <- match.arg(input_mode)
  paths <- c(final_matrices = final_path, clusters = cluster_path)
  if (input_mode == "training_fold") {
    if (is.null(normalized_path) || length(normalized_path) != 1L) {
      stop("training_fold requires normalized_path.")
    }
    paths <- c(paths, normalized = normalized_path)
  }
  if (any(!file.exists(paths))) {
    stop("Missing input: ", paste(paths[!file.exists(paths)], collapse = ", "))
  }
  spec <- pi_archive_spec()
  message("Checking input checksums...")
  hashes <- vapply(paths, function(p) unname(tools::md5sum(p)), "")
  if (hashes[["final_matrices"]] != spec$final_md5 ||
      (input_mode == "training_fold" && hashes[["normalized"]] != spec$normalized_md5)) {
    stop("The matrices do not match Zenodo record 20617051. Check their provenance ",
         "and layout before using this archive-specific adapter.")
  }
  cenv <- new.env(parent = emptyenv())
  load(cluster_path, envir = cenv)
  fenv <- new.env(parent = emptyenv())
  message(if (input_mode == "archived_final")
    "Loading final matrices; preserving all archived predictor values..." else
    "Loading final matrices; retaining clinical/mutation predictors only...")
  load(final_path, envir = fenv)
  if (!all(spec$objects %in% ls(fenv))) stop("Missing final matrix objects.")
  if (!length(frameworks) || any(!frameworks %in% names(spec$objects))) stop("Invalid frameworks.")
  # Free the matrices of frameworks not requested before building predictors:
  # each prepared matrix is several GB.
  rm(list = unname(spec$objects[!names(spec$objects) %in% frameworks]), envir = fenv)
  invisible(gc())
  extras <- audits <- vector("list", length(frameworks))
  names(extras) <- names(audits) <- frameworks
  for (framework in names(extras)) {
    object <- spec$objects[[framework]]
    final <- fenv[[object]]
    if (!identical(dim(final), spec$dimensions[[framework]])) {
      stop(framework, ": unexpected final matrix dimensions.")
    }
    ids <- rownames(final)
    pi_check_ids(ids, framework)
    id_col <- if (framework == "Drews") "patient" else "Tumor_Sample_Barcode"
    target <- spec$targets[[framework]]
    if (!setequal(names(final)[1:5], c(id_col, "cancer_type", target, "age", "purity")) ||
        !identical(as.character(final[[id_col]]), ids) || anyDuplicated(names(final))) {
      stop(framework, ": metadata layout or patient identity mismatch.")
    }
    first <- unname(spec$mutation_start[[framework]])
    n_exp <- unname(spec$final_expression[[framework]])
    n_meth <- unname(spec$final_methylation[[framework]])
    if (first != 6L + n_exp + n_meth) {
      stop(framework, ": inconsistent archived omics block boundaries.")
    }
    if (names(final)[first] != "cell_type" || ncol(final) - first + 1L != 13612L ||
        !identical(names(final)[first + c(1L, 3L, 4L)], c("HIPPO", "NOTCH", "NRF2"))) {
      stop(framework, ": mutation/pathway block does not match the archive.")
    }
    suffix <- spec$cluster_suffix[[framework]]
    labels <- cenv[[paste0("mc_", suffix, "_active_4")]]$classification
    active <- cenv[[paste0("active_", suffix)]]
    pi_check_ids(names(labels), paste(framework, "cluster labels"))
    if (!all(ids %in% names(labels)) || !all(ids %in% rownames(active))) {
      stop(framework, ": final patients not found in archived clustering.")
    }
    cancer <- as.character(final$cancer_type)
    if (anyNA(cancer) || any(cancer != as.character(active[ids, "cancer_type"]))) {
      stop(framework, ": cancer-type annotations disagree with the clustering.")
    }
    keep <- as.character(labels[ids]) %in% c("1", "2") & cancer != "LAML"
    cols <- if (input_mode == "archived_final") {
      names(final)[!names(final) %in% c(id_col, "cancer_type", target)]
    } else c("age", "purity", names(final)[seq.int(first, ncol(final))])
    if (!all(vapply(final[cols], is.numeric, logical(1)))) {
      stop(framework, ": nonnumeric archived predictor.")
    }
    other <- as.matrix(final[keep, cols, drop = FALSE])
    if (any(!is.finite(other))) stop(framework, ": incomplete archived predictors.")
    source_columns <- match(cols, names(final))
    types <- rep("other", length(cols))
    if (input_mode == "archived_final") {
      types[source_columns %in% seq.int(6L, 5L + n_exp)] <- "expression"
      types[source_columns %in% seq.int(6L + n_exp, first - 1L)] <- "methylation"
    }
    feature_map <- data.frame(
      feature = paste0(if (input_mode == "archived_final") "archived__" else "other__",
                       seq_along(cols)), source_name = cols,
      feature_type = types, source = "final_matrices.RData",
      source_column = source_columns, stringsAsFactors = FALSE
    )
    colnames(other) <- feature_map$feature
    patients <- data.frame(patient_id = ids[keep], cancer_type = cancer[keep],
                           truth = as.character(labels[ids[keep]]))
    # Complete archived row order and cancer types: needed to reconstruct the
    # historical seed-1234 partition, which was drawn on the whole matrix.
    archive_frame <- data.frame(patient_id = ids, cancer_type = cancer,
                                stringsAsFactors = FALSE)
    extras[[framework]] <- list(x = other, patients = patients, feature_map = feature_map,
                                archive_frame = archive_frame)
    audits[[framework]] <- data.frame(
      framework = framework, input_mode = input_mode, archived_n = nrow(final),
      selected_cluster12_n = sum(keep), excluded_other_clusters_or_LAML = sum(!keep),
      cluster1_n = sum(patients$truth == "1"), cluster2_n = sum(patients$truth == "2"),
      cancer_types = length(unique(patients$cancer_type)),
      clinical_mutation_features = sum(types == "other"),
      loaded_archived_features = ncol(other),
      missing_normalized_n = NA_integer_
    )
    rm(list = object, envir = fenv)
    rm(final, other)
    invisible(gc())
  }
  rm(fenv, cenv)
  invisible(gc())
  normalized <- normalized_map <- NULL
  if (input_mode == "training_fold") {
    message("Loading normalized expression and methylation before feature selection...")
    nenv <- new.env(parent = emptyenv())
    load(normalized_path, envir = nenv)
    if (!exists("ordinata", nenv, inherits = FALSE)) stop("Missing ordinata object.")
    normalized <- nenv$ordinata
    if (!identical(dim(normalized), c(8988L, 29650L)) ||
        !all(vapply(normalized, is.numeric, logical(1)))) {
      stop("Unexpected normalized matrix layout.")
    }
    pi_check_ids(rownames(normalized), "ordinata")
    # The block boundary is documented in notebook 00 and verified against this
    # immutable archive. Duplicate gene names across modalities are kept distinct.
    if (!identical(names(normalized)[c(1L, 15976L, 15977L)],
                   c("A1BG", "psiTPTE22", "A1CF"))) {
      stop("Unexpected expression/methylation boundary.")
    }
    types <- rep(c("expression", "methylation"),
                 c(spec$n_expression, spec$n_methylation))
    normalized_map <- data.frame(
      feature = paste0("omics__", seq_len(ncol(normalized))),
      source_name = names(normalized), feature_type = types,
      source = "exp_meth_post_normalization.RData",
      source_column = seq_len(ncol(normalized)), stringsAsFactors = FALSE
    )
    names(normalized) <- normalized_map$feature
    for (framework in names(extras)) {
      missing <- setdiff(extras[[framework]]$patients$patient_id, rownames(normalized))
      audits[[framework]]$missing_normalized_n <- length(missing)
      if (length(missing)) {
        stop(framework, ": ", length(missing),
             " eligible patients are absent from ordinata. Cohort must be reconciled ",
             "explicitly; this adapter will not silently drop them.")
      }
    }
  }
  audit <- do.call(rbind, audits)
  rownames(audit) <- NULL
  provenance <- data.frame(input = names(paths),
                           file = basename(paths), md5 = hashes, row.names = NULL)
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(audit, file.path(output_dir, "cohort_audit.csv"), row.names = FALSE)
    write.csv(provenance, file.path(output_dir, "input_provenance.csv"), row.names = FALSE)
  }
  list(input_mode = input_mode, normalized = normalized, normalized_map = normalized_map,
       extras = extras, audit = audit, provenance = provenance)
}

pi_prepare_framework <- function(inputs, framework) {
  if (!framework %in% names(inputs$extras)) stop("Unknown framework.")
  extra <- inputs$extras[[framework]]
  patients <- extra$patients[order(extra$patients$patient_id), , drop = FALSE]
  ids <- patients$patient_id
  if (identical(inputs$input_mode, "archived_final")) {
    x <- extra$x[ids, , drop = FALSE]
    feature_map <- extra$feature_map
  } else {
    x <- cbind(as.matrix(inputs$normalized[ids, , drop = FALSE]),
               extra$x[ids, , drop = FALSE])
    feature_map <- rbind(inputs$normalized_map, extra$feature_map)
  }
  rownames(patients) <- NULL
  list(x = x, patients = patients,
       feature_type = setNames(feature_map$feature_type, feature_map$feature),
       feature_map = feature_map, provenance = inputs$provenance,
       input_mode = inputs$input_mode, archive_frame = extra$archive_frame)
}

# Remove predictors by their source name (e.g. "purity" for R2).
pi_drop_features <- function(prepared, source_names) {
  drop <- prepared$feature_map$feature[prepared$feature_map$source_name %in% source_names]
  if (length(drop) != length(source_names)) {
    stop("Expected exactly one predictor for each of: ", paste(source_names, collapse = ", "))
  }
  keep <- !colnames(prepared$x) %in% drop
  prepared$x <- prepared$x[, keep, drop = FALSE]
  prepared$feature_type <- prepared$feature_type[colnames(prepared$x)]
  prepared$feature_map <- prepared$feature_map[!prepared$feature_map$feature %in% drop, ]
  prepared
}

# Restrict the predictors to the given feature types (e.g. "expression" for
# the expression-only models of the manuscript).
pi_keep_feature_types <- function(prepared, types) {
  keep <- prepared$feature_type[colnames(prepared$x)] %in% types
  if (!any(keep)) stop("No predictors of type: ", paste(types, collapse = ", "))
  prepared$x <- prepared$x[, keep, drop = FALSE]
  prepared$feature_type <- prepared$feature_type[colnames(prepared$x)]
  prepared$feature_map <- prepared$feature_map[prepared$feature_map$feature %in% colnames(prepared$x), ]
  prepared
}

# Scalar copy-number baseline (R2): ploidy, fraction of the genome with loss of
# heterozygosity (ASCAT frac_homo) and number of copy number alterations, from
# the ASCAT penalty-70 TCGA fits distributed with Drews et al. One primary
# tumour (sample type 01), representative ("rep") profile per patient. Patients
# without a profile keep missing values, which XGBoost handles natively, so the
# cohort and splits are identical to the multi-omic models.
pi_scalar_baseline <- function(prepared, ascat_path) {
  if (!file.exists(ascat_path)) stop("Missing ASCAT metadata: ", ascat_path)
  ascat <- as.data.frame(readRDS(ascat_path))
  needed <- c("patient", "barcodeTumour", "rep", "ploidy", "frac_homo", "CNAs")
  if (!all(needed %in% names(ascat))) stop("Unexpected ASCAT metadata layout.")
  ascat <- ascat[substr(ascat$barcodeTumour, 14, 15) == "01" & ascat$rep %in% TRUE, ]
  if (anyDuplicated(ascat$patient)) stop("More than one representative primary profile per patient.")
  ids <- prepared$patients$patient_id
  row <- match(ids, ascat$patient)
  scalars <- c("ploidy", "frac_homo", "CNAs")
  x <- vapply(scalars, function(v) as.numeric(ascat[[v]][row]), numeric(length(ids)))
  x <- matrix(x, nrow = length(ids), dimnames = list(ids, paste0("scalar__", scalars)))
  feature_map <- data.frame(feature = colnames(x), source_name = scalars,
                            feature_type = "other", source = basename(ascat_path),
                            source_column = match(scalars, names(ascat)), stringsAsFactors = FALSE)
  coverage <- data.frame(n = length(ids), with_profile = sum(!is.na(row)),
                         t(colSums(!is.na(x))))
  list(x = x, patients = prepared$patients,
       feature_type = setNames(feature_map$feature_type, feature_map$feature),
       feature_map = feature_map, provenance = prepared$provenance,
       input_mode = "scalar_baseline", archive_frame = prepared$archive_frame,
       coverage = coverage)
}
