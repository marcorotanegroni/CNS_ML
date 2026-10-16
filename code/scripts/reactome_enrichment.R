# Reactome over-representation of the high-importance predictors of the
# split-0 Cluster 1 versus Cluster 2 classifiers (Figure 8, Methods).
#
# STRING functional enrichment (category Reactome, whole-genome background,
# FDR < 0.05) through the STRING API with a fixed STRING version (12.0), on the
# gene symbols of the expression and methylation predictors accounting for up
# to 80% of the cumulative gain (suffixes ".1" of duplicated names removed;
# multi-gene methylation features split into their genes).
#
# Run from the repository root after code/scripts/figure8.R:
#   Rscript code/scripts/reactome_enrichment.R
# Inputs:  results/cluster_prediction/feature_importance/importance_split0_<framework>.csv
# Outputs: results/cluster_prediction/reactome/
# Only gene symbols are sent to STRING.

suppressMessages({library(httr); library(jsonlite)})
string_url <- "https://version-12-0.string-db.org"
frameworks <- c("Drews", "Steele", "Tao")
out <- "results/cluster_prediction/reactome"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

gene_list <- function(importance) {
  x <- importance[importance$in_80pct & importance$data_type %in% c("exp", "meth"), ]
  genes <- trimws(unlist(strsplit(x$variable, ";", fixed = TRUE)))
  unique(sub("\\.[0-9]+$", "", genes[genes != ""]))
}

string_post <- function(path, genes, ...) {
  res <- POST(paste0(string_url, path), encode = "form",
              body = list(identifiers = paste(genes, collapse = "\r"), species = 9606,
                          caller_identity = "cns_ml", ...))
  stop_for_status(res)
  res
}

summary <- do.call(rbind, lapply(frameworks, function(fw) {
  imp <- utils::read.csv(file.path("results/cluster_prediction/feature_importance",
                                   paste0("importance_split0_", fw, ".csv")))
  genes <- gene_list(imp)
  writeLines(genes, file.path(out, paste0("genes_", fw, ".txt")))
  enrichment <- fromJSON(content(string_post("/api/json/enrichment", genes), "text",
                                 encoding = "UTF-8"))
  reactome <- if (length(enrichment) && "category" %in% names(enrichment)) {
    enrichment[enrichment$category == "RCTM", , drop = FALSE]
  } else data.frame()
  if (nrow(reactome)) {
    reactome$inputGenes <- vapply(reactome$inputGenes, paste, "", collapse = ";")
    reactome$preferredNames <- vapply(reactome$preferredNames, paste, "", collapse = ";")
    utils::write.csv(reactome, file.path(out, paste0("reactome_", fw, ".csv")), row.names = FALSE)
    # STRING's own enrichment figure.
    fig <- string_post("/api/highres_image/enrichmentfigure", genes, category = "RCTM",
                       group_by_similarity = 0.8, color_palette = "mint_blue",
                       number_of_term_shown = 10, x_axis = "signal")
    writeBin(content(fig, "raw"), file.path(out, paste0("reactome_", fw, ".png")))
  }
  data.frame(framework = fw, n_genes = length(genes),
             n_reactome_fdr_005 = nrow(reactome),
             top_terms = if (nrow(reactome)) paste(head(reactome$description, 5), collapse = "; ") else "")
}))
version <- fromJSON(content(GET(paste0(string_url, "/api/json/version")), "text", encoding = "UTF-8"))
summary$string_version <- version$string_version[1]
summary$queried <- format(Sys.Date())
utils::write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
print(summary)
