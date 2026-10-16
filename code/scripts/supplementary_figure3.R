# Supplementary Figure 3: overlap of the signature-based clusters across
# compendia, from the saved cluster assignments (mc_*_active_4 in
# data/processed/clustering_data.RData). Each compendium is restricted to the
# tumours of its multi-omic predictor matrix (final_matrices.RData), and the
# UpSet panels show intersections of at least 50 tumours: (a) all clusters,
# (b) Clusters 1 and 2, (c) Venn diagram of Cluster 2.
#
# Run from the repository root:
#   Rscript code/scripts/supplementary_figure3.R
# Requires ComplexUpset and ggVennDiagram (CRAN) and final_matrices.RData
# (Zenodo; data/processed/zenodo/ or the repository root).
# Outputs: results/figures/Supplementary_Figure3.{png,pdf} and
#          results/clustering/supplementary_figure3_counts.csv

suppressMessages({library(ggplot2); library(ComplexUpset); library(ggVennDiagram); library(ggpubr)})
out <- "results/figures"; tables_out <- "results/clustering"
for (d in c(out, tables_out)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
final_path <- if (file.exists("data/processed/zenodo/final_matrices.RData")) {
  "data/processed/zenodo/final_matrices.RData"
} else "final_matrices.RData"

clusters <- new.env(); load("data/processed/clustering_data.RData", envir = clusters)
matrices <- new.env(); load(final_path, envir = matrices)
spec <- list(Drews = c("mc_drew_active_4", "mat_purity_mut_drews"),
             Steele = c("mc_steele_active_4", "mat_purity_mut_steele"),
             Tao = c("mc_tao_active_4", "mat_purity_mut_tao"))
labels <- lapply(spec, function(s) {
  cls <- clusters[[s[1]]]$classification
  cls[intersect(names(cls), rownames(matrices[[s[2]]]))]
})
rm(matrices); invisible(gc())

cluster_sets <- function(keep) {
  out <- list()
  for (fw in names(labels)) for (k in keep[[fw]]) {
    out[[paste0(fw, "_", k)]] <- names(labels[[fw]])[labels[[fw]] == k]
  }
  out
}
membership <- function(sets) {
  ids <- unique(unlist(sets))
  data.frame(ID = ids, sapply(sets, function(s) ids %in% s), check.names = FALSE)
}
upset_panel <- function(sets) {
  ComplexUpset::upset(
    membership(sets), intersect = names(sets), min_size = 50, width_ratio = 0.15,
    set_sizes = upset_set_size() + theme(axis.text.x = element_text(angle = 90)),
    base_annotations = list(`Intersection size` = intersection_size(
      text = list(colour = "black", angle = 90, hjust = -0.1, vjust = 0.5)) +
      scale_fill_manual(values = "black")))
}

all_sets <- cluster_sets(list(Drews = 1:4, Steele = 1:4, Tao = 1:3))
c12_sets <- cluster_sets(list(Drews = 1:2, Steele = 1:2, Tao = 1:2))
c2_sets <- cluster_sets(list(Drews = 2, Steele = 2, Tao = 2))

venn <- ggVennDiagram(c2_sets, label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "orange") +
  theme_void() + coord_cartesian(clip = "off")

figure <- ggarrange(plotlist = list(upset_panel(all_sets), upset_panel(c12_sets), venn),
                    labels = "auto", heights = c(1.75, 1), nrow = 2, ncol = 2, vjust = 3)
for (ext in c("png", "pdf")) {
  args <- list(filename = file.path(out, paste0("Supplementary_Figure3.", ext)), plot = figure,
               width = 13, height = 10, units = "in", bg = "white")
  if (ext == "png") args$dpi <- 300
  do.call(ggsave, args)
}

# Counts shown in the panels.
distinct_counts <- function(sets) {
  m <- membership(sets)
  key <- apply(m[, -1, drop = FALSE], 1, function(r) paste(names(sets)[r], collapse = " & "))
  t <- sort(table(key), decreasing = TRUE)
  data.frame(intersection = names(t), n = as.integer(t))
}
counts <- rbind(
  cbind(panel = "a", subset(distinct_counts(all_sets), n >= 50)),
  cbind(panel = "b", subset(distinct_counts(c12_sets), n >= 50)),
  cbind(panel = "c", distinct_counts(c2_sets)))
utils::write.csv(counts, file.path(tables_out, "supplementary_figure3_counts.csv"), row.names = FALSE)
print(table(counts$panel))
