# Figure 5 in its publication format: cluster heatmaps, signature activity by
# cluster and cluster composition by tumor type, from the saved cluster fits
# (mc_*_active_4 in data/processed/clustering_data.RData), the partitions used
# by every analysis (survival, prediction, FGA).
#
# Plotting code as in code/02_clustering.Rmd (pheatmap heatmaps, violins with
# mean and 95% CI, cluster composition by tumor type ordered by the Cluster 1
# share), assembled with ggarrange on an A4 landscape page (11.69 x 8.27 in)
# with lowercase panel letters, "Clusters" above the first row of violins and
# x-axis titles on the bottom row only. The page is written at 4080 px across
# and labelled 600 dpi. The ggplot panels use the Verdana font; where it is
# not installed, the graphics devices substitute another font.
#
# Run from the repository root:
#   Rscript code/scripts/figure5.R
# Outputs: results/figures/Figure5.{tiff,png,pdf} and
#          results/clustering/cluster_summary.csv (cluster sizes and mean
#          signature activities)

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(pheatmap); library(ggpubr)
})
load("data/processed/clustering_data.RData")
out <- "results/figures"; tables_out <- "results/clustering"
for (d in c(out, tables_out)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
# The ggplot panels and panel letters are in Verdana, the heatmaps in
# Helvetica (pheatmap default).
font <- "Verdana"
cluster_colours <- c("1" = "#72F59A", "2" = "#F57D72", "3" = "#F7C073", "4" = "#8672F5")

heatmap_panel <- function(active, cls) {
  ord <- order(cls)
  tmp <- active[ord, -1]
  ann <- data.frame(Cluster = as.factor(cls[ord]), row.names = rownames(tmp))
  p <- pheatmap(tmp, cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = ann,
           annotation_colors = list(Cluster = cluster_colours[levels(ann$Cluster)]),
           show_rownames = FALSE, angle_col = 90, cellwidth = 15, silent = TRUE)$gtable
  # Room on the top and left for the panel letter.
  gtable::gtable_add_padding(p, unit(c(16, 0, 0, 30), "pt"))
}

violin_panel <- function(active, cls, top, bottom) {
  sigs <- colnames(active)[-1]
  d <- data.frame(active[, -1], Group = cls, check.names = FALSE)
  d <- pivot_longer(d, all_of(sigs), names_to = "Variable", values_to = "Value")
  # One colour ramp per cluster, from white to the cluster colour.
  ramp <- unlist(lapply(sort(unique(cls)), function(k) {
    v <- colorRampPalette(c("white", cluster_colours[[as.character(k)]]))(length(sigs) + 1)[-1]
    setNames(v, paste(sigs, k, sep = "."))
  }))
  d$fill_var <- paste(d$Variable, d$Group, sep = ".")
  # A signature that is zero in every tumor of a cluster (Steele CN17 in
  # Clusters 2 and 4) has no density to draw: with trim = FALSE the default
  # bandwidth would spread it from -0.5 to 0.5. Only its mean is shown.
  varying <- d %>% group_by(Variable, Group) %>% filter(sd(Value) > 0) %>% ungroup()
  # No axis limits: the violins (trim = FALSE) and the means use all values.
  ggplot(d, aes(Variable, Value, fill = fill_var)) +
    geom_violin(data = varying, trim = FALSE, scale = "width", colour = "black", alpha = .8) +
    scale_fill_manual(values = ramp) +
    stat_summary(fun = mean, geom = "point", shape = 21, fill = "black", size = 2.5) +
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = .2, colour = "black") +
    facet_wrap(~Group, nrow = 1) +
    labs(x = if (bottom) "Signature" else NULL, title = if (top) "Clusters" else NULL) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = .5, size = 14), text = element_text(family = font))
}

composition_panel <- function(active, cls, n_clusters, bottom) {
  agg <- data.frame(cancer_type = as.character(active$cancer_type), cluster = cls) %>%
    count(cancer_type, cluster) %>%
    group_by(cancer_type) %>%
    mutate(relative_freq = n / sum(n)) %>%
    ungroup()
  first <- agg %>% filter(cluster == 1) %>% arrange(desc(relative_freq)) %>% pull(cancer_type)
  agg$cancer_type <- factor(agg$cancer_type, c(first, setdiff(unique(agg$cancer_type), first)))
  ggplot(agg, aes(cancer_type, relative_freq)) +
    geom_bar(aes(fill = as.factor(cluster)), stat = "identity", position = "fill") +
    labs(y = "Relative Frequency", x = if (bottom) "Tumor Type" else NULL, fill = "Cluster") +
    scale_fill_manual(values = cluster_colours[1:n_clusters]) +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.title = element_text(size = 14), plot.margin = margin(16, 5.5, 5.5, 24),
          text = element_text(family = font))
}

compendia <- list(Drews = list(active_drew, mc_drew_active_4),
                  Steele = list(active_steele, mc_steele_active_4),
                  Tao = list(active_tao, mc_tao_active_4))
plots <- list(); summary_rows <- list()
for (i in seq_along(compendia)) {
  fw <- names(compendia)[i]
  a <- compendia[[i]][[1]]; cls <- compendia[[i]][[2]]$classification
  stopifnot(length(cls) == nrow(a))
  plots <- c(plots, list(heatmap_panel(a, cls),
                         violin_panel(a, cls, top = i == 1, bottom = i == 3),
                         composition_panel(a, cls, length(unique(cls)), bottom = i == 3)))
  means <- aggregate(a[, -1], list(cluster = cls), mean)
  tab <- table(cls)[as.character(means$cluster)]
  summary_rows[[fw]] <- pivot_longer(
    data.frame(compendium = fw, cluster = means$cluster, n = as.vector(tab),
               share = as.vector(tab) / length(cls), round(means[, -1], 3), check.names = FALSE),
    -(1:4), names_to = "signature", values_to = "mean_activity")
}
summary_table <- as.data.frame(do.call(rbind, summary_rows))
utils::write.csv(summary_table, file.path(tables_out, "cluster_summary.csv"), row.names = FALSE)

# ggarrange (through cowplot) measures the text on a null device when it builds
# the panels: use ragg, which knows Verdana, instead of R's PDF device.
cowplot::set_null_device("agg")
page <- c(width = 841 / 72, height = 595 / 72)
figure <- ggarrange(plotlist = plots, labels = letters[1:9], nrow = 3, ncol = 3,
                    widths = c(.68, 1.37, 1.15), heights = c(.9, 1, 1.1),
                    font.label = list(size = 16, face = "plain", family = font)) +
  theme(plot.margin = margin(18, 14, 10, 14))
# Quartz (macOS) for the PDF, so that Verdana is embedded (cairo elsewhere);
# ragg for the bitmaps.
pdf_device <- if (capabilities("aqua")) function(filename, width, height, ...)
  grDevices::quartz(type = "pdf", file = filename, width = width, height = height) else
  function(filename, width, height, ...) grDevices::cairo_pdf(filename, width = width, height = height)
ggsave(file.path(out, "Figure5.pdf"), figure, width = page[["width"]], height = page[["height"]],
       bg = "white", device = pdf_device)
# 4080 px across the page; labelled 600 dpi (6.8 x 4.81 in).
png_file <- file.path(out, "Figure5.png")
ggsave(png_file, figure, width = page[["width"]], height = page[["height"]],
       dpi = 4080 / page[["width"]], bg = "white", device = ragg::agg_png)
img <- magick::image_read(png_file)
magick::image_write(img, png_file, density = "600x600")
magick::image_write(img, file.path(out, "Figure5.tiff"), format = "tiff",
                    compression = "LZW", density = "600x600")
print(magick::image_info(img))
