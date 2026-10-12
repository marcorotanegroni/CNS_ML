# Supplementary Figure 9: fraction genome altered (FGA) as a scalar copy
# number baseline (Reviewer 2).
#   a  FGA by cluster in each compendium, with the share of FGA variance
#      explained by the clusters (R2) and by the clusters beyond cancer type.
#   b  Cluster 2 vs Cluster 1 hazard ratios of the Figure 7 Cox models, on the
#      patients with an FGA value, without FGA and with FGA (penalized spline),
#      with their 95% CIs. The Holm-adjusted p-values (as in Figure 7) are in
#      results/revision/fga_baseline/survival_cox.csv and in the text.
# Panel c (classifiers) is added once the classification without purity
# (Reviewer 2, comment 1) is available.
#
# Run from the repository root after 15_fga_survival.R and 16_fga_clustering.R:
#   Rscript code/revision/17_supplementary_figure9.R
# Output: results/revision/supplementary_figure9/

suppressMessages({library(ggplot2); library(patchwork)})
frameworks <- c("Drews", "Steele", "Tao")
colours <- c(Drews = "#1f78b4", Steele = "#33a02c", Tao = "#e31a1c")
fig_theme <- theme_minimal(base_size = 7) +
  theme(axis.text = element_text(size = 6, colour = "black"),
        axis.title = element_text(size = 7), legend.text = element_text(size = 6),
        legend.title = element_text(size = 6.5), strip.text = element_text(size = 6.5, face = "bold"),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 9, face = "bold"))
out <- "results/revision/supplementary_figure9"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fga <- utils::read.csv("data/processed/fga_pancanatlas.csv")
fga <- setNames(fga$fga, fga$patient)

# Panel a: FGA by cluster.
cl <- new.env(); load("data/processed/clustering_data.RData", envir = cl)
fits <- c(Drews = "mc_drew_active_4", Steele = "mc_steele_active_4", Tao = "mc_tao_active_4")
long <- do.call(rbind, lapply(frameworks, function(fw) {
  cls <- cl[[fits[fw]]]$classification
  data.frame(framework = fw, cluster = paste("Cluster", cls), fga = unname(fga[names(cls)]))
}))
long <- long[!is.na(long$fga), ]
long$framework <- factor(long$framework, frameworks)
r2 <- utils::read.csv("results/revision/fga_baseline/clustering_fga_variance.csv")
r2$label1 <- sprintf("italic(R)^2 == '%.2f'", r2$r2_cluster)
r2$label2 <- sprintf("beyond cancer type: %.2f", r2$partial_r2_cluster_given_cancer_type)
r2$framework <- factor(r2$framework, frameworks)
panel_a <- ggplot(long, aes(cluster, fga, fill = framework)) +
  geom_boxplot(outlier.size = .2, outlier.alpha = .4, linewidth = .3, alpha = .7) +
  geom_text(data = r2, aes(x = Inf, y = 1.15, label = label1), inherit.aes = FALSE,
            parse = TRUE, hjust = 1.02, vjust = 1, size = 2) +
  geom_text(data = r2, aes(x = Inf, y = 1.06, label = label2), inherit.aes = FALSE,
            hjust = 1.02, vjust = 1, size = 2) +
  facet_wrap(~framework, scales = "free_x") +
  scale_fill_manual(values = colours, guide = "none") +
  scale_y_continuous(limits = c(0, 1.16), breaks = seq(0, 1, .25)) +
  labs(x = NULL, y = "Fraction genome altered") +
  fig_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Panel b: Cluster 2 vs Cluster 1 without and with FGA, same patients.
cox <- utils::read.csv("results/revision/fga_baseline/survival_cox.csv")
models <- c(figure7_fga_n = "Without FGA", with_fga_spline = "With FGA")
f <- cox[cox$term == "cluster2" & cox$model %in% names(models), ]
f$model <- factor(models[f$model], models)
f$framework <- factor(f$framework, frameworks)
f$row <- factor(paste(f$framework, f$endpoint), rev(paste(rep(frameworks, each = 2), c("OS", "PFI"))))
f$hr_text <- sprintf("%.2f (%.2f-%.2f)", f$hr, f$lower, f$upper)
dodge <- c(`Without FGA` = .17, `With FGA` = -.17)
f$y <- as.numeric(f$row) + dodge[as.character(f$model)]
forest <- ggplot(f, aes(hr, y)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50", linewidth = .3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, colour = framework), height = .15, linewidth = .4) +
  geom_point(aes(shape = model, colour = framework), fill = "white", size = 1.8, stroke = .5) +
  scale_shape_manual(values = c(`Without FGA` = 16, `With FGA` = 21), name = NULL) +
  scale_colour_manual(values = colours, guide = "none") +
  scale_y_continuous(breaks = seq_along(levels(f$row)), labels = levels(f$row),
                     expand = expansion(add = .5)) +
  scale_x_log10(breaks = c(.5, .6, .7, .8, 1, 1.2)) +
  labs(x = "Hazard ratio, Cluster 2 vs Cluster 1 (95% CI)", y = NULL) +
  fig_theme + theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0),
                    legend.box.spacing = unit(0, "pt"))
table_b <- ggplot(f, aes(y = y)) +
  geom_text(aes(x = 0, label = hr_text), hjust = 0, size = 2) +
  annotate("text", x = 0, y = length(levels(f$row)) + .65, label = "HR (95% CI)",
           hjust = 0, size = 2, fontface = "bold") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = range(f$y) + c(-.5, .9), expand = c(0, 0)) +
  theme_void()
forest <- forest + coord_cartesian(xlim = c(.5, 1.2), ylim = range(f$y) + c(-.5, .9), expand = FALSE)

figure <- panel_a + labs(tag = "a") + forest + labs(tag = "b") + table_b +
  plot_layout(design = "AA\nBC", heights = c(1, 1), widths = c(.78, .22))
size <- c(width = 6.8, height = 5.2)
ggsave(file.path(out, "Supplementary_Figure9.png"), figure, width = size[["width"]],
       height = size[["height"]], dpi = 600, bg = "white")
ggsave(file.path(out, "Supplementary_Figure9.pdf"), figure, width = size[["width"]],
       height = size[["height"]])
ggsave(file.path(out, "Supplementary_Figure9.tiff"), figure, width = size[["width"]],
       height = size[["height"]], dpi = 600, compression = "lzw", bg = "white")
message("Supplementary Figure 9 written to ", out)
