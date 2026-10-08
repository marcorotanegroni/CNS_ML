# Figure 8 from the revised cluster classifiers. Panels a-e: the models fitted
# on split 0 (cancer-type-stratified 80/20 partition, seed 1234) with the
# saved cluster assignments. Panel f: Cluster 1 test-set F1 on split 0 with
# its 95% bootstrap interval and on the 20 repeated splits. The submitted
# panel f (Reactome enrichment of the Drews predictors) is no longer shown: no
# Reactome term is significant for the revised classifiers
# (code/revision/10_reactome_enrichment.R).
#
# Run from the repository root after the prediction workflow:
#   Rscript code/revision/08_figure8.R
# Inputs:  results/revision/prediction/fixed/archived/<framework>/predictions.csv
#          results/revision/prediction/fixed/all_split_metrics.csv
#          results/revision/prediction/fixed/archived/<framework>/checkpoints/split_000_seed_1234.rds
#          final_matrices.RData (feature names; data/processed/zenodo/ or root)
# Outputs: results/revision/prediction/figure8/

suppressMessages({library(ggplot2); library(patchwork); library(xgboost)})
for (f in c("02_prediction_resampling.R", "03_prediction_training.R", "04_prediction_inputs.R")) {
  source(file.path("code/revision", f))
}
frameworks <- c("Drews", "Steele", "Tao")
colours <- c(Drews = "#1f78b4", Steele = "#33a02c", Tao = "#e31a1c")
fixed <- "results/revision/prediction/fixed/archived"
out <- "results/revision/prediction/figure8"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
final_path <- if (file.exists("data/processed/zenodo/final_matrices.RData")) {
  "data/processed/zenodo/final_matrices.RData"
} else "final_matrices.RData"

# Split-0 held-out predictions.
predictions <- do.call(rbind, lapply(frameworks, function(fw) {
  p <- read.csv(file.path(fixed, fw, "predictions.csv"), colClasses = c(split_id = "character"))
  p[p$split_id == "0", ]
}))

# Panels a-c: overall and cancer-type-specific test accuracy (types with at
# least four test tumours), as in the original figure.
accuracy_panel <- function(fw) {
  p <- predictions[predictions$framework == fw, ]
  correct <- p$truth == p$predicted
  by_type <- aggregate(correct, list(cancer_type = p$cancer_type), mean)
  n <- table(p$cancer_type)
  by_type <- by_type[n[by_type$cancer_type] >= 4, ]
  by_type <- by_type[order(-by_type$x), ]
  d <- rbind(data.frame(cancer_type = "Overall", x = mean(correct), overall = TRUE),
             data.frame(cancer_type = by_type$cancer_type, x = by_type$x, overall = FALSE))
  d$cancer_type <- factor(d$cancer_type, levels = d$cancer_type)
  ggplot(d, aes(cancer_type, x, fill = overall)) +
    geom_col(width = .85) +
    scale_fill_manual(values = c(`TRUE` = "#1f78b4", `FALSE` = "#8fb8d9"), guide = "none") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25), expand = c(0, 0)) +
    labs(x = NULL, y = "Accuracy") +
    theme_minimal(base_size = 5.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
          axis.text.y = element_text(colour = "black"))
}

# Panel d: test-set sample counts by cancer type and framework.
counts <- as.data.frame(table(framework = predictions$framework, cancer_type = predictions$cancer_type))
order_types <- names(sort(tapply(counts$Freq, counts$cancer_type, sum)))
counts$cancer_type <- factor(counts$cancer_type, levels = order_types)
counts$framework <- factor(counts$framework, levels = frameworks)
panel_d <- ggplot(counts, aes(Freq, cancer_type, fill = framework)) +
  geom_col(alpha = .6) + facet_wrap(~framework, nrow = 1) +
  scale_fill_manual(values = colours, name = "Model") +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_continuous(breaks = seq(0, 150, 30)) +
  labs(x = "n\u00b0 samples", y = NULL) +
  theme_minimal(base_size = 5.5) +
  theme(strip.text = element_text(face = "bold"),
        axis.text = element_text(colour = "black"))

# Panel e: predictors accounting for 80% of cumulative gain in each split-0
# model, by data type. Feature importance tables are written for reference.
importance <- do.call(rbind, lapply(frameworks, function(fw) {
  inputs <- pi_load_inputs(final_path, cluster_path = "data/processed/clustering_data.RData",
                           frameworks = fw)
  map <- inputs$extras[[fw]]$feature_map
  rm(inputs); invisible(gc())
  ck <- readRDS(file.path(fixed, fw, "checkpoints", "split_000_seed_1234.rds"))
  booster <- tryCatch(xgb.load.raw(ck$model_raw, as_booster = TRUE),
                      error = function(e) xgb.load.raw(ck$model_raw))
  imp <- as.data.frame(xgb.importance(model = booster))
  row <- match(imp$Feature, map$feature)
  imp$variable <- map$source_name[row]
  # Categories as in the original panel: clinical (age, purity), expression,
  # methylation, pathway-level mutations (the first 272 columns of the
  # mutation block) and gene-level mutations.
  type <- c(expression = "exp", methylation = "meth", other = "other")[map$feature_type[row]]
  first_mut <- pi_archive_spec()$mutation_start[[fw]]
  column <- map$source_column[row]
  type[type == "other" & imp$variable %in% c("age", "purity")] <- "clin"
  type[type == "other" & column < first_mut + 272L] <- "path_mut"
  type[type == "other"] <- "mut"
  imp$data_type <- unname(type)
  imp <- imp[order(-imp$Gain), ]
  imp$cumulative_gain <- cumsum(imp$Gain) / sum(imp$Gain)
  imp$in_80pct <- imp$cumulative_gain <= .80
  imp$rank <- seq_len(nrow(imp))
  utils::write.csv(imp[, c("rank", "variable", "data_type", "Gain", "cumulative_gain", "in_80pct")],
                   file.path(out, paste0("importance_split0_", fw, ".csv")), row.names = FALSE)
  cbind(framework = fw, imp)
}))
top <- importance[importance$in_80pct, ]
type_levels <- c("clin", "exp", "meth", "path_mut", "mut")
panel_e_data <- as.data.frame(table(framework = factor(top$framework, frameworks),
                                    data_type = factor(top$data_type, type_levels)))
panel_e_data <- panel_e_data[panel_e_data$Freq > 0, ]
panel_e <- ggplot(panel_e_data, aes(framework, data_type, size = Freq, fill = data_type)) +
  geom_point(shape = 21, colour = "black", stroke = .3) +
  scale_size_continuous(range = c(1.5, 7.5), breaks = c(100, 500, 1000, 1500),
                        limits = c(1, NA), name = "Count") +
  # Original colours per data type; categories without retained predictors
  # (e.g. gene-level mutations) are not shown.
  scale_fill_manual(values = setNames(RColorBrewer::brewer.pal(5, "Set1"), type_levels),
                    guide = "none") +
  scale_y_discrete(drop = TRUE) +
  labs(x = "Signature", y = "Omic") +
  theme_minimal(base_size = 5.5) +
  theme(axis.text = element_text(colour = "black"), axis.title = element_text(size = 7.5),
        legend.justification = "top")
utils::write.csv(panel_e_data, file.path(out, "panel_e_counts.csv"), row.names = FALSE)

# Panel f: Cluster 1 test-set F1. Diamonds: reconstructed original split with
# its 95% bootstrap interval; points: the 20 repeated splits (vertical jitter
# only, so F1 values are exact). Labels give the Cluster 1 proportion, since
# the three classifiers predict different targets.
metrics <- read.csv(file.path(dirname(fixed), "all_split_metrics.csv"),
                    colClasses = c(split_id = "character"))
metrics <- metrics[metrics$variant == "archived", ]
lanes <- rev(frameworks)
metrics$row <- match(metrics$framework, lanes)
original <- metrics[metrics$split_id == "0", ]
repeated <- metrics[metrics$split_id != "0", ]
prevalence <- round(100 * tapply(metrics$prevalence_cluster1, metrics$framework, mean))
point_lab <- "Repeated splits (n = 20)"; diamond_lab <- "Original split, 95% CI"
panel_f <- ggplot() +
  geom_point(data = repeated, aes(f1, row - .13, colour = framework, shape = point_lab),
             position = position_jitter(width = 0, height = .065, seed = 42),
             size = .85, alpha = .65) +
  geom_errorbar(data = original, aes(xmin = f1_lower_95, xmax = f1_upper_95, y = row + .13,
                                     colour = framework),
                orientation = "y", width = .075, linewidth = .35) +
  geom_point(data = original, aes(f1, row + .13, colour = framework, shape = diamond_lab),
             fill = "white", size = 1.65, stroke = .4) +
  scale_colour_manual(values = colours, guide = "none") +
  scale_shape_manual(name = NULL, breaks = c(point_lab, diamond_lab),
                     values = setNames(c(16, 23), c(point_lab, diamond_lab))) +
  scale_y_continuous(breaks = seq_along(lanes),
                     labels = sprintf("%s\nCluster 1: %d%%", lanes, prevalence[lanes])) +
  scale_x_continuous(breaks = seq(.5, 1, .1)) +
  coord_cartesian(xlim = c(.5, 1), ylim = c(.65, 3.35)) +
  labs(x = "Test-set F1 (Cluster 1)", y = NULL) +
  theme_minimal(base_size = 5.5) +
  theme(axis.text = element_text(colour = "black"), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), legend.position = "top",
        legend.key.width = grid::unit(8, "pt"), legend.key.height = grid::unit(8, "pt"),
        legend.margin = margin(0, 0, 2, 0))

# Same format as the submitted Figure 8: 6.8 x 7.35 in at 600 dpi
# (4080 x 4411 px), LZW-compressed TIFF, a/b, c/d, e/f layout.
figure <- (accuracy_panel("Drews") | accuracy_panel("Steele")) /
  (accuracy_panel("Tao") | panel_d) /
  (panel_e + panel_f + plot_layout(widths = c(.5, .5))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9))
size <- c(width = 4080 / 600, height = 4411 / 600)
ggsave(file.path(out, "Figure8.tiff"), figure, width = size[["width"]], height = size[["height"]],
       dpi = 600, compression = "lzw", bg = "white")
ggsave(file.path(out, "Figure8.png"), figure, width = size[["width"]], height = size[["height"]],
       dpi = 300, bg = "white")
ggsave(file.path(out, "Figure8.pdf"), figure, width = size[["width"]], height = size[["height"]])

summary <- do.call(rbind, lapply(frameworks, function(fw) {
  p <- predictions[predictions$framework == fw, ]
  i <- importance[importance$framework == fw, ]
  data.frame(framework = fw, n_test = nrow(p), accuracy = mean(p$truth == p$predicted),
             predictors_80pct = sum(i$in_80pct),
             purity_rank = if (any(i$variable == "purity")) i$rank[i$variable == "purity"] else NA,
             purity_in_80pct = any(i$variable == "purity" & i$in_80pct))
}))
utils::write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
print(summary)
