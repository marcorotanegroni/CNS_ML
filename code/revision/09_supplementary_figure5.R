# Supplementary Figure 5 (Kaplan-Meier risk tables of Figure 6) regenerated
# from the saved cluster assignments, with the same steps as
# code/03_survival_analysis.Rmd. In the submitted version the labels of
# Clusters 1 and 3 were swapped in the Tao overall-survival table (panel e);
# the curves in Figure 6e and all other tables are unchanged.
#
# Run from the repository root:
#   Rscript code/revision/09_supplementary_figure5.R
# Output: results/revision/survival/

suppressMessages({library(survival); library(survminer); library(ggplot2); library(patchwork)})
load("data/processed/surv_data.RData")
out <- "results/revision/survival"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

palette <- c("#72F59A", "#F57D72", "#F7C073", "#8672F5")
spec <- list(Drews = list(drew_f, c("CX1", "CX2", "CX3", "CX5"), mc_drew_active_4),
             Steele = list(steele_f, c("CN1", "CN2", "CN9", "CN17"), mc_steele_active_4),
             Tao = list(tao_f, c("Sig1", "Sig2", "Sig5", "Sig7"), mc_tao_active_4))

risk_table <- function(framework, endpoint) {
  s <- spec[[framework]]
  cut <- surv_cutpoint(s[[1]], time = paste0(endpoint, ".time"), event = endpoint,
                       variables = s[[2]], progressbar = FALSE)
  res.cat <- surv_categorize(cut)
  res.cat$cluster <- s[[3]]$classification[rownames(res.cat)]
  fit <- surv_fit(as.formula(paste0("Surv(", endpoint, ".time, ", endpoint, ") ~ cluster")),
                  data = res.cat)
  k <- length(unique(res.cat$cluster))
  plot <- ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = FALSE, pval = TRUE,
                     palette = palette[seq_len(k)], censor = FALSE,
                     fontsize = 2.6,
                     tables.theme = theme_survminer(font.main = 9, font.tickslab = c(7, "bold")))
  at <- summary(fit, times = ggplot_build(plot$table)$layout$panel_params[[1]]$x$breaks, extend = TRUE)
  numbers <- data.frame(framework = framework, endpoint = endpoint,
                        cluster = sub("cluster=", "", at$strata), time = at$time, n_risk = at$n.risk)
  list(table = plot$table, numbers = numbers)
}

panels <- list(c("Drews", "OS"), c("Drews", "PFI"), c("Steele", "OS"),
               c("Steele", "PFI"), c("Tao", "OS"), c("Tao", "PFI"))
tables <- lapply(panels, function(p) risk_table(p[1], p[2]))
plots <- lapply(seq_along(tables), function(i) {
  g <- tables[[i]]$table +
    theme(axis.text.x = element_text(size = 6.5, face = "plain", colour = "black"),
          axis.title.x = element_text(size = 7.5, colour = "black"),
          axis.title.y = element_text(size = 7.5, colour = "black"))
  g <- g + labs(title = if (i <= 2) "Number at risk" else NULL,
                x = if (i >= 5) "Time in days" else NULL,
                y = if (i %% 2 == 1) "Strata" else NULL)
  g
})
figure <- wrap_plots(plots, ncol = 2) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 12))
ggsave(file.path(out, "Supplementary_Figure5.png"), figure, width = 6.8, height = 4.6, dpi = 600, bg = "white")
ggsave(file.path(out, "Supplementary_Figure5.pdf"), figure, width = 6.8, height = 4.6)
numbers <- do.call(rbind, lapply(tables, `[[`, "numbers"))
utils::write.csv(numbers, file.path(out, "supplementary_figure5_numbers_at_risk.csv"), row.names = FALSE)
print(numbers[numbers$framework == "Tao" & numbers$endpoint == "OS", ])
