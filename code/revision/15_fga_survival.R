# Scalar copy number baseline in the survival analysis (Reviewer 2): the
# cancer-type-stratified Cox models of Figure 7 (cluster + age + purity,
# strata(cancer_type); code/03_survival_analysis.Rmd), refitted with fraction
# genome altered (FGA) as an additional covariate.
#
# For each compendium and endpoint:
#   figure7       the Figure 7 model on all its patients (reproduces Figure 7);
#   figure7_fga_n the same model restricted to patients with an FGA value;
#   with_fga      that model plus FGA (same patients as figure7_fga_n);
#   fga_only      FGA + age + purity, without clusters (same patients).
# FGA hazard ratios are per 0.1 increase. Holm adjustment of the cluster terms
# across the six models follows Figure 7, separately for each model set.
#
# Run from the repository root after 13_fetch_fga.R:
#   Rscript code/revision/15_fga_survival.R
# Output: results/revision/fga_baseline/survival_cox.csv

suppressMessages({library(survival); library(survminer)})
load("data/processed/surv_data.RData")
fga <- utils::read.csv("data/processed/fga_pancanatlas.csv")
fga <- setNames(fga$fga, fga$patient)
out <- "results/revision/fga_baseline"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Same data preparation as run_stratified_cox() in 03_survival_analysis.Rmd.
cox_data <- function(df, cluster_vec, sig_vars, time_var, event_var) {
  a <- surv_cutpoint(df, time = time_var, event = event_var,
                     variables = sig_vars, progressbar = FALSE)
  res.cat <- surv_categorize(a)
  res.cat$cluster <- cluster_vec[rownames(res.cat)]
  res.cat$cancer_type <- df$cancer_type
  tt <- intersect(rownames(res.cat), rownames(purezza))
  d <- data.frame(time = res.cat[tt, time_var], event = res.cat[tt, event_var],
                  cluster = factor(res.cat[tt, "cluster"]), cancer_type = res.cat[tt, "cancer_type"],
                  purity = purezza[tt, "purity"], age = df[tt, 1], row.names = tt)
  d$fga10 <- unname(fga[tt]) * 10
  d
}

terms <- function(fit, model, framework, endpoint) {
  s <- summary(fit)
  co <- s$coefficients; ci <- s$conf.int
  keep <- grepl("^cluster|^fga10$", rownames(co))
  data.frame(framework = framework, endpoint = endpoint, model = model,
             n = fit$n, events = fit$nevent,
             term = sub("^cluster", "cluster", rownames(co)[keep]),
             hr = ci[keep, "exp(coef)"], lower = ci[keep, "lower .95"],
             upper = ci[keep, "upper .95"], p = co[keep, "Pr(>|z|)"], row.names = NULL)
}

spec <- list(Drews = list(drew_f, mc_drew_active_4, c("CX1", "CX2", "CX3", "CX5")),
             Steele = list(steele_f, mc_steele_active_4, c("CN1", "CN2", "CN9", "CN17")),
             Tao = list(tao_f, mc_tao_active_4, c("Sig1", "Sig2", "Sig5", "Sig7")))
results <- list()
for (fw in names(spec)) for (ep in c("OS", "PFI")) {
  s <- spec[[fw]]
  d <- cox_data(s[[1]], s[[2]]$classification, s[[3]], paste0(ep, ".time"), ep)
  f0 <- coxph(Surv(time, event) ~ cluster + age + purity + strata(cancer_type), data = d)
  dn <- d[!is.na(d$fga10), ]
  f1 <- coxph(Surv(time, event) ~ cluster + age + purity + strata(cancer_type), data = dn)
  f2 <- coxph(Surv(time, event) ~ cluster + fga10 + age + purity + strata(cancer_type), data = dn)
  f3 <- coxph(Surv(time, event) ~ fga10 + age + purity + strata(cancer_type), data = dn)
  results[[paste(fw, ep)]] <- rbind(terms(f0, "figure7", fw, ep), terms(f1, "figure7_fga_n", fw, ep),
                                    terms(f2, "with_fga", fw, ep), terms(f3, "fga_only", fw, ep))
}
res <- do.call(rbind, results)
res$p_holm <- NA_real_
for (m in unique(res$model)) {
  i <- res$model == m & res$term != "fga10"
  res$p_holm[i] <- p.adjust(res$p[i], method = "holm")
}
utils::write.csv(res, file.path(out, "survival_cox.csv"), row.names = FALSE)
show <- res[res$term %in% c("cluster2", "fga10"), ]
show[, c("hr", "lower", "upper")] <- round(show[, c("hr", "lower", "upper")], 2)
show$p_holm <- signif(show$p_holm, 2)
print(show[, c("framework", "endpoint", "model", "n", "term", "hr", "lower", "upper", "p_holm")], row.names = FALSE)
