# Fraction genome altered (FGA) across the signature-based clusters (notebook
# 06, Supplementary Figure 9a): how much of the variation in fraction genome altered (FGA) the signature-based
# clusters capture, on all tumours of each compendium with an FGA value.
#
# For each compendium: FGA median and quartiles by cluster, and the share of
# FGA variance explained by cluster membership (R2 of FGA ~ cluster). Because
# cancer types differ in FGA, the share explained by clusters beyond cancer type
# is also given (partial R2 of cluster in FGA ~ cancer_type + cluster).
#
# Run from the repository root:
#   Rscript code/scripts/fga_clustering.R
# Outputs: results/fraction_genome_altered/clustering_fga_by_cluster.csv
#          results/fraction_genome_altered/clustering_fga_variance.csv

cl <- new.env(); load("data/processed/clustering_data.RData", envir = cl)
sv <- new.env(); load("data/processed/surv_data.RData", envir = sv)
fga <- utils::read.csv("data/processed/fga_pancanatlas.csv")
fga <- setNames(fga$fga, fga$patient)
out <- "results/fraction_genome_altered"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

fits <- c(Drews = "mc_drew_active_4", Steele = "mc_steele_active_4", Tao = "mc_tao_active_4")
actives <- c(Drews = "active_drew", Steele = "active_steele", Tao = "active_tao")
by_cluster <- list(); variance <- list()
for (fw in names(fits)) {
  cls <- cl[[fits[fw]]]$classification
  cancer <- as.character(cl[[actives[fw]]][names(cls), "cancer_type"])
  d <- data.frame(patient = names(cls), cluster = factor(cls), cancer_type = cancer,
                  fga = unname(fga[names(cls)]))
  n_all <- nrow(d)
  d <- d[!is.na(d$fga), ]
  q <- do.call(rbind, lapply(split(d$fga, d$cluster), function(x)
    data.frame(n = length(x), median = median(x), q1 = unname(quantile(x, .25)),
               q3 = unname(quantile(x, .75)))))
  by_cluster[[fw]] <- data.frame(framework = fw, cluster = rownames(q), q, row.names = NULL)
  r2 <- summary(lm(fga ~ cluster, data = d))$r.squared
  base <- lm(fga ~ cancer_type, data = d)
  full <- lm(fga ~ cancer_type + cluster, data = d)
  rss <- function(m) sum(residuals(m)^2)
  variance[[fw]] <- data.frame(framework = fw, n_with_fga = nrow(d), n_clustered = n_all,
    r2_cluster = r2, r2_cancer_type = summary(base)$r.squared,
    r2_cancer_type_plus_cluster = summary(full)$r.squared,
    partial_r2_cluster_given_cancer_type = (rss(base) - rss(full)) / rss(base))
}
by_cluster <- do.call(rbind, by_cluster); variance <- do.call(rbind, variance)
utils::write.csv(by_cluster, file.path(out, "clustering_fga_by_cluster.csv"), row.names = FALSE)
utils::write.csv(variance, file.path(out, "clustering_fga_variance.csv"), row.names = FALSE)
print(by_cluster, digits = 2, row.names = FALSE)
print(variance, digits = 3, row.names = FALSE)
