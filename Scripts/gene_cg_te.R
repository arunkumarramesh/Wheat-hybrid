#!/usr/bin/env Rscript

library(data.table)

dt <- fread("merged_CG_symmetric_te.txt.gz")
out_suffix <- "_cg"

cols_no_v2 <- c("chr","start","end","te_id","te_strand","te_class","te_family","te_consensus","te_consensus_pct","te_status","te_copie","te_compo","te_post","te_length","subgenome","nearest_gene_id","distance_to_gene","centromere_status","centromere_interval","n_sites","pct_CS","cov_CS","pct_CSxP","cov_CSxP","pct_P","cov_P")
cols_with_v2 <- c("chr","start","end","te_id","te_strand","te_class","te_family","te_consensus","te_consensus_pct","te_status","te_copie","te_compo","te_post","te_length","subgenome","nearest_gene_id_v1.1","nearest_gene_id_v2.1","distance_to_gene","centromere_status","centromere_interval","n_sites","pct_CS","cov_CS","pct_CSxP","cov_CSxP","pct_P","cov_P")

if (ncol(dt) == length(cols_no_v2)) setnames(dt, cols_no_v2)
if (ncol(dt) == length(cols_with_v2)) setnames(dt, cols_with_v2)
if (!ncol(dt) %in% c(length(cols_no_v2), length(cols_with_v2))) stop("Unexpected number of columns: ", ncol(dt))

dt <- dt[cov_CS > 10 & cov_CSxP > 10]
dt <- dt[!is.na(pct_CS) & !is.na(pct_CSxP)]
dt <- dt[subgenome %in% c("A","B","D")]
dt$subgenome <- factor(dt$subgenome, levels = c("A","B","D"))

dt$distance_to_gene <- as.numeric(dt$centromere_status)
dt[, distance_to_gene_class := ifelse(distance_to_gene < 1000, "<1000", ifelse(distance_to_gene > 5000, ">5000", "1000-5000"))]
dt$distance_to_gene_class <- factor(dt$distance_to_gene_class, levels = c("<1000","1000-5000",">5000"))

dt[, diff_CSxP_minus_CS := pct_CSxP - pct_CS]

ci_summary <- function(x, by_cols, value_col) {
  out <- x[, list(n_TEs = uniqueN(te_id), total_sites = sum(n_sites, na.rm = TRUE), N = sum(!is.na(get(value_col))), mean_value = mean(get(value_col), na.rm = TRUE), sd_value = sd(get(value_col), na.rm = TRUE)), by = by_cols]
  out$se <- out$sd_value / sqrt(out$N)
  out$ci <- qt(0.975, df = out$N - 1) * out$se
  out$ci_low <- out$mean_value - out$ci
  out$ci_high <- out$mean_value + out$ci
  out
}

make_meth_summary <- function(dt, group_col, grouping_type) {
  meth_long <- rbind(data.table(grouping_type = grouping_type, group = dt[[group_col]], sample = "CS", pct = dt$pct_CS, te_id = dt$te_id, n_sites = dt$n_sites), data.table(grouping_type = grouping_type, group = dt[[group_col]], sample = "CSxP", pct = dt$pct_CSxP, te_id = dt$te_id, n_sites = dt$n_sites))
  meth_summary <- ci_summary(meth_long, c("grouping_type","group","sample"), "pct")
  setnames(meth_summary, "mean_value", "mean_pct")
  meth_summary
}

make_diff_summary <- function(dt, group_col, grouping_type) {
  diff_summary <- ci_summary(dt[, .(grouping_type = grouping_type, group = get(group_col), diff_CSxP_minus_CS, te_id, n_sites)], c("grouping_type","group"), "diff_CSxP_minus_CS")
  setnames(diff_summary, "mean_value", "mean_diff_CSxP_minus_CS")
  diff_summary
}

meth_summary_all <- rbindlist(list(make_meth_summary(dt, "subgenome", "subgenome"), make_meth_summary(dt, "te_class", "te_class"), make_meth_summary(dt, "te_status", "te_status"), make_meth_summary(dt, "centromere_interval", "centromere_interval"), make_meth_summary(dt, "distance_to_gene_class", "distance_to_gene_class")), use.names = TRUE, fill = TRUE)

diff_summary_all <- rbindlist(list(make_diff_summary(dt, "subgenome", "subgenome"), make_diff_summary(dt, "te_class", "te_class"), make_diff_summary(dt, "te_status", "te_status"), make_diff_summary(dt, "centromere_interval", "centromere_interval"), make_diff_summary(dt, "distance_to_gene_class", "distance_to_gene_class")), use.names = TRUE, fill = TRUE)

write.table(meth_summary_all, file = paste0("met_classify", out_suffix, "_all_groups.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(diff_summary_all, file = paste0("met_diff_", out_suffix, "_all_groups.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)


dt <- fread("merged_CG_symmetric_te.txt.gz")
out_suffix <- "_cg"

cols_no_v2 <- c("chr","start","end","te_id","te_strand","te_class","te_family","te_consensus","te_consensus_pct","te_status","te_copie","te_compo","te_post","te_length","subgenome","nearest_gene_id","distance_to_gene","centromere_status","centromere_interval","n_sites","pct_CS","cov_CS","pct_CSxP","cov_CSxP","pct_P","cov_P")
cols_with_v2 <- c("chr","start","end","te_id","te_strand","te_class","te_family","te_consensus","te_consensus_pct","te_status","te_copie","te_compo","te_post","te_length","subgenome","nearest_gene_id_v1.1","nearest_gene_id_v2.1","distance_to_gene","centromere_status","centromere_interval","n_sites","pct_CS","cov_CS","pct_CSxP","cov_CSxP","pct_P","cov_P")

if (ncol(dt) == length(cols_no_v2)) setnames(dt, cols_no_v2)
if (ncol(dt) == length(cols_with_v2)) setnames(dt, cols_with_v2)
if (!ncol(dt) %in% c(length(cols_no_v2), length(cols_with_v2))) stop("Unexpected number of columns: ", ncol(dt))

dt <- dt[centromere_status < 1000]

dt[, meth_CS := pct_CS * cov_CS / 100]
dt[, meth_CSxP := pct_CSxP * cov_CSxP / 100]
dt[, meth_P := pct_P * cov_P / 100]

gene_meth <- dt[, .(meth_CS = sum(meth_CS, na.rm = TRUE), cov_CS = sum(cov_CS, na.rm = TRUE), meth_CSxP = sum(meth_CSxP, na.rm = TRUE), cov_CSxP = sum(cov_CSxP, na.rm = TRUE), meth_P = sum(meth_P, na.rm = TRUE), cov_P = sum(cov_P, na.rm = TRUE), n_TEs = uniqueN(te_id), n_sites = sum(n_sites, na.rm = TRUE)), by = distance_to_gene]

gene_meth[, pct_CS := 100 * meth_CS / cov_CS]
gene_meth[, pct_CSxP := 100 * meth_CSxP / cov_CSxP]
gene_meth[, pct_P := 100 * meth_P / cov_P]

gene_meth <- gene_meth[, .(gene_id = distance_to_gene, n_TEs, n_sites, pct_CS, cov_CS, pct_CSxP, cov_CSxP, pct_P, cov_P)]
gene_meth <- gene_meth[cov_CS > 10 & cov_CSxP > 10]

fwrite(gene_meth, "gene_level_TE_CG.tsv", sep = "\t")
