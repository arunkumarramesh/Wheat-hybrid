#!/usr/bin/env Rscript

library(data.table)

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
