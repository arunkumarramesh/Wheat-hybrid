#!/usr/bin/env Rscript

library(data.table)

dt <- fread("merged_CG_symmetric_all.txt.gz")
out_suffix <- "_cg"

dt <- dt[cov_CS > 10 & cov_CSxP > 10 & cov_P > 10]
dt[, subgenome := factor(sub("^chr[0-9]+([ABD])$", "\\1", chr),levels = c("A","B","D"))]
dt <- dt[!is.na(subgenome)]

summarize_ci <- function(x,value_col,mean_col){
  out <- x[,.(N = sum(!is.na(get(value_col))),mean_value = mean(get(value_col),na.rm = TRUE),sd_value = sd(get(value_col),na.rm = TRUE)),by = .(subgenome,sample)]
  out[,se := sd_value / sqrt(N)]
  out[,ci := qt(0.975,df = N - 1) * se]
  out[,`:=`(ci_low = mean_value - ci,ci_high = mean_value + ci)]
  setnames(out,"mean_value",mean_col)
  out[]
}

meth_long <- melt(dt,id.vars = "subgenome",measure.vars = c("pct_CS","pct_CSxP","pct_P"),variable.name = "sample",value.name = "pct")
meth_long[,sample := factor(sub("^pct_","",sample),levels = c("CS","CSxP","P"))]

cov_long <- melt(dt,id.vars = "subgenome",measure.vars = c("cov_CS","cov_CSxP","cov_P"),variable.name = "sample",value.name = "coverage")
cov_long[,sample := factor(sub("^cov_","",sample),levels = c("CS","CSxP","P"))]

meth_summary <- summarize_ci(meth_long,"pct","mean_pct")
cov_summary <- summarize_ci(cov_long,"coverage","mean_cov")

fwrite(meth_summary,paste0("met_classify",out_suffix,".tsv"),sep = "\t")
fwrite(cov_summary,paste0("mean_coverage",out_suffix,".tsv"),sep = "\t")
