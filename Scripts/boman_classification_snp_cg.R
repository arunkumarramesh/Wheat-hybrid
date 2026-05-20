#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(scales)

dt <- fread("merged_CG_symmetric_all.txt.gz")
snp_counts <- fread("cs_par_snps_50kb_counts.tsv",header = FALSE)
out_suffix <- "_cg"

colnames(snp_counts) <- c("chr","window_start","window_end","SNP_count_50kb")
snp_counts[, chr := sub("^triticum_aestivum\\.([0-9]+[ABD])$","chr\\1",chr)]

dt <- dt[cov_CS > 10 & cov_CSxP > 10 & cov_P > 10]

dt$subgenome <- sub("^chr[0-9]+([ABD])$","\\1",dt$chr)
dt <- dt[dt$subgenome %in% c("A","B","D"),]
dt$subgenome <- factor(dt$subgenome,levels = c("A","B","D"))

dt[, window_start := floor(pos / 50000) * 50000]
dt[, window_end := window_start + 50000]

dt <- merge(dt,snp_counts,by = c("chr","window_start","window_end"),all.x = FALSE)

dt[, SNP_density_group := fifelse(SNP_count_50kb > 1000,"High SNP density",
                            fifelse(SNP_count_50kb < 10,"Low SNP density",NA_character_))]

dt <- dt[!is.na(SNP_density_group)]
dt$SNP_density_group <- factor(dt$SNP_density_group,levels = c("Low SNP density","High SNP density"))

dt$A <- dt$pct_CS / 100
dt$H <- dt$pct_CSxP / 100
dt$B <- dt$pct_P / 100

meth_long <- rbind(
  data.table(subgenome = dt$subgenome,SNP_density_group = dt$SNP_density_group,sample = "CS",pct = dt$pct_CS),
  data.table(subgenome = dt$subgenome,SNP_density_group = dt$SNP_density_group,sample = "CSxP",pct = dt$pct_CSxP),
  data.table(subgenome = dt$subgenome,SNP_density_group = dt$SNP_density_group,sample = "P",pct = dt$pct_P)
)

ci_summary <- function(x,by_cols,value_col) {
  out <- x[, list(N = sum(!is.na(get(value_col))),mean_value = mean(get(value_col),na.rm = TRUE),sd_value = sd(get(value_col),na.rm = TRUE)), by = by_cols]
  out$se <- out$sd_value / sqrt(out$N)
  out$ci <- qt(0.975,df = out$N - 1) * out$se
  out$ci_low <- out$mean_value - out$ci
  out$ci_high <- out$mean_value + out$ci
  return(out)
}

meth_summary <- ci_summary(meth_long,c("subgenome","SNP_density_group","sample"),"pct")
setnames(meth_summary,"mean_value","mean_pct")

dt$x <- dt$H - dt$A
dt$y <- dt$H - dt$B
dt$radius <- sqrt(dt$x^2 + dt$y^2)
dt$angle_deg <- (atan2(dt$y,dt$x) * 180 / pi) %% 360

circ_dist <- function(a,b) {
  d <- abs(a - b)
  pmin(d,360 - d)
}

nearest_sector_center <- function(angle_deg) {
  centers <- c(0,45,90,135,180,225,270,315)
  centers[which.min(circ_dist(angle_deg,centers))]
}

sector_to_class <- function(sector_center) {
  if (sector_center %in% c(0,180)) {
    "P_dominant"
  } else if (sector_center %in% c(90,270)) {
    "CS_dominant"
  } else if (sector_center %in% c(135,315)) {
    "additive"
  } else if (sector_center == 45) {
    "overdominant"
  } else if (sector_center == 225) {
    "underdominant"
  } else {
    NA_character_
  }
}

dt$sector_center <- vapply(dt$angle_deg,nearest_sector_center,numeric(1))
dt$max_mC <- pmax(dt$A,dt$H,dt$B,na.rm = TRUE)
dt <- dt[max_mC >= 0.01]

dt$category <- ifelse(dt$radius < 0.1,"conserved_mC",vapply(dt$sector_center,sector_to_class,character(1)))
dt$category <- factor(dt$category,levels = c("conserved_mC","additive","CS_dominant","P_dominant","overdominant","underdominant"))

summary_dt <- as.data.table(table(dt$subgenome,dt$SNP_density_group,dt$category))
colnames(summary_dt) <- c("subgenome","SNP_density_group","category","N")
summary_dt$subgenome <- factor(summary_dt$subgenome,levels = c("A","B","D"))
summary_dt$SNP_density_group <- factor(summary_dt$SNP_density_group,levels = c("Low SNP density","High SNP density"))
summary_dt$category <- factor(summary_dt$category,levels = c("conserved_mC","additive","CS_dominant","P_dominant","overdominant","underdominant"))
summary_dt$proportion <- ave(summary_dt$N,summary_dt$subgenome,summary_dt$SNP_density_group,FUN = function(x) x / sum(x))

summary_dt_all <- dt[, .N, by = .(SNP_density_group,category)]
summary_dt_all[, proportion := N / sum(N), by = SNP_density_group]
summary_dt_all$SNP_density_group <- factor(summary_dt_all$SNP_density_group,levels = c("Low SNP density","High SNP density"))
summary_dt_all$category <- factor(summary_dt_all$category,levels = c("conserved_mC","additive","CS_dominant","P_dominant","overdominant","underdominant"))

write.table(meth_summary,file = paste0("met_classify_snp_density",out_suffix,".tsv"),sep = "\t",quote = FALSE,row.names = FALSE)
write.table(summary_dt,file = paste0("percent_met_snp_density",out_suffix,".tsv"),sep = "\t",quote = FALSE,row.names = FALSE)
write.table(summary_dt_all,file = paste0("percent_met_snp_density_all",out_suffix,".tsv"),sep = "\t",quote = FALSE,row.names = FALSE)