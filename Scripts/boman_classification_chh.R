#!/usr/bin/env Rscript

# Classification from: https://pmc.ncbi.nlm.nih.gov/articles/PMC11134931/

library(data.table)
library(ggplot2)
library(scales)

dt <- fread("merged_CHH_all.txt.gz")
dt <- dt[, -3]
out_suffix <- "_chh"

# Keep only sites with per-sample coverage > 10
dt <- dt[cov_CS > 10 & cov_CSxP > 10 & cov_P > 10]

dt$subgenome <- sub("^chr[0-9]+([ABD])$", "\\1", dt$chr)
dt <- dt[dt$subgenome %in% c("A", "B", "D"), ]
dt$subgenome <- factor(dt$subgenome, levels = c("A", "B", "D"))

dt$A <- dt$pct_CS / 100
dt$H <- dt$pct_CSxP / 100
dt$B <- dt$pct_P / 100

meth_long <- rbind(
  data.table(subgenome = dt$subgenome, sample = "CS", pct = dt$pct_CS),
  data.table(subgenome = dt$subgenome, sample = "CSxP", pct = dt$pct_CSxP),
  data.table(subgenome = dt$subgenome, sample = "P", pct = dt$pct_P)
)

cov_long <- rbind(
  data.table(subgenome = dt$subgenome, sample = "CS", coverage = dt$cov_CS),
  data.table(subgenome = dt$subgenome, sample = "CSxP", coverage = dt$cov_CSxP),
  data.table(subgenome = dt$subgenome, sample = "P", coverage = dt$cov_P)
)

ci_summary <- function(x, by_cols, value_col) {
  out <- x[, list(N = sum(!is.na(get(value_col))), mean_value = mean(get(value_col), na.rm = TRUE), sd_value = sd(get(value_col), na.rm = TRUE)), by = by_cols]
  
  out$se <- out$sd_value / sqrt(out$N)
  out$ci <- qt(0.975, df = out$N - 1) * out$se
  out$ci_low <- out$mean_value - out$ci
  out$ci_high <- out$mean_value + out$ci
  
  return(out)
}

meth_summary <- ci_summary(meth_long, c("subgenome", "sample"), "pct")
setnames(meth_summary, "mean_value", "mean_pct")

cov_summary <- ci_summary(cov_long, c("subgenome", "sample"), "coverage")
setnames(cov_summary, "mean_value", "mean_cov")

# Coordinates in the Boman plane

dt$x <- dt$H - dt$A
dt$y <- dt$H - dt$B
dt$radius <- sqrt(dt$x^2 + dt$y^2)
dt$angle_deg <- (atan2(dt$y, dt$x) * 180 / pi) %% 360

circ_dist <- function(a, b) {
  d <- abs(a - b)
  pmin(d, 360 - d)
}

nearest_sector_center <- function(angle_deg) {
  centers <- c(0, 45, 90, 135, 180, 225, 270, 315)
  centers[which.min(circ_dist(angle_deg, centers))]
}

sector_to_class <- function(sector_center) {
  if (sector_center %in% c(0, 180)) {
    "P_dominant"
  } else if (sector_center %in% c(90, 270)) {
    "CS_dominant"
  } else if (sector_center %in% c(135, 315)) {
    "additive"
  } else if (sector_center == 45) {
    "overdominant"
  } else if (sector_center == 225) {
    "underdominant"
  } else {
    NA_character_
  }
}
dt$sector_center <- vapply(dt$angle_deg, nearest_sector_center, numeric(1))

dt$max_mC <- pmax(dt$A,dt$H,dt$B,na.rm = TRUE)
dt <- dt[max_mC >= 0.01]

# between 0 to 1
dt$category <- ifelse(dt$radius < 0.1, "conserved_mC",vapply(dt$sector_center, sector_to_class, character(1)))

dt$category <- factor(dt$category,levels = c("conserved_mC","additive","CS_dominant","P_dominant","overdominant","underdominant"))

summary_dt <- as.data.table(table(dt$subgenome, dt$category))
colnames(summary_dt) <- c("subgenome", "category", "N")
summary_dt$subgenome <- factor(summary_dt$subgenome, levels = c("A", "B", "D"))
summary_dt$category <- factor(summary_dt$category,levels = c("conserved_mC","additive","CS_dominant","P_dominant","overdominant","underdominant"))

summary_dt$proportion <- ave(summary_dt$N,summary_dt$subgenome,FUN = function(x) x / sum(x))

panel_counts_all <- as.data.table(table(dt$subgenome))
colnames(panel_counts_all) <- c("subgenome", "N")
panel_counts_all$subgenome <- factor(panel_counts_all$subgenome, levels = c("A", "B", "D"))
panel_counts_all <- panel_counts_all[order(panel_counts_all$subgenome), ]

count_map_all <- setNames(paste0(panel_counts_all$subgenome, " (n=", panel_counts_all$N, ")"),as.character(panel_counts_all$subgenome))

write.table(meth_summary,file = paste0("met_classify", out_suffix, ".tsv"),sep = "\t",quote = FALSE,row.names = FALSE)
write.table(summary_dt,file = paste0("percent_met", out_suffix, ".tsv"),sep = "\t",quote = FALSE,row.names = FALSE)
write.table(cov_summary,file = paste0("mean_coverage", out_suffix, ".tsv"),sep = "\t",quote = FALSE,row.names = FALSE)