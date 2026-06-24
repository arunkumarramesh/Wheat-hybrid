library(data.table)
library(dplyr)
library(ggplot2)
library(scales)

## first just by gene or promoter

cat_cols <- c(conserved_mC = "#1A1A1A", additive = "#D9D9D9", CS_dominant = "#0072B2", P_dominant = "#CC79A7", overdominant = "#E69F00", underdominant = "#1B7837")
sample_cols <- c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7")

ci_summary <- function(x, by_cols, value_col) {
  out <- x[, list(N = sum(!is.na(get(value_col))), mean_value = mean(get(value_col), na.rm = TRUE), sd_value = sd(get(value_col), na.rm = TRUE)), by = by_cols]
  out$se <- out$sd_value / sqrt(out$N)
  out$ci <- qt(0.975, df = out$N - 1) * out$se
  out$ci_low <- out$mean_value - out$ci
  out$ci_high <- out$mean_value + out$ci
  return(out)
}

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

process_methylation_file <- function(infile, out_suffix, context, region, remove_third_col = FALSE) {
  dt <- fread(infile)
  
  if (remove_third_col) {
    dt <- dt[, -3, with = FALSE]
  }
  
  ## temp
  #n <- min(10000, nrow(dt))
  #idx <- sort(sample.int(nrow(dt), n))
  #dt <- dt[idx]
  ## temp end
  
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
  
  meth_summary <- ci_summary(meth_long, c("subgenome", "sample"), "pct")
  setnames(meth_summary, "mean_value", "mean_pct")
  
  meth_summary$context <- context
  meth_summary$region <- region
  
  # Coordinates in the Boman plane
  
  dt$x <- dt$H - dt$A
  dt$y <- dt$H - dt$B
  dt$radius <- sqrt(dt$x^2 + dt$y^2)
  dt$angle_deg <- (atan2(dt$y, dt$x) * 180 / pi) %% 360
  
  dt$sector_center <- vapply(dt$angle_deg, nearest_sector_center, numeric(1))
  dt$max_mC <- pmax(dt$A,dt$H,dt$B,na.rm = TRUE)
  dt <- dt[max_mC >= 0.01] # only keep sites where at last one sample at least 1% methylation
  
  # between 0 to 1
  dt$category <- ifelse(dt$radius < 0.1, "conserved_mC",vapply(dt$sector_center, sector_to_class, character(1)))
  
  dt$category <- factor(dt$category,levels = c("conserved_mC","additive","CS_dominant","P_dominant","overdominant","underdominant"))
  
  summary_dt <- as.data.table(table(dt$subgenome, dt$category))
  colnames(summary_dt) <- c("subgenome", "category", "N")
  summary_dt$subgenome <- factor(summary_dt$subgenome, levels = c("A", "B", "D"))
  summary_dt$category <- factor(summary_dt$category,levels = c("conserved_mC","additive","CS_dominant","P_dominant","overdominant","underdominant"))
  summary_dt$proportion <- ave(summary_dt$N,summary_dt$subgenome,FUN = function(x) x / sum(x))
  summary_dt_all <- dt[, .N, by = category]
  summary_dt_all[, proportion := N / sum(N)]
  summary_dt_all$category <- factor(summary_dt_all$category,levels = c("conserved_mC","additive","CS_dominant","P_dominant","overdominant","underdominant"))
  
  panel_counts_all <- as.data.table(table(dt$subgenome))
  colnames(panel_counts_all) <- c("subgenome", "N")
  panel_counts_all$subgenome <- factor(panel_counts_all$subgenome, levels = c("A", "B", "D"))
  panel_counts_all <- panel_counts_all[order(panel_counts_all$subgenome), ]
  
  count_map_all <- setNames(paste0(panel_counts_all$subgenome, " (n=", panel_counts_all$N, ")"),as.character(panel_counts_all$subgenome))
  
  p1 <- ggplot(summary_dt_all,aes(x = category,y = proportion,fill = category)) +
    geom_col(width = 0.8,linewidth = 0.2) +
    geom_text(aes(label = percent(proportion,accuracy = 0.1)),vjust = -0.35,size = 3.2,fontface = "bold") +
    scale_fill_manual(values = cat_cols) +
    scale_y_continuous(labels = percent_format(accuracy = 1),expand = expansion(mult = c(0,0.12))) +
    labs(title = paste0(region," ",context," (n=",sum(summary_dt_all$N),")"),x = NULL,y = "Percentage of sites") +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(linewidth = 0.2,colour = "grey90"),plot.title = element_text(face = "bold",hjust = 0.5),axis.title = element_text(face = "bold"),axis.text.x = element_text(angle = 35,hjust = 1),legend.position = "none")
  
  p2 <- ggplot(summary_dt, aes(x = category, y = proportion, fill = category)) +
    geom_col(width = 0.8, fill = rep(cat_cols,3), linewidth = 0.2) +
    geom_text(aes(label = percent(proportion, accuracy = 0.1)), vjust = -0.35, size = 3.2, fontface = "bold")    +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12))) +
    facet_wrap(~ subgenome, nrow = 1, labeller = labeller(subgenome = count_map_all)) +
    labs(x = NULL, y = "Percentage of sites") +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), legend.position = "right", plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
  
  pdf(paste0("met_classify", out_suffix, ".pdf"), width = 7, height = 3)
  print(p2)
  dev.off()
  
  pdf(paste0("met_classify_whole", out_suffix, ".pdf"), width = 3.1, height = 3)
  print(p1)
  dev.off()
  
  return(meth_summary)
}

all_meth_summary <- list()
all_meth_summary[["cg_gene"]] <- process_methylation_file(infile = "merged_CG_symmetric_CDS.txt.gz",out_suffix = "_cg_gene",context = "CG",region = "CDS")
all_meth_summary[["cg_promoter"]] <- process_methylation_file(infile = "merged_CG_symmetric_promoter1kb.txt.gz",out_suffix = "_cg_promoter",context = "CG",region = "Promoter")
all_meth_summary[["chg_gene"]] <- process_methylation_file(infile = "merged_CHG_symmetric_CDS.txt.gz",out_suffix = "_chg_gene",context = "CHG",region = "CDS")
all_meth_summary[["chg_promoter"]] <- process_methylation_file(infile = "merged_CHG_symmetric_promoter1kb.txt.gz",out_suffix = "_chg_promoter",context = "CHG",region = "Promoter")
all_meth_summary[["chh_gene"]] <- process_methylation_file(infile = "merged_CHH_all_CDS.txt.gz",out_suffix = "_chh_gene",context = "CHH",region = "CDS",remove_third_col = TRUE)
all_meth_summary[["chh_promoter"]] <- process_methylation_file(infile = "merged_CHH_promoter1kb.txt.gz",out_suffix = "_chh_promoter",context = "CHH",region = "Promoter",remove_third_col = TRUE)

meth_summary_all <- rbindlist(all_meth_summary)
meth_summary_all$context <- factor(meth_summary_all$context, levels = c("CG", "CHG", "CHH"))
meth_summary_all$region <- factor(meth_summary_all$region, levels = c("CDS", "Promoter"))
meth_summary_all$subgenome <- factor(meth_summary_all$subgenome, levels = c("A", "B", "D"))
meth_summary_all$sample <- factor(meth_summary_all$sample, levels = c("CS", "CSxP", "P"))

meth_summary_all$x_group <- paste(meth_summary_all$region, meth_summary_all$subgenome, sep = ".")

x_levels <- c("CDS.A","CDS.B","CDS.D","Promoter.A","Promoter.B","Promoter.D")
meth_summary_all$x_group <- factor(meth_summary_all$x_group, levels = x_levels)

p_percent_all <- ggplot(meth_summary_all, aes(x = x_group, y = mean_pct, fill = sample)) +
  geom_vline(xintercept = 3.5, linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_col(aes(group = sample),position = position_dodge(width = 0.75),width = 0.65,colour = "black",linewidth = 0.25) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high, group = sample),position = position_dodge(width = 0.75),width = 0.2,colour = "black",linewidth = 0.4) +
  scale_fill_manual(values = sample_cols) +
  scale_x_discrete(labels = c("CDS.A" = "A","CDS.B" = "B","CDS.D" = "D","Promoter.A" = "A","Promoter.B" = "B","Promoter.D" = "D")) +
  facet_wrap(~ context, nrow = 1) +
  labs(x = NULL, y = "Mean Methylation (%)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), legend.position = "right", plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(size = 10),legend.position = "right",plot.margin = margin(5.5, 5.5, 28, 5.5)) +
  annotate("text", x = 2, y = -Inf, label = "CDS", vjust = 3.2, size = 3.5, fontface = "bold") +
  annotate("text", x = 5, y = -Inf, label = "Promoter", vjust = 3.2, size = 3.5, fontface = "bold") +
  coord_cartesian(clip = "off")

pdf("percent_met_combined.pdf", width = 6.5, height = 3)
p_percent_all
dev.off()


## across all sites

cat_cols <- c(conserved_mC = "#1A1A1A", additive = "#D9D9D9", CS_dominant = "#0072B2", P_dominant = "#CC79A7", overdominant = "#E69F00", underdominant = "#1B7837")

## methylation classification according to Boman et al. 2024

summary_dt <- read.table(file="percent_met_cg.tsv",header=T)
summary_dt$category <- factor(summary_dt$category, levels = names(cat_cols))
count_map_all <- aggregate(N ~ subgenome, data = summary_dt, FUN = sum, na.rm = TRUE)
count_map_all$label <- paste0(count_map_all$subgenome, " (", format(count_map_all$N, big.mark = ","), ")")
count_map_all <- setNames(count_map_all$label, count_map_all$subgenome)

pdf(paste0("met_classify_cg.pdf"), width = 7, height = 3)
ggplot(summary_dt, aes(x = category, y = proportion, fill = category)) +
  geom_col(width = 0.8, colour = "black", linewidth = 0.2) +
  geom_text(aes(label = percent(proportion, accuracy = 0.1)), vjust = -0.35, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = cat_cols, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12))) +
  facet_wrap(~ subgenome, nrow = 1, labeller = labeller(subgenome = count_map_all)) +
  labs(x = NULL, y = "Percentage of sites") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), legend.position = "right", plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
dev.off()

summary_dt <- read.table(file="percent_met_chg.tsv",header=T)
summary_dt$category <- factor(summary_dt$category, levels = names(cat_cols))
count_map_all <- aggregate(N ~ subgenome, data = summary_dt, FUN = sum, na.rm = TRUE)
count_map_all$label <- paste0(count_map_all$subgenome, " (", format(count_map_all$N, big.mark = ","), ")")
count_map_all <- setNames(count_map_all$label, count_map_all$subgenome)

pdf(paste0("met_classify_chg.pdf"), width = 7, height = 3)
ggplot(summary_dt, aes(x = category, y = proportion, fill = category)) +
  geom_col(width = 0.8, colour = "black", linewidth = 0.2) +
  geom_text(aes(label = percent(proportion, accuracy = 0.1)), vjust = -0.35, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = cat_cols, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12))) +
  facet_wrap(~ subgenome, nrow = 1, labeller = labeller(subgenome = count_map_all)) +
  labs(x = NULL, y = "Percentage of sites") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), legend.position = "right", plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
dev.off()

summary_dt <- read.table(file="percent_met_chh.tsv",header=T)
summary_dt$category <- factor(summary_dt$category, levels = names(cat_cols))
count_map_all <- aggregate(N ~ subgenome, data = summary_dt, FUN = sum, na.rm = TRUE)
count_map_all$label <- paste0(count_map_all$subgenome, " (", format(count_map_all$N, big.mark = ","), ")")
count_map_all <- setNames(count_map_all$label, count_map_all$subgenome)

pdf(paste0("met_classify_chh.pdf"), width = 7, height = 3)
ggplot(summary_dt, aes(x = category, y = proportion, fill = category)) +
  geom_col(width = 0.8, colour = "black", linewidth = 0.2) +
  geom_text(aes(label = percent(proportion, accuracy = 0.1)), vjust = -0.35, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = cat_cols, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12))) +
  facet_wrap(~ subgenome, nrow = 1, labeller = labeller(subgenome = count_map_all)) +
  labs(x = NULL, y = "Percentage of sites") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), legend.position = "right", plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
dev.off()

## percent methylation
meth_cg <- read.table(file = "met_classify_cg.tsv", header = TRUE)
meth_chg <- read.table(file = "met_classify_chg.tsv", header = TRUE)
meth_chh <- read.table(file = "met_classify_chh.tsv", header = TRUE)

meth_cg$context <- "CG"
meth_chg$context <- "CHG"
meth_chh$context <- "CHH"

meth_summary <- rbind(meth_cg, meth_chg, meth_chh)

meth_summary$context <- factor(meth_summary$context, levels = c("CG", "CHG", "CHH"))
meth_summary$sample <- factor(meth_summary$sample, levels = c("CS", "CSxP", "P"))
meth_summary$subgenome <- factor(meth_summary$subgenome, levels = c("A", "B", "D"))

sample_cols <- c(CS="#0072B2",CSxP="#E69F00",P="#CC79A7")

pdf(paste0("percent_meth.pdf"), width = 5, height = 2.5)
ggplot(meth_summary, aes(x = subgenome, y = mean_pct, fill = sample, group = sample)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65, colour = "black", linewidth = 0.25) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(width = 0.7), width = 0.2, colour = "black", linewidth = 0.4) +
  scale_fill_manual(values = sample_cols) +
  facet_wrap(~ context, nrow = 1, scales="free_y") +
  labs(x = NULL, y = "Mean Methylation (%)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), legend.position = "right", plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "right")
dev.off()

## methylation coverage
coverage_cg <- read.table(file = "mean_coverage_cg.tsv", header = TRUE)
coverage_chg <- read.table(file = "mean_coverage_chg.tsv", header = TRUE)
coverage_chh <- read.table(file = "mean_coverage_chh.tsv", header = TRUE)

coverage_cg$context <- "CG"
coverage_chg$context <- "CHG"
coverage_chh$context <- "CHH"

coverage_summary <- rbind(coverage_cg, coverage_chg, coverage_chh)

coverage_summary$context <- factor(coverage_summary$context, levels = c("CG", "CHG", "CHH"))
coverage_summary$sample <- factor(coverage_summary$sample, levels = c("CS", "CSxP", "P"))
coverage_summary$subgenome <- factor(coverage_summary$subgenome, levels = c("A", "B", "D"))

pdf(paste0("coverage_meth.pdf"), width = 4.5, height = 2.5)
ggplot(coverage_summary, aes(x = subgenome, y = mean_cov, fill = sample, group = sample)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65, colour = "black", linewidth = 0.25) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(width = 0.7), width = 0.2, colour = "black", linewidth = 0.4) +
  scale_fill_manual(values = sample_cols) +
  facet_wrap(~ context, nrow = 1) +
  labs(x = NULL, y = "Coverage", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), legend.position = "right", plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "right")
dev.off()

## using CDS, plot what the different boman methylation categories mean 

library(data.table)
library(ggplot2)
library(cowplot)
library(ragg)

dt <- fread("merged_CG_symmetric_CDS.txt.gz")

## temp
#n <- min(10000, nrow(dt))
#idx <- sort(sample.int(nrow(dt), n))
#dt <- dt[idx]
## temp end

dt <- dt[cov_CS > 10 & cov_CSxP > 10 & cov_P > 10]

dt$subgenome <- sub("^chr[0-9]+([ABD])$", "\\1", dt$chr)
dt <- dt[subgenome %in% c("A", "B", "D")]
dt$subgenome <- factor(dt$subgenome, levels = c("A", "B", "D"))

dt$A <- dt$pct_CS / 100
dt$H <- dt$pct_CSxP / 100
dt$B <- dt$pct_P / 100

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
dt <- dt[max_mC >= 0.1]
dt$category <- ifelse(dt$radius < 0.1,"conserved_mC",vapply(dt$sector_center, sector_to_class, character(1)))
category_levels <- c("conserved_mC","additive","CS_dominant","P_dominant","overdominant","underdominant")
dt$category <- factor(dt$category, levels = category_levels)


cat_cols <- c(conserved_mC = "#1A1A1A", additive = "#D9D9D9", CS_dominant = "#0072B2", P_dominant = "#CC79A7", overdominant = "#E69F00", underdominant = "#1B7837")

category_counts <- dt %>% 
  count(category) %>% 
  mutate(percentage = n / sum(n) * 100) %>% 
  mutate(label = paste0(category, " (", round(percentage, 2), "%)"))

legend_labels <- setNames(category_counts$label, category_counts$category)

a_plot <- ggplot(dt, aes(x = x*100, y = y*100)) +
  geom_point(aes(color = category), size = 1, alpha = 1) +
  geom_abline(intercept = 0, slope = tan(c(1, 3, 5, 7) * pi / 8), color = "gray70", linetype = "dotted") +
  geom_hline(yintercept = 0, color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = cat_cols) + 
  labs(x = "Mean CSxP - P methylation %", y = "Mean CSxP - CS methylation %",color = "") +
  theme_minimal(base_size = 12) +
  theme( plot.title = element_text(face = "bold", hjust = 0.5),plot.subtitle = element_text(size = 10, color = "gray30", hjust = 0.5),legend.position = "right",panel.grid.minor = element_blank())

b_plot <- ggplot(dt, aes(x = (pct_CS-pct_P), y = pct_CSxP - ((pct_CS+pct_P)/2)) ) +
  geom_point(aes(color = category), size = 1, alpha = 1) +
  geom_hline(yintercept = 0, color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = cat_cols) + 
  labs(x = "Mean CS - P methylation %", y = "Mean CSxP - MPV methylation %",color = "") +
  theme_minimal(base_size = 12) +
  theme( plot.title = element_text(face = "bold", hjust = 0.5),plot.subtitle = element_text(size = 10, color = "gray30", hjust = 0.5),legend.position = "right",panel.grid.minor = element_blank())

ggsave(filename = "boman_class_geometric.png",plot = plot_grid(a_plot, b_plot, ncol=2),width = 10,height = 3,units = "in",dpi = 300,device = ragg::agg_png,bg = "white")

agg_png( filename = "boman_additive.png", width = 4, height = 4, units = "in", res = 300, background = "white" )
heatscatter(dt[dt$category %in% "additive",]$pct_CS, dt[dt$category %in% "additive",]$pct_P,xlab = "CS Methylation %",ylab = "Paragon Methylation %", main = "")
dev.off()

ci_summary <- function(x, by_cols, value_col) {
  out <- x[,list(N = sum(!is.na(get(value_col))),mean_value = mean(get(value_col), na.rm = TRUE),sd_value = sd(get(value_col), na.rm = TRUE)),by = by_cols]
  
  out$se <- out$sd_value / sqrt(out$N)
  out$ci <- qt(0.975, df = out$N - 1) * out$se
  out$ci_low <- out$mean_value - out$ci
  out$ci_high <- out$mean_value + out$ci
  
  out
}

meth_long_cat <- rbind(
  data.table(category = dt$category, sample = "CS", pct = dt$pct_CS),
  data.table(category = dt$category, sample = "CSxP", pct = dt$pct_CSxP),
  data.table(category = dt$category, sample = "P", pct = dt$pct_P)
)

meth_long_cat$category <- factor(meth_long_cat$category, levels = category_levels)
meth_long_cat$sample <- factor(meth_long_cat$sample, levels = c("CS", "CSxP", "P"))

meth_cat_summary <- ci_summary(meth_long_cat, c("category", "sample"), "pct")
setnames(meth_cat_summary, "mean_value", "mean_pct")

meth_cat_summary$category <- factor(meth_cat_summary$category, levels = category_levels)
meth_cat_summary$sample <- factor(meth_cat_summary$sample, levels = c("CS", "CSxP", "P"))

sample_cols <- c(CS="#0072B2",CSxP="#E69F00",P="#CC79A7")

p_meth_cat <- ggplot(meth_cat_summary, aes(x = sample, y = mean_pct, fill = sample)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, linewidth = 0.4) +
  facet_wrap(~category, ncol = 3, nrow = 2) +
  scale_fill_manual(values = sample_cols) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Mean methylation level (%)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),strip.background = element_rect(fill = "white", colour = "black"),legend.position = "none",plot.title = element_text(face = "bold"),axis.title = element_text(face = "bold"))

pdf(file = "mean_methylation_by_category_cg_cds.pdf", width = 4, height = 4)
print(p_meth_cat)
dev.off()


dt <- fread("merged_CHG_symmetric_CDS.txt.gz")
dt <- dt[cov_CS > 10 & cov_CSxP > 10 & cov_P > 10]
dt$subgenome <- sub("^chr[0-9]+([ABD])$", "\\1", dt$chr)
dt <- dt[subgenome %in% c("A", "B", "D")]
dt$A <- dt$pct_CS / 100
dt$H <- dt$pct_CSxP / 100
dt$B <- dt$pct_P / 100
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
dt$max_mC <- pmax(dt$A, dt$H, dt$B, na.rm=TRUE)
dt <- dt[max_mC >= 0.1]
dt$category <- ifelse(dt$radius < 0.1, "conserved_mC", vapply(dt$sector_center, sector_to_class, character(1)))
category_levels <- c("conserved_mC", "additive", "CS_dominant", "P_dominant", "overdominant", "underdominant")
dt$category <- factor(dt$category, levels=category_levels)
ci_summary <- function(x, by_cols, value_col) {
  out <- x[, list(N=sum(!is.na(get(value_col))), mean_value=mean(get(value_col), na.rm=TRUE), sd_value=sd(get(value_col), na.rm=TRUE)), by=by_cols]
  out$se <- out$sd_value / sqrt(out$N)
  out$ci <- qt(0.975, df=out$N - 1) * out$se
  out$ci_low <- out$mean_value - out$ci
  out$ci_high <- out$mean_value + out$ci
  out
}
meth_long_cat <- rbind(
  data.table(category=dt$category, sample="CS", pct=dt$pct_CS),
  data.table(category=dt$category, sample="CSxP", pct=dt$pct_CSxP),
  data.table(category=dt$category, sample="P", pct=dt$pct_P)
)
meth_long_cat$category <- factor(meth_long_cat$category, levels=category_levels)
meth_long_cat$sample <- factor(meth_long_cat$sample, levels=c("CS", "CSxP", "P"))
meth_cat_summary <- ci_summary(meth_long_cat, c("category", "sample"), "pct")
setnames(meth_cat_summary, "mean_value", "mean_pct")
meth_cat_summary$category <- factor(meth_cat_summary$category, levels=category_levels)
meth_cat_summary$sample <- factor(meth_cat_summary$sample, levels=c("CS", "CSxP", "P"))
sample_cols <- c(CS="#0072B2", CSxP="#E69F00", P="#CC79A7")
p_meth_cat <- ggplot(meth_cat_summary, aes(x=sample, y=mean_pct, fill=sample)) +
  geom_col(width=0.7, colour="black", linewidth=0.3) +
  geom_errorbar(aes(ymin=ci_low, ymax=ci_high), width=0.2, linewidth=0.4) +
  facet_wrap(~category, ncol=3, nrow=2) +
  scale_fill_manual(values=sample_cols) +
  scale_y_continuous(limits=c(0, 100), expand=expansion(mult=c(0, 0.05))) +
  labs(x=NULL, y="Mean methylation level (%)", fill=NULL) +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linewidth=0.2, colour="grey90"), strip.background=element_rect(fill="white", colour="black"), legend.position="none", plot.title=element_text(face="bold"), axis.title=element_text(face="bold"))
pdf(file="mean_methylation_by_category_chg_cds.pdf", width=4, height=4)
print(p_meth_cat)
dev.off()


dt <- fread("merged_CHH_all_CDS.txt.gz")
dt <- dt[cov_CS > 10 & cov_CSxP > 10 & cov_P > 10]
dt$subgenome <- sub("^chr[0-9]+([ABD])$", "\\1", dt$chr)
dt <- dt[subgenome %in% c("A", "B", "D")]
dt$A <- dt$pct_CS / 100
dt$H <- dt$pct_CSxP / 100
dt$B <- dt$pct_P / 100
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
dt$max_mC <- pmax(dt$A, dt$H, dt$B, na.rm=TRUE)
dt <- dt[max_mC >= 0.1]
dt$category <- ifelse(dt$radius < 0.1, "conserved_mC", vapply(dt$sector_center, sector_to_class, character(1)))
category_levels <- c("conserved_mC", "additive", "CS_dominant", "P_dominant", "overdominant", "underdominant")
dt$category <- factor(dt$category, levels=category_levels)
ci_summary <- function(x, by_cols, value_col) {
  out <- x[, list(N=sum(!is.na(get(value_col))), mean_value=mean(get(value_col), na.rm=TRUE), sd_value=sd(get(value_col), na.rm=TRUE)), by=by_cols]
  out$se <- out$sd_value / sqrt(out$N)
  out$ci <- qt(0.975, df=out$N - 1) * out$se
  out$ci_low <- out$mean_value - out$ci
  out$ci_high <- out$mean_value + out$ci
  out
}
meth_long_cat <- rbind(
  data.table(category=dt$category, sample="CS", pct=dt$pct_CS),
  data.table(category=dt$category, sample="CSxP", pct=dt$pct_CSxP),
  data.table(category=dt$category, sample="P", pct=dt$pct_P)
)
meth_long_cat$category <- factor(meth_long_cat$category, levels=category_levels)
meth_long_cat$sample <- factor(meth_long_cat$sample, levels=c("CS", "CSxP", "P"))
meth_cat_summary <- ci_summary(meth_long_cat, c("category", "sample"), "pct")
setnames(meth_cat_summary, "mean_value", "mean_pct")
meth_cat_summary$category <- factor(meth_cat_summary$category, levels=category_levels)
meth_cat_summary$sample <- factor(meth_cat_summary$sample, levels=c("CS", "CSxP", "P"))
sample_cols <- c(CS="#0072B2", CSxP="#E69F00", P="#CC79A7")
p_meth_cat <- ggplot(meth_cat_summary, aes(x=sample, y=mean_pct, fill=sample)) +
  geom_col(width=0.7, colour="black", linewidth=0.3) +
  geom_errorbar(aes(ymin=ci_low, ymax=ci_high), width=0.2, linewidth=0.4) +
  facet_wrap(~category, ncol=3, nrow=2) +
  scale_fill_manual(values=sample_cols) +
  scale_y_continuous(limits=c(0, 100), expand=expansion(mult=c(0, 0.05))) +
  labs(x=NULL, y="Mean methylation level (%)", fill=NULL) +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linewidth=0.2, colour="grey90"), strip.background=element_rect(fill="white", colour="black"), legend.position="none", plot.title=element_text(face="bold"), axis.title=element_text(face="bold"))
pdf(file="mean_methylation_by_category_chh_cds.pdf", width=4, height=4)
print(p_meth_cat)
dev.off()


library(data.table)
library(ggplot2)
library(cowplot)
library(ragg)

dt <- fread("merged_CHG_symmetric_CDS.txt.gz")

dt <- dt[cov_CS > 10 & cov_CSxP > 10 & cov_P > 10]

dt$subgenome <- sub("^chr[0-9]+([ABD])$", "\\1", dt$chr)
dt <- dt[subgenome %in% c("A", "B", "D")]
dt$subgenome <- factor(dt$subgenome, levels = c("A", "B", "D"))

dt$A <- dt$pct_CS / 100
dt$H <- dt$pct_CSxP / 100
dt$B <- dt$pct_P / 100

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
dt <- dt[max_mC >= 0.1]
dt$category <- ifelse(dt$radius < 0.1,"conserved_mC",vapply(dt$sector_center, sector_to_class, character(1)))
category_levels <- c("conserved_mC","additive","CS_dominant","P_dominant","overdominant","underdominant")
dt$category <- factor(dt$category, levels = category_levels)


cat_cols <- c(conserved_mC = "#1A1A1A", additive = "#D9D9D9", CS_dominant = "#0072B2", P_dominant = "#CC79A7", overdominant = "#E69F00", underdominant = "#1B7837")

category_counts <- dt %>% 
  count(category) %>% 
  mutate(percentage = n / sum(n) * 100) %>% 
  mutate(label = paste0(category, " (", round(percentage, 2), "%)"))

legend_labels <- setNames(category_counts$label, category_counts$category)

a_plot <- ggplot(dt, aes(x = x*100, y = y*100)) +
  geom_point(aes(color = category), size = 1, alpha = 1) +
  geom_abline(intercept = 0, slope = tan(c(1, 3, 5, 7) * pi / 8), color = "gray70", linetype = "dotted") +
  geom_hline(yintercept = 0, color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = cat_cols) + 
  labs(x = "Mean CSxP - P methylation %", y = "Mean CSxP - CS methylation %",color = "") +
  theme_minimal(base_size = 12) +
  theme( plot.title = element_text(face = "bold", hjust = 0.5),plot.subtitle = element_text(size = 10, color = "gray30", hjust = 0.5),legend.position = "right",panel.grid.minor = element_blank())

b_plot <- ggplot(dt, aes(x = (pct_CS-pct_P), y = pct_CSxP - ((pct_CS+pct_P)/2)) ) +
  geom_point(aes(color = category), size = 1, alpha = 1) +
  geom_hline(yintercept = 0, color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = cat_cols) + 
  labs(x = "Mean CS - P methylation %", y = "Mean CSxP - MPV methylation %",color = "") +
  theme_minimal(base_size = 12) +
  theme( plot.title = element_text(face = "bold", hjust = 0.5),plot.subtitle = element_text(size = 10, color = "gray30", hjust = 0.5),legend.position = "right",panel.grid.minor = element_blank())

ggsave(filename = "boman_class_geometric_chg.png",plot = plot_grid(a_plot, b_plot, ncol=2),width = 10,height = 3,units = "in",dpi = 300,device = ragg::agg_png,bg = "white")



dt <- fread("merged_CHH_all_CDS.txt.gz")
dt <- dt[cov_CS > 10 & cov_CSxP > 10 & cov_P > 10]

dt$subgenome <- sub("^chr[0-9]+([ABD])$", "\\1", dt$chr)
dt <- dt[subgenome %in% c("A", "B", "D")]
dt$subgenome <- factor(dt$subgenome, levels = c("A", "B", "D"))

dt$A <- dt$pct_CS / 100
dt$H <- dt$pct_CSxP / 100
dt$B <- dt$pct_P / 100

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
dt <- dt[max_mC >= 0.1]
dt$category <- ifelse(dt$radius < 0.1,"conserved_mC",vapply(dt$sector_center, sector_to_class, character(1)))
category_levels <- c("conserved_mC","additive","CS_dominant","P_dominant","overdominant","underdominant")
dt$category <- factor(dt$category, levels = category_levels)


cat_cols <- c(conserved_mC = "#1A1A1A", additive = "#D9D9D9", CS_dominant = "#0072B2", P_dominant = "#CC79A7", overdominant = "#E69F00", underdominant = "#1B7837")

category_counts <- dt %>% 
  count(category) %>% 
  mutate(percentage = n / sum(n) * 100) %>% 
  mutate(label = paste0(category, " (", round(percentage, 2), "%)"))

legend_labels <- setNames(category_counts$label, category_counts$category)

a_plot <- ggplot(dt, aes(x = x*100, y = y*100)) +
  geom_point(aes(color = category), size = 1, alpha = 1) +
  geom_abline(intercept = 0, slope = tan(c(1, 3, 5, 7) * pi / 8), color = "gray70", linetype = "dotted") +
  geom_hline(yintercept = 0, color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = cat_cols) + 
  labs(x = "Mean CSxP - P methylation %", y = "Mean CSxP - CS methylation %",color = "") +
  theme_minimal(base_size = 12) +
  theme( plot.title = element_text(face = "bold", hjust = 0.5),plot.subtitle = element_text(size = 10, color = "gray30", hjust = 0.5),legend.position = "right",panel.grid.minor = element_blank())

b_plot <- ggplot(dt, aes(x = (pct_CS-pct_P), y = pct_CSxP - ((pct_CS+pct_P)/2)) ) +
  geom_point(aes(color = category), size = 1, alpha = 1) +
  geom_hline(yintercept = 0, color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = cat_cols) + 
  labs(x = "Mean CS - P methylation %", y = "Mean CSxP - MPV methylation %",color = "") +
  theme_minimal(base_size = 12) +
  theme( plot.title = element_text(face = "bold", hjust = 0.5),plot.subtitle = element_text(size = 10, color = "gray30", hjust = 0.5),legend.position = "right",panel.grid.minor = element_blank())

ggsave(filename = "boman_class_geometric_chh.png",plot = plot_grid(a_plot, b_plot, ncol=2),width = 10,height = 3,units = "in",dpi = 300,device = ragg::agg_png,bg = "white")

library(data.table)
library(ggplot2)
library(cowplot)
library(ragg)

dt <- fread("merged_CHG_symmetric_CDS.txt.gz")

dt <- dt[cov_CS > 10 & cov_CSxP > 10 & cov_P > 10]

dt$subgenome <- sub("^chr[0-9]+([ABD])$", "\\1", dt$chr)
dt <- dt[subgenome %in% c("A", "B", "D")]
dt$subgenome <- factor(dt$subgenome, levels = c("A", "B", "D"))

dt$A <- dt$pct_CS / 100
dt$H <- dt$pct_CSxP / 100
dt$B <- dt$pct_P / 100

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
dt <- dt[max_mC >= 0.1]
dt$category <- ifelse(dt$radius < 0.1,"conserved_mC",vapply(dt$sector_center, sector_to_class, character(1)))
category_levels <- c("conserved_mC","additive","CS_dominant","P_dominant","overdominant","underdominant")
dt$category <- factor(dt$category, levels = category_levels)


cat_cols <- c(conserved_mC = "#1A1A1A", additive = "#D9D9D9", CS_dominant = "#0072B2", P_dominant = "#CC79A7", overdominant = "#E69F00", underdominant = "#1B7837")

category_counts <- dt %>% 
  count(category) %>% 
  mutate(percentage = n / sum(n) * 100) %>% 
  mutate(label = paste0(category, " (", round(percentage, 2), "%)"))

legend_labels <- setNames(category_counts$label, category_counts$category)

a_plot <- ggplot(dt, aes(x = x*100, y = y*100)) +
  geom_point(aes(color = category), size = 1, alpha = 1) +
  geom_abline(intercept = 0, slope = tan(c(1, 3, 5, 7) * pi / 8), color = "gray70", linetype = "dotted") +
  geom_hline(yintercept = 0, color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = cat_cols) + 
  labs(x = "Mean CSxP - P methylation %", y = "Mean CSxP - CS methylation %",color = "") +
  theme_minimal(base_size = 12) +
  theme( plot.title = element_text(face = "bold", hjust = 0.5),plot.subtitle = element_text(size = 10, color = "gray30", hjust = 0.5),legend.position = "right",panel.grid.minor = element_blank())

b_plot <- ggplot(dt, aes(x = (pct_CS-pct_P), y = pct_CSxP - ((pct_CS+pct_P)/2)) ) +
  geom_point(aes(color = category), size = 1, alpha = 1) +
  geom_hline(yintercept = 0, color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = cat_cols) + 
  labs(x = "Mean CS - P methylation %", y = "Mean CSxP - MPV methylation %",color = "") +
  theme_minimal(base_size = 12) +
  theme( plot.title = element_text(face = "bold", hjust = 0.5),plot.subtitle = element_text(size = 10, color = "gray30", hjust = 0.5),legend.position = "right",panel.grid.minor = element_blank())

ggsave(filename = "boman_class_geometric_chg.png",plot = plot_grid(a_plot, b_plot, ncol=2),width = 10,height = 3,units = "in",dpi = 300,device = ragg::agg_png,bg = "white")



dt <- fread("merged_CHH_all_CDS.txt.gz")
dt <- dt[cov_CS > 10 & cov_CSxP > 10 & cov_P > 10]

dt$subgenome <- sub("^chr[0-9]+([ABD])$", "\\1", dt$chr)
dt <- dt[subgenome %in% c("A", "B", "D")]
dt$subgenome <- factor(dt$subgenome, levels = c("A", "B", "D"))

dt$A <- dt$pct_CS / 100
dt$H <- dt$pct_CSxP / 100
dt$B <- dt$pct_P / 100

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
dt <- dt[max_mC >= 0.1]
dt$category <- ifelse(dt$radius < 0.1,"conserved_mC",vapply(dt$sector_center, sector_to_class, character(1)))
category_levels <- c("conserved_mC","additive","CS_dominant","P_dominant","overdominant","underdominant")
dt$category <- factor(dt$category, levels = category_levels)


cat_cols <- c(conserved_mC = "#1A1A1A", additive = "#D9D9D9", CS_dominant = "#0072B2", P_dominant = "#CC79A7", overdominant = "#E69F00", underdominant = "#1B7837")

category_counts <- dt %>% 
  count(category) %>% 
  mutate(percentage = n / sum(n) * 100) %>% 
  mutate(label = paste0(category, " (", round(percentage, 2), "%)"))

legend_labels <- setNames(category_counts$label, category_counts$category)

a_plot <- ggplot(dt, aes(x = x*100, y = y*100)) +
  geom_point(aes(color = category), size = 1, alpha = 1) +
  geom_abline(intercept = 0, slope = tan(c(1, 3, 5, 7) * pi / 8), color = "gray70", linetype = "dotted") +
  geom_hline(yintercept = 0, color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = cat_cols) + 
  labs(x = "Mean CSxP - P methylation %", y = "Mean CSxP - CS methylation %",color = "") +
  theme_minimal(base_size = 12) +
  theme( plot.title = element_text(face = "bold", hjust = 0.5),plot.subtitle = element_text(size = 10, color = "gray30", hjust = 0.5),legend.position = "right",panel.grid.minor = element_blank())

b_plot <- ggplot(dt, aes(x = (pct_CS-pct_P), y = pct_CSxP - ((pct_CS+pct_P)/2)) ) +
  geom_point(aes(color = category), size = 1, alpha = 1) +
  geom_hline(yintercept = 0, color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = cat_cols) + 
  labs(x = "Mean CS - P methylation %", y = "Mean CSxP - MPV methylation %",color = "") +
  theme_minimal(base_size = 12) +
  theme( plot.title = element_text(face = "bold", hjust = 0.5),plot.subtitle = element_text(size = 10, color = "gray30", hjust = 0.5),legend.position = "right",panel.grid.minor = element_blank())

ggsave(filename = "boman_class_geometric_chh.png",plot = plot_grid(a_plot, b_plot, ncol=2),width = 10,height = 3,units = "in",dpi = 300,device = ragg::agg_png,bg = "white")

## expression level by cytosine context

library(data.table)
library(ggplot2)

cs_tpm <- read.table(file = "cs_tpm.tsv")
cs_tpm <- as.data.table(cs_tpm, keep.rownames = "gene_id")
cs_tpm <- cs_tpm[!grepl("LC$", gene_id)]
colnames(cs_tpm) <- c("gene_id", sub("^PXCS", "PxCS", sub("_.*", "", colnames(cs_tpm)[-1]))); sample_cols <- colnames(cs_tpm)[-1]
cs_tpm <- cs_tpm[rowSums(as.data.frame(cs_tpm)[, sample_cols] >= 1) >= 3]
cs_tpm <- cs_tpm[, 1:10]

tpm_gene <- data.table(
  gene_id = cs_tpm$gene_id,
  CS = rowMeans(cs_tpm[, c("CS1", "CS2", "CS3"), with = FALSE], na.rm = TRUE),
  CSxP = rowMeans(cs_tpm[, c("CSxP1", "CSxP2", "CSxP3"), with = FALSE], na.rm = TRUE),
  P = rowMeans(cs_tpm[, c("P1", "P2", "P3"), with = FALSE], na.rm = TRUE)
)

tpm_long <- melt(tpm_gene, id.vars = "gene_id", variable.name = "sample", value.name = "TPM")

process_context <- function(infile, context_name, remove_col3 = FALSE) {
  dt <- fread(infile)
  if (remove_col3) {
    dt <- dt[, -3, with = FALSE]
  }
  dt <- dt[!is.na(gene_id)]
  
  dt$meth_CS <- round((dt$pct_CS / 100) * dt$cov_CS)
  dt$meth_CSxP <- round((dt$pct_CSxP / 100) * dt$cov_CSxP)
  dt$meth_P <- round((dt$pct_P / 100) * dt$cov_P)
  dt$site_CS <- dt$cov_CS > 0
  dt$site_CSxP <- dt$cov_CSxP > 0
  dt$site_P <- dt$cov_P > 0
  meth_gene <- aggregate(cbind(meth_CS, cov_CS, site_CS, meth_CSxP, cov_CSxP, site_CSxP, meth_P, cov_P, site_P) ~ gene_id, data = dt, FUN = sum, na.rm = TRUE)
  meth_gene <- meth_gene[meth_gene$cov_CS >= 10 & meth_gene$cov_CSxP >= 10 & meth_gene$cov_P >= 10 & meth_gene$site_CS > 5 & meth_gene$site_CSxP > 5 & meth_gene$site_P > 5, ]
  meth_gene$pct_CS <- 100 * meth_gene$meth_CS / meth_gene$cov_CS
  meth_gene$pct_CSxP <- 100 * meth_gene$meth_CSxP / meth_gene$cov_CSxP
  meth_gene$pct_P <- 100 * meth_gene$meth_P / meth_gene$cov_P
  meth_gene <- meth_gene[, c("gene_id", "pct_CS", "pct_CSxP", "pct_P")]
  meth_gene <- as.data.table(meth_gene)
  meth_long <- melt(meth_gene, id.vars = "gene_id", variable.name = "sample", value.name = "methylation_pct")
  meth_long$sample <- sub("^pct_", "", meth_long$sample)
  meth_long$context <- context_name
  return(meth_long)
}

cg_long <- process_context("merged_CG_symmetric_CDS.txt.gz", "CG", remove_col3 = FALSE)
chg_long <- process_context("merged_CHG_symmetric_CDS.txt.gz", "CHG", remove_col3 = FALSE)
chh_long <- process_context("merged_CHH_all_CDS.txt.gz", "CHH", remove_col3 = TRUE)

meth_long <- rbind(cg_long, chg_long, chh_long)
write.table(meth_long,file="CDS_meth_pct.txt",row.names=F,sep="\t")

meth_long <- read.table(file="CDS_meth_pct.txt",header=T)

plot_dt <- merge(meth_long, tpm_long, by = c("gene_id", "sample"))
plot_dt$log2_TPM <- log2(plot_dt$TPM + 1)
plot_dt$context <- factor(plot_dt$context, levels = c("CG", "CHG", "CHH"))
cor_results <- data.table(context = c("CG", "CHG", "CHH"), r = NA_real_, p = NA_real_)
for (i in 1:nrow(cor_results)) {
  x <- plot_dt[plot_dt$context == cor_results$context[i], ]
  test <- cor.test(x$TPM, x$methylation_pct, method = "pearson")
  cor_results$r[i] <- unname(test$estimate)
  cor_results$p[i] <- test$p.value
}
cor_results$label <- paste0(cor_results$context, ": r = ", sprintf("%.3f", cor_results$r), ", P = ", format.pval(cor_results$p, digits = 2, eps = 1e-300))

pdf(file="tpm_meth_cds.pdf",height=2.5,width=2.5)
ggplot(plot_dt, aes(x = log2_TPM, y = methylation_pct, colour = context)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  geom_text(data = cor_results, aes(x = 13, y = c(10,15,20), label = label, colour = context), inherit.aes = FALSE, hjust = 1, vjust = 1, size = 3) +
  scale_colour_manual(values = c(CG = "#D55E00",CHG = "#56B4E9",CHH = "#999999")) +
  labs(x = "log2(TPM)", y = "Mean Methylation (%)", colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), axis.title = element_text(face = "bold"))
dev.off()

meth_long <- read.table(file="CDS_meth_pct.txt",header=T)


homologies <- read.csv(file="homologies.csv")
homologies <- homologies[1:4]
homologies <- as.data.table(homologies)
hom_long <- melt(homologies, id.vars = "group_id", measure.vars = c("A","B","D"), variable.name = "homoeolog", value.name = "gene_id")
meth_long <- meth_long[meth_long$gene_id %in% hom_long$gene_id,]
meth_long <- left_join(meth_long,hom_long,by="gene_id")
colnames(meth_long)[2] <- "genotype"

bias_categories <- read.csv(file="bias_category_all_samples_inc_orig_expr.csv")
bias_categories$sample <- gsub("_.*","",bias_categories$sample)
bias_categories <- bias_categories[!bias_categories$sample %in% "PxCS3",]
bias_categories$sample <- gsub("PXCS2","PxCS2",bias_categories$sample)
bias_categories <- bias_categories[!bias_categories$sample %in% c("PxCS1", "PxCS2"),]
bias_categories <- bias_categories %>%
  mutate(genotype = sub("[0-9]+$","",as.character(sample))) %>%
  group_by(group_id,genotype) %>%
  summarise(A_tpm = mean(A_tpm,na.rm = TRUE),B_tpm = mean(B_tpm,na.rm = TRUE),D_tpm = mean(D_tpm,na.rm = TRUE),.groups = "drop") %>%
  mutate(triad_tpm = A_tpm + B_tpm + D_tpm,A = A_tpm / triad_tpm,B = B_tpm / triad_tpm,D = D_tpm / triad_tpm)
bias_categories$CV <- apply(bias_categories[7:9],1,sd)/apply(bias_categories[7:9],1,mean)
bias_categories$group_id <- as.numeric(gsub("X","",bias_categories$group_id))
bias_categories_part <- bias_categories[c(1,2,6,10)]

joined_df <- full_join(bias_categories_part,meth_long,by = c("group_id","genotype"))
joined_df <- joined_df[complete.cases(joined_df),]

cor_df <- joined_df %>%
  group_by(context) %>%
  summarise(r = cor(methylation_pct,CV,use = "complete.obs",method = "pearson"),p = cor.test(methylation_pct,CV,method = "pearson")$p.value,.groups = "drop") %>%
  mutate(label = paste0(context,": R = ",round(r,2),", ",ifelse(p < 0.001,"P < 0.001",paste0("P = ",signif(p,2)))))

pdf(file="heb_cds_meth.pdf",height=2.5,width=2.5)
ggplot(data = joined_df,aes(y = methylation_pct,x = CV,color = context)) +
  geom_smooth(method = "lm",se = T,linewidth = 1) +
  geom_text(data = cor_df,aes(y = c(30,35,40),x = 1.8,label = label,color = context),hjust = 1.05,vjust = 1.2,size = 2,inherit.aes = FALSE,show.legend = FALSE) +
  scale_colour_manual(values = c(CG = "#D55E00",CHG = "#56B4E9",CHH = "#999999")) +
  labs(y = "Mean Methylation (%)",x = "HEB",color=NULL) +
  theme_bw(base_size = 12) +
  theme(axis.title = element_text(face = "bold"),panel.grid.minor = element_blank(),legend.position = "none")
dev.off()

cg_long <- process_context("merged_CG_symmetric_promoter1kb.txt.gz", "CG", remove_col3 = FALSE)
chg_long <- process_context("merged_CHG_symmetric_promoter1kb.txt.gz", "CHG", remove_col3 = FALSE)
chh_long <- process_context("merged_CHH_promoter1kb.txt.gz", "CHH", remove_col3 = TRUE)

meth_long <- rbind(cg_long, chg_long, chh_long)
write.table(meth_long,file="promoter_meth_pct.txt",row.names=F,sep="\t")

meth_long <- read.table(file="promoter_meth_pct.txt",header=T)

plot_dt <- merge(meth_long, tpm_long, by = c("gene_id", "sample"))
plot_dt$log2_TPM <- log2(plot_dt$TPM + 1)
plot_dt$context <- factor(plot_dt$context, levels = c("CG", "CHG", "CHH"))
cor_results <- data.table(context = c("CG", "CHG", "CHH"), r = NA_real_, p = NA_real_)
for (i in 1:nrow(cor_results)) {
  x <- plot_dt[plot_dt$context == cor_results$context[i], ]
  test <- cor.test(x$TPM, x$methylation_pct, method = "pearson")
  cor_results$r[i] <- unname(test$estimate)
  cor_results$p[i] <- test$p.value
}
cor_results$label <- paste0(cor_results$context, ": r = ", sprintf("%.3f", cor_results$r), ", P = ", format.pval(cor_results$p, digits = 2, eps = 1e-300))

pdf(file="tpm_meth_promoter.pdf",height=2.5,width=2.5)
ggplot(plot_dt, aes(x = log2_TPM, y = methylation_pct, colour = context)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  geom_text(data = cor_results, aes(x = 13, y = c(8,13,18), label = label, colour = context), inherit.aes = FALSE, hjust = 1, vjust = 1, size = 3) +
  scale_colour_manual(values = c(CG = "#D55E00",CHG = "#56B4E9",CHH = "#999999")) +
  labs(x = "log2(TPM)", y = "Mean Methylation (%)", colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), axis.title = element_text(face = "bold"))
dev.off()

meth_long <- read.table(file="promoter_meth_pct.txt",header=T)

homologies <- read.csv(file="homologies.csv")
homologies <- homologies[1:4]
homologies <- as.data.table(homologies)
hom_long <- melt(homologies, id.vars = "group_id", measure.vars = c("A","B","D"), variable.name = "homoeolog", value.name = "gene_id")
meth_long <- meth_long[meth_long$gene_id %in% hom_long$gene_id,]
meth_long <- left_join(meth_long,hom_long,by="gene_id")
colnames(meth_long)[2] <- "genotype"

bias_categories <- read.csv(file="bias_category_all_samples_inc_orig_expr.csv")
bias_categories$sample <- gsub("_.*","",bias_categories$sample)
bias_categories <- bias_categories[!bias_categories$sample %in% "PxCS3",]
bias_categories$sample <- gsub("PXCS2","PxCS2",bias_categories$sample)
bias_categories <- bias_categories[!bias_categories$sample %in% c("PxCS1", "PxCS2"),]
bias_categories <- bias_categories %>%
  mutate(genotype = sub("[0-9]+$","",as.character(sample))) %>%
  group_by(group_id,genotype) %>%
  summarise(A_tpm = mean(A_tpm,na.rm = TRUE),B_tpm = mean(B_tpm,na.rm = TRUE),D_tpm = mean(D_tpm,na.rm = TRUE),.groups = "drop") %>%
  mutate(triad_tpm = A_tpm + B_tpm + D_tpm,A = A_tpm / triad_tpm,B = B_tpm / triad_tpm,D = D_tpm / triad_tpm)
bias_categories$CV <- apply(bias_categories[7:9],1,sd)/apply(bias_categories[7:9],1,mean)
bias_categories$group_id <- as.numeric(gsub("X","",bias_categories$group_id))
bias_categories_part <- bias_categories[c(1,2,6,10)]

joined_df <- full_join(bias_categories_part,meth_long,by = c("group_id","genotype"))
joined_df <- joined_df[complete.cases(joined_df),]

cor_df <- joined_df %>%
  group_by(context) %>%
  summarise(r = cor(methylation_pct,CV,use = "complete.obs",method = "pearson"),p = cor.test(methylation_pct,CV,method = "pearson")$p.value,.groups = "drop") %>%
  mutate(label = paste0(context,": R = ",round(r,2),", ",ifelse(p < 0.001,"P < 0.001",paste0("P = ",signif(p,2)))))

pdf(file="heb_promoter_meth.pdf",height=2.5,width=2.5)
ggplot(data = joined_df,aes(y = methylation_pct,x = CV,color = context)) +
  geom_smooth(method = "lm",se = T,linewidth = 1) +
  geom_text(data = cor_df,aes(y = c(10,15,20),x = 1,label = label,color = context),hjust = 1.05,vjust = 1.2,size = 2,inherit.aes = FALSE,show.legend = FALSE) +
  scale_colour_manual(values = c(CG = "#D55E00",CHG = "#56B4E9",CHH = "#999999")) +
  labs(y = "Mean Methylation (%)",x = "HEB",color=NULL) +
  theme_bw(base_size = 12) +
  theme(axis.title = element_text(face = "bold"),panel.grid.minor = element_blank(),legend.position = "none")
dev.off()

## compare methylation distance from parents

e3 <- function(a1,b1,d1, a2,b2,d2) sqrt((a1 - a2)^2 + (b1 - b2)^2 + (d1 - d2)^2)

meth_bias_df <- function(infile, context_label, remove_col3 = FALSE) {
  dt <- fread(infile)

  if (remove_col3) {
    dt <- dt[, -3, with = FALSE]
  }
  
  dt$meth_CS <- round((dt$pct_CS / 100) * dt$cov_CS)
  dt$meth_CSxP <- round((dt$pct_CSxP / 100) * dt$cov_CSxP)
  dt$meth_P <- round((dt$pct_P / 100) * dt$cov_P)
  dt$site_CS <- dt$cov_CS > 0
  dt$site_CSxP <- dt$cov_CSxP > 0
  dt$site_P <- dt$cov_P > 0
  gene_counts <- aggregate(cbind(meth_CS, cov_CS, site_CS, meth_CSxP, cov_CSxP, site_CSxP, meth_P, cov_P, site_P) ~ gene_id, data = dt, FUN = sum, na.rm = TRUE)
  gene_counts <- gene_counts[gene_counts$cov_CS >= 10 & gene_counts$cov_CSxP >= 10 & gene_counts$cov_P >= 10 & gene_counts$site_CS > 5 & gene_counts$site_CSxP > 5 & gene_counts$site_P > 5, ]
  gene_counts$pct_CS <- gene_counts$meth_CS / gene_counts$cov_CS
  gene_counts$pct_CSxP <- gene_counts$meth_CSxP / gene_counts$cov_CSxP
  gene_counts$pct_P <- gene_counts$meth_P / gene_counts$cov_P
  
  homologies <- read.csv(file="homologies.csv")
  homologies$A_pct_CS <- gene_counts$pct_CS[match(homologies$A, gene_counts$gene_id)]
  homologies$B_pct_CS <- gene_counts$pct_CS[match(homologies$B, gene_counts$gene_id)]
  homologies$D_pct_CS <- gene_counts$pct_CS[match(homologies$D, gene_counts$gene_id)]
  homologies$A_pct_P <- gene_counts$pct_P[match(homologies$A, gene_counts$gene_id)]
  homologies$B_pct_P <- gene_counts$pct_P[match(homologies$B, gene_counts$gene_id)]
  homologies$D_pct_P <- gene_counts$pct_P[match(homologies$D, gene_counts$gene_id)]
  homologies$A_pct_CSxP <- gene_counts$pct_CSxP[match(homologies$A, gene_counts$gene_id)]
  homologies$B_pct_CSxP <- gene_counts$pct_CSxP[match(homologies$B, gene_counts$gene_id)]
  homologies$D_pct_CSxP <- gene_counts$pct_CSxP[match(homologies$D, gene_counts$gene_id)]
  homologies <- homologies[complete.cases(homologies[, c("A_pct_CS", "B_pct_CS", "D_pct_CS","A_pct_P", "B_pct_P", "D_pct_P","A_pct_CSxP", "B_pct_CSxP", "D_pct_CSxP")]), ]
  homologies$meth_CV_CS <- apply(homologies[, c("A_pct_CS", "B_pct_CS", "D_pct_CS")], 1, sd, na.rm = TRUE) / apply(homologies[, c("A_pct_CS", "B_pct_CS", "D_pct_CS")], 1, mean, na.rm = TRUE)
  homologies$meth_CV_P <- apply(homologies[, c("A_pct_P", "B_pct_P", "D_pct_P")], 1, sd, na.rm = TRUE) / apply(homologies[, c("A_pct_P", "B_pct_P", "D_pct_P")], 1, mean, na.rm = TRUE)
  homologies$meth_CV_CSxP <- apply(homologies[, c("A_pct_CSxP", "B_pct_CSxP", "D_pct_CSxP")], 1, sd, na.rm = TRUE) / apply(homologies[, c("A_pct_CSxP", "B_pct_CSxP", "D_pct_CSxP")], 1, mean, na.rm = TRUE)
  homologies <- homologies %>%
    mutate(meth_dist_CS_P = e3(A_pct_CS, B_pct_CS, D_pct_CS, A_pct_P, B_pct_P, D_pct_P))
  homologies <- homologies %>%
    mutate(meth_dist_CS_hybrid = e3(A_pct_CS, B_pct_CS, D_pct_CS, A_pct_CSxP, B_pct_CSxP, D_pct_CSxP))
  homologies <- homologies %>%
    mutate(meth_dist_P_hybrid = e3(A_pct_P, B_pct_P, D_pct_P, A_pct_CSxP, B_pct_CSxP, D_pct_CSxP))
  
  homologies
}

meth_bias_cg <- meth_bias_df("merged_CG_symmetric_CDS.txt.gz", "CG")
meth_bias_cg$context <- "CG"
meth_bias_chg <- meth_bias_df("merged_CHG_symmetric_CDS.txt.gz", "CHG")
meth_bias_chg$context <- "CHG"
meth_bias_chh <- meth_bias_df("merged_CHH_all_CDS.txt.gz", "CHH", remove_col3 = TRUE)
meth_bias_chh$context <- "CHH"
meth_bias_all <- rbind(meth_bias_cg, meth_bias_chg, meth_bias_chh)
meth_bias_all$context <- factor(meth_bias_all$context, levels = c("CG", "CHG", "CHH"))

cor_labs <- meth_bias_all %>%
  filter(!is.na(meth_dist_CS_hybrid), !is.na(meth_dist_P_hybrid)) %>%
  group_by(context) %>%
  summarise(r = cor(meth_dist_CS_hybrid, meth_dist_P_hybrid),p = cor.test(meth_dist_CS_hybrid, meth_dist_P_hybrid)$p.value,.groups = "drop") %>%
  mutate(label = paste0("R = ", sprintf("%.2f", r), "\nP = ", format.pval(p, digits = 2, eps = 1e-3)))

pdf(file="methylation_dist_parents_cds.pdf",height=2,width=4.8)
ggplot(data = meth_bias_all, aes(x = meth_dist_CS_hybrid, y = meth_dist_P_hybrid, color = context)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey50") +
  geom_point(shape = 1, size = 1.2, stroke = 0.4) +
  geom_text(data = cor_labs, aes(x = 0.1, y = 1.2, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3) +
  facet_grid(~context) +
  scale_colour_manual(values = c(CG = "#D55E00", CHG = "#56B4E9", CHH = "#999999")) +
  labs(x = "Methylation Bias Distance to Chinese Spring", y = "Methylation Bias\nDistance to Paragon") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), plot.title = element_text(face = "bold"), legend.position = "none")
dev.off()

e3 <- function(a1,b1,d1, a2,b2,d2) sqrt((a1 - a2)^2 + (b1 - b2)^2 + (d1 - d2)^2)

meth_bias_df <- function(infile, context_label, remove_col3 = FALSE) {
  dt <- fread(infile)
  
  if (remove_col3) {
    dt <- dt[, -3, with = FALSE]
  }
  
  dt$meth_CS <- round((dt$pct_CS / 100) * dt$cov_CS)
  dt$meth_CSxP <- round((dt$pct_CSxP / 100) * dt$cov_CSxP)
  dt$meth_P <- round((dt$pct_P / 100) * dt$cov_P)
  dt$site_CS <- dt$cov_CS > 0
  dt$site_CSxP <- dt$cov_CSxP > 0
  dt$site_P <- dt$cov_P > 0
  gene_counts <- aggregate(cbind(meth_CS, cov_CS, site_CS, meth_CSxP, cov_CSxP, site_CSxP, meth_P, cov_P, site_P) ~ gene_id, data = dt, FUN = sum, na.rm = TRUE)
  gene_counts <- gene_counts[gene_counts$cov_CS >= 10 & gene_counts$cov_CSxP >= 10 & gene_counts$cov_P >= 10 & gene_counts$site_CS > 5 & gene_counts$site_CSxP > 5 & gene_counts$site_P > 5, ]
  gene_counts$pct_CS <- gene_counts$meth_CS / gene_counts$cov_CS
  gene_counts$pct_CSxP <- gene_counts$meth_CSxP / gene_counts$cov_CSxP
  gene_counts$pct_P <- gene_counts$meth_P / gene_counts$cov_P
  
  homologies <- read.csv(file="homologies.csv")
  homologies$A_pct_CS <- gene_counts$pct_CS[match(homologies$A, gene_counts$gene_id)]
  homologies$B_pct_CS <- gene_counts$pct_CS[match(homologies$B, gene_counts$gene_id)]
  homologies$D_pct_CS <- gene_counts$pct_CS[match(homologies$D, gene_counts$gene_id)]
  homologies$A_pct_P <- gene_counts$pct_P[match(homologies$A, gene_counts$gene_id)]
  homologies$B_pct_P <- gene_counts$pct_P[match(homologies$B, gene_counts$gene_id)]
  homologies$D_pct_P <- gene_counts$pct_P[match(homologies$D, gene_counts$gene_id)]
  homologies$A_pct_CSxP <- gene_counts$pct_CSxP[match(homologies$A, gene_counts$gene_id)]
  homologies$B_pct_CSxP <- gene_counts$pct_CSxP[match(homologies$B, gene_counts$gene_id)]
  homologies$D_pct_CSxP <- gene_counts$pct_CSxP[match(homologies$D, gene_counts$gene_id)]
  homologies <- homologies[complete.cases(homologies[, c("A_pct_CS", "B_pct_CS", "D_pct_CS","A_pct_P", "B_pct_P", "D_pct_P","A_pct_CSxP", "B_pct_CSxP", "D_pct_CSxP")]), ]
  homologies$meth_CV_CS <- apply(homologies[, c("A_pct_CS", "B_pct_CS", "D_pct_CS")], 1, sd, na.rm = TRUE) / apply(homologies[, c("A_pct_CS", "B_pct_CS", "D_pct_CS")], 1, mean, na.rm = TRUE)
  homologies$meth_CV_P <- apply(homologies[, c("A_pct_P", "B_pct_P", "D_pct_P")], 1, sd, na.rm = TRUE) / apply(homologies[, c("A_pct_P", "B_pct_P", "D_pct_P")], 1, mean, na.rm = TRUE)
  homologies$meth_CV_CSxP <- apply(homologies[, c("A_pct_CSxP", "B_pct_CSxP", "D_pct_CSxP")], 1, sd, na.rm = TRUE) / apply(homologies[, c("A_pct_CSxP", "B_pct_CSxP", "D_pct_CSxP")], 1, mean, na.rm = TRUE)
  homologies <- homologies %>%
    mutate(meth_dist_CS_P = e3(A_pct_CS, B_pct_CS, D_pct_CS, A_pct_P, B_pct_P, D_pct_P))
  homologies <- homologies %>%
    mutate(meth_dist_CS_hybrid = e3(A_pct_CS, B_pct_CS, D_pct_CS, A_pct_CSxP, B_pct_CSxP, D_pct_CSxP))
  homologies <- homologies %>%
    mutate(meth_dist_P_hybrid = e3(A_pct_P, B_pct_P, D_pct_P, A_pct_CSxP, B_pct_CSxP, D_pct_CSxP))
  
  homologies
}

meth_bias_cg <- meth_bias_df("merged_CG_symmetric_promoter1kb.txt.gz", "CG")
meth_bias_cg$context <- "CG"
meth_bias_chg <- meth_bias_df("merged_CHG_symmetric_promoter1kb.txt.gz", "CHG")
meth_bias_chg$context <- "CHG"
meth_bias_chh <- meth_bias_df("merged_CHH_promoter1kb.txt.gz", "CHH", remove_col3 = TRUE)
meth_bias_chh$context <- "CHH"
meth_bias_all <- rbind(meth_bias_cg, meth_bias_chg, meth_bias_chh)
meth_bias_all$context <- factor(meth_bias_all$context, levels = c("CG", "CHG", "CHH"))

cor_labs <- meth_bias_all %>%
  filter(!is.na(meth_dist_CS_hybrid), !is.na(meth_dist_P_hybrid)) %>%
  group_by(context) %>%
  summarise(r = cor(meth_dist_CS_hybrid, meth_dist_P_hybrid),p = cor.test(meth_dist_CS_hybrid, meth_dist_P_hybrid)$p.value,.groups = "drop") %>%
  mutate(label = paste0("R = ", sprintf("%.2f", r), "\nP = ", format.pval(p, digits = 2, eps = 1e-3)))

pdf(file="methylation_dist_parents_promoter.pdf",height=2,width=4.8)
ggplot(data = meth_bias_all, aes(x = meth_dist_CS_hybrid, y = meth_dist_P_hybrid, color = context)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey50") +
  geom_point(shape = 1, size = 1.2, stroke = 0.4) +
  geom_text(data = cor_labs, aes(x = 0.1, y = 1.2, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3) +
  facet_grid(~context) +
  scale_colour_manual(values = c(CG = "#D55E00", CHG = "#56B4E9", CHH = "#999999")) +
  labs(x = "Methylation Bias Distance to Chinese Spring", y = "Methylation Bias\nDistance to Paragon") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), plot.title = element_text(face = "bold"), legend.position = "none")
dev.off()


## methylation inheritance classification by chromatin states

cat_cols <- c(conserved_mC = "#1A1A1A", additive = "#D9D9D9", CS_dominant = "#0072B2", P_dominant = "#CC79A7", overdominant = "#E69F00", underdominant = "#1B7837")

summary_dt <- read.table(file="percent_met_chromatin_state_cg.tsv",header=T)
summary_dt$category <- factor(summary_dt$category, levels = names(cat_cols))
count_map_all <- aggregate(N ~ chromatin_state, data = summary_dt, FUN = sum, na.rm = TRUE)
count_map_all$label <- paste0(count_map_all$chromatin_state, " (", format(count_map_all$N, big.mark = ","), ")")
count_map_all <- setNames(count_map_all$label, count_map_all$chromatin_state)

pdf(paste0("met_classify_chromatin_cg.pdf"), width = 5, height = 3)
ggplot(summary_dt, aes(x = category, y = proportion, fill = category)) +
  geom_col(width = 0.8, colour = "black", linewidth = 0.2) +
  geom_text(aes(label = percent(proportion, accuracy = 0.1)), vjust = -0.35, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = cat_cols, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12))) +
  facet_wrap(~ chromatin_state, nrow = 1, labeller = labeller(chromatin_state = count_map_all)) +
  labs(x = NULL, y = "Percentage of sites") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), legend.position = "right", plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
dev.off()


summary_dt <- read.table(file="percent_met_chromatin_state_chg.tsv",header=T)
summary_dt$category <- factor(summary_dt$category, levels = names(cat_cols))
count_map_all <- aggregate(N ~ chromatin_state, data = summary_dt, FUN = sum, na.rm = TRUE)
count_map_all$label <- paste0(count_map_all$chromatin_state, " (", format(count_map_all$N, big.mark = ","), ")")
count_map_all <- setNames(count_map_all$label, count_map_all$chromatin_state)

pdf(paste0("met_classify_chromatin_chg.pdf"), width = 5, height = 3)
ggplot(summary_dt, aes(x = category, y = proportion, fill = category)) +
  geom_col(width = 0.8, colour = "black", linewidth = 0.2) +
  geom_text(aes(label = percent(proportion, accuracy = 0.1)), vjust = -0.35, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = cat_cols, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12))) +
  facet_wrap(~ chromatin_state, nrow = 1, labeller = labeller(chromatin_state = count_map_all)) +
  labs(x = NULL, y = "Percentage of sites") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), legend.position = "right", plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
dev.off()


summary_dt <- read.table(file="percent_met_chromatin_state_chh.tsv",header=T)
summary_dt$category <- factor(summary_dt$category, levels = names(cat_cols))
count_map_all <- aggregate(N ~ chromatin_state, data = summary_dt, FUN = sum, na.rm = TRUE)
count_map_all$label <- paste0(count_map_all$chromatin_state, " (", format(count_map_all$N, big.mark = ","), ")")
count_map_all <- setNames(count_map_all$label, count_map_all$chromatin_state)

pdf(paste0("met_classify_chromatin_chh.pdf"), width = 5, height = 3)
ggplot(summary_dt, aes(x = category, y = proportion, fill = category)) +
  geom_col(width = 0.8, colour = "black", linewidth = 0.2) +
  geom_text(aes(label = percent(proportion, accuracy = 0.1)), vjust = -0.35, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = cat_cols, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12))) +
  facet_wrap(~ chromatin_state, nrow = 1, labeller = labeller(chromatin_state = count_map_all)) +
  labs(x = NULL, y = "Percentage of sites") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), legend.position = "right", plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
dev.off()

