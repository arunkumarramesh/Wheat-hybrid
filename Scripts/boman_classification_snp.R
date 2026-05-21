library(data.table)
library(ggplot2)
library(scales)

snp_counts <- read.table(file="cs_par_snps_50kb_counts.tsv")
colnames(snp_counts) <- c("chr","start","end","SNP_count")
snp_counts$subgenome <- sub("^triticum_aestivum\\.[0-9]+([ABD])$","\\1",snp_counts$chr)
snp_counts$subgenome <- factor(snp_counts$subgenome,levels = c("A","B","D"))

pdf(file="SNP_50kb.pdf",height=3,width=7)
ggplot(snp_counts,aes(x = SNP_count)) +
  geom_histogram(bins = 80,colour = "black",fill = "grey70") +
  facet_wrap(~subgenome) +
  scale_x_log10() +
  labs(x = "Number of SNVs per 50kb",y = "Number of 50kb intervals") +
  theme_bw(base_size = 12) +
  theme(axis.title = element_text(face = "bold"),panel.grid.minor = element_blank())
dev.off()

process_methylation_file_snp <- function(infile,snp_count_file,out_suffix,context,region,remove_third_col = FALSE) {
  dt <- fread(infile)
  
  if (remove_third_col) {
    dt <- dt[, -3, with = FALSE]
  }
  
  snp_counts <- as.data.table(read.table(file="cs_par_snps_50kb_counts.tsv"))
  colnames(snp_counts) <- c("chr","window_start","window_end","SNP_count_50kb")
  
  snp_counts[, chr := sub("^triticum_aestivum\\.([0-9]+[ABD])$","chr\\1",chr)]
  
  dt <- dt[cov_CS > 10 & cov_CSxP > 10 & cov_P > 10]
  
  dt$subgenome <- sub("^chr[0-9]+([ABD])$","\\1",dt$chr)
  dt <- dt[dt$subgenome %in% c("A","B","D"),]
  dt$subgenome <- factor(dt$subgenome,levels = c("A","B","D"))
  
  dt[, window_start := floor(pos / 50000) * 50000]
  dt[, window_end := window_start + 50000]
  
  dt <- merge(dt,snp_counts,by = c("chr","window_start","window_end"),all.x = FALSE)
  
  dt[, SNP_density_group := fifelse(SNP_count_50kb > 1000,"High SNV density",
                                    fifelse(SNP_count_50kb < 10,"Low SNV density",NA_character_))]
  
  dt <- dt[!is.na(SNP_density_group)]
  dt$SNP_density_group <- factor(dt$SNP_density_group,levels = c("Low SNV density","High SNV density"))
  
  dt$A <- dt$pct_CS / 100
  dt$H <- dt$pct_CSxP / 100
  dt$B <- dt$pct_P / 100
  
  meth_long <- rbind(
    data.table(subgenome = dt$subgenome,SNP_density_group = dt$SNP_density_group,sample = "CS",pct = dt$pct_CS),
    data.table(subgenome = dt$subgenome,SNP_density_group = dt$SNP_density_group,sample = "CSxP",pct = dt$pct_CSxP),
    data.table(subgenome = dt$subgenome,SNP_density_group = dt$SNP_density_group,sample = "P",pct = dt$pct_P)
  )
  
  meth_summary <- ci_summary(meth_long,c("subgenome","SNP_density_group","sample"),"pct")
  setnames(meth_summary,"mean_value","mean_pct")
  meth_summary$context <- context
  meth_summary$region <- region
  
  dt$x <- dt$H - dt$A
  dt$y <- dt$H - dt$B
  dt$radius <- sqrt(dt$x^2 + dt$y^2)
  dt$angle_deg <- (atan2(dt$y,dt$x) * 180 / pi) %% 360
  
  dt$sector_center <- vapply(dt$angle_deg,nearest_sector_center,numeric(1))
  dt$max_mC <- pmax(dt$A,dt$H,dt$B,na.rm = TRUE)
  dt <- dt[max_mC >= 0.01]
  
  dt$category <- ifelse(dt$radius < 0.1,"conserved_mC",vapply(dt$sector_center,sector_to_class,character(1)))
  dt$category <- factor(dt$category,levels = c("conserved_mC","additive","CS_dominant","P_dominant","overdominant","underdominant"))
  
  summary_dt <- as.data.table(table(dt$subgenome,dt$SNP_density_group,dt$category))
  colnames(summary_dt) <- c("subgenome","SNP_density_group","category","N")
  summary_dt$subgenome <- factor(summary_dt$subgenome,levels = c("A","B","D"))
  summary_dt$SNP_density_group <- factor(summary_dt$SNP_density_group,levels = c("Low SNV density","High SNV density"))
  summary_dt$category <- factor(summary_dt$category,levels = c("conserved_mC","additive","CS_dominant","P_dominant","overdominant","underdominant"))
  summary_dt$proportion <- ave(summary_dt$N,summary_dt$subgenome,summary_dt$SNP_density_group,FUN = function(x) x / sum(x))
  
  summary_dt_all <- dt[, .N, by = .(SNP_density_group,category)]
  summary_dt_all[, proportion := N / sum(N), by = SNP_density_group]
  summary_dt_all$SNP_density_group <- factor(summary_dt_all$SNP_density_group,levels = c("Low SNV density","High SNV density"))
  summary_dt_all$category <- factor(summary_dt_all$category,levels = c("conserved_mC","additive","CS_dominant","P_dominant","overdominant","underdominant"))
  
  panel_counts_all <- dt[, .N, by = SNP_density_group]
  panel_counts_all$SNP_density_group <- factor(panel_counts_all$SNP_density_group,levels = c("Low SNV density","High SNV density"))
  panel_counts_all <- panel_counts_all[order(panel_counts_all$SNP_density_group),]
  count_map_all <- setNames(paste0(panel_counts_all$SNP_density_group," (n=",panel_counts_all$N,")"),as.character(panel_counts_all$SNP_density_group))
  
  p1 <- ggplot(summary_dt_all,aes(x = category,y = proportion,fill = category)) +
    geom_col(width = 0.8,linewidth = 0.2) +
    geom_text(aes(label = percent(proportion,accuracy = 0.1)),vjust = -0.35,size = 3.2,fontface = "bold") +
    scale_fill_manual(values = cat_cols) +
    scale_y_continuous(labels = percent_format(accuracy = 1),expand = expansion(mult = c(0,0.12))) +
    facet_wrap(~ SNP_density_group,nrow = 1,labeller = labeller(SNP_density_group = count_map_all)) +
    labs(title = paste0(region," ",context," (n=",sum(summary_dt_all$N),")"),x = NULL,y = "Percentage of sites") +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(linewidth = 0.2,colour = "grey90"),strip.background = element_rect(fill = "white",colour = "black"),plot.title = element_text(face = "bold",hjust = 0.5),axis.title = element_text(face = "bold"),axis.text.x = element_text(angle = 35,hjust = 1),legend.position = "none")
  
  pdf(paste0("met_classify_snp",out_suffix,".pdf"),width = 5.8,height = 3)
  print(p1)
  dev.off()
  
  return(meth_summary)
}

all_meth_summary <- list()
all_meth_summary[["cg_gene"]] <- process_methylation_file_snp(infile = "merged_CG_symmetric_CDS.txt.gz",out_suffix = "_cg_gene",context = "CG",region = "CDS")
all_meth_summary[["cg_promoter"]] <- process_methylation_file_snp(infile = "merged_CG_symmetric_promoter1kb.txt.gz",out_suffix = "_cg_promoter",context = "CG",region = "Promoter")
all_meth_summary[["chg_gene"]] <- process_methylation_file_snp(infile = "merged_CHG_symmetric_CDS.txt.gz",out_suffix = "_chg_gene",context = "CHG",region = "CDS")
all_meth_summary[["chg_promoter"]] <- process_methylation_file_snp(infile = "merged_CHG_symmetric_promoter1kb.txt.gz",out_suffix = "_chg_promoter",context = "CHG",region = "Promoter")
all_meth_summary[["chh_gene"]] <- process_methylation_file_snp(infile = "merged_CHH_all_CDS.txt.gz",out_suffix = "_chh_gene",context = "CHH",region = "CDS",remove_third_col = TRUE)
all_meth_summary[["chh_promoter"]] <- process_methylation_file_snp(infile = "merged_CHH_promoter1kb.txt.gz",out_suffix = "_chh_promoter",context = "CHH",region = "Promoter",remove_third_col = TRUE)

meth_summary_all <- rbindlist(all_meth_summary)
meth_summary_all$context <- factor(meth_summary_all$context,levels = c("CG","CHG","CHH"))
meth_summary_all$region <- factor(meth_summary_all$region,levels = c("CDS","Promoter"))
meth_summary_all$subgenome <- factor(meth_summary_all$subgenome,levels = c("A","B","D"))
meth_summary_all$sample <- factor(meth_summary_all$sample,levels = c("CS","CSxP","P"))
meth_summary_all$SNP_density_group <- factor(meth_summary_all$SNP_density_group,levels = c("Low SNV density","High SNV density"))

meth_summary_all$x_group <- paste(meth_summary_all$region,meth_summary_all$subgenome,sep = ".")
x_levels <- c("CDS.A","CDS.B","CDS.D","Promoter.A","Promoter.B","Promoter.D")
meth_summary_all$x_group <- factor(meth_summary_all$x_group,levels = x_levels)

p_percent_all <- ggplot(meth_summary_all,aes(x = x_group,y = mean_pct,fill = sample)) +
  geom_vline(xintercept = 3.5,linetype = "dashed",linewidth = 0.3,colour = "grey50") +
  geom_col(aes(group = sample),position = position_dodge(width = 0.75),width = 0.65,colour = "black",linewidth = 0.25) +
  geom_errorbar(aes(ymin = ci_low,ymax = ci_high,group = sample),position = position_dodge(width = 0.75),width = 0.2,colour = "black",linewidth = 0.4) +
  scale_fill_manual(values = sample_cols) +
  scale_x_discrete(labels = c("CDS.A" = "A","CDS.B" = "B","CDS.D" = "D","Promoter.A" = "A","Promoter.B" = "B","Promoter.D" = "D")) +
  facet_grid(SNP_density_group ~ context) +
  labs(x = NULL,y = "Mean Methylation (%)",fill = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(linewidth = 0.2,colour = "grey90"),strip.background = element_rect(fill = "white",colour = "black"),legend.position = "right",plot.title = element_text(face = "bold"),axis.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(size = 10),legend.position = "right",plot.margin = margin(5.5,5.5,28,5.5)) +
  annotate("text",x = 2,y = -Inf,label = "CDS",vjust = 3.2,size = 3.5,fontface = "bold") +
  annotate("text",x = 5,y = -Inf,label = "Promoter",vjust = 3.2,size = 3.5,fontface = "bold") +
  coord_cartesian(clip = "off")

pdf("percent_met_snp.pdf", width = 6.5, height = 3.4)
p_percent_all
dev.off()



## across all sites

cat_cols <- c(conserved_mC = "#1A1A1A", additive = "#D9D9D9", CS_dominant = "#0072B2", P_dominant = "#CC79A7", overdominant = "#E69F00", underdominant = "#1B7837")

## methylation classification according to Boman et al. 2024

summary_dt <- read.table(file="percent_met_snp_density_all_cg.tsv",header=T,sep="\t")
summary_dt$SNP_density_group <- gsub("High SNP density","High SNV density",summary_dt$SNP_density_group)
summary_dt$SNP_density_group <- gsub("Low SNP density","Low SNV density",summary_dt$SNP_density_group)
summary_dt$category <- factor(summary_dt$category, levels = names(cat_cols))
count_map_all <- aggregate(N ~ SNP_density_group, data = summary_dt, FUN = sum, na.rm = TRUE)
count_map_all$label <- paste0(count_map_all$SNP_density_group, " (", format(count_map_all$N, big.mark = ","), ")")
count_map_all <- setNames(count_map_all$label, count_map_all$SNP_density_group)

pdf(paste0("met_classify_cg_snp.pdf"), width = 5, height = 3)
ggplot(summary_dt, aes(x = category, y = proportion, fill = category)) +
  geom_col(width = 0.8, colour = "black", linewidth = 0.2) +
  geom_text(aes(label = percent(proportion, accuracy = 0.1)), vjust = -0.35, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = cat_cols, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12))) +
  facet_wrap(~ SNP_density_group, nrow = 1, labeller = labeller(SNP_density_group = count_map_all)) +
  labs(x = NULL, y = "Percentage of sites") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), legend.position = "right", plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
dev.off()

summary_dt <- read.table(file="percent_met_snp_density_all_chg.tsv",header=T,sep="\t")
summary_dt$SNP_density_group <- gsub("High SNP density","High SNV density",summary_dt$SNP_density_group)
summary_dt$SNP_density_group <- gsub("Low SNP density","Low SNV density",summary_dt$SNP_density_group)
summary_dt$category <- factor(summary_dt$category, levels = names(cat_cols))
count_map_all <- aggregate(N ~ SNP_density_group, data = summary_dt, FUN = sum, na.rm = TRUE)
count_map_all$label <- paste0(count_map_all$SNP_density_group, " (", format(count_map_all$N, big.mark = ","), ")")
count_map_all <- setNames(count_map_all$label, count_map_all$SNP_density_group)

pdf(paste0("met_classify_chg_snp.pdf"), width = 5, height = 3)
ggplot(summary_dt, aes(x = category, y = proportion, fill = category)) +
  geom_col(width = 0.8, colour = "black", linewidth = 0.2) +
  geom_text(aes(label = percent(proportion, accuracy = 0.1)), vjust = -0.35, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = cat_cols, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12))) +
  facet_wrap(~ SNP_density_group, nrow = 1, labeller = labeller(SNP_density_group = count_map_all)) +
  labs(x = NULL, y = "Percentage of sites") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), legend.position = "right", plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
dev.off()

summary_dt <- read.table(file="percent_met_snp_density_all_chh.tsv",header=T,sep="\t")
summary_dt$SNP_density_group <- gsub("High SNP density","High SNV density",summary_dt$SNP_density_group)
summary_dt$SNP_density_group <- gsub("Low SNP density","Low SNV density",summary_dt$SNP_density_group)
summary_dt$category <- factor(summary_dt$category, levels = names(cat_cols))
count_map_all <- aggregate(N ~ SNP_density_group, data = summary_dt, FUN = sum, na.rm = TRUE)
count_map_all$label <- paste0(count_map_all$SNP_density_group, " (", format(count_map_all$N, big.mark = ","), ")")
count_map_all <- setNames(count_map_all$label, count_map_all$SNP_density_group)

pdf(paste0("met_classify_chh_snp.pdf"), width = 5, height = 3)
ggplot(summary_dt, aes(x = category, y = proportion, fill = category)) +
  geom_col(width = 0.8, colour = "black", linewidth = 0.2) +
  geom_text(aes(label = percent(proportion, accuracy = 0.1)), vjust = -0.35, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = cat_cols, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12))) +
  facet_wrap(~ SNP_density_group, nrow = 1, labeller = labeller(SNP_density_group = count_map_all)) +
  labs(x = NULL, y = "Percentage of sites") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), legend.position = "right", plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
dev.off()



## percent methylation
meth_cg <- read.table(file = "met_classify_snp_density_cg.tsv", header = TRUE,sep="\t")
meth_chg <- read.table(file = "met_classify_snp_density_chg.tsv", header = TRUE,sep="\t")
meth_chh <- read.table(file = "met_classify_snp_density_chh.tsv", header = TRUE,sep="\t")

meth_cg$context <- "CG"
meth_chg$context <- "CHG"
meth_chh$context <- "CHH"

meth_summary <- rbind(meth_cg, meth_chg, meth_chh)
meth_summary$SNP_density_group <- gsub("High SNP density","High SNV density",meth_summary$SNP_density_group)
meth_summary$SNP_density_group <- gsub("Low SNP density","Low SNV density",meth_summary$SNP_density_group)

meth_summary$context <- factor(meth_summary$context, levels = c("CG", "CHG", "CHH"))
meth_summary$sample <- factor(meth_summary$sample, levels = c("CS", "CSxP", "P"))
meth_summary$subgenome <- factor(meth_summary$subgenome, levels = c("A", "B", "D"))

sample_cols <- c(CS="#0072B2",CSxP="#E69F00",P="#CC79A7")

pdf(paste0("percent_meth_snp.pdf"), width = 5, height = 4.5)
ggplot(meth_summary, aes(x = SNP_density_group, y = mean_pct, fill = sample, group = sample)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65, colour = "black", linewidth = 0.25) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(width = 0.7), width = 0.2, colour = "black", linewidth = 0.4) +
  scale_fill_manual(values = sample_cols) +
  facet_grid(context ~ subgenome,scales = "free_y") +
  labs(x = NULL, y = "Mean Methylation (%)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), legend.position = "right", plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
dev.off()
