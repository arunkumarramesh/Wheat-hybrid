library(ggplot2)
library(scales)

## percent methylation
meth_cg <- read.table(file = "met_classify_cg_all_groups.tsv", header = TRUE)
meth_chg <- read.table(file = "met_classify_chg_all_groups.tsv", header = TRUE)
meth_chh <- read.table(file = "met_classify_chh_all_groups.tsv", header = TRUE)

meth_cg$context <- "CG"
meth_chg$context <- "CHG"
meth_chh$context <- "CHH"

meth_summary <- rbind(meth_cg, meth_chg, meth_chh)
meth_summary <- meth_summary[meth_summary$grouping_type %in% "subgenome",]

meth_summary$context <- factor(meth_summary$context, levels = c("CG", "CHG", "CHH"))
meth_summary$sample <- factor(meth_summary$sample, levels = c("CS", "CSxP", "P"))

sample_cols <- c(CS="#0072B2",CSxP="#E69F00",P="#CC79A7")

pdf("te_subgenome.pdf",height=2.5,width=4)
ggplot(meth_summary, aes(x = group, y = mean_pct, fill = sample, group = sample)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65, colour = "black", linewidth = 0.25) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(width = 0.7), width = 0.2, colour = "black", linewidth = 0.4) +
  scale_fill_manual(values = sample_cols) +
  facet_wrap(~ context, nrow = 1, scales="free_y") +
  labs(x = "Subgenome", y = "Mean Methylation (%)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), legend.position = "right", plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "right")
dev.off()


meth_cg <- read.table(file = "met_classify_cg_all_groups.tsv", header = TRUE)
meth_chg <- read.table(file = "met_classify_chg_all_groups.tsv", header = TRUE)
meth_chh <- read.table(file = "met_classify_chh_all_groups.tsv", header = TRUE)

meth_cg$context <- "CG"
meth_chg$context <- "CHG"
meth_chh$context <- "CHH"

meth_summary <- rbind(meth_cg, meth_chg, meth_chh)
meth_summary <- meth_summary[meth_summary$grouping_type %in% "distance_to_gene_class",]
meth_summary$group <- factor(meth_summary$group, levels = c("<1000","1000-5000",">5000"))
meth_summary$context <- factor(meth_summary$context, levels = c("CG", "CHG", "CHH"))
meth_summary$sample <- factor(meth_summary$sample, levels = c("CS", "CSxP", "P"))

sample_cols <- c(CS="#0072B2",CSxP="#E69F00",P="#CC79A7")

pdf("te_gene_dist.pdf",height=3,width=5)
ggplot(meth_summary, aes(x = group, y = mean_pct, fill = sample, group = sample)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65, colour = "black", linewidth = 0.25) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(width = 0.7), width = 0.2, colour = "black", linewidth = 0.4) +
  scale_fill_manual(values = sample_cols) +
  facet_wrap(~ context, nrow = 1, scales="free_y") +
  labs(x = "Distance to gene (bp)", y = "Mean Methylation (%)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), legend.position = "right", plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 0.5), legend.position = "right")
dev.off()

meth_cg <- read.table(file = "met_classify_cg_all_groups.tsv", header = TRUE)
meth_chg <- read.table(file = "met_classify_chg_all_groups.tsv", header = TRUE)
meth_chh <- read.table(file = "met_classify_chh_all_groups.tsv", header = TRUE)

meth_cg$context <- "CG"
meth_chg$context <- "CHG"
meth_chh$context <- "CHH"

meth_summary <- rbind(meth_cg, meth_chg, meth_chh)
meth_summary <- meth_summary[meth_summary$grouping_type %in% "te_class",]
meth_summary <- meth_summary[!meth_summary$group %in% c("XXX","no","exon%2C"),]
te_names <- c(DHH = "Helitron", DTA = "hAT", DTC = "CACTA", DTH = "PIF-Harbinger", DTM = "Mutator", DTT = "Tc1-Mariner", RLC = "Copia", RLG = "Gypsy")
meth_summary$group <- ifelse(meth_summary$group %in% names(te_names), te_names[meth_summary$group], meth_summary$group)
meth_summary$context <- factor(meth_summary$context, levels = c("CG", "CHG", "CHH"))
meth_summary$sample <- factor(meth_summary$sample, levels = c("CS", "CSxP", "P"))

sample_cols <- c(CS="#0072B2",CSxP="#E69F00",P="#CC79A7")

pdf("te_class.pdf",height=3,width=7)
ggplot(meth_summary, aes(x = group, y = mean_pct, fill = sample, group = sample)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65, colour = "black", linewidth = 0.25) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(width = 0.7), width = 0.2, colour = "black", linewidth = 0.4) +
  scale_fill_manual(values = sample_cols) +
  facet_wrap(~ context, nrow = 1, scales="free_y") +
  labs(x = "TE class", y = "Mean Methylation (%)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), legend.position = "right", plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5), legend.position = "right")
dev.off()

## CS vs hybrid methylation difference

library(ggplot2)
library(scales)

cols <- c(CG = "#D55E00", CHG = "#56B4E9", CHH = "#999999")

## percent methylation
meth_cg <- read.table(file = "met_diff__cg_all_groups.tsv", header = TRUE)
meth_chg <- read.table(file = "met_diff__chg_all_groups.tsv", header = TRUE)
meth_chh <- read.table(file = "met_diff__chh_all_groups.tsv", header = TRUE)

meth_cg$context <- "CG"
meth_chg$context <- "CHG"
meth_chh$context <- "CHH"

meth_summary <- rbind(meth_cg, meth_chg, meth_chh)
meth_summary <- meth_summary[meth_summary$grouping_type %in% "subgenome",]

meth_summary$context <- factor(meth_summary$context, levels = c("CG", "CHG", "CHH"))

pdf("te_diff.pdf",height=2,width=3.5)
ggplot(meth_summary, aes(x = group, y = mean_diff_CSxP_minus_CS, fill=context)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65, colour = "black", linewidth = 0.25) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(width = 0.7), width = 0.2, colour = "black", linewidth = 0.4) +
  scale_fill_manual(values = cols) +
  facet_wrap(~ context, nrow = 1) +
  labs(x = "Subgenome", y = "CSxP- CS\nMethylation (%)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), legend.position = "right", plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "right")
dev.off()

library(data.table)
library(ggplot2)
sample_cols <- c(CS = "#0072B2", CSxP = "#E69F00")

## TE methylation bias

gene_TE_CG <- fread(file="gene_level_TE_CG.tsv")
gene_TE_CG <- gene_TE_CG[,1:7, with = FALSE]
gene_TE_CG$context <- "CG"
gene_TE_CHG <- fread(file="gene_level_TE_CHG.tsv")
gene_TE_CHG <- gene_TE_CHG[,1:7, with = FALSE]
gene_TE_CHG$context <- "CHG"
gene_TE_CHH <- fread(file="gene_level_TE_CHH.tsv")
gene_TE_CHH <- gene_TE_CHH[,1:7, with = FALSE]
gene_TE_CHH$context <- "CHH"
gene_TE <- rbind(gene_TE_CG,gene_TE_CHG,gene_TE_CHH)
gene_TE <- gene_TE[gene_TE$n_sites > 5,]

homologies <- read.csv(file="homologies.csv")
homologies <- homologies[1:4]

library(data.table)

homologies <- as.data.table(homologies)
gene_TE <- as.data.table(gene_TE)

hom_long <- melt(homologies, id.vars = "group_id", measure.vars = c("A","B","D"), variable.name = "homoeolog", value.name = "gene_id")
gene_TE_sub <- gene_TE[gene_id %in% hom_long$gene_id, .(gene_id, context, pct_CS, pct_CSxP)]
merged <- merge(hom_long, gene_TE_sub, by = "gene_id", allow.cartesian = TRUE)
triad_meth <- dcast(merged, group_id + context ~ homoeolog, value.var = c("gene_id","pct_CS","pct_CSxP"))
triad_meth <- triad_meth[!is.na(pct_CS_A) & !is.na(pct_CS_B) & !is.na(pct_CS_D) & !is.na(pct_CSxP_A) & !is.na(pct_CSxP_B) & !is.na(pct_CSxP_D)]
triad_meth <- triad_meth[, .(A = gene_id_A, B = gene_id_B, D = gene_id_D, group_id, context, CS_A = pct_CS_A / 100, CS_B = pct_CS_B / 100, CS_D = pct_CS_D / 100, CSxP_A = pct_CSxP_A / 100, CSxP_B = pct_CSxP_B / 100, CSxP_D = pct_CSxP_D / 100)]
triad_meth[, CS_cv := apply(.SD, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)), .SDcols = c("CS_A","CS_B","CS_D")]
triad_meth[, CSxP_cv := apply(.SD, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)), .SDcols = c("CSxP_A","CSxP_B","CSxP_D")]

triad_meth <- as.data.table(triad_meth)

cols <- c(CG = "#D55E00", CHG = "#56B4E9", CHH = "#999999")

cor_lab <- triad_meth[, {
  ct <- cor.test(CS_cv, CSxP_cv)
  p_lab <- ifelse(ct$p.value < 0.001, "P < 0.001", sprintf("P = %.3f", ct$p.value))
  .(R = unname(ct$estimate), P = ct$p.value, label = sprintf("R = %.2f\n%s", unname(ct$estimate), p_lab))
}, by = context]

pdf("te_bias.pdf",height=2,width=4.5)
ggplot(data = triad_meth, aes(x = CS_cv, y = CSxP_cv, colour = context)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey50") +
  geom_point(shape = 1, size = 1.4, stroke = 0.4, alpha = 0.8) +
  geom_text(data = cor_lab, aes(x = 1.7, y = 0.1, label = label), hjust = 1, vjust = 0, size = 3.2, inherit.aes = FALSE) +
  scale_colour_manual(values = cols) +
  facet_wrap(~context, scales = "free") +
  labs(x = "CS TE methylation bias", y = "CSxP TE\nmethylation bias", colour = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"), strip.background = element_rect(fill = "white", colour = "black"), axis.title = element_text(face = "bold"), legend.position = "none")
dev.off()
