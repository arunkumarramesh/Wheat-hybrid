
library(data.table)
library(dplyr)
library(ggplot2)
library(scales)

## identify gbM genes

test_gbM <- function(file, context_name, remove_col3 = FALSE) {
  dt <- fread(file)
  if (remove_col3) {
    dt <- dt[, -3, with = FALSE]
  }
  dt <- dt[!is.na(gene_id)]
  samples <- c("CS", "CSxP", "P")
  context_results <- data.table()
  for (sample_name in samples) {
    pct_col <- paste0("pct_", sample_name)
    cov_col <- paste0("cov_", sample_name)
    x <- dt
    x$meth_count <- round((x[[pct_col]] / 100) * x[[cov_col]])
    x$total_count <- x[[cov_col]]
    genome_p <- sum(x$meth_count, na.rm = TRUE) / sum(x$total_count, na.rm = TRUE)
    gene_summary <- x[,list(sites = .N,nC = sum(total_count, na.rm = TRUE),mC = sum(meth_count, na.rm = TRUE)),by = gene_id]
    gene_summary <- gene_summary[nC >= 10]
    gene_summary$sample <- sample_name
    gene_summary$context <- context_name
    gene_summary$genome_p <- genome_p
    gene_summary$gene_methylation_level <- gene_summary$mC / gene_summary$nC
    gene_summary$p_value <- pbinom(q = gene_summary$mC - 1,size = gene_summary$nC,prob = gene_summary$genome_p,lower.tail = FALSE)
    gene_summary$padj <- p.adjust(gene_summary$p_value, method = "BH")
    gene_summary$context_candidate <- gene_summary$padj < 0.05 & gene_summary$gene_methylation_level > gene_summary$genome_p
    context_results <- rbind(context_results, gene_summary)
  }
  return(context_results)
}

cg_results <- test_gbM(file = "merged_CG_symmetric_CDS.txt.gz",context_name = "CG",remove_col3 = FALSE)
chg_results <- test_gbM(file = "merged_CHG_symmetric_CDS.txt.gz",context_name = "CHG",remove_col3 = FALSE)
chh_results <- test_gbM(file = "merged_CHH_all_CDS.txt.gz",context_name = "CHH",remove_col3 = TRUE)

all_results <- rbind(cg_results, chg_results, chh_results)
all_results <- dcast(all_results, gene_id + sample ~ context, value.var = c("sites", "nC", "mC", "gene_methylation_level", "genome_p", "p_value", "padj", "context_candidate"))
all_results <- all_results %>%
  group_by(gene_id) %>%
  filter(all(sites_CG > 5), all(sites_CHG > 5),all(sites_CHH > 5)) %>%
  ungroup()
all_results$gbM_candidate <- all_results$context_candidate_CG == TRUE &all_results$context_candidate_CHG == FALSE &all_results$context_candidate_CHH == FALSE
all_results <- all_results[, -c(18:20, 24:26), with = FALSE]
fwrite(all_results, "gene_body_methylation.tsv", sep = "\t")

library(data.table)
library(ggplot2)

gbm_status <- fread("gene_body_methylation.tsv")
gbm_status <- gbm_status[, c("gene_id", "sample", "gbM_candidate")]
gbm_status <- dcast(gbm_status, gene_id ~ sample, value.var = "gbM_candidate")
gbm_status$parent_status <- ifelse(gbm_status$CS == TRUE & gbm_status$P == TRUE,"gbM in both parents",
                                   ifelse(gbm_status$CS == TRUE & gbm_status$P == FALSE,"CS-specific gbM",
                                          ifelse(gbm_status$CS == FALSE & gbm_status$P == TRUE,"P-specific gbM","non-gbM in both parents"
                                          )
                                   )
)
gbm_status$hybrid_status <- ifelse(gbm_status$CSxP == TRUE, "gbM", "non-gbM")
gbm_status$parent_status <- factor(gbm_status$parent_status,levels = c("gbM in both parents","CS-specific gbM","P-specific gbM","non-gbM in both parents"))
gbm_status$hybrid_status <- factor(gbm_status$hybrid_status,levels = c("gbM", "non-gbM"))
gbm_status$subgenome <- sub("^TraesCS[0-9]+([ABD]).*", "\\1", gbm_status$gene_id)
gbm_status <- gbm_status[gbm_status$subgenome %in% c("A", "B", "D"), ]
gbm_status$subgenome <- factor(gbm_status$subgenome, levels = c("A", "B", "D"))
gbm_status <- gbm_status[!is.na(gbm_status$parent_status),]
gbm_status <- gbm_status[!is.na(gbm_status$hybrid_status),]

# classify gbM by subgenome
plot_dt <- gbm_status %>%
  group_by(subgenome, parent_status, hybrid_status) %>%
  summarise(N = n(), .groups = "drop")
plot_dt$subgenome <- factor(plot_dt$subgenome, levels = c("A", "B", "D"))
plot_dt$parent_status <- factor(plot_dt$parent_status,levels = c("gbM in both parents", "CS-specific gbM", "P-specific gbM", "non-gbM in both parents"))
plot_dt$hybrid_status <- factor(plot_dt$hybrid_status,levels = c("gbM", "non-gbM"))

pdf(file="gbM_parents_hybrid.pdf",height=3,width=6.2)
ggplot(plot_dt, aes(x = hybrid_status, y = parent_status, fill = N)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = comma(N)), size = 3) +
  facet_wrap(~subgenome, nrow = 1) +
  scale_fill_gradient(low = "white", high = "#0072B2", labels = comma) +
  labs(x = "Hybrid", y = NULL, fill = "Number of Genes") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),strip.background = element_rect(fill = "white", colour = "black"),axis.text.x = element_text(angle = 45, hjust = 1),axis.title = element_text(face = "bold"))
dev.off()

# classify gbM by gene regulatory classes
Reg_Class <- read.csv("classified_all.csv", header = TRUE)
Reg_Class <- Reg_Class[!Reg_Class$category %in% "Ambiguous", ]
colnames(Reg_Class)[1] <- "gene_id"

gbm_status_reg <- inner_join(gbm_status,Reg_Class,by="gene_id")
write.csv(gbm_status_reg,file="gbM_by_regulatory_state.csv",row.names = F)

plot_dt <- gbm_status_reg %>%
  group_by(category, parent_status, hybrid_status) %>%
  summarise(N = n(), .groups = "drop")
plot_dt$subgenome <- factor(plot_dt$category)
plot_dt$parent_status <- factor(plot_dt$parent_status,levels = c("gbM in both parents", "CS-specific gbM", "P-specific gbM", "non-gbM in both parents"))
plot_dt$hybrid_status <- factor(plot_dt$hybrid_status,levels = c("gbM", "non-gbM"))

plot_dt$category_label <- plot_dt$category

plot_dt$category_label <- gsub("Cis", "italic('Cis')", plot_dt$category_label)
plot_dt$category_label <- gsub("Trans", "italic('Trans')", plot_dt$category_label)
plot_dt$category_label <- gsub(" only", "~only", plot_dt$category_label)
plot_dt$category_label <- gsub(" \\+ ", "~'+'~", plot_dt$category_label)

pdf(file="gbM_parents_hybrid_reg_class.pdf",height=3,width=7.5)
ggplot(plot_dt, aes(x = hybrid_status, y = parent_status, fill = N)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = comma(N)), size = 3) +
  facet_wrap(~category_label, nrow = 1, labeller = label_parsed) +
  scale_fill_gradient(low = "white", high = "#0072B2", labels = comma) +
  labs(x = "Hybrid", y = NULL, fill = "Number of Genes") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),strip.background = element_rect(fill = "white", colour = "black"),axis.text.x = element_text(angle = 45, hjust = 1),axis.title = element_text(face = "bold"))
dev.off()

## gbM by expression level
library(data.table)
library(ggplot2)
cs_tpm <- read.table(file="cs_tpm.tsv")
cs_tpm <- as.data.table(cs_tpm, keep.rownames = "gene_id")
cs_tpm <- cs_tpm[!grepl("LC$", gene_id)]
colnames(cs_tpm) <- c("gene_id", sub("^PXCS", "PxCS", sub("_.*", "", colnames(cs_tpm)[-1])))
sample_cols <- colnames(cs_tpm)[-1]
cs_tpm <- cs_tpm[rowSums(as.data.frame(cs_tpm)[, sample_cols] >= 1) >= 3]
cs_tpm <- cs_tpm[,1:10]
gbm_status <- read.table("gene_body_methylation.tsv",header=T)
gbm_status <- gbm_status[c(1,2,21)]

expr_long <- melt(cs_tpm, id.vars = "gene_id", variable.name = "replicate", value.name = "TPM")
expr_long$sample <- sub("[0-9]+$", "", expr_long$replicate)
expr_long <- expr_long[expr_long$sample %in% c("CS", "CSxP", "P"), ]

plot_dt <- merge(expr_long, gbm_status, by = c("gene_id", "sample"))
plot_dt$gbM_status <- ifelse(plot_dt$gbM_candidate, "gbM", "non-gbM")
plot_dt$gbM_status <- factor(plot_dt$gbM_status, levels = c("non-gbM", "gbM"))
plot_dt$sample <- factor(plot_dt$sample, levels = c("CS", "CSxP", "P"))
plot_dt$log2_TPM <- log2(plot_dt$TPM + 1)
plot_dt <- plot_dt[!is.na(plot_dt$gbM_status),]

median_labs <- plot_dt %>%
  group_by(gbM_status) %>%
  summarise(median_log2_TPM = median(log2_TPM, na.rm = TRUE), .groups = "drop") %>%
  mutate(label = round(median_log2_TPM, 2))

p_lab <- data.frame(label = paste0("italic(P)~'='~",signif(wilcox.test(log2_TPM ~ gbM_status, data = plot_dt)$p.value, 3)))

pdf("expression_gbM_vs_non_gbM_by_genotype.pdf", width = 2, height = 3)
ggplot(plot_dt, aes(x = gbM_status, y = log2_TPM, fill = gbM_status)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.3) +
  geom_text(data = median_labs, aes(x = gbM_status, y = median_log2_TPM, label = label), inherit.aes = FALSE, size = 3, vjust = -0.4) +
  geom_text(data = p_lab, aes(x = 1.5, y = 15, label = paste0("italic(P)~'='~", signif(wilcox.test(plot_dt$log2_TPM ~ plot_dt$gbM_status)$p.value, 3))), inherit.aes = FALSE, size = 4, parse = TRUE) +
  scale_fill_manual(values = c("non-gbM" = "grey70", "gbM" = "#0072B2")) +
  labs(x = NULL, y = "Mean log2(TPM + 1)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),strip.background = element_rect(fill = "white", colour = "black"),legend.position = "none",axis.title = element_text(face = "bold"))
dev.off()
