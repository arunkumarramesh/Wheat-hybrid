library(data.table)
library(ggplot2)
library(dplyr)

choulet <- fread("choulet_URGI_tpm.tsv.gz")
expr_cols <- setdiff(names(choulet), "gene")
long <- melt(choulet, id.vars="gene", measure.vars=expr_cols, variable.name="sample", value.name="TPM")
long[, condition := sub("_rep[0-9]+$", "", sample)]
cond <- long[, .(mean_TPM=mean(TPM, na.rm=TRUE), n_reps_TPM_ge_1=sum(TPM >= 1, na.rm=TRUE)), by=.(gene, condition)]
wide <- dcast(cond, gene ~ condition, value.var="mean_TPM")
leaf_filter <- cond[condition == "leaf_Z10", .(gene, leaf_Z10_TPM=mean_TPM, leaf_Z10_reps_TPM_ge_1=n_reps_TPM_ge_1)]
wide <- merge(wide, leaf_filter, by="gene")
condition_cols <- setdiff(names(wide), c("gene", "leaf_Z10_TPM", "leaf_Z10_reps_TPM_ge_1"))

calc_tau <- function(x) {
  x <- as.numeric(x)
  x_non_na <- x[!is.na(x)]
  if (length(x_non_na) <= 1 || max(x_non_na) == 0) {
    return(NA_real_)
  }
  x_scaled <- x_non_na / max(x_non_na)
  sum(1 - x_scaled) / (length(x_non_na) - 1)
}

wide[, tau := apply(.SD, 1, calc_tau), .SDcols=condition_cols]
wide[, max_TPM := apply(.SD, 1, max, na.rm=TRUE), .SDcols=condition_cols]
wide[, leaf_Z10_fraction_of_max := leaf_Z10_TPM / max_TPM]
wide[, n_conditions_TPM_ge_1 := rowSums(.SD >= 1, na.rm=TRUE), .SDcols=condition_cols]
wide[, frac_conditions_TPM_ge_1 := n_conditions_TPM_ge_1 / length(condition_cols)]

gene_labels <- rbind(
  wide[leaf_Z10_TPM >= 1 & leaf_Z10_reps_TPM_ge_1 >= 1 & tau >= 0.8 & leaf_Z10_fraction_of_max >= 0.8, .(label="seedling_specific", gene)],
  wide[leaf_Z10_TPM >= 1 & leaf_Z10_reps_TPM_ge_1 >= 1 & tau <= 0.4 & frac_conditions_TPM_ge_1 >= 0.7, .(label="broad", gene)]
)

iwgsc_refseq_all_correspondances <- read.table("iwgsc_refseq_all_correspondances.csv", header=T)

## convert v1.1 to v2.1
map_unique <- iwgsc_refseq_all_correspondances %>%
  transmute(v11 = `v1.1`, v21 = `v2.1`) %>%
  filter(!is.na(v11), !is.na(v21), v11 != "-", v21 != "-") %>%
  distinct(v11, v21) %>%  
  group_by(v11) %>%
  filter(n_distinct(v21) == 1) %>% 
  dplyr::slice(1) %>% 
  ungroup()

gene_labels[, gene := map_unique$v21[match(gene, map_unique$v11)]]
gene_labels <- gene_labels[!is.na(gene)]
gene_labels <- gene_labels[!grepl("LC$", gene)]

fwrite(gene_labels, "choulet_gene_tissue.tsv", sep="\t")

gbm_status <- fread("gene_body_methylation.tsv")
gbm_status <- gbm_status[, c("gene_id", "sample", "gbM_candidate")]
gbm_status <- dcast(gbm_status, gene_id ~ sample, value.var="gbM_candidate")
gbm_status$parent_status <- ifelse(gbm_status$CS == TRUE & gbm_status$P == TRUE, "gbM in both parents",
                                   ifelse(gbm_status$CS == TRUE & gbm_status$P == FALSE, "CS-specific gbM",
                                          ifelse(gbm_status$CS == FALSE & gbm_status$P == TRUE, "P-specific gbM", "non-gbM in both parents"
                                          )
                                   )
)
gbm_status$hybrid_status <- ifelse(gbm_status$CSxP == TRUE, "gbM", "non-gbM")
gbm_status$parent_status <- factor(gbm_status$parent_status, levels=c("gbM in both parents", "CS-specific gbM", "P-specific gbM", "non-gbM in both parents"))
gbm_status$hybrid_status <- factor(gbm_status$hybrid_status, levels=c("gbM", "non-gbM"))
gbm_status <- gbm_status[CS == CSxP & CSxP == P]

choulet_tissue <- read.table(file="choulet_gene_tissue.tsv", header=T)
colnames(choulet_tissue)[2] <- "gene_id"
gbm_status <- inner_join(gbm_status, choulet_tissue, by="gene_id")
gbm_status$label <- gsub("seedling_specific", "Seedling Specific", gbm_status$label)
gbm_status$label <- gsub("broad", "Broad", gbm_status$label)

tab <- table(gbm_status$hybrid_status, gbm_status$label)
n_dt <- as.data.table(colSums(tab), keep.rownames="label")
setnames(n_dt, c("label", "n"))
n_dt[, n_label := paste0("n=", n)]

chisq.test(tab)

gbm_status2 <- as.data.table(gbm_status)[, .N, by=.(label, hybrid_status)]
gbm_status2[, prop := N / sum(N), by=label]

pdf(file="gbm_tissue.pdf", height=2.5, width=4)
ggplot(gbm_status2, aes(x=label, y=prop, fill=hybrid_status)) +
  geom_col(width=0.7) +
  geom_text(data=n_dt, aes(x=label, y=1.05, label=n_label), inherit.aes=FALSE, size=4) +
  scale_fill_manual(values = c("non-gbM" = "grey70", "gbM" = "#0072B2")) +
  scale_y_continuous(labels=scales::percent, limits=c(0, 1.15)) +
  labs(x="", y="Proportion of genes", fill="gbM status") +
  annotate("text", x=1.5, y=1.1, label="***", size=4) +
  theme_classic()
dev.off()

