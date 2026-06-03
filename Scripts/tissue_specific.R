library(data.table)
library(ggplot2)

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
  if (length(x_non_na) == 0 || max(x_non_na) == 0) {
    return(NA)
  }
  x <- x / max(x, na.rm=TRUE)
  sum(1 - x, na.rm=TRUE) / (length(x) - 1)
}

wide[,tau := apply(.SD, 1, calc_tau), .SDcols=condition_cols]
wide[,max_TPM := apply(.SD, 1, max, na.rm=TRUE), .SDcols=condition_cols]
wide[,leaf_Z10_fraction_of_max := leaf_Z10_TPM / max_TPM]
wide[,n_conditions_TPM_ge_1 := rowSums(.SD >= 1, na.rm=TRUE), .SDcols=condition_cols]
wide[,frac_conditions_TPM_ge_1 := n_conditions_TPM_ge_1 / length(condition_cols)]

gene_labels <- rbind(
  wide[leaf_Z10_TPM >= 1 & leaf_Z10_reps_TPM_ge_1 >= 1 & tau >= 0.8 & leaf_Z10_fraction_of_max >= 0.8, .(label="seedling_specific", gene)],
  wide[leaf_Z10_TPM >= 1 & leaf_Z10_reps_TPM_ge_1 >= 1 & tau <= 0.4 & frac_conditions_TPM_ge_1 >= 0.7, .(label="broad", gene)]
)

iwgsc_refseq_all_correspondances <- read.table("iwgsc_refseq_all_correspondances.csv",header=T)

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


choulet_tissue <- read.table(file="choulet_gene_tissue.tsv",header=T)

classified_all <- read.csv(file="classified_all.csv")
classified_all_broad <- classified_all[classified_all$gene %in% choulet_tissue[choulet_tissue$label %in% "broad",]$gene,]

counts_all_broad <- as.data.frame(table(classified_all_broad$category))
names(counts_all_broad) <- c("category", "count")
props_all_broad <- prop.table(table(classified_all_broad$category))
df_all_broad <- counts_all_broad %>%
  mutate(prop = as.numeric(props_all_broad[as.character(category)]), category = factor(category, levels = category))

pdf(file="all_classification_broad.pdf",height=3.5,width=4)
ggplot(df_all_broad, aes(x = category, y = prop, fill = category)) +
  geom_col() +
  geom_text(aes(label = percent(prop, accuracy = 0.1)),vjust = -0.3, size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = c("#66C2A5","#E78AC3","#A6D854","#FFD92F","#E5C494"), guide = "none") +
  scale_x_discrete(labels = function(x) {
    x <- gsub("\\b[Cc]is\\b", "<i>cis</i>", x, perl = TRUE)
    x <- gsub("\\b[Tt]rans\\b", "<i>trans</i>", x, perl = TRUE)
    x
  }) +
  labs(x = NULL, y = "Proportion of genes", title = paste("Broad:",nrow(classified_all_broad))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +
  coord_cartesian(clip = "off")
dev.off()


classified_all_seedling_specific <- classified_all[classified_all$gene %in% choulet_tissue[choulet_tissue$label %in% "seedling_specific",]$gene,]

counts_all_seedling_specific <- as.data.frame(table(classified_all_seedling_specific$category))
names(counts_all_seedling_specific) <- c("category", "count")
props_all_seedling_specific  <- prop.table(table(classified_all_seedling_specific$category))
df_all_seedling_specific <- counts_all_seedling_specific %>%
  mutate(prop = as.numeric(props_all_seedling_specific[as.character(category)]), category = factor(category, levels = category))

pdf(file="all_classification_seedling_specific.pdf",height=3.5,width=4)
ggplot(df_all_seedling_specific, aes(x = category, y = prop, fill = category)) +
  geom_col() +
  geom_text(aes(label = percent(prop, accuracy = 0.1)),vjust = -0.3, size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = c("#66C2A5","#E78AC3","#A6D854","#FFD92F","#E5C494"), guide = "none") +
  scale_x_discrete(labels = function(x) {
    x <- gsub("\\b[Cc]is\\b", "<i>cis</i>", x, perl = TRUE)
    x <- gsub("\\b[Tt]rans\\b", "<i>trans</i>", x, perl = TRUE)
    x
  }) +
  labs(x = NULL, y = "Proportion of genes", title = paste("Seedling Specific:",nrow(classified_all_seedling_specific))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +
  coord_cartesian(clip = "off")
dev.off()


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
gbm_status <- gbm_status[CS == CSxP & CSxP == P]
choulet_tissue <- read.table(file="choulet_gene_tissue.tsv",header=T)
colnames(choulet_tissue)[2] <- "gene_id"
gbm_status <- inner_join(gbm_status,choulet_tissue,by="gene_id")
gbm_status$label <- gsub("seedling_specific","Seedling Specific",gbm_status$label)
gbm_status$label <- gsub("broad","Broad",gbm_status$label)
tab <- table(gbm_status$hybrid_status, gbm_status$label)
n_dt <- as.data.table(colSums(tab), keep.rownames="label")
n_dt[, n_label := paste0("n=", V2)]

chisq.test(tab)

pdf(file="gbm_tissue.pdf",height=2.5,width=4)
ggplot(gbm_status2, aes(x=label, y=prop, fill=hybrid_status)) +
  geom_col(width=0.7) +
  geom_text(data=n_dt, aes(x=V1, y=1.05, label=n_label), inherit.aes=FALSE, size=4) +
  scale_fill_manual(values = c("non-gbM" = "grey70", "gbM" = "#0072B2")) +
  scale_y_continuous(labels=scales::percent, limits=c(0, 1.15)) +
  labs(x="", y="Proportion of genes", fill="gbM status") +
  annotate("text", x=1.5, y=1.1, label="***", size=4) +
  theme_classic()
dev.off()

