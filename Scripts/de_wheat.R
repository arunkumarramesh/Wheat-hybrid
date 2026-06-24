setwd("~/Downloads/wheat/")

## get list of core genes
library(dplyr)
library(tidyr)
core <- read.table(file="core_genes.txt",header = T)
core <- core %>%
  pivot_longer(cols = everything(),values_to = "gene") %>%
  separate_rows(gene, sep = ",\\s*") %>%
  mutate(gene = sub("\\.[0-9]+$", "", gene)) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(gene)

iwgsc_refseq_all_correspondances <- read.table("iwgsc_refseq_all_correspondances.csv",header=T)
# convert v1.1 to v2.1
map_unique <- iwgsc_refseq_all_correspondances %>%
  transmute(v11 = `v1.1`, v21 = `v2.1`) %>%
  filter(!is.na(v11), !is.na(v21), v11 != "-", v21 != "-") %>%
  distinct(v11, v21) %>%  
  group_by(v11) %>%
  filter(n_distinct(v21) == 1) %>% 
  dplyr::slice(1) %>% 
  ungroup()

core <- core %>%
  left_join(map_unique, by = c("gene" = "v11"))

## different expression analysis for CS and PAR and hybrids

## first pass, justification for removing PxCS3
library(edgeR)
library(factoextra)
library(pheatmap)

read.counts <- read.table("cs_count.tsv", header = TRUE) 
read.counts <- read.counts[!grepl("LC$", rownames(read.counts)), , drop = FALSE]
read.counts <- read.counts[rownames(read.counts) %in% core$v21,]
sample_info.edger <- factor(c( rep("CS", 3), rep("CSxP", 3), rep("P", 3), rep("PxCS", 3)))
edgeR.DGElist <- DGEList(counts = read.counts, group = sample_info.edger) 
keep <- rowSums( cpm(edgeR.DGElist) >= 2) >= 4
edgeR.DGElist <- edgeR.DGElist[keep,]
edgeR.DGElist$samples$lib.size <- colSums(edgeR.DGElist$counts)
edgeR.DGElist <- calcNormFactors(edgeR.DGElist, method = "TMM")

cpm_log <- cpm(edgeR.DGElist, log = TRUE)
cpm_nolog <- cpm(edgeR.DGElist, log = FALSE)
colnames(cpm_log) <- sub("_.*","",colnames(cpm_log))
colnames(cpm_log) <- sub("_.*","",colnames(cpm_log))
cpm_nolog_relative <- cpm_nolog/rowMeans(cpm_nolog)
colnames(cpm_nolog_relative) <- sub("_.*","",colnames(cpm_nolog_relative))
colnames(cpm_nolog_relative) <- sub("_.*","",colnames(cpm_nolog_relative))
colnames(cpm_nolog) <- sub("_.*","",colnames(cpm_nolog))
colnames(cpm_nolog) <- sub("_.*","",colnames(cpm_nolog))

Group <- edgeR.DGElist$samples[1]
Group <- as.factor(unlist(Group))
cpm_log_forpca <- cpm_log
pca <- prcomp(t(cpm_log_forpca), scale. = TRUE) ## do pca
## plot of pca with groups in ellipses
CSpca <- fviz_pca_ind(pca,
                      col.ind = Group,
                      palette = c("#0072B2", "#E69F00", "#009E73", "#CC79A7"),
                      legend.title = "Genotypes",
                      repel = TRUE,
                      pointshape = 16,
                      pointsize  = 3,
                      mean.point = FALSE, 
                      title = ""
) + guides(color = guide_legend(override.aes = list(shape = 16, size = 3)))

pdf("first_pass_PCA.pdf",height=3.5,width=4.5)
CSpca
dev.off()

## second pass, justification for group PXCS and CSxP, only 10 genes are diff expressed between these two groups
library(edgeR)
library(factoextra)
library(pheatmap)

read.counts <- read.table("cs_count.tsv", header = TRUE) 
read.counts <- read.counts[!grepl("LC$", rownames(read.counts)), , drop = FALSE]
read.counts <- read.counts[rownames(read.counts) %in% core$v21,]
read.counts <- read.counts[1:11]
sample_info.edger <- factor(c( rep("CS", 3), rep("CSxP", 3), rep("P", 3), rep("PxCS", 2)))
edgeR.DGElist <- DGEList(counts = read.counts, group = sample_info.edger) 
keep <- rowSums( cpm(edgeR.DGElist) >= 2) >= 4
edgeR.DGElist <- edgeR.DGElist[keep,]
edgeR.DGElist$samples$lib.size <- colSums(edgeR.DGElist$counts)
edgeR.DGElist <- calcNormFactors(edgeR.DGElist, method = "TMM")

cpm_log <- cpm(edgeR.DGElist, log = TRUE)
cpm_nolog <- cpm(edgeR.DGElist, log = FALSE)
colnames(cpm_log) <- sub("_.*","",colnames(cpm_log))
colnames(cpm_log) <- sub("_.*","",colnames(cpm_log))
cpm_nolog_relative <- cpm_nolog/rowMeans(cpm_nolog)
colnames(cpm_nolog_relative) <- sub("_.*","",colnames(cpm_nolog_relative))
colnames(cpm_nolog_relative) <- sub("_.*","",colnames(cpm_nolog_relative))
colnames(cpm_nolog) <- sub("_.*","",colnames(cpm_nolog))
colnames(cpm_nolog) <- sub("_.*","",colnames(cpm_nolog))

mm <- model.matrix(~0+edgeR.DGElist$samples$group, data = edgeR.DGElist$samples)
colnames(mm) <- levels(edgeR.DGElist$samples$group)
y <- voom(edgeR.DGElist, mm, plot = F)
fit <- lmFit(y, mm)

CSxPvPxCS <- eBayes(contrasts.fit(fit, contrast = c(0, 1, 0, -1)))
top.table <- topTable(CSxPvPxCS, sort.by = "P", n = Inf)
all.CSxPvPxCS <- topTable(CSxPvPxCS, sort.by = "none", n = Inf,p.value=1,lfc=0) 
length(which(top.table$adj.P.Val < 0.05))
sig_genes <- subset(top.table, top.table$adj.P.Val < 0.05) 
write.csv(sig_genes, file = "CSxPvPxCS sig genes.csv")
DGEgenes <- rownames(subset(top.table, top.table$adj.P.Val < 0.05))
mat_DGEgenes <- cpm_nolog_relative[rownames(sig_genes), ]

pdf("poi_genes_CS.pdf",height=3,width=5)
pheatmap(mat_DGEgenes)
dev.off()

library(edgeR)
library(factoextra)
library(ComplexHeatmap)
library(magick)
library(dplyr)
library(ggplot2)
library(forcats)
library(cowplot)
library(goseq)
library(forcats)
library(eulerr)
library(ComplexUpset)

## now read counts

read.counts <- read.table("cs_count.tsv", header = TRUE) ## read counts from feature counts following STAR mapping
read.counts <- read.counts[!grepl("LC$", rownames(read.counts)), , drop = FALSE]
read.counts <- read.counts[rownames(read.counts) %in% core$v21,]
read.counts <- read.counts[1:11]
CSxPvPxCS_sig_genes <- read.csv(file = "CSxPvPxCS sig genes.csv")
read.counts <- read.counts[!rownames(read.counts) %in% CSxPvPxCS_sig_genes$X,]
sample_info.edger <- factor(c( rep("CS", 3), rep("CSxP", 3), rep("P", 3), rep("CSxP", 2))) ### treatment as grouping variables
edgeR.DGElist <- DGEList(counts = read.counts, group = sample_info.edger) ### group read counts by treatment
keep <- rowSums( cpm(edgeR.DGElist) >= 2) >= 4
edgeR.DGElist <- edgeR.DGElist[keep,]
edgeR.DGElist$samples$lib.size <- colSums(edgeR.DGElist$counts)
edgeR.DGElist <- calcNormFactors(edgeR.DGElist, method = "TMM")
mm <- model.matrix(~0+edgeR.DGElist$samples$group, data = edgeR.DGElist$samples)
colnames(mm) <- levels(edgeR.DGElist$samples$group)
y <- voom(edgeR.DGElist, mm, plot = F)
fit <- lmFit(y, mm)

top.table <- topTable(eBayes(contrasts.fit(fit, contrast = c(-0.5, 1, -0.5))), sort.by = "P", n = Inf)
DGEgenes_CS_PvCSxP <- rownames(top.table[top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > 0.58, ])

top.table <- topTable(eBayes(contrasts.fit(fit, contrast = c(-1, 1, 0))), sort.by = "P", n = Inf)
DGEgenes_CS_PvCS <- rownames(top.table[top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > 0.58, ])

top.table <- topTable(eBayes(contrasts.fit(fit, contrast = c(0, 1, -1))), sort.by = "P", n = Inf)
DGEgenes_CS_PvP <- rownames(top.table[top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > 0.58, ])

upset_df <- data.frame(gene = unique(c(unique(na.omit(DGEgenes_CS_PvCSxP)), unique(na.omit(DGEgenes_CS_PvCS)), unique(na.omit(DGEgenes_CS_PvP)))), `Hybrid v. Midparent` = unique(c(unique(na.omit(DGEgenes_CS_PvCSxP)), unique(na.omit(DGEgenes_CS_PvCS)), unique(na.omit(DGEgenes_CS_PvP)))) %in% unique(na.omit(DGEgenes_CS_PvCSxP)), `Hybrid v. CS` = unique(c(unique(na.omit(DGEgenes_CS_PvCSxP)), unique(na.omit(DGEgenes_CS_PvCS)), unique(na.omit(DGEgenes_CS_PvP)))) %in% unique(na.omit(DGEgenes_CS_PvCS)), `Hybrid v. Paragon` = unique(c(unique(na.omit(DGEgenes_CS_PvCSxP)), unique(na.omit(DGEgenes_CS_PvCS)), unique(na.omit(DGEgenes_CS_PvP)))) %in% unique(na.omit(DGEgenes_CS_PvP)), check.names = FALSE)

pdf("Hybrid_vParent_diff.pdf",height=3,width=5.3)
upset(upset_df, c("Hybrid v. Midparent", "Hybrid v. CS", "Hybrid v. Paragon"), set_sizes = FALSE, name = "DE Comparison", base_annotations = list("Intersection size" = intersection_size(counts = TRUE, text_mapping = aes(label = !!upset_text_percentage(digits = 1)))  + labs(y = "Number of genes", title = NULL)))
dev.off()

cpm_log <- cpm(edgeR.DGElist, log = TRUE)
cpm_nolog <- cpm(edgeR.DGElist, log = FALSE)
colnames(cpm_log) <- sub("_.*","",colnames(cpm_log))
colnames(cpm_log) <- sub("_.*","",colnames(cpm_log))
cpm_nolog_relative <- cpm_nolog/rowMeans(cpm_nolog)
colnames(cpm_nolog_relative) <- sub("_.*","",colnames(cpm_nolog_relative))
colnames(cpm_nolog_relative) <- sub("_.*","",colnames(cpm_nolog_relative))
colnames(cpm_nolog) <- sub("_.*","",colnames(cpm_nolog))
colnames(cpm_nolog) <- sub("_.*","",colnames(cpm_nolog))

Group <- edgeR.DGElist$samples[1]
Group <- as.factor(unlist(Group))
cpm_log_forpca <- cpm_log
pca <- prcomp(t(cpm_log_forpca), scale. = TRUE) ## do pca
## plot of pca with groups in ellipses
CSpca <- fviz_pca_ind(pca,
                      col.ind = Group,
                      palette = c("#0072B2", "#E69F00", "#CC79A7"),
                      legend.title = "Genotypes",
                      repel = TRUE,
                      pointshape = 16,
                      pointsize  = 3,
                      mean.point = FALSE, 
                      title = paste("n=",nrow(cpm_log)," genes",sep="")
) + guides(color = guide_legend(override.aes = list(shape = 16, size = 3)))

pdf("CS_PCA.pdf",height=3.5,width=4.5)
CSpca
dev.off()

csscree <- fviz_screeplot(pca, ncp=10,title = "")
csscree

CSvP <- eBayes(contrasts.fit(fit, contrast = c(1, 0, -1))) ## CS upregulated, P downregulated
top.table <- topTable(CSvP, sort.by = "P", n = Inf) ## sort by most significantly DE genes
all.CSvP <- topTable(CSvP, sort.by = "none", n = Inf,p.value=1,lfc=0) ## get all genes with logFC and pvalues
length(which(top.table$adj.P.Val < 0.05)) ## how many significantly DE genes
prop_sig_CSvP <- length(which(top.table$adj.P.Val < 0.05))/nrow(top.table)
length_CSvP <- nrow(top.table)
sig_genes <- subset(top.table, top.table$adj.P.Val < 0.05) ##  significantly DE genes only
write.csv(all.CSvP, file = "CSvP all genes.csv")
DGEgenes <- rownames(subset(top.table, top.table$adj.P.Val < 0.05))
mat_DGEgenes <- cpm_nolog_relative[DGEgenes, ]

pdf("CSvsP_heatmap.pdf",height=3.5,width=4)
Heatmap(mat_DGEgenes, name = "Scaled CPM", show_row_names = FALSE, use_raster = F)
dev.off()

## check if hybrid differ from midparent value
CS_PvCSxP <- eBayes(contrasts.fit(fit, contrast = c(-0.5, 1, -0.5)))
top.table <- topTable(CS_PvCSxP, sort.by = "P", n = Inf) ## sort by most significantly DE genes
all.CS_PvCSxP <- topTable(CS_PvCSxP, sort.by = "none", n = Inf,p.value=1,lfc=0) ## get all genes with logFC and pvalues
length(which(top.table$adj.P.Val < 0.05)) ## how many significantly DE genes
prop_sig_CS_PvCSxP <- length(which(top.table$adj.P.Val < 0.05))/nrow(top.table)
length_CS_PvCSxP <- nrow(top.table)
sig_genes <- subset(top.table, top.table$adj.P.Val < 0.05) ##  significantly DE genes only
DGEgenes <- rownames(subset(top.table, top.table$adj.P.Val < 0.05))
mat_DGEgenes <- cpm_nolog_relative[DGEgenes, ]
write.csv(all.CS_PvCSxP, file = "CS_PvCSxP all genes.csv")

pdf("CS_PvCSxP_heatmap.pdf",height=3.5,width=4)
Heatmap(mat_DGEgenes, name = "Scaled CPM", show_row_names = FALSE, use_raster = F)
dev.off()

nrow(all.CSvP[all.CSvP$adj.P.Val < 0.05 & all.CSvP$logFC > 0.58 & all.CS_PvCSxP$adj.P.Val < 0.05 & all.CS_PvCSxP$logFC > 0.58,])
nrow(all.CSvP[all.CSvP$adj.P.Val < 0.05 & all.CSvP$logFC < -0.58 & all.CS_PvCSxP$adj.P.Val < 0.05 & all.CS_PvCSxP$logFC < -0.58,])
nrow(all.CSvP[all.CSvP$adj.P.Val < 0.05 & all.CSvP$logFC > 0.58,])
nrow(all.CSvP[all.CSvP$adj.P.Val < 0.05 & all.CSvP$logFC < -0.58,])
nrow(all.CS_PvCSxP[all.CS_PvCSxP$adj.P.Val < 0.05 & all.CS_PvCSxP$logFC > 0.58,])
nrow(all.CS_PvCSxP[all.CS_PvCSxP$adj.P.Val < 0.05 & all.CS_PvCSxP$logFC < -0.58,])
nrow(all.CSvP)

## for venn diagram
keep_A_CSvP <- grepl("^TraesCS[0-9]+A", rownames(all.CSvP))
keep_A_PvH  <- grepl("^TraesCS[0-9]+A", rownames(all.CS_PvCSxP))
all.CSvP_A <- all.CSvP[keep_A_CSvP, ]
all.CS_PvCSxP_A <- all.CS_PvCSxP[keep_A_PvH, ]
sig_CSvP_A <- rownames(all.CSvP_A)[all.CSvP_A$adj.P.Val < 0.05 & abs(all.CSvP_A$logFC) > 0.58]
sig_PvH_A <- rownames(all.CS_PvCSxP_A)[all.CS_PvCSxP_A$adj.P.Val < 0.05 & abs(all.CS_PvCSxP_A$logFC) > 0.58]
total_A_genes <- length(intersect(rownames(all.CSvP_A), rownames(all.CS_PvCSxP_A)))

vennfit <- euler(c("CS vs Paragon" = length(sig_CSvP_A),"Parents vs Hybrids" = length(sig_PvH_A),"CS vs Paragon&Parents vs Hybrids" = length(intersect(sig_CSvP_A, sig_PvH_A))))

pdf("venn_A_subgenome.pdf", height = 3, width = 4)
plot(vennfit,quantities = TRUE,legend = TRUE,main = paste0("A=", total_A_genes,sep=""),fills = c("#0072B2", "#E69F00"))
dev.off()

keep_B_CSvP <- grepl("^TraesCS[0-9]+B", rownames(all.CSvP))
keep_B_PvH  <- grepl("^TraesCS[0-9]+B", rownames(all.CS_PvCSxP))
all.CSvP_B <- all.CSvP[keep_B_CSvP, ]
all.CS_PvCSxP_B <- all.CS_PvCSxP[keep_B_PvH, ]
sig_CSvP_B <- rownames(all.CSvP_B)[all.CSvP_B$adj.P.Val < 0.05 & abs(all.CSvP_B$logFC) > 0.58]
sig_PvH_B <- rownames(all.CS_PvCSxP_B)[all.CS_PvCSxP_B$adj.P.Val < 0.05 & abs(all.CS_PvCSxP_B$logFC) > 0.58]
total_B_genes <- length(intersect(rownames(all.CSvP_B), rownames(all.CS_PvCSxP_B)))

vennfit <- euler(c("CS vs Paragon" = length(sig_CSvP_B),"Parents vs Hybrids" = length(sig_PvH_B),"CS vs Paragon&Parents vs Hybrids" = length(intersect(sig_CSvP_B, sig_PvH_B))))

pdf("venn_B_subgenome.pdf", height = 3, width = 4)
plot(vennfit,quantities = TRUE,legend = TRUE,main = paste0("B=", total_B_genes,sep=""),fills = c("#0072B2", "#E69F00"))
dev.off()

keep_D_CSvP <- grepl("^TraesCS[0-9]+D", rownames(all.CSvP))
keep_D_PvH  <- grepl("^TraesCS[0-9]+D", rownames(all.CS_PvCSxP))
all.CSvP_D <- all.CSvP[keep_D_CSvP, ]
all.CS_PvCSxP_D <- all.CS_PvCSxP[keep_D_PvH, ]
sig_CSvP_D <- rownames(all.CSvP_D)[all.CSvP_D$adj.P.Val < 0.05 & abs(all.CSvP_D$logFC) > 0.58]
sig_PvH_D <- rownames(all.CS_PvCSxP_D)[all.CS_PvCSxP_D$adj.P.Val < 0.05 & abs(all.CS_PvCSxP_D$logFC) > 0.58]
total_D_genes <- length(intersect(rownames(all.CSvP_D), rownames(all.CS_PvCSxP_D)))

vennfit <- euler(c("CS vs Paragon" = length(sig_CSvP_D),"Parents vs Hybrids" = length(sig_PvH_D),"CS vs Paragon&Parents vs Hybrids" = length(intersect(sig_CSvP_D, sig_PvH_D))))

pdf("venn_D_subgenome.pdf", height = 3, width = 4)
plot(vennfit,quantities = TRUE,legend = TRUE,main = paste0("D=", total_D_genes,sep=""),fills = c("#0072B2", "#E69F00"))
dev.off()


## venn diagram, tissue specific
choulet_tissue <- read.table(file="choulet_gene_tissue.tsv", header=T)
colnames(choulet_tissue)[2] <- "gene_id"
all.CSvP_seedling <- all.CSvP[rownames(all.CSvP) %in% choulet_tissue[choulet_tissue$label %in% "seedling_specific",]$gene_id,]
all.CS_PvCSxP_seedling <- all.CS_PvCSxP[rownames(all.CS_PvCSxP) %in% choulet_tissue[choulet_tissue$label %in% "seedling_specific",]$gene_id,]

sig_CSvP <- rownames(all.CSvP_seedling)[all.CSvP_seedling$adj.P.Val < 0.05 & abs(all.CSvP_seedling$logFC) > 0.58]
sig_PvH <- rownames(all.CS_PvCSxP_seedling)[all.CS_PvCSxP_seedling$adj.P.Val < 0.05 & abs(all.CS_PvCSxP_seedling$logFC) > 0.58]
total_genes <- length(intersect(rownames(all.CSvP_seedling), rownames(all.CS_PvCSxP_seedling)))
vennfit <- euler(c("CS vs Paragon" = length(sig_CSvP), "Parents vs Hybrids" = length(sig_PvH), "CS vs Paragon&Parents vs Hybrids" = length(intersect(sig_CSvP, sig_PvH))))

pdf("venn_seedling_subgenome.pdf", height = 3, width = 4)
plot(vennfit,quantities = TRUE,legend = TRUE,main = paste0("Seedling\nspecific:", total_genes,sep=""),fills = c("#0072B2", "#E69F00"))
dev.off()


all.CSvP_broad <- all.CSvP[rownames(all.CSvP) %in% choulet_tissue[choulet_tissue$label %in% "broad",]$gene_id,]
all.CS_PvCSxP_broad <- all.CS_PvCSxP[rownames(all.CS_PvCSxP) %in% choulet_tissue[choulet_tissue$label %in% "broad",]$gene_id,]

sig_CSvP <- rownames(all.CSvP_broad)[all.CSvP_broad$adj.P.Val < 0.05 & abs(all.CSvP_broad$logFC) > 0.58]
sig_PvH <- rownames(all.CS_PvCSxP_broad)[all.CS_PvCSxP_broad$adj.P.Val < 0.05 & abs(all.CS_PvCSxP_broad$logFC) > 0.58]
total_genes <- length(intersect(rownames(all.CSvP_broad), rownames(all.CS_PvCSxP_broad)))
vennfit <- euler(c("CS vs Paragon" = length(sig_CSvP), "Parents vs Hybrids" = length(sig_PvH), "CS vs Paragon&Parents vs Hybrids" = length(intersect(sig_CSvP, sig_PvH))))

pdf("venn_broad_subgenome.pdf", height = 3, width = 4)
plot(vennfit,quantities = TRUE,legend = TRUE,main = paste0("Broad:", total_genes,sep=""),fills = c("#0072B2", "#E69F00"))
dev.off()


