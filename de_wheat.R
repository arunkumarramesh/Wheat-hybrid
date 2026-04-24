setwd("~/Downloads/wheat/")

## different expression analysis for CS and PAR and hybrids

## first pass, justification for removing PxCS3
library(edgeR)
library(factoextra)
library(pheatmap)
read.counts <- read.table("cs_count.tsv", header = TRUE) 
read.counts <- read.counts[!grepl("LC$", rownames(read.counts)), , drop = FALSE]
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
                      col.ind = Group, # color by groups
                      palette = c("#0072B2", "#E69F00", "#009E73", "#CC79A7"),
                      legend.title = "Genotypes",
                      repel = TRUE,
                      pointshape = 16,                                  # filled circles
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
library(EnhancedVolcano)
library(goseq)
library(forcats)
library(eulerr)

## file for GO analysis
gene_length <- read.table(file="gene.gff3")
gene_length$V6 <- abs(gene_length$V4 - gene_length$V3)
gene_length <- gene_length[5:6]

bp_go <- read.csv('BP.csv',header = F)
mf_go <- read.csv('MF.csv',header = F)
cc_go <- read.csv('CC.csv',header = F)
all_go <- rbind(bp_go,mf_go,cc_go)
colnames(all_go) <- c('Gene','GO_term','Term')

## now read counts

read.counts <- read.table("cs_count.tsv", header = TRUE) ## read counts from feature counts following STAR mapping
read.counts <- read.counts[!grepl("LC$", rownames(read.counts)), , drop = FALSE]
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
                      col.ind = Group, # color by groups
                      palette = c("#0072B2", "#E69F00", "#CC79A7"),
                      legend.title = "Genotypes",
                      repel = TRUE,
                      pointshape = 16,                                  # filled circles
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

all.CSvP_volcano <- EnhancedVolcano(all.CSvP, lab = "", x = 'logFC', y = 'adj.P.Val', pCutoff = 0.05, FCcutoff = 1, title = 'CS v P', subtitle = "",  ylim = c(0,6), xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(-Log[10]~ 'FDR'))

all_go_subset <- subset(all_go, Gene %in% rownames(all.CSvP))
degenes <- rownames(all.CSvP)
names(degenes) <- degenes
degenes <- setNames(rep(0L, length(degenes)), names(degenes))
degenes[names(degenes) %in% rownames(all.CSvP[(all.CSvP$adj.P.Val < 0.05) & (all.CSvP$logFC > 0.58),])] <- 1
gene_length_subsest <- gene_length[gene_length$V5 %in% names(degenes),]
len_lookup <- setNames(as.numeric(gene_length$V6), gene_length$V5)
degenes_len <- len_lookup[names(degenes)]
pwf = nullp(degenes, bias.data = degenes_len, plot.fit = TRUE)
GO.wall_CS = goseq(pwf, gene2cat = all_go_subset[1:2])
GO.wall_CS$over_rep_padj=p.adjust(GO.wall_CS$over_represented_pvalue, method="BH")
GO.wall_CS$under_rep_padj=p.adjust(GO.wall_CS$under_represented_pvalue, method="BH")
GO.wall_CS_over <- GO.wall_CS[GO.wall_CS$over_rep_padj < 0.05,]
GO.wall_CS_under <- GO.wall_CS[GO.wall_CS$under_rep_padj < 0.05,]

all_go_subset <- subset(all_go, Gene %in% rownames(all.CSvP))
degenes <- rownames(all.CSvP)
names(degenes) <- degenes
degenes <- setNames(rep(0L, length(degenes)), names(degenes))
degenes[names(degenes) %in% rownames(all.CSvP[(all.CSvP$adj.P.Val < 0.05) & (all.CSvP$logFC < -0.58),])] <- 1
gene_length_subsest <- gene_length[gene_length$V5 %in% names(degenes),]
len_lookup <- setNames(as.numeric(gene_length$V6), gene_length$V5)
degenes_len <- len_lookup[names(degenes)]
pwf = nullp(degenes, bias.data = degenes_len, plot.fit = TRUE)
GO.wall_PAR = goseq(pwf, gene2cat = all_go_subset[1:2])
GO.wall_PAR$over_rep_padj=p.adjust(GO.wall_PAR$over_represented_pvalue, method="BH")
GO.wall_PAR$under_rep_padj=p.adjust(GO.wall_PAR$under_represented_pvalue, method="BH")
GO.wall_PAR_over <- GO.wall_PAR[GO.wall_PAR$over_rep_padj < 0.05,]
GO.wall_PAR_under <- GO.wall_PAR[GO.wall_PAR$under_rep_padj < 0.05,]

# Over-represented in CS
go_over_cs <- GO.wall_CS_over %>%
  transmute(term_label = paste0(term, " (", category, ")"), ontology, numDEInCat, numInCat, padj = over_rep_padj) %>%
  filter(!is.na(padj)) %>%
  mutate(log10p = -log10(pmax(padj, .Machine$double.xmin)), gene_ratio = numDEInCat / numInCat) %>%
  arrange(desc(log10p)) %>%
  slice_head(n = 12)

p_cs_over <- ggplot(go_over_cs,aes(x = log10p, y = fct_reorder(term_label, log10p), size = numDEInCat, color = ontology)) +
  geom_point() +
  labs(title = "Over-represented in Chinese Spring", x = expression(-log[10]("adjusted p-value")), y = NULL, size = "# DE genes", color = "Ontology") +
  theme_minimal(base_size = 12)

# Over-represented in PAR
go_over_par <- GO.wall_PAR_over %>%
  transmute(term_label = paste0(term, " (", category, ")"), ontology, numDEInCat, numInCat, padj = over_rep_padj) %>%
  filter(!is.na(padj)) %>%
  mutate(log10p = -log10(pmax(padj, .Machine$double.xmin)),  gene_ratio = numDEInCat / numInCat) %>%
  arrange(desc(log10p)) %>%
  slice_head(n = 12)

p_par_over <- ggplot(go_over_par, aes(x = log10p, y = fct_reorder(term_label, log10p), size = numDEInCat, color = ontology)) +
  geom_point() +
  labs(title = "Over-represented in Paragon", x = expression(-log[10]("adjusted p-value")), y = NULL, size = "# DE genes", color = "Ontology") +
  theme_minimal(base_size = 12)

pdf("CSvP_GO.pdf",height=5,width=6.5)
plot_grid(p_cs_over,p_par_over,ncol=1,rel_heights = c(0.6,1),labels="AUTO")
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

#pdf("CS_PvCSxP_heatmap.pdf",height=3.5,width=4)
#Heatmap(mat_DGEgenes, name = "Scaled CPM", show_row_names = FALSE, use_raster = F)
#dev.off()

nrow(all.CSvP[all.CSvP$adj.P.Val < 0.05 & all.CSvP > 0.58 & all.CS_PvCSxP$adj.P.Val < 0.05 & all.CS_PvCSxP$logFC > 0.58,])
nrow(all.CSvP[all.CSvP$adj.P.Val < 0.05 & all.CSvP < -0.58 & all.CS_PvCSxP$adj.P.Val < 0.05 & all.CS_PvCSxP$logFC < -0.58,])
nrow(all.CSvP[all.CSvP$adj.P.Val < 0.05 & all.CSvP > 0.58,])
nrow(all.CSvP[all.CSvP$adj.P.Val < 0.05 & all.CSvP < -0.58,])
nrow(all.CS_PvCSxP[all.CS_PvCSxP$adj.P.Val < 0.05 & all.CS_PvCSxP$logFC > 0.58,])
nrow(all.CS_PvCSxP[all.CS_PvCSxP$adj.P.Val < 0.05 & all.CS_PvCSxP$logFC < -0.58,])
nrow(all.CSvP)

## for venn diagram
vennfit <- euler(c(
  "CS vs Paragon" = nrow(all.CSvP[all.CSvP$adj.P.Val < 0.05 & abs(all.CSvP$logFC) > 0.58,]),
  "Parents vs Hybrids" = nrow(all.CS_PvCSxP[all.CS_PvCSxP$adj.P.Val < 0.05 & abs(all.CS_PvCSxP$logFC) > 0.58,]),
  "CS vs Paragon&Parents vs Hybrids" = nrow(all.CSvP[all.CSvP$adj.P.Val < 0.05 & abs(all.CSvP$logFC) > 0.58 & rownames(all.CSvP) %in% rownames(all.CS_PvCSxP[all.CS_PvCSxP$adj.P.Val < 0.05 & abs(all.CS_PvCSxP$logFC) > 0.58,]),])
))
pdf("venn.pdf",height=3,width = 4)
plot(vennfit, quantities = TRUE, legend = TRUE, main = "",fills = c("#0072B2", "#E69F00"))
dev.off()

all.CS_PvCSxP_volcano <- EnhancedVolcano(all.CS_PvCSxP, lab = "", x = 'logFC', y = 'adj.P.Val', pCutoff = 0.05, FCcutoff = 1, title = 'Hybrids v Parents', subtitle = "", ylim = c(0,6), xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(-Log[10]~ 'FDR'))

pdf("volcano_cs.pdf",height=10,width = 7)
plot_grid(all.CSvP_volcano,all.CS_PvCSxP_volcano,ncol=1,labels="AUTO")
dev.off()


### now using par reference

library(edgeR)
library(factoextra)
library(pheatmap)
library(ComplexHeatmap)
library(magick)
library(EnhancedVolcano)
library(cowplot)
library(eulerr)

## first pass to get genes that different in direction in hybrids

read.counts <- read.table("par_count.tsv", header = TRUE) ## read counts from feature counts following STAR mapping
read.counts <- read.counts[1:11]
sample_info.edger <- factor(c( rep("CS", 3), rep("CSxP", 3), rep("P", 3), rep("PxCS", 2))) ### treatment as grouping variables
edgeR.DGElist <- DGEList(counts = read.counts, group = sample_info.edger) ### group read counts by treatment
keep <- rowSums( cpm(edgeR.DGElist) >= 2) >= 4
edgeR.DGElist <- edgeR.DGElist[keep,]
edgeR.DGElist$samples$lib.size <- colSums(edgeR.DGElist$counts)
edgeR.DGElist <- calcNormFactors(edgeR.DGElist, method = "TMM")
mm <- model.matrix(~0+edgeR.DGElist$samples$group, data = edgeR.DGElist$samples)
colnames(mm) <- levels(edgeR.DGElist$samples$group)
y <- voom(edgeR.DGElist, mm, plot = T)
fit <- lmFit(y, mm)

CSxPvPxCS <- eBayes(contrasts.fit(fit, contrast = c(0, 1, 0, -1)))
top.table <- topTable(CSxPvPxCS, sort.by = "P", n = Inf)
all.CSxPvPxCS <- topTable(CSxPvPxCS, sort.by = "none", n = Inf,p.value=1,lfc=0)
length(which(top.table$adj.P.Val < 0.05)) 
length(rownames(top.table) %in% SingleCopyOrthologues$V2)
sig_genes <- subset(top.table, top.table$adj.P.Val < 0.05) 
write.csv(sig_genes, file = "CSxPvPxCS sig genes PAR.csv")

cpm_log <- cpm(edgeR.DGElist, log = TRUE)
cpm_nolog <- cpm(edgeR.DGElist, log = FALSE)
colnames(cpm_log) <- sub("_.*","",colnames(cpm_log))
colnames(cpm_log) <- sub("_.*","",colnames(cpm_log))
cpm_nolog_relative <- cpm_nolog/rowMeans(cpm_nolog)
colnames(cpm_nolog_relative) <- sub("_.*","",colnames(cpm_nolog_relative))
colnames(cpm_nolog_relative) <- sub("_.*","",colnames(cpm_nolog_relative))
colnames(cpm_nolog) <- sub("_.*","",colnames(cpm_nolog))
colnames(cpm_nolog) <- sub("_.*","",colnames(cpm_nolog))

DGEgenes <- rownames(subset(top.table, top.table$adj.P.Val < 0.05))
mat_DGEgenes <- cpm_nolog_relative[DGEgenes, ]

pdf("poi_genes_PAR.pdf",height=6,width=5)
Heatmap(mat_DGEgenes, name = "CPM", use_raster = F,column_labels = DGEgenes)
dev.off()

## second pass, now running analyses

read.counts <- read.table("par_count.tsv", header = TRUE)
read.counts <- read.counts[1:11]
CSxPvPxCS_sig_genes <- read.csv(file = "CSxPvPxCS sig genes PAR.csv")
read.counts <- read.counts[!rownames(read.counts) %in% CSxPvPxCS_sig_genes$X,]
sample_info.edger <- factor(c( rep("CS", 3), rep("CSxP", 3), rep("P", 3), rep("CSxP", 2)))
edgeR.DGElist <- DGEList(counts = read.counts, group = sample_info.edger)
keep <- rowSums( cpm(edgeR.DGElist) >= 2) >= 4
edgeR.DGElist <- edgeR.DGElist[keep,]
edgeR.DGElist$samples$lib.size <- colSums(edgeR.DGElist$counts)
edgeR.DGElist <- calcNormFactors(edgeR.DGElist, method = "TMM")
mm <- model.matrix(~0+edgeR.DGElist$samples$group, data = edgeR.DGElist$samples)
colnames(mm) <- levels(edgeR.DGElist$samples$group)
y <- voom(edgeR.DGElist, mm, plot = T)
fit <- lmFit(y, mm)
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
PARpca <- fviz_pca_ind(pca,
                       col.ind = Group, # color by groups
                       palette = c("#0072B2","#E69F00","#CC79A7"),
                       legend.title = "Genotype",
                       repel = TRUE,
                       pointshape = 16,                                  # filled circles
                       pointsize  = 3,
                       mean.point = FALSE, 
                       title = paste("n=",nrow(cpm_log)," genes",sep="")
) + guides(color = guide_legend(override.aes = list(shape = 16, size = 3)))

pdf("PAR_PCA.pdf",height=3.5,width=4.5)
PARpca
dev.off()

parscree <- fviz_screeplot(pca, ncp=10,title = "CS")
parscree

CSvP <- eBayes(contrasts.fit(fit, contrast = c(1, 0, -1)))
top.table <- topTable(CSvP, sort.by = "P", n = Inf)
all.CSvP <- topTable(CSvP, sort.by = "none", n = Inf,p.value=1,lfc=0)
length(which(top.table$adj.P.Val < 0.05))
prop_sig_CSvP_PAR <- length(which(top.table$adj.P.Val < 0.05))/nrow(top.table)
length_CSvP_PAR <- nrow(top.table)
sig_genes <- subset(top.table, top.table$adj.P.Val < 0.05)
write.csv(all.CSvP, file = "CSvP all genes_PAR.csv")

all.CSvP_volcano <- EnhancedVolcano(all.CSvP, lab = "", x = 'logFC', y = 'adj.P.Val', pCutoff = 0.05, FCcutoff = 1, title = 'CS v P', subtitle = "",  ylim = c(0,6), xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(-Log[10]~ 'FDR'))

DGEgenes <- rownames(subset(top.table, top.table$adj.P.Val < 0.05))
mat_DGEgenes <- cpm_nolog_relative[DGEgenes, ]

pdf("CSvsP_heatmap_PAR.pdf",height=3.5,width=4)
Heatmap(mat_DGEgenes, name = "Scaled CPM", show_row_names = FALSE, use_raster = F)
dev.off()

## check if hybrid differ from midparent value
CS_PvCSxP <- eBayes(contrasts.fit(fit, contrast = c(-0.5, 1, -0.5)))
top.table <- topTable(CS_PvCSxP, sort.by = "P", n = Inf)
all.CS_PvCSxP <- topTable(CS_PvCSxP, sort.by = "none", n = Inf,p.value=1,lfc=0) 
length(which(top.table$adj.P.Val < 0.05))
prop_sig_CS_PvCSxP_PAR <- length(which(top.table$adj.P.Val < 0.05))/nrow(top.table)
length_CS_PvCSxP_PAR <- nrow(top.table)
sig_genes <- subset(top.table, top.table$adj.P.Val < 0.05)
write.csv(all.CS_PvCSxP, file = "CS_PvCSxP all genes_PAR.csv")

all.CS_PvCSxP_volcano <- EnhancedVolcano(all.CS_PvCSxP, lab = "", x = 'logFC', y = 'adj.P.Val', pCutoff = 0.05, FCcutoff = 1, title = 'Hybrids v Parents', subtitle = "", ylim = c(0,6), xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(-Log[10]~ 'FDR'))

DGEgenes <- rownames(subset(top.table, top.table$adj.P.Val < 0.05))
mat_DGEgenes <- cpm_nolog_relative[DGEgenes, ]

#pdf("CS_PvCSxP_heatmap_par.pdf",height=3.5,width=4)
#Heatmap(mat_DGEgenes, name = "Scaled CPM", show_row_names = FALSE, use_raster = F)
#dev.off()

## for venn diagram

nrow(all.CSvP[all.CSvP$adj.P.Val < 0.05 & all.CSvP > 0.58 & all.CS_PvCSxP$adj.P.Val < 0.05 & all.CS_PvCSxP$logFC > 0.58,])
nrow(all.CSvP[all.CSvP$adj.P.Val < 0.05 & all.CSvP < -0.58 & all.CS_PvCSxP$adj.P.Val < 0.05 & all.CS_PvCSxP$logFC < -0.58,])
nrow(all.CSvP[all.CSvP$adj.P.Val < 0.05 & all.CSvP > 0.58,])
nrow(all.CSvP[all.CSvP$adj.P.Val < 0.05 & all.CSvP < -0.58,])
nrow(all.CS_PvCSxP[all.CS_PvCSxP$adj.P.Val < 0.05 & all.CS_PvCSxP$logFC > 0.58,])
nrow(all.CS_PvCSxP[all.CS_PvCSxP$adj.P.Val < 0.05 & all.CS_PvCSxP$logFC < -0.58,])
nrow(all.CSvP)

vennfit <- euler(c(
  "CS vs Paragon" = nrow(all.CSvP[all.CSvP$adj.P.Val < 0.05 & abs(all.CSvP$logFC) > 0.58,]),
  "Parents vs Hybrids" = nrow(all.CS_PvCSxP[all.CS_PvCSxP$adj.P.Val < 0.05 & abs(all.CS_PvCSxP$logFC) > 0.58,]),
  "CS vs Paragon&Parents vs Hybrids" = nrow(all.CSvP[all.CSvP$adj.P.Val < 0.05 & abs(all.CSvP$logFC) > 0.58 & rownames(all.CSvP) %in% rownames(all.CS_PvCSxP[all.CS_PvCSxP$adj.P.Val < 0.05 & abs(all.CS_PvCSxP$logFC) > 0.58,]),])
))
pdf("venn_par.pdf",height=3,width = 4)
plot(vennfit, quantities = TRUE, legend = TRUE, main = "",fills = c("#0072B2", "#E69F00"))
dev.off()

pdf("volcano_par.pdf",height=10,width = 7)
plot_grid(all.CSvP_volcano,all.CS_PvCSxP_volcano,ncol=1,labels="AUTO")
dev.off()

## compare DE from two ref genomes
library(ggpubr)
library(scales)
library(dplyr)
library(edgeR)

## files for 1:1 mapping
SingleCopyOrthologues <- read.table(file="SingleCopyOrthologues_matrix.tsv")
SingleCopyOrthologues$V3 <- gsub("\\..*","",SingleCopyOrthologues$V3 )
SingleCopyOrthologues$V2 <- sub("^([^.]*\\.[^.]*)\\..*$", "\\1", SingleCopyOrthologues$V2)

SingleCopyOrthologues_mod <- SingleCopyOrthologues[2:3]
SingleCopyOrthologues_mod <- SingleCopyOrthologues_mod %>%
  add_count(V2, name = "nV2") %>%
  add_count(V3, name = "nV3") %>%
  filter(nV2 == 1, nV3 == 1) %>%
  select(-nV2, -nV3)
colnames(SingleCopyOrthologues_mod) <- c("PAR_id","CS_id")

## with cs reference

read.counts <- read.table("cs_count.tsv", header = TRUE)
read.counts <- read.counts[!grepl("LC$", rownames(read.counts)), , drop = FALSE]
read.counts <- read.counts[1:11]
CSxPvPxCS_sig_genes <- read.csv(file = "CSxPvPxCS sig genes.csv")
read.counts <- read.counts[!rownames(read.counts) %in% CSxPvPxCS_sig_genes$X,]
read.counts <- read.counts[rownames(read.counts) %in% SingleCopyOrthologues_mod$CS_id,]
sample_info.edger <- factor(c( rep("CS", 3), rep("CSxP", 3), rep("P", 3), rep("CSxP", 2)))
edgeR.DGElist <- DGEList(counts = read.counts, group = sample_info.edger)
keep <- rowSums( cpm(edgeR.DGElist) >= 2) >= 4
edgeR.DGElist <- edgeR.DGElist[keep,]
edgeR.DGElist$samples$lib.size <- colSums(edgeR.DGElist$counts)
edgeR.DGElist <- calcNormFactors(edgeR.DGElist, method = "TMM")
mm <- model.matrix(~0+edgeR.DGElist$samples$group, data = edgeR.DGElist$samples)
colnames(mm) <- levels(edgeR.DGElist$samples$group)
y <- voom(edgeR.DGElist, mm, plot = F)
fit <- lmFit(y, mm)

## compare parents
CSvP <- eBayes(contrasts.fit(fit, contrast = c(1, 0, -1)))
top.table <- topTable(CSvP, sort.by = "P", n = Inf) 
all.CSvP <- topTable(CSvP, sort.by = "none", n = Inf,p.value=1,lfc=0) 
length(which(top.table$adj.P.Val < 0.05))
colnames(all.CSvP) <- gsub("$","_CS_ref", colnames(all.CSvP))
all.CSvP$CS_id <- rownames(all.CSvP)
all.CSvP <- inner_join(all.CSvP,SingleCopyOrthologues_mod,by="CS_id")

## check if hybrid differ from midparent value
CS_PvCSxP <- eBayes(contrasts.fit(fit, contrast = c(-0.5, 1, -0.5)))
top.table <- topTable(CS_PvCSxP, sort.by = "P", n = Inf) ## sort by most significantly DE genes
all.CS_PvCSxP <- topTable(CS_PvCSxP, sort.by = "none", n = Inf,p.value=1,lfc=0) ## get all genes with logFC and pvalues
length(which(top.table$adj.P.Val < 0.05)) ## how many significantly DE genes
colnames(all.CS_PvCSxP) <- gsub("$","_CS_ref", colnames(all.CS_PvCSxP))
all.CS_PvCSxP$CS_id <- rownames(all.CS_PvCSxP)
all.CS_PvCSxP <- inner_join(all.CS_PvCSxP,SingleCopyOrthologues_mod,by="CS_id")

## with par reference
read.counts <- read.table("par_count.tsv", header = TRUE) ## read counts from feature counts following STAR mapping
read.counts <- read.counts[1:11]
CSxPvPxCS_sig_genes <- read.csv(file = "CSxPvPxCS sig genes PAR.csv")
read.counts <- read.counts[!rownames(read.counts) %in% CSxPvPxCS_sig_genes$X,]
read.counts <- read.counts[rownames(read.counts) %in% SingleCopyOrthologues_mod$PAR_id,]
sample_info.edger <- factor(c( rep("CS", 3), rep("CSxP", 3), rep("P", 3), rep("CSxP", 2))) ### treatment as grouping variables
edgeR.DGElist <- DGEList(counts = read.counts, group = sample_info.edger) ### group read counts by treatment
keep <- rowSums( cpm(edgeR.DGElist) >= 2) >= 4
edgeR.DGElist <- edgeR.DGElist[keep,]
edgeR.DGElist$samples$lib.size <- colSums(edgeR.DGElist$counts)
edgeR.DGElist <- calcNormFactors(edgeR.DGElist, method = "TMM")
mm <- model.matrix(~0+edgeR.DGElist$samples$group, data = edgeR.DGElist$samples)
colnames(mm) <- levels(edgeR.DGElist$samples$group)
y <- voom(edgeR.DGElist, mm, plot = T)
fit <- lmFit(y, mm)

## compare parents
CSvP <- eBayes(contrasts.fit(fit, contrast = c(1, 0, -1)))
top.table <- topTable(CSvP, sort.by = "P", n = Inf) ## sort by most significantly DE genes
all.CSvP_par <- topTable(CSvP, sort.by = "none", n = Inf,p.value=1,lfc=0) ## get all genes with logFC and pvalues
length(which(top.table$adj.P.Val < 0.05)) ## how many significantly DE genes
colnames(all.CSvP_par) <- gsub("$","_PAR_ref", colnames(all.CSvP_par))
all.CSvP_par$PAR_id <- rownames(all.CSvP_par)
CSvP_all_both <- inner_join(all.CSvP,all.CSvP_par,by="PAR_id")
write.csv(CSvP_all_both,file="CSvP_all_both.csv",row.names = F)

## check if hybrid differ from midparent value
CS_PvCSxP <- eBayes(contrasts.fit(fit, contrast = c(-0.5, 1, -0.5)))
top.table <- topTable(CS_PvCSxP, sort.by = "P", n = Inf) ## sort by most significantly DE genes
all.CS_PvCSxP_par <- topTable(CS_PvCSxP, sort.by = "none", n = Inf,p.value=1,lfc=0) ## get all genes with logFC and pvalues
colnames(all.CS_PvCSxP_par) <- gsub("$","_PAR_ref", colnames(all.CS_PvCSxP_par))
length(which(top.table$adj.P.Val < 0.05)) ## how many significantly DE genes
all.CS_PvCSxP_par$PAR_id <- rownames(all.CS_PvCSxP_par)
CS_PvCSxP_all_both <- inner_join(all.CS_PvCSxP,all.CS_PvCSxP_par,by="PAR_id")
write.csv(CS_PvCSxP_all_both,file="CS_PvCSxP_all_both.csv",row.names = F)

## now perform comparisons for parents

p <- c(
  PAR_ref_1to1 = sum(CSvP_all_both$adj.P.Val_PAR_ref < 0.05, na.rm = TRUE) / nrow(CSvP_all_both),
  CS_ref_1to1  = sum(CSvP_all_both$adj.P.Val_CS_ref  < 0.05, na.rm = TRUE) / nrow(CSvP_all_both),
  Both_ref_1to1 = sum((CSvP_all_both$adj.P.Val_CS_ref < 0.05) & (CSvP_all_both$adj.P.Val_PAR_ref < 0.05), na.rm = TRUE) / nrow(CSvP_all_both),
  PAR_ref = prop_sig_CSvP_PAR,
  CS_ref = prop_sig_CSvP
)

df <- data.frame(group = names(p), prop = as.numeric(p))
df$total <- c(nrow(CSvP_all_both),nrow(CSvP_all_both),nrow(CSvP_all_both),nrow(read.csv(file="CSvP all genes_PAR.csv")),nrow(read.csv(file="CSvP all genes.csv")))

bias_plot_a <- ggplot(df, aes(x = group, y = prop)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = total), vjust = -0.3, size = 3.5) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, max(df$prop) * 1.15)) +
  labs(x = NULL, y = "Proportion of genes with\nFDR adjusted P < 0.05", title = "") +
  theme_minimal()

bias_plot_b <- ggscatter(data=CSvP_all_both,x="logFC_CS_ref",y="logFC_PAR_ref", add = "reg.line", conf.int = F, cor.coef = TRUE, cor.method = "pearson",shape = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(x = "log2FC using CS reference", y = "log2FC using\nParagon reference")
bias_plot_c <- ggscatter(data=CSvP_all_both,x="AveExpr_CS_ref",y="AveExpr_PAR_ref", add = "reg.line", conf.int = F, cor.coef = TRUE, cor.method = "pearson",shape = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(x = "Mean Expression using CS reference", y = "Mean Expression\nusing Paragon reference")
bias_plot_d <- ggplot(CSvP_all_both %>% filter(P.Value_CS_ref > 0, P.Value_PAR_ref > 0) %>% transmute(x = -log10(P.Value_CS_ref), y = -log10(P.Value_PAR_ref)), aes(x, y)) +
  geom_point(shape = 1) +
  stat_cor(method = "pearson",label.x.npc = "right", label.y.npc = "top",hjust = 1, vjust = 1, size = 4) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(x = "-log10(CS reference p-value)", y = "-log10 (Paragon\nreference p-value)")

pdf(file="CSvP_ref_bias.pdf",height = 6,width = 9)
plot_grid(bias_plot_a,bias_plot_c,bias_plot_b,bias_plot_d,ncol=2)
dev.off()

## now doing it for hybrids vs parents

p <- c(
  PAR_ref_1to1 = sum(CS_PvCSxP_all_both$adj.P.Val_PAR_ref < 0.05, na.rm = TRUE) / nrow(CS_PvCSxP_all_both),
  CS_ref_1to1  = sum(CS_PvCSxP_all_both$adj.P.Val_CS_ref  < 0.05, na.rm = TRUE) / nrow(CS_PvCSxP_all_both),
  Both_ref_1to1 = sum((CS_PvCSxP_all_both$adj.P.Val_CS_ref < 0.05) & (CS_PvCSxP_all_both$adj.P.Val_PAR_ref < 0.05), na.rm = TRUE) / nrow(CS_PvCSxP_all_both),
  PAR_ref = prop_sig_CS_PvCSxP_PAR,
  CS_ref = prop_sig_CS_PvCSxP
)

df <- data.frame(group = names(p), prop = as.numeric(p))
df$total <- c(nrow(CS_PvCSxP_all_both),nrow(CS_PvCSxP_all_both),nrow(CS_PvCSxP_all_both),nrow(read.csv(file="CS_PvCSxP all genes_PAR.csv")),nrow(read.csv(file="CS_PvCSxP all genes.csv")))

bias_plot_a <- ggplot(df, aes(x = group, y = prop)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = total), vjust = -0.3, size = 3.5) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, max(df$prop) * 1.15)) +
  labs(x = NULL, y = "Proportion of genes with\nFDR adjusted P < 0.05", title = "") +
  theme_minimal()

bias_plot_b <- ggscatter(data=CS_PvCSxP_all_both,x="logFC_CS_ref",y="logFC_PAR_ref", add = "reg.line", conf.int = F, cor.coef = TRUE, cor.method = "pearson",shape = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(x = "log2FC using CS reference", y = "log2FC using\nParagon reference")
bias_plot_c <- ggscatter(data=CS_PvCSxP_all_both,x="AveExpr_CS_ref",y="AveExpr_PAR_ref", add = "reg.line", conf.int = F, cor.coef = TRUE, cor.method = "pearson",shape = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(x = "Mean Expression using CS reference", y = "Mean Expression\nusing Paragon reference")
bias_plot_d <- ggplot(CS_PvCSxP_all_both %>% filter(P.Value_CS_ref > 0, P.Value_PAR_ref > 0) %>% transmute(x = -log10(P.Value_CS_ref), y = -log10(P.Value_PAR_ref)), aes(x, y)) +
  geom_point(shape = 1) +
  stat_cor(method = "pearson",label.x.npc = "right", label.y.npc = "top",hjust = 1, vjust = 1, size = 4) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(x = "-log10(CS reference p-value)", y = "-log10 (Paragon\nreference p-value)")

pdf(file="CS_PvCSxP_ref_bias.pdf",height = 6,width = 9)
plot_grid(bias_plot_a,bias_plot_b,bias_plot_d,ncol=2)
dev.off()
