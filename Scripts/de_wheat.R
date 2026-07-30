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
cpm_nolog_relative <- cpm_nolog/rowMeans(cpm_nolog)
colnames(cpm_nolog_relative) <- sub("_.*","",colnames(cpm_nolog_relative))
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
read.counts <- read.table("cs_count.tsv",header=TRUE)
read.counts <- read.counts[!grepl("LC$",rownames(read.counts)),,drop=FALSE]
read.counts <- read.counts[rownames(read.counts) %in% core$v21,]
read.counts <- read.counts[1:11]

CSxPvPxCS_sig_genes <- read.csv(file="CSxPvPxCS sig genes.csv")
read.counts <- read.counts[!rownames(read.counts) %in% CSxPvPxCS_sig_genes$X,]

sample_info.edger <- factor(c(rep("CS",3),rep("CSxP",3),rep("P",3),rep("CSxP",2)),levels=c("CS","CSxP","P"))

edgeR.DGElist <- DGEList(counts=read.counts,group=sample_info.edger)
keep <- rowSums(cpm(edgeR.DGElist)>=2)>=4
edgeR.DGElist <- edgeR.DGElist[keep,]
edgeR.DGElist$samples$lib.size <- colSums(edgeR.DGElist$counts)
edgeR.DGElist <- calcNormFactors(edgeR.DGElist,method="TMM")

mm <- model.matrix(~0+edgeR.DGElist$samples$group,data=edgeR.DGElist$samples)
colnames(mm) <- levels(edgeR.DGElist$samples$group)

y <- voom(edgeR.DGElist,mm,plot=FALSE)
fit <- lmFit(y,mm)

effective_lib_sizes <- edgeR.DGElist$samples$lib.size*edgeR.DGElist$samples$norm.factors
target_lib_size <- median(effective_lib_sizes)
normalized_counts <- sweep(edgeR.DGElist$counts,2,effective_lib_sizes,"/")*target_lib_size

CS_cols <- which(sample_info.edger=="CS")
CSxP_cols <- which(sample_info.edger=="CSxP")
P_cols <- which(sample_info.edger=="P")

MP_counts <- (normalized_counts[,CS_cols,drop=FALSE]+normalized_counts[,P_cols,drop=FALSE])/2
colnames(MP_counts) <- paste0("MP",1:ncol(MP_counts))

MP_test_counts <- cbind(normalized_counts[,CSxP_cols,drop=FALSE],MP_counts)
MP_group <- factor(c(rep("CSxP",length(CSxP_cols)),rep("MP",ncol(MP_counts))),levels=c("MP","CSxP"))
mm_MP <- model.matrix(~0+MP_group)
colnames(mm_MP) <- levels(MP_group)

y_MP <- voom(MP_test_counts,mm_MP,lib.size=rep(target_lib_size,ncol(MP_test_counts)),plot=FALSE)
fit_MP <- lmFit(y_MP,mm_MP)

top_MP <- topTable(eBayes(contrasts.fit(fit_MP,contrast=c(-1,1))),sort.by="none",n=Inf)
top_HCS <- topTable(eBayes(contrasts.fit(fit,contrast=c(-1,1,0))),sort.by="none",n=Inf)
top_HP <- topTable(eBayes(contrasts.fit(fit,contrast=c(0,1,-1))),sort.by="none",n=Inf)
top_PCS <- topTable(eBayes(contrasts.fit(fit,contrast=c(-1,0,1))),sort.by="none",n=Inf)

genes <- rownames(fit)

top_MP <- top_MP[genes,,drop=FALSE]
top_HCS <- top_HCS[genes,,drop=FALSE]
top_HP <- top_HP[genes,,drop=FALSE]
top_PCS <- top_PCS[genes,,drop=FALSE]

MP_diff <- top_MP$adj.P.Val<0.05 & abs(top_MP$logFC)>1
HCS_diff <- top_HCS$adj.P.Val<0.05 & abs(top_HCS$logFC)>1
HP_diff <- top_HP$adj.P.Val<0.05 & abs(top_HP$logFC)>1
PCS_diff <- top_PCS$adj.P.Val<0.05 & abs(top_PCS$logFC)>1

MP_same <- top_MP$adj.P.Val>=0.05
HCS_same <- top_HCS$adj.P.Val>=0.05
HP_same <- top_HP$adj.P.Val>=0.05
PCS_same <- top_PCS$adj.P.Val>=0.05

DGEgenes_CS_PvCSxP_up <- genes[MP_diff & top_MP$logFC>0]
DGEgenes_CS_PvCS_up <- genes[HCS_diff & top_HCS$logFC>0]
DGEgenes_CS_PvP_up <- genes[HP_diff & top_HP$logFC>0]

upset_genes <- unique(c(DGEgenes_CS_PvCSxP_up,DGEgenes_CS_PvCS_up,DGEgenes_CS_PvP_up))

upset_df <- data.frame(gene=upset_genes,`Hybrid v. Midparent`=upset_genes %in% DGEgenes_CS_PvCSxP_up,`Hybrid v. CS`=upset_genes %in% DGEgenes_CS_PvCS_up,`Hybrid v. Paragon`=upset_genes %in% DGEgenes_CS_PvP_up,check.names=FALSE)

pdf("Hybrid_vParent_diff_up.pdf",height=3,width=5.3)
upset(upset_df,c("Hybrid v. Midparent","Hybrid v. CS","Hybrid v. Paragon"),set_sizes=FALSE,name="DE Comparison",base_annotations=list("Intersection size"=intersection_size(counts=TRUE,text_mapping=aes(label=!!upset_text_percentage(digits=1)))+labs(y="Number of genes",title="Upregulated")))
dev.off()

DGEgenes_CS_PvCSxP_down <- genes[MP_diff & top_MP$logFC<0]
DGEgenes_CS_PvCS_down <- genes[HCS_diff & top_HCS$logFC<0]
DGEgenes_CS_PvP_down <- genes[HP_diff & top_HP$logFC<0]

upset_genes2 <- unique(c(DGEgenes_CS_PvCSxP_down,DGEgenes_CS_PvCS_down,DGEgenes_CS_PvP_down))

upset_df2 <- data.frame(gene=upset_genes2,`Hybrid v. Midparent`=upset_genes2 %in% DGEgenes_CS_PvCSxP_down,`Hybrid v. CS`=upset_genes2 %in% DGEgenes_CS_PvCS_down,`Hybrid v. Paragon`=upset_genes2 %in% DGEgenes_CS_PvP_down,check.names=FALSE)

pdf("Hybrid_vParent_diff_down.pdf",height=3,width=5.3)
upset(upset_df2,c("Hybrid v. Midparent","Hybrid v. CS","Hybrid v. Paragon"),set_sizes=FALSE,name="DE Comparison",base_annotations=list("Intersection size"=intersection_size(counts=TRUE,text_mapping=aes(label=!!upset_text_percentage(digits=1)))+labs(y="Number of genes",title="Downregulated")))
dev.off()

transgressive_high <- MP_diff & top_MP$logFC>0 & HCS_diff & top_HCS$logFC>0 & HP_diff & top_HP$logFC>0
transgressive_low <- MP_diff & top_MP$logFC<0 & HCS_diff & top_HCS$logFC<0 & HP_diff & top_HP$logFC<0
transgressive <- genes[transgressive_high | transgressive_low]

CS_dominant <- PCS_diff & MP_diff & HCS_same & HP_diff & top_HP$logFC*top_PCS$logFC<0 & top_MP$logFC*top_PCS$logFC<0
P_dominant <- PCS_diff & MP_diff & HP_same & HCS_diff & top_HCS$logFC*top_PCS$logFC>0 & top_MP$logFC*top_PCS$logFC>0
dominant <- genes[CS_dominant | P_dominant]

MP_expression <- rowMeans(MP_counts)
CSxP_expression <- rowMeans(normalized_counts[,CSxP_cols,drop=FALSE])
CS_expression <- rowMeans(normalized_counts[,CS_cols,drop=FALSE])
P_expression <- rowMeans(normalized_counts[,P_cols,drop=FALSE])
within_parent_range <- CSxP_expression>=pmin(CS_expression,P_expression) & CSxP_expression<=pmax(CS_expression,P_expression)

additive <- genes[PCS_diff & MP_same & within_parent_range]

nonDE <- genes[MP_same & HCS_same & HP_same & PCS_same]

write.table(transgressive,file="transgressive_genes.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(dominant,file="dominant_genes.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(additive,file="additive_genes.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(nonDE,file="nonDE_genes.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

c(transgressive=length(transgressive),dominant=length(dominant),additive=length(additive),nonDE=length(nonDE))

cpm_log <- cpm(edgeR.DGElist,log=TRUE)
cpm_nolog <- cpm(edgeR.DGElist,log=FALSE)
colnames(cpm_log) <- sub("_.*","",colnames(cpm_log))
cpm_nolog_relative <- cpm_nolog/rowMeans(cpm_nolog)
colnames(cpm_nolog_relative) <- sub("_.*","",colnames(cpm_nolog_relative))
colnames(cpm_nolog) <- sub("_.*","",colnames(cpm_nolog))

Group <- edgeR.DGElist$samples[1]
Group <- as.factor(unlist(Group))
cpm_log_forpca <- cpm_log
pca <- prcomp(t(cpm_log_forpca),scale.=TRUE)

CSpca <- fviz_pca_ind(pca,col.ind=Group,palette=c("#0072B2","#E69F00","#CC79A7"),legend.title="Genotypes",repel=TRUE,pointshape=16,pointsize=3,mean.point=FALSE,title=paste("n=",nrow(cpm_log)," genes",sep=""))+guides(color=guide_legend(override.aes=list(shape=16,size=3)))

pdf("CS_PCA.pdf",height=3.5,width=4.5)
CSpca
dev.off()

CSvP <- eBayes(contrasts.fit(fit,contrast=c(1,0,-1)))
top.table <- topTable(CSvP,sort.by="P",n=Inf)
all.CSvP <- topTable(CSvP,sort.by="none",n=Inf,p.value=1,lfc=0)
length(which(top.table$adj.P.Val<0.05))
prop_sig_CSvP <- length(which(top.table$adj.P.Val<0.05))/nrow(top.table)
length_CSvP <- nrow(top.table)
sig_genes <- subset(top.table,top.table$adj.P.Val<0.05)
write.csv(all.CSvP,file="CSvP all genes.csv")
DGEgenes <- rownames(top.table)[top.table$adj.P.Val<0.05 & abs(top.table$logFC)>1]
mat_DGEgenes <- cpm_nolog_relative[DGEgenes,]

pdf("CSvsP_heatmap.pdf",height=3.5,width=4)
Heatmap(mat_DGEgenes,name="Scaled CPM",show_row_names=FALSE,use_raster=FALSE)
dev.off()

## check if hybrid differ from midparent value
CS_PvCSxP <- eBayes(contrasts.fit(fit_MP,contrast=c(-1,1)))
top.table <- topTable(CS_PvCSxP,sort.by="P",n=Inf)
all.CS_PvCSxP <- topTable(CS_PvCSxP,sort.by="none",n=Inf,p.value=1,lfc=0)
length(which(top.table$adj.P.Val<0.05))
prop_sig_CS_PvCSxP <- length(which(top.table$adj.P.Val<0.05))/nrow(top.table)
length_CS_PvCSxP <- nrow(top.table)
sig_genes <- subset(top.table,top.table$adj.P.Val<0.05)
DGEgenes <- rownames(top.table)[top.table$adj.P.Val<0.05 & abs(top.table$logFC)>1]
mat_DGEgenes <- cpm_nolog_relative[DGEgenes,]
write.csv(all.CS_PvCSxP,file="CS_PvCSxP all genes.csv")

pdf("CS_PvCSxP_heatmap.pdf",height=3.5,width=4)
Heatmap(mat_DGEgenes,name="Scaled CPM",show_row_names=FALSE,use_raster=FALSE)
dev.off()

nrow(all.CSvP[all.CSvP$adj.P.Val<0.05 & all.CSvP$logFC>1 & all.CS_PvCSxP$adj.P.Val<0.05 & all.CS_PvCSxP$logFC>1,])
nrow(all.CSvP[all.CSvP$adj.P.Val<0.05 & all.CSvP$logFC< -1 & all.CS_PvCSxP$adj.P.Val<0.05 & all.CS_PvCSxP$logFC< -1,])
nrow(all.CSvP[all.CSvP$adj.P.Val<0.05 & all.CSvP$logFC>1,])
nrow(all.CSvP[all.CSvP$adj.P.Val<0.05 & all.CSvP$logFC< -1,])
nrow(all.CS_PvCSxP[all.CS_PvCSxP$adj.P.Val<0.05 & all.CS_PvCSxP$logFC>1,])
nrow(all.CS_PvCSxP[all.CS_PvCSxP$adj.P.Val<0.05 & all.CS_PvCSxP$logFC< -1,])
nrow(all.CSvP)

## for venn diagram
keep_A_CSvP <- grepl("^TraesCS[0-9]+A",rownames(all.CSvP))
keep_A_PvH <- grepl("^TraesCS[0-9]+A",rownames(all.CS_PvCSxP))
all.CSvP_A <- all.CSvP[keep_A_CSvP,]
all.CS_PvCSxP_A <- all.CS_PvCSxP[keep_A_PvH,]
sig_CSvP_A <- rownames(all.CSvP_A)[all.CSvP_A$adj.P.Val<0.05 & abs(all.CSvP_A$logFC)>1]
sig_PvH_A <- rownames(all.CS_PvCSxP_A)[all.CS_PvCSxP_A$adj.P.Val<0.05 & abs(all.CS_PvCSxP_A$logFC)>1]
total_A_genes <- length(intersect(rownames(all.CSvP_A),rownames(all.CS_PvCSxP_A)))

vennfit <- euler(list( `CS vs Paragon`=sig_CSvP_A, `Parents vs Hybrids`=sig_PvH_A ))

pdf("venn_A_subgenome.pdf",height=3,width=4)
plot(vennfit,quantities=TRUE,legend=TRUE,main=paste0("A=",total_A_genes,sep=""),fills=c("#0072B2","#E69F00"))
dev.off()

keep_B_CSvP <- grepl("^TraesCS[0-9]+B",rownames(all.CSvP))
keep_B_PvH <- grepl("^TraesCS[0-9]+B",rownames(all.CS_PvCSxP))
all.CSvP_B <- all.CSvP[keep_B_CSvP,]
all.CS_PvCSxP_B <- all.CS_PvCSxP[keep_B_PvH,]
sig_CSvP_B <- rownames(all.CSvP_B)[all.CSvP_B$adj.P.Val<0.05 & abs(all.CSvP_B$logFC)>1]
sig_PvH_B <- rownames(all.CS_PvCSxP_B)[all.CS_PvCSxP_B$adj.P.Val<0.05 & abs(all.CS_PvCSxP_B$logFC)>1]
total_B_genes <- length(intersect(rownames(all.CSvP_B),rownames(all.CS_PvCSxP_B)))

vennfit <- euler(list( `CS vs Paragon`=sig_CSvP_B, `Parents vs Hybrids`=sig_PvH_B ))

pdf("venn_B_subgenome.pdf",height=3,width=4)
plot(vennfit,quantities=TRUE,legend=TRUE,main=paste0("B=",total_B_genes,sep=""),fills=c("#0072B2","#E69F00"))
dev.off()

keep_D_CSvP <- grepl("^TraesCS[0-9]+D",rownames(all.CSvP))
keep_D_PvH <- grepl("^TraesCS[0-9]+D",rownames(all.CS_PvCSxP))
all.CSvP_D <- all.CSvP[keep_D_CSvP,]
all.CS_PvCSxP_D <- all.CS_PvCSxP[keep_D_PvH,]
sig_CSvP_D <- rownames(all.CSvP_D)[all.CSvP_D$adj.P.Val<0.05 & abs(all.CSvP_D$logFC)>1]
sig_PvH_D <- rownames(all.CS_PvCSxP_D)[all.CS_PvCSxP_D$adj.P.Val<0.05 & abs(all.CS_PvCSxP_D$logFC)>1]
total_D_genes <- length(intersect(rownames(all.CSvP_D),rownames(all.CS_PvCSxP_D)))

vennfit <- euler(list( `CS vs Paragon`=sig_CSvP_D, `Parents vs Hybrids`=sig_PvH_D ))

pdf("venn_D_subgenome.pdf",height=3,width=4)
plot(vennfit,quantities=TRUE,legend=TRUE,main=paste0("D=",total_D_genes,sep=""),fills=c("#0072B2","#E69F00"))
dev.off()

tab <- rbind(A=c(DE=length(sig_PvH_A),not_DE=total_A_genes-length(sig_PvH_A)),B=c(DE=length(sig_PvH_B),not_DE=total_B_genes-length(sig_PvH_B)),D=c(DE=length(sig_PvH_D),not_DE=total_D_genes-length(sig_PvH_D)))
chisq.test(tab)
prop.table(tab,margin=1)

tab2 <- rbind(A=c(DE=length(sig_CSvP_A),not_DE=total_A_genes-length(sig_CSvP_A)),B=c(DE=length(sig_CSvP_B),not_DE=total_B_genes-length(sig_CSvP_B)),D=c(DE=length(sig_CSvP_D),not_DE=total_D_genes-length(sig_CSvP_D)))
chisq.test(tab2)
prop.table(tab2,margin=1)

pct_DE_CSvP <- 100*sum(all.CSvP$adj.P.Val<0.05 & abs(all.CSvP$logFC)>1)/nrow(all.CSvP)
pct_DE_HybridvMPV <- 100*sum(all.CS_PvCSxP$adj.P.Val<0.05 & abs(all.CS_PvCSxP$logFC)>1)/nrow(all.CS_PvCSxP)
c(`CS vs Paragon (%)`=pct_DE_CSvP,`Hybrid vs MPV (%)`=pct_DE_HybridvMPV)

