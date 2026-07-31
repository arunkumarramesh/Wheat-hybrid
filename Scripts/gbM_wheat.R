library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(cowplot)

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
    x <- dt[!is.na(get(cov_col)) & get(cov_col) > 2]
    x$meth_count <- round((x[[pct_col]] / 100) * x[[cov_col]])
    x$total_count <- x[[cov_col]]
    x$site_p_value <- pbinom(q = x$meth_count - 1,size = x$total_count,prob = 0.002,lower.tail = FALSE)
    x$site_padj <- p.adjust(x$site_p_value,method="BH")
    x$site_methylated <- x$site_padj < 0.001
    genome_p <- mean(x$site_methylated)
    gene_summary <- x[,list(sites = .N,nC = .N,mC = sum(site_methylated)),by = gene_id]
    gene_summary <- gene_summary[sites >= 20]
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
all_results <- all_results[!is.na(padj_CG) & !is.na(padj_CHG) & !is.na(padj_CHH)]
all_results$gbM_candidate <- all_results$padj_CG <= 0.05 & all_results$padj_CHG > 0.05 & all_results$padj_CHH > 0.05
all_results$gbM_class <- ifelse(all_results$gbM_candidate == TRUE,"gbM",ifelse(all_results$padj_CG > 0.05 & all_results$padj_CHG > 0.05 & all_results$padj_CHH > 0.05,"UM",NA))
all_results <- all_results[!is.na(all_results$gbM_class),]
all_results <- all_results[, -c(18:20, 24:26), with = FALSE]
fwrite(all_results, "gene_body_methylation.tsv", sep = "\t")


## gbM by expression level

library(data.table)
library(ggplot2)

cs_tpm <- read.table(file="cs_tpm.tsv")
cs_tpm <- as.data.table(cs_tpm,keep.rownames="gene_id")
cs_tpm <- cs_tpm[!grepl("LC$",gene_id)]
colnames(cs_tpm) <- c("gene_id",sub("^PXCS","PxCS",sub("_.*","",colnames(cs_tpm)[-1])))
cs_tpm <- cs_tpm[,c("gene_id","CS1","CS2","CS3","CSxP1","CSxP2","CSxP3","P1","P2","P3"),with=FALSE]
sample_cols <- colnames(cs_tpm)[-1]
cs_tpm <- cs_tpm[rowSums(as.data.frame(cs_tpm)[,sample_cols]>=1)>=3]
tpm_gene <- data.table(gene_id=cs_tpm$gene_id,CS=rowMeans(cs_tpm[,c("CS1","CS2","CS3"),with=FALSE],na.rm=TRUE),CSxP=rowMeans(cs_tpm[,c("CSxP1","CSxP2","CSxP3"),with=FALSE],na.rm=TRUE),P=rowMeans(cs_tpm[,c("P1","P2","P3"),with=FALSE],na.rm=TRUE))
tpm_long <- melt(tpm_gene,id.vars="gene_id",variable.name="genotype",value.name="TPM")
gbm_status <- fread("gene_body_methylation.tsv")[,.(gene_id,sample,gbM_class)]
gbm_status <- dcast(gbm_status,gene_id~sample,value.var="gbM_class")
gbm_status <- gbm_status[!is.na(CS) & !is.na(CSxP) & !is.na(P)]
gbm_status <- gbm_status[CS==CSxP & CSxP==P]
gbm_status[,gbM_status:=factor(CS,levels=c("UM","gbM"))]
gbm_status <- gbm_status[,.(gene_id,gbM_status)]
gbm_subgenome <- gbm_status
gbm_subgenome[,subgenome:=sub("^TraesCS[0-9]+([ABD]).*$","\\1",gene_id)]
gbm_subgenome[,subgenome:=factor(subgenome,levels=c("A","B","D"))]

tab_subgenome <- table(gbm_subgenome$subgenome,gbm_subgenome$gbM_status)
print(tab_subgenome)
print(prop.table(tab_subgenome,margin=1))
print(chisq.test(tab_subgenome))

gbm_subgenome_summary <- gbm_subgenome[,.(n_genes=.N,n_gbM=sum(gbM_status=="gbM"),prop_gbM=mean(gbM_status=="gbM")),by=subgenome]
print(gbm_subgenome_summary)

plot_dt <- merge(tpm_long,gbm_status,by="gene_id")
plot_dt[,genotype:=factor(genotype,levels=c("CS","CSxP","P"))]
plot_dt[,log2_TPM:=log2(TPM+1)]
median_labs <- plot_dt[,.(median_log2_TPM=median(log2_TPM,na.rm=TRUE)),by=gbM_status]
median_labs[,label:=round(median_log2_TPM,2)]
wilcox_p <- wilcox.test(log2_TPM~gbM_status,data=plot_dt)$p.value

plot_dt_plot <- ggplot(plot_dt,aes(x=gbM_status,y=log2_TPM,fill=gbM_status)) +
  geom_boxplot(outlier.size=0.2,linewidth=0.3) +
  geom_text(data=median_labs,aes(x=gbM_status,y=median_log2_TPM,label=label),inherit.aes=FALSE,size=3,vjust=-0.4) +
  geom_text(data=data.frame(x=1.5,y=13.5),aes(x=x,y=y,label=paste0("italic(P)~'='~",signif(wilcox_p,3))),inherit.aes=FALSE,size=2.5,parse=TRUE) +
  scale_fill_manual(values=c("UM"="grey70","gbM"="#0072B2")) +
  labs(x=NULL,y=expression(log[2]*"(TPM+1)"),fill=NULL) +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_line(linewidth=0.2,colour="grey90"),strip.background=element_rect(fill="white",colour="black"),legend.position="none",axis.title=element_text(face="bold"))


## gbM and bias
library(data.table)
library(dplyr)
library(ggplot2)
library(agricolae)

bias_categories <- fread("bias_category_all_samples_inc_orig_expr.csv")
bias_categories[,sample:=sub("^PXCS","PxCS",sub("_.*","",sample))]
bias_categories <- bias_categories[!sample %in% c("PxCS1","PxCS2","PxCS3")]
bias_categories[,group_id:=as.numeric(gsub("^X","",group_id))]
bias_categories <- bias_categories %>% mutate(genotype=sub("[0-9]+$","",as.character(sample))) %>% group_by(group_id,genotype) %>% summarise(A_tpm=mean(A_tpm,na.rm=TRUE),B_tpm=mean(B_tpm,na.rm=TRUE),D_tpm=mean(D_tpm,na.rm=TRUE),.groups="drop") %>% mutate(triad_tpm=A_tpm+B_tpm+D_tpm,A=A_tpm/triad_tpm,B=B_tpm/triad_tpm,D=D_tpm/triad_tpm) %>% as.data.table()
bias_categories[,CV:=apply(.SD,1,sd)/apply(.SD,1,mean),.SDcols=c("A","B","D")]
bias_categories <- bias_categories[is.finite(CV)]

homologies <- fread("homologies.csv")
homologies[,group_id:=as.numeric(gsub("^X","",group_id))]
hom_long <- melt(homologies,id.vars=c("group_id","cardinality","synteny","group","chrs","origin"),measure.vars=c("A","B","D"),variable.name="subgenome",value.name="gene_id")

gbm_status <- fread("gene_body_methylation.tsv")[,.(gene_id,sample,gbM_class)]
gbm_status <- dcast(gbm_status,gene_id~sample,value.var="gbM_class")
gbm_status <- gbm_status[!is.na(CS) & !is.na(CSxP) & !is.na(P)]
gbm_status <- gbm_status[CS==CSxP & CSxP==P]
gbm_status <- gbm_status[,.(gene_id,gbM_state=CS)]

triad_gbm_long <- merge(hom_long[,.(group_id,subgenome,gene_id)],gbm_status,by="gene_id")
triad_gbm_wide <- dcast(triad_gbm_long,group_id~subgenome,value.var="gbM_state")
triad_gbm_wide <- triad_gbm_wide[!is.na(A) & !is.na(B) & !is.na(D)]
setnames(triad_gbm_wide,c("A","B","D"),c("A_gbM_state","B_gbM_state","D_gbM_state"))
triad_gbm_wide[,n_gbM:=rowSums(.SD=="gbM"),.SDcols=c("A_gbM_state","B_gbM_state","D_gbM_state")]
triad_gbm_wide[,triad_gbM_state:=fcase(n_gbM==0,"All UM",n_gbM==1,"1 gbM + 2 UM",n_gbM==2,"2 gbM + 1 UM",n_gbM==3,"All gbM")]
triad_gbm_wide[,triad_gbM_state:=factor(triad_gbM_state,levels=c("All UM","1 gbM + 2 UM","2 gbM + 1 UM","All gbM"))]

bias_categories_gbm <- merge(bias_categories,triad_gbm_wide[,.(group_id,triad_gbM_state)],by="group_id")
bias_categories_gbm <- bias_categories_gbm[,.(CV=mean(CV,na.rm=TRUE),n_genotypes=uniqueN(genotype)),by=.(group_id,triad_gbM_state)]
bias_categories_gbm <- bias_categories_gbm[n_genotypes==3]

tuk <- HSD.test(aov(CV~triad_gbM_state,data=bias_categories_gbm),"triad_gbM_state",group=TRUE)
letters_df <- data.frame(triad_gbM_state=rownames(tuk$groups),letters=tuk$groups$groups)
label_df <- bias_categories_gbm %>% group_by(triad_gbM_state) %>% summarise(median_CV=median(CV,na.rm=TRUE),n=n_distinct(group_id),.groups="drop") %>% left_join(letters_df,by="triad_gbM_state")

bias_categories_gbm_plot <- ggplot(bias_categories_gbm,aes(x=triad_gbM_state,y=CV)) +
  geom_boxplot(outlier.shape=1,outlier.size=1) +
  geom_text(data=label_df,aes(x=triad_gbM_state,y=median_CV+0.09,label=round(median_CV,2)),size=2) +
  geom_text(data=label_df,aes(x=triad_gbM_state,y=1.8,label=letters),size=4,fontface="bold") +
  geom_text(data=label_df,aes(x=triad_gbM_state,y=2.03,label=paste0("n=",n)),size=2.5) +
  labs(x="",y="HEB") +
  coord_flip(ylim=c(0,2.1),clip="off") +
  theme_bw(base_size=12) +
  theme(axis.title=element_text(face="bold"),panel.grid.minor=element_blank(),plot.margin=margin(5.5,30,5.5,5.5))


## transgressive methylation

library(data.table)
library(ggplot2)

gbm_status <- fread("gene_body_methylation.tsv")[,.(gene_id,sample,gbM_class)]
stopifnot(all(gbm_status[, .N,by=.(gene_id,sample)]$N==1))
gbm_status <- dcast(gbm_status,gene_id~sample,value.var="gbM_class")
gbm_status <- gbm_status[!is.na(CS) & !is.na(CSxP) & !is.na(P)]
gbm_status <- gbm_status[CS==P & P==CSxP]
gbm_status[,gbM_status:=factor(CS,levels=c("gbM","UM"))]
transgressive_genes <- fread("transgressive_genes.txt",header=FALSE)[,.(gene_id=V1,label="Transgressive")]
dominant_genes <- fread("dominant_genes.txt",header=FALSE)[,.(gene_id=V1,label="Dominant")]
additive_genes <- fread("additive_genes.txt",header=FALSE)[,.(gene_id=V1,label="Additive")]
nonDE_genes <- fread("nonDE_genes.txt",header=FALSE)[,.(gene_id=V1,label="non-DE")]
gene_class <- unique(rbindlist(list(transgressive_genes,dominant_genes,additive_genes,nonDE_genes)))
gbm_status <- merge(gbm_status,gene_class,by="gene_id")
gbm_status[,label:=factor(label,levels=c("Transgressive","Dominant","Additive","non-DE"))]
tab <- table(gbm_status$gbM_status,gbm_status$label)
print(tab)
print(chisq.test(tab))
gbm_status2 <- gbm_status[, .N,by=.(label,gbM_status)]
gbm_status2[,prop:=N/sum(N),by=label]
n_dt <- gbm_status[,.(n=.N),by=label]
n_dt[,n_label:=paste0("n=",n)]
print(gbm_status2)

gbm_status2_plot <- ggplot(gbm_status2,aes(x=label,y=prop,fill=gbM_status)) +
  geom_col(width=0.7) +
  geom_text(data=n_dt,aes(x=label,y=0.9,label=n_label),inherit.aes=FALSE,size=3) +
  scale_fill_manual(values=c("UM"="grey70","gbM"="#0072B2")) +
  labs(x="",y="Proportion of genes",fill="") +
  theme_classic() +
  theme(legend.position = "top")

pdf(file="gbM_DE.pdf",height=3.5,width=5)
plot_grid(plot_grid(plot_dt_plot,gbm_status2_plot,ncol=2,labels=c("A","C"),rel_widths=c(0.65,1.7)),bias_categories_gbm_plot,ncol=1,labels=c("","B"),rel_heights = c(1.7,1))
dev.off()

