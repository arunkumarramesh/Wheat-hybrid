
library(edgeR)
library(limma)
library(ComplexHeatmap)
library(ggplot2)
library(scales)
library(cowplot)

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

## files for 1:1 mapping
SingleCopyOrthologues <- read.table(file="SingleCopyOrthologues_matrix.tsv")
SingleCopyOrthologues$V3 <- gsub("\\..*","",SingleCopyOrthologues$V3 )
SingleCopyOrthologues$V2 <- sub("^([^.]*\\.[^.]*)\\..*$", "\\1", SingleCopyOrthologues$V2)
SingleCopyOrthologues_unique <- SingleCopyOrthologues %>%
  add_count(V2, name = "nV2") %>%
  add_count(V3, name = "nV3") %>%
  filter(nV2 == 1, nV3 == 1) %>%
  select(-nV2, -nV3)
SingleCopyOrthologues_unique <- SingleCopyOrthologues_unique[2:3]
colnames(SingleCopyOrthologues_unique) <- c("Par_gene","CS_gene")

cs_longest <- read.table(file="iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna_longest.fasta.fai")[1:2]
cs_longest <- cs_longest[!grepl("LC", cs_longest$V1), ]
cs_longest$V1 <- sub("\\.[0-9]+$", "", cs_longest$V1)
colnames(cs_longest) <- c("CS_gene","CS_length")
par_longest <- read.table(file="Triticum_aestivum_paragon.GCA949126075v1.cdna.longest.fa.fai")[1:2]
par_longest$V1 <- sub("\\.[0-9]+$", "", par_longest$V1)
colnames(par_longest) <- c("Par_gene","Par_length")

SingleCopyOrthologues_unique <- inner_join(SingleCopyOrthologues_unique,cs_longest,by="CS_gene")
SingleCopyOrthologues_unique <- inner_join(SingleCopyOrthologues_unique,par_longest,by="Par_gene")
hist(SingleCopyOrthologues_unique$CS_length - SingleCopyOrthologues_unique$Par_length)
SingleCopyOrthologues_unique <- SingleCopyOrthologues_unique[abs(SingleCopyOrthologues_unique$CS_length - SingleCopyOrthologues_unique$Par_length) < 100,]
SingleCopyOrthologues_unique <- SingleCopyOrthologues_unique[SingleCopyOrthologues_unique$CS_gene %in% core$v21,]

files <- c("CS1_RNA_MKRN250026357-1A_22VTNMLT4_L3.tsv", "CS2_RNA_MKRN250026358-1A_22VTNMLT4_L4.tsv", "CS3_RNA_MKRN250026359-1A_22VTNMLT4_L3.tsv", "P1_RNA_MKRN250026354-1A_22VTNMLT4_L3.tsv", "P2_RNA_MKRN250026355-1A_22VTNMLT4_L3.tsv", "P3_RNA_MKRN250026356-1A_22VTNMLT4_L4.tsv")

process_file <- function(file) {
  x <- read.table(file, sep = "\t")
  colnames(x) <- c("gene_id", "read_count")
  x <- x %>%
    filter(!grepl("LC", gene_id)) %>%
    mutate(gene_id = sub("\\.[0-9]+$", "", gene_id))
  SingleCopyOrthologues_unique %>%
    left_join(x %>% select(gene_id, CS_count = read_count), by = c("CS_gene" = "gene_id")) %>%
    left_join(x %>% select(gene_id, Par_count = read_count), by = c("Par_gene" = "gene_id")) %>%
    mutate(CS_count = replace_na(CS_count, 0), Par_count = replace_na(Par_count, 0)) %>%
    filter(CS_count > 3, Par_count < 1)
}

CS <- lapply(files[1:3], process_file)
names(CS) <- c("CS1", "CS2", "CS3")
CS_common_genes <- data.frame(CS_gene = Reduce(intersect, lapply(CS, function(x) x$CS_gene)))

process_par_file <- function(file) {
  x <- read.table(file, sep = "\t")
  colnames(x) <- c("gene_id", "read_count")
  x <- x %>%
    filter(!grepl("LC", gene_id)) %>%
    mutate(gene_id = sub("\\.[0-9]+$", "", gene_id))
  SingleCopyOrthologues_unique %>%
    left_join(x %>% select(gene_id, CS_count = read_count), by = c("CS_gene" = "gene_id")) %>%
    left_join(x %>% select(gene_id, Par_count = read_count), by = c("Par_gene" = "gene_id")) %>%
    mutate(CS_count = replace_na(CS_count, 0), Par_count = replace_na(Par_count, 0)) %>%
    filter(Par_count > 3, CS_count < 1)
}

Par <- lapply(files[4:6], process_par_file)
names(Par) <- c("P1", "P2", "P3")
Par_common_genes <- data.frame(CS_gene = Reduce(intersect, lapply(Par, function(x) x$CS_gene)))

CS_common_genes <- data.frame(CS_gene = intersect(CS_common_genes$CS_gene, Par_common_genes$CS_gene))

files <- c("CSxP1_MKRN250026363-1A_22VTNMLT4_L3.tsv", "CSxP2_MKRN250026364-1A_22VTNMLT4_L3.tsv", "CSxP3_RNA_MKRN250033262-1A_22VTNMLT4_L3.tsv", "PxCS1_MKRN250026360-1A_22VTNMLT4_L3.tsv", "PXCS2_RNA_MKRN250033261-1A_22VTNMLT4_L4.tsv")

process_file <- function(file, sample) {
  x <- read.table(file, sep = "\t")
  colnames(x) <- c("gene_id", "read_count")
  x <- x %>%
    filter(!grepl("LC", gene_id)) %>%
    mutate(gene_id = sub("\\.[0-9]+$", "", gene_id))
  SingleCopyOrthologues_unique %>%
    left_join(x %>% select(gene_id, ref = read_count), by = c("CS_gene" = "gene_id")) %>%
    left_join(x %>% select(gene_id, alt = read_count), by = c("Par_gene" = "gene_id")) %>%
    mutate(ref = replace_na(ref, 0), alt = replace_na(alt, 0)) %>%
    filter(CS_gene %in% CS_common_genes$CS_gene) %>%
    select(gene = CS_gene, ref, alt) %>%
    dplyr::rename(!!paste0(sample, "_ref") := ref, !!paste0(sample, "_alt") := alt)
}

sample_names <- sub("(_RNA)?_MKRN.*", "", files)
x <- mapply(process_file, files, sample_names, SIMPLIFY = FALSE)
hybrid_counts <- Reduce(function(a, b) inner_join(a, b, by = "gene"), x)
colnames(hybrid_counts)[10:11] <- c("PxCS2_ref","PxCS2_alt")
hybrid_counts <- hybrid_counts %>%
  select(gene, ends_with("_ref"), ends_with("_alt"))
rownames(hybrid_counts) <- hybrid_counts$gene
hybrid_counts <- hybrid_counts[-c(1)]

ratio_df <- hybrid_counts %>%
  mutate(gene = rownames(.)) %>%
  pivot_longer(-gene, names_to = c("sample", "allele"), names_pattern = "(.*)_(ref|alt)") %>%
  pivot_wider(names_from = allele, values_from = value) %>%
  mutate(ratio = ref / (ref + alt)) %>%
  filter(ref + alt > 0)

med_df <- ratio_df %>%
  group_by(sample) %>%
  summarise(med = median(ratio, na.rm = TRUE))

pdf(file="ref_prop.pdf",height=2.5,width=9)
ggplot(ratio_df, aes(ratio)) +
  geom_histogram(bins = 50,fill="#E69F00") +
  geom_vline(data = med_df, aes(xintercept = med), linetype = "dashed") +
  geom_text(data = med_df, aes(x = med, y = Inf, label = round(med, 3)), vjust = 1.5, hjust = 1.5) +
  facet_wrap(~sample,ncol=5) +
  xlab("Proportion of Chinese Spring reads") +
  ylab("Number of genes") +
  theme_minimal()
dev.off()

## first testing if ref counts for reciprocal hybrids differ

total_counts <- hybrid_counts[,1:5] + hybrid_counts[,6:10]
keep <- rowSums(cpm(DGEList(total_counts)) >= 2) >= 4
dge <- DGEList(counts = hybrid_counts[keep,])
nf <- calcNormFactors(DGEList(counts = total_counts[keep,]), method = "TMM")
dge$samples$lib.size <- rep(nf$samples$lib.size, 2)
dge$samples$norm.factors <- rep(nf$samples$norm.factors, 2)
cpm_log <- cpm(dge, log = TRUE)
cpm_nolog <- cpm(dge, log = FALSE)
cpm_nolog_relative <- cpm_nolog/rowMeans(cpm_nolog)
allele_alt <- c(rep(0, 5), rep(1, 5))
cross_pxcs <- rep(c(0, 0, 0, 1, 1), 2)
sample <- factor(rep(1:5, 2))
design <- model.matrix(~ sample + allele_alt + I(allele_alt * cross_pxcs))
colnames(design)[ncol(design)] <- "allele_by_cross"
v <- voom(dge, design, plot = TRUE)
fit <- eBayes(lmFit(v, design))
reciprocal_ase <- topTable(fit, coef = "allele_by_cross", n = Inf, sort.by = "P")
length(which(reciprocal_ase$adj.P.Val < 0.05)) ## how many significantly DE genes
ref_reciprocal_sig <- reciprocal_ase %>%
  filter(adj.P.Val < 0.05)
mat_DGEgenes <- cpm_nolog_relative[rownames(ref_reciprocal_sig), ]
pdf("poi_genes_ase.pdf",height=4,width=6)
Heatmap(mat_DGEgenes, name = "Scaled CPM", show_row_names = T, use_raster = F)
dev.off()

# now testing of CS and Paragon alleles differ after removing genes that show allele differences in reciprocal hybrids
sample_info.edger <- factor(c(rep("ref", 5), rep("alt", 5)))
hybrid_counts <- hybrid_counts[!rownames(hybrid_counts) %in% rownames(ref_reciprocal_sig),]
edgeR.DGElist <- DGEList(counts = hybrid_counts, group = sample_info.edger)
total_counts <- hybrid_counts[, 1:5] + hybrid_counts[, 6:10]
keep <- rowSums( cpm(DGEList(total_counts)) >= 2 ) >= 4
edgeR.DGElist <- edgeR.DGElist[keep,]

total_counts <- hybrid_counts[1:5] + hybrid_counts[6:10]
total_counts <- total_counts[keep,]
edgeR.DGElist$samples$lib.size <- rep(calcNormFactors(DGEList(counts = total_counts), method = "TMM")$samples[,2], 2)
edgeR.DGElist$samples$norm.factors <- rep(calcNormFactors(DGEList(counts = total_counts), method = "TMM")$samples[,3], 2)
allele <- factor(c(rep("ref", 5), rep("alt", 5)), levels = c("ref", "alt"))
sample <- factor(rep(1:5, 2))
mm <- model.matrix(~ sample + allele)
y <- voom(edgeR.DGElist, mm, plot = TRUE)
fit <- lmFit(y, mm)
asetest <- eBayes(fit)
top.table <- topTable(asetest, coef = "allelealt", sort.by = "P", n = Inf)
asetest_pvals <- topTable(asetest, coef = "allelealt", sort.by = "none", n = Inf, p.value = 1, lfc = 0)
dim(asetest_pvals[abs(asetest_pvals$logFC) > 1 & asetest_pvals$adj.P.Val < 0.05,])
write.csv(asetest_pvals,file="Ref_vs_Alt.csv")

DGEgenes <- rownames(subset(top.table, top.table$adj.P.Val < 0.05))
mat_DGEgenes <- cpm_nolog_relative[DGEgenes, ]
pdf("ase_heatmap.pdf",height=3.5,width=4)
Heatmap(mat_DGEgenes, name = "Scaled CPM", show_row_names = FALSE, use_raster = F)
dev.off()


category_levels <- c("Conserved","Cis only","Trans only","Cis + trans","Cis × trans","Compensatory","Ambiguous")
category_colors <- c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494")

asetest_pvals_sub <- asetest_pvals[c(1,5)]
asetest_pvals_sub$H <- abs(asetest_pvals_sub$logFC) > 1 & asetest_pvals_sub$adj.P.Val < 0.05
asetest_pvals_sub <- asetest_pvals_sub[-2]
colnames(asetest_pvals_sub) <- c("H_FC","H")
asetest_pvals_sub$gene <- rownames(asetest_pvals_sub)

all.CSvP <- read.csv("CSvP all genes.csv",row.names=1)
all.CSvP$P <- abs(all.CSvP$logFC) > 1 & all.CSvP$adj.P.Val < 0.05
all.CSvP <- all.CSvP[c("logFC","P")]
colnames(all.CSvP) <- c("P_FC","P")
all.CSvP$gene <- rownames(all.CSvP)

## Mcmanus approach

parent_counts_base <- data.frame(gene=CS_common_genes$CS_gene,CS1=CS[[1]]$CS_count[match(CS_common_genes$CS_gene,CS[[1]]$CS_gene)],CS2=CS[[2]]$CS_count[match(CS_common_genes$CS_gene,CS[[2]]$CS_gene)],CS3=CS[[3]]$CS_count[match(CS_common_genes$CS_gene,CS[[3]]$CS_gene)],P1=Par[[1]]$Par_count[match(CS_common_genes$CS_gene,Par[[1]]$CS_gene)],P2=Par[[2]]$Par_count[match(CS_common_genes$CS_gene,Par[[2]]$CS_gene)],P3=Par[[3]]$Par_count[match(CS_common_genes$CS_gene,Par[[3]]$CS_gene)])

fisher_parent_counts <- parent_counts_base %>%
  mutate(parent_CS=CS1+CS2+CS3,parent_Par=P1+P2+P3) %>%
  select(gene,parent_CS,parent_Par)

fisher_ref_cols <- grep("_ref$",colnames(hybrid_counts),value=TRUE)
fisher_alt_cols <- grep("_alt$",colnames(hybrid_counts),value=TRUE)

fisher_hybrid_counts <- data.frame(gene=rownames(hybrid_counts),hybrid_CS=rowSums(hybrid_counts[,fisher_ref_cols,drop=FALSE]),hybrid_Par=rowSums(hybrid_counts[,fisher_alt_cols,drop=FALSE]))

fisher_T_results <- inner_join(fisher_parent_counts,fisher_hybrid_counts,by="gene")

fisher_T_results$T_p <- mapply(function(parent_CS,parent_Par,hybrid_CS,hybrid_Par) fisher.test(matrix(c(parent_CS,parent_Par,hybrid_CS,hybrid_Par),nrow=2,byrow=TRUE))$p.value,fisher_T_results$parent_CS,fisher_T_results$parent_Par,fisher_T_results$hybrid_CS,fisher_T_results$hybrid_Par)

fisher_T_results$T <- p.adjust(fisher_T_results$T_p,method="BH") < 0.05

fisher_all_genes <- asetest_pvals_sub %>%
  inner_join(all.CSvP,by="gene") %>%
  inner_join(fisher_T_results[c("gene","T")],by="gene")

fisher_classified <- fisher_all_genes %>%
  mutate(T_FC=P_FC-H_FC,h_sign=sign(H_FC),t_sign=sign(T_FC),category=case_when(P & H & !T ~ "Cis only",P & !H & T ~ "Trans only",P & H & T & h_sign == t_sign ~ "Cis + trans",P & H & T & h_sign != t_sign ~ "Cis × trans",!P & H & T ~ "Compensatory",!P & !H & !T ~ "Conserved",TRUE ~ "Ambiguous"))

fisher_classified$category <- factor(fisher_classified$category,levels=category_levels)

fisher_df <- fisher_classified %>%
  count(category,.drop=FALSE,name="count") %>%
  mutate(prop=count/sum(count))

fisher_plot <- ggplot(fisher_df,aes(x=category,y=prop,fill=category)) +
  geom_col() +
  geom_text(aes(label=percent(prop,accuracy=0.1)),vjust=-0.3,size=3) +
  scale_y_continuous(labels=percent_format(accuracy=0.1),expand=expansion(mult=c(0,0.08))) +
  scale_fill_manual(values=category_colors,guide="none") +
  scale_x_discrete(labels=function(x) {
    x <- gsub("\\b[Cc]is\\b","<i>cis</i>",x,perl=TRUE)
    x <- gsub("\\b[Tt]rans\\b","<i>trans</i>",x,perl=TRUE)
    x
  }) +
  labs(x=NULL,y="Proportion of genes",title=paste0("Fisher test (n=",nrow(fisher_classified),")")) +
  theme_minimal(base_size=12) +
  theme(axis.text.x=ggtext::element_markdown(angle=90,hjust=1)) +
  coord_cartesian(clip="off")

## Haas approach

barley_ref_cols <- grep("_ref$",colnames(hybrid_counts),value=TRUE)
barley_hybrid_samples <- sub("_ref$","",barley_ref_cols)
barley_alt_cols <- paste0(barley_hybrid_samples,"_alt")

barley_genes <- Reduce(intersect,list(asetest_pvals_sub$gene,all.CSvP$gene,parent_counts_base$gene,rownames(hybrid_counts)))

barley_parent_counts <- parent_counts_base[match(barley_genes,parent_counts_base$gene),]
barley_hybrid_counts <- hybrid_counts[barley_genes,,drop=FALSE]

barley_parent_matrix <- as.matrix(barley_parent_counts[,c("CS1","CS2","CS3","P1","P2","P3"),drop=FALSE])
rownames(barley_parent_matrix) <- barley_parent_counts$gene
storage.mode(barley_parent_matrix) <- "numeric"

barley_hybrid_CS <- as.matrix(barley_hybrid_counts[,barley_ref_cols,drop=FALSE])
barley_hybrid_Par <- as.matrix(barley_hybrid_counts[,barley_alt_cols,drop=FALSE])
colnames(barley_hybrid_CS) <- barley_hybrid_samples
colnames(barley_hybrid_Par) <- barley_hybrid_samples
storage.mode(barley_hybrid_CS) <- "numeric"
storage.mode(barley_hybrid_Par) <- "numeric"

barley_hybrid_total <- barley_hybrid_CS + barley_hybrid_Par
barley_biological_counts <- cbind(barley_parent_matrix,barley_hybrid_total)
colnames(barley_biological_counts) <- c("CS1","CS2","CS3","P1","P2","P3",barley_hybrid_samples)

barley_joint_counts <- cbind(barley_parent_matrix,barley_hybrid_CS,barley_hybrid_Par)
colnames(barley_joint_counts) <- c("CS1_parent","CS2_parent","CS3_parent","P1_parent","P2_parent","P3_parent",paste0(barley_hybrid_samples,"_CS"),paste0(barley_hybrid_samples,"_Par"))

barley_nf <- calcNormFactors(DGEList(counts=barley_biological_counts),method="TMM")
barley_dge <- DGEList(counts=barley_joint_counts)
barley_dge$samples$lib.size <- c(barley_nf$samples$lib.size[1:6],barley_nf$samples$lib.size[7:11],barley_nf$samples$lib.size[7:11])
barley_dge$samples$norm.factors <- c(barley_nf$samples$norm.factors[1:6],barley_nf$samples$norm.factors[7:11],barley_nf$samples$norm.factors[7:11])

barley_group <- factor(c(rep("Parent_CS",3),rep("Parent_Par",3),rep("Hybrid_CS",5),rep("Hybrid_Par",5)),levels=c("Parent_CS","Parent_Par","Hybrid_CS","Hybrid_Par"))
barley_design <- model.matrix(~0+barley_group)
colnames(barley_design) <- levels(barley_group)

barley_pair <- c(rep("parent",6),barley_hybrid_samples,barley_hybrid_samples)

for (sample_name in barley_hybrid_samples[-1]) {
  barley_design <- cbind(barley_design,as.integer(barley_pair == sample_name))
  colnames(barley_design)[ncol(barley_design)] <- paste0("pair_",sample_name)
}

barley_voom <- voom(barley_dge,barley_design,plot=FALSE)
barley_fit <- lmFit(barley_voom,barley_design)
barley_fit <- contrasts.fit(barley_fit,makeContrasts(T=(Parent_Par-Parent_CS)-(Hybrid_Par-Hybrid_CS),levels=barley_design))
barley_fit <- eBayes(barley_fit)

barley_T_results <- topTable(barley_fit,coef="T",sort.by="none",n=Inf,p.value=1,lfc=0)
barley_T_results$gene <- rownames(barley_T_results)
barley_T_results <- barley_T_results %>% transmute(gene,T_FC=logFC,T=adj.P.Val < 0.05)

barley_all_genes <- asetest_pvals_sub %>%
  inner_join(all.CSvP,by="gene") %>%
  inner_join(barley_T_results,by="gene")

barley_classified <- barley_all_genes %>%
  mutate(h_sign=sign(H_FC),t_sign=sign(T_FC),category=case_when(P & H & !T ~ "Cis only",P & !H & T ~ "Trans only",P & H & T & h_sign == t_sign ~ "Cis + trans",P & H & T & h_sign != t_sign ~ "Cis × trans",!P & H & T ~ "Compensatory",!P & !H & !T ~ "Conserved",TRUE ~ "Ambiguous"))

barley_classified$category <- factor(barley_classified$category,levels=category_levels)

barley_df <- barley_classified %>%
  count(category,.drop=FALSE,name="count") %>%
  mutate(prop=count/sum(count))

barley_plot <- ggplot(barley_df,aes(x=category,y=prop,fill=category)) +
  geom_col() +
  geom_text(aes(label=percent(prop,accuracy=0.1)),vjust=-0.3,size=3) +
  scale_y_continuous(labels=percent_format(accuracy=0.1),expand=expansion(mult=c(0,0.08))) +
  scale_fill_manual(values=category_colors,guide="none") +
  scale_x_discrete(labels=function(x) {
    x <- gsub("\\b[Cc]is\\b","<i>cis</i>",x,perl=TRUE)
    x <- gsub("\\b[Tt]rans\\b","<i>trans</i>",x,perl=TRUE)
    x
  }) +
  labs(x=NULL,y="Proportion of genes",title=paste0("Replicate-aware test (n=",nrow(barley_classified),")")) +
  theme_minimal(base_size=12) +
  theme(axis.text.x=ggtext::element_markdown(angle=90,hjust=1)) +
  coord_cartesian(clip="off")

combined_plot <- plot_grid(fisher_plot,barley_plot,nrow=1,labels=c("A","B"),align="hv")

pdf("limma_classification_combined.pdf",height=4,width=8)
print(combined_plot)
dev.off()