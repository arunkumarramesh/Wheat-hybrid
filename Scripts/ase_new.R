setwd("~/Downloads/wheat")
library(edgeR)
library(limma)
library(ComplexHeatmap)
library(ggplot2)
library(scales)

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
hybrid_counts <- hybrid_counts %>%
  filter(
    CSxP1_ref + CSxP1_alt > 0,
    CSxP2_ref + CSxP2_alt > 0,
    CSxP3_ref + CSxP3_alt > 0,
    PxCS1_ref + PxCS1_alt > 0,
    PxCS2_ref + PxCS2_alt > 0
  )

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
keep <- rowSums(cpm(DGEList(total_counts)) >= 2) >= 2
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
keep <- rowSums(cpm(edgeR.DGElist) >= 2) >= 4
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
dim(asetest_pvals[abs(asetest_pvals$logFC) > 0.58 & asetest_pvals$adj.P.Val < 0.05,])
write.csv(asetest_pvals,file="Ref_vs_Alt.csv")

DGEgenes <- rownames(subset(top.table, top.table$adj.P.Val < 0.05))
mat_DGEgenes <- cpm_nolog_relative[DGEgenes, ]
pdf("ase_heatmap.pdf",height=3.5,width=4)
Heatmap(mat_DGEgenes, name = "Scaled CPM", show_row_names = FALSE, use_raster = F)
dev.off()

asetest_pvals_sub <- asetest_pvals[c(1,5)]
asetest_pvals_sub$Sig <- F
asetest_pvals_sub$Sig[abs(asetest_pvals_sub$logFC) > 0.58 & asetest_pvals_sub$adj.P.Val < 0.05] <- TRUE
asetest_pvals_sub <- asetest_pvals_sub[c(-2)]
colnames(asetest_pvals_sub) <- c("H_FC", "H")
asetest_pvals_sub$gene <- rownames(asetest_pvals_sub)
all.CSvP <- read.csv(file = "CSvP all genes.csv",row.names = 1)
all.CSvP <- all.CSvP[c(1,5)]
all.CSvP$Sig <- F
all.CSvP$Sig[abs(all.CSvP$logFC) > 0.58 & all.CSvP$adj.P.Val < 0.05] <- TRUE
all.CSvP <- all.CSvP[c(-2)]
colnames(all.CSvP) <- c("P_FC", "P")
all.CSvP$gene <- rownames(all.CSvP)
all.CS_PvCSxP <- read.csv(file = "CS_PvCSxP all genes.csv",row.names = 1)
all.CS_PvCSxP <- all.CS_PvCSxP[c(1,5)]
all.CS_PvCSxP$Sig <- F
all.CS_PvCSxP$Sig[abs(all.CS_PvCSxP$logFC) > 0.58 & all.CS_PvCSxP$adj.P.Val < 0.05] <- TRUE
all.CS_PvCSxP <- all.CS_PvCSxP[c(-2)]
colnames(all.CS_PvCSxP) <- c("T_FC", "T")
all.CS_PvCSxP$gene <- rownames(all.CS_PvCSxP)
all_genes <- inner_join(asetest_pvals_sub,all.CSvP,by="gene")
all_genes <- inner_join(all_genes,all.CS_PvCSxP,by="gene")

classified <- all_genes %>%
  # convert to comparable signs
  mutate(p_sign = case_when(is.na(P_FC) ~ NA_integer_, P_FC > 0 ~  1L, P_FC < 0  ~ -1L, TRUE   ~  0L), 
         h_sign = case_when(is.na(H_FC) ~ NA_integer_, H_FC > 0 ~ 1L, H_FC < 0 ~ -1L, TRUE ~ 0L) ) %>%
  mutate( category = case_when(
    # Cis only: sig in P and H, NOT sig in T
    P & H & !T ~ "Cis only",
    # Trans only: sig in P, NOT H, but sig in T
    P & !H & T ~ "Trans only",
    # Cis + trans: sig in P, H, T; same sign
    P & H & T & !is.na(p_sign) & !is.na(h_sign) & (p_sign == h_sign) ~ "Cis + trans",
    # Cis × trans: sig in P, H, T; opposite sign
    P & H & T & !is.na(p_sign) & !is.na(h_sign) & (p_sign != h_sign) ~ "Cis × trans",
    # Compensatory: sig in H, NOT P, and sig in T
    !P & H & T ~ "Compensatory",
    # Conserved: none are significant
    !P & !H & !T ~ "Conserved",
    # Everything else
    TRUE ~ "Ambiguous")) %>%
  select(gene, category)

classified_limma <- classified
write.csv(classified_limma,file="classified_limma.csv",row.names = F)

counts_limma <- as.data.frame(table(classified_limma$category))
names(counts_limma) <- c("category", "count")
props_limma <- prop.table(table(classified_limma$category))
df_limma <- counts_limma %>%
  mutate(prop = as.numeric(props_limma[as.character(category)]), category = factor(category, levels = category))

pdf(file="limma_classification_subgenome.pdf",height=3.5,width=4)
ggplot(df_limma, aes(x = category, y = prop, fill = category)) +
  geom_col() +
  geom_text(aes(label = percent(prop, accuracy = 0.1)),vjust = -0.3, size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494"), guide = "none") +
  scale_x_discrete(labels = function(x) {
    x <- gsub("\\b[Cc]is\\b", "<i>cis</i>", x, perl = TRUE)
    x <- gsub("\\b[Tt]rans\\b", "<i>trans</i>", x, perl = TRUE)
    x
  }) +
  labs(x = NULL, y = "Proportion of genes", title = paste("n=",nrow(classified_limma))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +
  coord_cartesian(clip = "off")
dev.off()
