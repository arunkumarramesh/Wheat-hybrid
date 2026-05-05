
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(GenomicRanges)
library(IRanges)
library(rtracklayer)
library(SummarizedExperiment)
library(MBASED)
library(metap)
library(multtest)
library(purrr)
library(tibble)
library(ComplexHeatmap)
library(purrr)
library(cowplot)
library(ggtext)
library(edgeR)
library(factoextra)
library(goseq)
library(ggpubr)
library(scales)

## file for GO analysis
gene_length <- read.table(file="gene.gff3")
gene_length$V6 <- abs(gene_length$V4 - gene_length$V3)
gene_length <- gene_length[5:6]
bp_go <- read.csv('BP.csv',header = F)
mf_go <- read.csv('MF.csv',header = F)
cc_go <- read.csv('CC.csv',header = F)
all_go <- rbind(bp_go,mf_go,cc_go)
colnames(all_go) <- c('Gene','GO_term','Term')

## get reference allele proportions  with CS reference
files <- c("P1.wasp.ase.tsv", "P2.wasp.ase.tsv", "P3.wasp.ase.tsv", "CS1.wasp.ase.tsv", "CS2.wasp.ase.tsv", "CS3.wasp.ase.tsv")

all_ase_parents <- map_dfr(files, function(f) {
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(sample = sub(".wasp.ase.tsv", "", basename(f), perl = TRUE))
})

toremove <- all_ase_parents[(all_ase_parents$sample %in% c("P1","P2","P3") ) & (all_ase_parents$refCount/all_ase_parents$totalCount > 0.6),]

all_ase_parents_filtered <- all_ase_parents %>%
  anti_join(toremove %>% distinct(contig, position), by = c("contig", "position"))

ref_prop_parents_plot <- ggplot(data=all_ase_parents_filtered,aes(x=refCount/totalCount)) +
  geom_histogram(fill="#0072B2") +
  facet_grid(~sample) +
  xlab("Proportion of Chinese Spring allele") +
  ylab("Number of sites")+
  theme_bw()

files <- c("CSxP1.wasp.ase.tsv", "CSxP2.wasp.ase.tsv", "CSxP3.wasp.ase.tsv", "PxCS1.wasp.ase.tsv", "PxCS2.wasp.ase.tsv")

all_ase <- map_dfr(files, function(f) {
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(sample = sub(".wasp.ase.tsv", "", basename(f), perl = TRUE))
})

all_ase_filtered <- all_ase %>%
  anti_join(toremove %>% distinct(contig, position), by = c("contig", "position"))

all_ase_filtered %>% group_by(sample) %>% summarise(median_prop = median(refCount/totalCount))
meds <- all_ase_filtered %>% group_by(sample) %>% summarise(median_prop = median(refCount/totalCount))

ref_prop_plot <- ggplot(data=all_ase_filtered,aes(x=refCount/totalCount)) +
  geom_histogram(fill="#E69F00") +
  facet_grid(~sample) +
  geom_text(data = meds,aes(x = 0.25, y = 2000,label = percent(median_prop, accuracy = 0.01)), vjust = 1.1, position = position_nudge(x = 0.01), size = 3, inherit.aes = FALSE) +
  geom_vline(data = all_ase_filtered %>% group_by(sample) %>% summarise(median_prop = median(refCount/totalCount)), aes(xintercept = median_prop), linetype = "dashed") +
  xlab("Proportion of Chinese Spring allele") +
  ylab("Number of sites") +
  theme_bw()

pdf(file="ref_prop_parents.pdf",height=2.5,width=11)
ref_prop_parents_plot
dev.off()

pdf(file="ref_prop.pdf",height=2.5,width=9)
ref_prop_plot
dev.off()

## now get reference allele proportions using Paragon reference

files <- c("P1.par.ase.tsv", "P2.par.ase.tsv", "P3.par.ase.tsv", "CS1.par.ase.tsv", "CS2.par.ase.tsv", "CS3.par.ase.tsv")

all_ase_parents <- map_dfr(files, function(f) {
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(sample = sub(".par.ase.tsv", "", basename(f), perl = TRUE))
})

toremove_par <- all_ase_parents[(all_ase_parents$sample %in% c("CS1","CS2","CS3") ) & (all_ase_parents$refCount/all_ase_parents$totalCount > 0.6),]

all_ase_parents_filtered <- all_ase_parents %>%
  anti_join(toremove_par %>% distinct(contig, position), by = c("contig", "position"))

ref_prop_parents_plot <- ggplot(data=all_ase_parents_filtered,aes(x=refCount/totalCount)) +
  geom_histogram(fill="#CC79A7") +
  facet_grid(~sample) +
  xlab("Proportion of Paragon allele") +
  ylab("Number of sites")

files <- c("CSxP1.par.ase.tsv", "CSxP2.par.ase.tsv", "CSxP3.par.ase.tsv", "PxCS1.par.ase.tsv", "PxCS2.par.ase.tsv")

all_ase <- map_dfr(files, function(f) {
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(sample = sub(".par.ase.tsv", "", basename(f), perl = TRUE))
})

all_ase_filtered <- all_ase %>%
  anti_join(toremove_par %>% distinct(contig, position), by = c("contig", "position"))

all_ase_filtered %>% group_by(sample) %>% summarise(median_prop = median(refCount/totalCount))
meds_par <- all_ase_filtered %>% group_by(sample) %>% summarise(median_prop = median(refCount/totalCount))

ref_prop_plot <- ggplot(data=all_ase_filtered,aes(x=refCount/totalCount)) +
  geom_histogram(fill="#E69F00") +
  facet_grid(~sample) +
  geom_text(data = meds_par,aes(x = 0.25, y = 2000,label = percent(median_prop, accuracy = 0.01)), vjust = 1.1, position = position_nudge(x = 0.01), size = 3, inherit.aes = FALSE) +
  geom_vline(data = all_ase_filtered %>% group_by(sample) %>% summarise(median_prop = median(refCount/totalCount)), aes(xintercept = median_prop), linetype = "dashed") +
  xlab("Proportion of Paragon allele") +
  ylab("Number of sites")

pdf(file="ref_prop_parents_PAR.pdf",height=2.5,width=11)
ref_prop_parents_plot
dev.off()

pdf(file="ref_prop.pdf_PAR.pdf",height=2.5,width=9)
ref_prop_plot
dev.off()

## Now doing ASE analyses with CS reference

ase_files <- c("CSxP1.wasp.ase.tsv", "CSxP2.wasp.ase.tsv", "CSxP3.wasp.ase.tsv", "PxCS1.wasp.ase.tsv", "PxCS2.wasp.ase.tsv")

read_one_ase <- function(fp) {
  smp <- sub(".wasp.ase.tsv", "", basename(fp), ignore.case = TRUE)
  tb <- readr::read_tsv(fp, col_types = cols(contig = col_character(),position = col_integer(),refCount = col_integer(),altCount = col_integer()),progress = FALSE) %>%
    mutate( locus = paste0(contig, "_", position),!!smp := paste0(pmax(refCount, 0L), ",",pmax(altCount, 0L))) %>%
    select(locus, !!smp)
  tb
}
ase_list <- lapply(ase_files, read_one_ase)

datum2c <- purrr::reduce(ase_list, full_join, by = "locus") %>%
  arrange(locus)

# Replace NAs with "0,0"
for (j in setdiff(colnames(datum2c), "locus")) {
  datum2c[[j]][is.na(datum2c[[j]])] <- "0,0"
}

datum2c <- as.data.frame(datum2c,check.names = FALSE)
rownames(datum2c) <- datum2c$locus
datum2c$locus <- NULL

datum2c <- datum2c %>%
  rownames_to_column("site") %>%
  anti_join(toremove %>% transmute(site = paste0(contig, "_", as.integer(position))) %>% distinct(), by = "site") %>%
  column_to_rownames("site")

## convert part names in variant calls back to whole for gtf compatability
part_sizes <- read.table("part_sizes.txt", header = FALSE, col.names = c("part", "size"))
part_sizes <- part_sizes[-c(nrow(part_sizes)),]
ps2 <- part_sizes %>%
  filter(grepl("_part[12]$", part)) %>%
  mutate(base = sub("_part[12]$", "", part), part_num = ifelse(grepl("_part2$", part), 2L, 1L))
p1 <- ps2 %>% filter(part_num == 1L) %>% select(base, size1 = size) %>% distinct()
out <- ps2 %>% left_join(p1, by = "base") %>%
  mutate(offset = ifelse(part_num == 1L, 0L, size1))
offset_lookup <- setNames(out$offset, out$part)
df <- datum2c %>% tibble::rownames_to_column("locus")
m <- stringr::str_match(df$locus, "^(Chr[^_]+)_part([12])_(\\d+)$")
hit <- !is.na(m[,1])
seqname <- sub("_(\\d+)$", "", df$locus)
pos <- as.integer(sub("^.*_(\\d+)$", "\\1", df$locus))
seqname[hit] <- m[hit, 2]
part_tag <- paste0(m[hit, 2], "_part", m[hit, 3])
off <- unname(offset_lookup[part_tag]); off[is.na(off)] <- 0L
pos[hit] <- pos[hit] + off
coord <- tibble(seqname, pos)
rowRanges <- GRanges(seqnames = coord$seqname, ranges = IRanges(coord$pos, width = 1), aseID = df$locus, allele1  = rep("REF", nrow(df)), allele2 = rep("ALT", nrow(df)))
names(rowRanges) <- paste0(coord$seqname, ":", coord$pos)

#  Map SNPs to genes
genes_gr <- rtracklayer::import("genes_refseqv2_HC.gff3")
genes_gr <- genes_gr[genes_gr$type == "gene"]
id_field <- intersect(c("ID","gene_id","Name"), colnames(mcols(genes_gr)))[1]
mcols(genes_gr)$gid <- as.character(mcols(genes_gr)[[id_field]])
hits <- findOverlaps(rowRanges, genes_gr, select = "first")
mcols(rowRanges)$aseID <- mcols(genes_gr)$gid[hits]
keep_overlap <- !is.na(mcols(rowRanges)$aseID)
rowRanges <- rowRanges[keep_overlap]
df <- df[keep_overlap, , drop = FALSE]
gene_ids  <- mcols(rowRanges)$aseID
feature_ids <- names(rowRanges)

## MBASED
run_one_sample <- function(sample_col, serial = TRUE) {
  x <- df[[sample_col]]
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "0,0"
  x <- gsub("\\s+", "", x)
  parts <- strsplit(x, ",", fixed = TRUE)
  ref <- as.integer(vapply(parts, function(v) if (length(v) >= 1) v[1] else "0", FUN.VALUE = character(1)))
  alt <- as.integer(vapply(parts, function(v) if (length(v) >= 2) v[2] else "0", FUN.VALUE = character(1)))
  ref[is.na(ref)] <- 0L; alt[is.na(alt)] <- 0L
  rna <- data.frame(ref = ref, alt = alt, stringsAsFactors = FALSE)
  tot <- rna$ref + rna$alt
  ok <- is.finite(tot) & tot > 0
  rr   <- rowRanges[ok]
  rref <- rna$ref[ok]
  ralt <- rna$alt[ok]
  feats <- names(rr)
  gids  <- mcols(rr)$aseID
  # comparing against median values
  pr <- rep(meds[meds$sample %in% sample_col,]$median_prop, length(rr))
  #pr <- rep(0.5, length(rr))
  SE <- SummarizedExperiment(
    assays = list(
      lociAllele1Counts = matrix(rref, ncol = 1, dimnames = list(feats, sample_col)),
      lociAllele2Counts = matrix(ralt, ncol = 1, dimnames = list(feats, sample_col)),
      lociAllele1CountsNoASEProbs = matrix(pr,   ncol = 1, dimnames = list(feats, sample_col))
    ),
    rowRanges = rr
  )
  
  set.seed(1)
  bp <- if (serial) BiocParallel::SerialParam() else BiocParallel::bpparam()
  # numSim=0 for speed; increase for full inference
  mb <- runMBASED(SE, isPhased = FALSE, numSim = 10^6, BPPARAM = bp)
  res <- tibble::tibble(gene = rownames(assays(mb)$pValueASE), sample = sample_col, pval = assays(mb)$pValueASE[, 1], adjpval = p.adjust(assays(mb)$pValueASE[, 1], method = "BH"), MAF = assays(mb)$majorAlleleFrequency[, 1])
  dir_df <- tibble::tibble( snv = feats, gene = gids, rref = rref, ralt = ralt) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarize(delta = log2(((sum(rref)+0.0000001) / (sum(ralt)+0.0000001) )), rref = sum(rref) , ralt = sum(ralt), .groups = "drop") %>%
    dplyr::mutate(sample = sample_col)
  out <- dplyr::left_join(res, dplyr::select(dir_df, gene, sample, delta, rref, ralt),by = c("gene", "sample"))
  return(out)
}

hyb_cols <-  c("CSxP1","CSxP2","CSxP3","PxCS1","PxCS2")
#all_res <- lapply(hyb_cols, run_one_sample) %>%
#  dplyr::bind_rows()

all_res <- as_tibble(read.csv(file="all_res_CS.csv"))

length(unique(all_res$gene)) ## number of testable genes

pval_mat <- all_res %>%
  select(gene, sample, pval) %>%
  group_by(gene, sample) %>%
  summarise(pval = dplyr::first(pval), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = pval) %>% 
  column_to_rownames("gene") %>% 
  as.matrix()

counts_mat <- all_res %>% 
  select(gene, sample, c(rref, ralt)) %>% 
  pivot_longer(cols = c(rref, ralt), names_to = "allele", values_to = "count") %>% 
  pivot_wider(id_cols = gene, names_from = c(sample, allele), values_from = count, names_glue = "{sample}_{allele}")%>% 
  column_to_rownames("gene") %>% 
  as.matrix()

testdat <- pval_mat[complete.cases(pval_mat),]
testdat <- testdat[apply(testdat,1,sd)>0,]

## Some genes vary depending on direction of the cross

##  delta 
delta_mat <- all_res %>%
  select(gene, sample, delta) %>%
  dplyr::group_by(gene, sample) %>%
  dplyr::summarise(delta = dplyr::first(delta), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = delta) %>% 
  column_to_rownames("gene") %>% 
  as.matrix()
testdat2 <- delta_mat[complete.cases(delta_mat),]
testdat2 <- testdat2[(apply(testdat2,1,sd)>0),]

ase_poi_genes2 <- as.data.frame(testdat2)
ase_poi_genes2 <- rbind(ase_poi_genes2[ase_poi_genes2$PxCS1 < -20 & ase_poi_genes2$PxCS2 < -20 & ase_poi_genes2$CSxP1 > -3 & ase_poi_genes2$CSxP2 > -3 & ase_poi_genes2$CSxP3 > -3,],
                        ase_poi_genes2[ase_poi_genes2$PxCS1 > 20 & ase_poi_genes2$PxCS2 > 20 & ase_poi_genes2$CSxP1 < 3 & ase_poi_genes2$CSxP2 < 3 & ase_poi_genes2$CSxP3 < 3,],
                        ase_poi_genes2[ase_poi_genes2$PxCS1 > -3 & ase_poi_genes2$PxCS2 > -3 & ase_poi_genes2$CSxP1 < -20 & ase_poi_genes2$CSxP2 < -20 & ase_poi_genes2$CSxP3 < -20,],
                        ase_poi_genes2[ase_poi_genes2$PxCS1 < 3 & ase_poi_genes2$PxCS2 < 3 & ase_poi_genes2$CSxP1 > 20 & ase_poi_genes2$CSxP2 > 20 & ase_poi_genes2$CSxP3 > 20,])
all_res <- all_res[!(all_res$gene %in% rownames(ase_poi_genes2)),] ## remove those genes from further analyses
csxp <- c("CSxP1","CSxP2","CSxP3")
pxcs <- c("PxCS1","PxCS2")
poi_ase_genes <- all_res %>%
  arrange(gene, sample) %>%
  distinct(gene, sample, .keep_all = TRUE) %>%
  filter(sample %in% c(csxp, pxcs)) %>%
  group_by(gene) %>%
  summarise(
    csxp_n    = sum(sample %in% csxp & !is.na(delta) & !is.na(adjpval)),
    pxcs_n    = sum(sample %in% pxcs & !is.na(delta) & !is.na(adjpval)),
    csxp_sig  = all(adjpval[sample %in% csxp] < 0.05, na.rm = TRUE),
    pxcs_sig  = all(adjpval[sample %in% pxcs] < 0.05, na.rm = TRUE),
    csxp_pos  = all(delta[sample %in% csxp] >  0, na.rm = TRUE),
    csxp_neg  = all(delta[sample %in% csxp] < 0, na.rm = TRUE),
    pxcs_pos  = all(delta[sample %in% pxcs] >  0, na.rm = TRUE),
    pxcs_neg  = all(delta[sample %in% pxcs] < 0, na.rm = TRUE),
    .groups = "drop") %>%
  filter(pxcs_n == 2,csxp_n >= 2,csxp_sig, pxcs_sig,
         (csxp_pos & pxcs_neg) | (csxp_neg & pxcs_pos)) %>%
  mutate(pattern = dplyr::case_when(csxp_pos & pxcs_neg ~ "CSxP +  /  PxCS -", csxp_neg & pxcs_pos ~ "CSxP -  /  PxCS +")) %>%
  select(gene)
ase_poi_genes2 <- rbind(ase_poi_genes2,delta_mat[rownames(delta_mat) %in% poi_ase_genes$gene,])
all_res <- all_res[!(all_res$gene %in% rownames(ase_poi_genes2)),] ## remove those genes from further analyses

write.csv(ase_poi_genes2,file="ase_poi_genes2.csv")
pdf("ase_cross_direction.pdf",height=8,width=5.5)
Heatmap(log2(pmax(counts_mat[rownames(counts_mat) %in% rownames(ase_poi_genes2), , drop = FALSE][complete.cases(counts_mat[rownames(counts_mat) %in% rownames(ase_poi_genes2), ]), , drop = FALSE], 1)),cluster_columns = F,cluster_rows = F, name = "Log2 (Allele counts)")
dev.off()

##  delta plot

testdat2 <- testdat2[!rownames(testdat2) %in% rownames(ase_poi_genes2),]

aseout <- all_res %>%
  reframe(n_pos = sum(adjpval < 0.05 & delta >=  0.58, na.rm = TRUE), n_neg = sum(adjpval < 0.05 & delta <= -0.58, na.rm = TRUE), .by = gene) %>%
  mutate(hit = (n_pos >= 4 & n_neg == 0) | (n_neg >= 4 & n_pos == 0), direction = case_when(hit & n_pos >= 4 ~ "positive", hit & n_neg >= 4 ~ "negative", TRUE ~ NA_character_)) %>%
  select(gene, hit, direction)

colnames(aseout)[2:3] <- c("H","H_FC")

all.CSvP <- read.csv(file = "CSvP all genes.csv",row.names = 1)
all.CSvP <- all.CSvP[c(1,5)]
all.CSvP$Sig <- F
all.CSvP[(abs(all.CSvP$logFC)> 0.58) & (all.CSvP$adj.P.Val < 0.05), ]$Sig  <- T
all.CSvP <- all.CSvP[c(-2)]
colnames(all.CSvP) <- c("P_FC", "P")
all.CSvP$gene <- rownames(all.CSvP)

all.CS_PvCSxP <- read.csv(file = "CS_PvCSxP all genes.csv",row.names = 1)
all.CS_PvCSxP <- all.CS_PvCSxP[c(1,5)]
all.CS_PvCSxP$Sig <- F
all.CS_PvCSxP[(abs(all.CS_PvCSxP$logFC)> 0.58) & (all.CS_PvCSxP$adj.P.Val < 0.05), ]$Sig  <- T
all.CS_PvCSxP <- all.CS_PvCSxP[c(-2)]
colnames(all.CS_PvCSxP) <- c("T_FC", "T")
all.CS_PvCSxP$gene <- rownames(all.CS_PvCSxP)

all_genes <- inner_join(aseout,all.CSvP,by="gene")
all_genes <- inner_join(all_genes,all.CS_PvCSxP,by="gene")

counts_HPT <- rbind(
  H = table(all_genes$H),
  P = table(all_genes$P),
  T = table(all_genes$T)
)
counts_HPT <- as.data.frame(as.table(counts_HPT))

pdf("Overall_patterns.pdf",height=3.5,width=5)
ggplot(counts_HPT, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col(position = "stack") +
  geom_text(aes(label = paste0(round(100 * Freq / ave(Freq, Var1, FUN = sum), 1), "%")),position = position_stack(vjust = 0.5),size = 3) +
  labs(x = NULL, y = "Number of Genes", fill = "") +
  scale_x_discrete(breaks = c("H", "P", "T"),labels = c(H = "ASE",P = "DE between\nparents",T = "DE in parents\n& hybrids")) +
  theme_minimal(base_size = 12)
dev.off()

classified <- all_genes %>%
  mutate(
    # convert to comparable signs
    p_sign = case_when(is.na(P_FC) ~ NA_integer_, P_FC > 0 ~  1L, P_FC < 0  ~ -1L, TRUE   ~  0L), 
    h_sign = case_when(H_FC == "positive" ~  1L, H_FC == "negative" ~ -1L, TRUE ~ NA_integer_)) %>%
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

classified_McManus <- classified
write.csv(classified_McManus,file="classified_McManus.csv",row.names = F)
classified_A <- classified %>% filter(grepl("^TraesCS[0-9]+A", gene))

counts_McManus_A <- as.data.frame(table(classified_A$category))
names(counts_McManus_A) <- c("category", "count")
props_McManus_A <- prop.table(table(classified_A$category))
df_McManus_A <- counts_McManus_A %>%
  mutate(prop = as.numeric(props_McManus_A[as.character(category)]), category = factor(category, levels = category))

pdf(file="McManus_classification_A_subgenome.pdf",height=3.5,width=4)
ggplot(df_McManus_A, aes(x = category, y = prop, fill = category)) +
  geom_col() +
  geom_text(aes(label = percent(prop, accuracy = 0.1)),vjust = -0.3, size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494"), guide = "none") +
  scale_x_discrete(labels = function(x) {
    x <- gsub("\\b[Cc]is\\b", "<i>cis</i>", x, perl = TRUE)
    x <- gsub("\\b[Tt]rans\\b", "<i>trans</i>", x, perl = TRUE)
    x
  }) +
  labs(x = NULL, y = "Proportion of genes", title = paste("A:",nrow(classified_A))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +
  coord_cartesian(clip = "off")
dev.off()

classified_B <- classified %>% filter(grepl("^TraesCS[0-9]+B", gene))

counts_McManus_B <- as.data.frame(table(classified_B$category))
names(counts_McManus_B) <- c("category", "count")
props_McManus_B <- prop.table(table(classified_B$category))
df_McManus_B <- counts_McManus_B %>%
  mutate(prop = as.numeric(props_McManus_B[as.character(category)]), category = factor(category, levels = category))

pdf(file="McManus_classification_B_subgenome.pdf",height=3.5,width=4)
ggplot(df_McManus_B, aes(x = category, y = prop, fill = category)) +
  geom_col() +
  geom_text(aes(label = percent(prop, accuracy = 0.1)),vjust = -0.3, size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494"), guide = "none") +
  scale_x_discrete(labels = function(x) {
    x <- gsub("\\b[Cc]is\\b", "<i>cis</i>", x, perl = TRUE)
    x <- gsub("\\b[Tt]rans\\b", "<i>trans</i>", x, perl = TRUE)
    x
  }) +
  labs(x = NULL, y = "Proportion of genes", title = paste("B:",nrow(classified_B))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +
  coord_cartesian(clip = "off")
dev.off()

classified_D <- classified %>% filter(grepl("^TraesCS[0-9]+D", gene))

counts_McManus_D <- as.data.frame(table(classified_D$category))
names(counts_McManus_D) <- c("category", "count")
props_McManus_D <- prop.table(table(classified_D$category))
df_McManus_D <- counts_McManus_D %>%
  mutate(prop = as.numeric(props_McManus_D[as.character(category)]), category = factor(category, levels = category))

pdf(file="McManus_classification_D_subgenome.pdf",height=3.5,width=4)
ggplot(df_McManus_D, aes(x = category, y = prop, fill = category)) +
  geom_col() +
  geom_text(aes(label = percent(prop, accuracy = 0.1)),vjust = -0.3, size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494"), guide = "none") +
  scale_x_discrete(labels = function(x) {
    x <- gsub("\\b[Cc]is\\b", "<i>cis</i>", x, perl = TRUE)
    x <- gsub("\\b[Tt]rans\\b", "<i>trans</i>", x, perl = TRUE)
    x
  }) +
  labs(x = NULL, y = "Proportion of genes", title = paste("D:",nrow(classified_D))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +
  coord_cartesian(clip = "off")
dev.off()

## now using limma for detecting ASE, comparing ref and alt counts

gene_id_vec <- mcols(rowRanges)$aseID

get_gene_counts_for_sample <- function(sname) {
  x <- df[[sname]]
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "0,0"
  x <- gsub("\\s+", "", x)
  parts <- strsplit(x, ",", fixed = TRUE)
  ref <- as.integer(vapply(parts, function(v) if (length(v) >= 1) v[1] else "0", FUN.VALUE = character(1)))
  alt <- as.integer(vapply(parts, function(v) if (length(v) >= 2) v[2] else "0", FUN.VALUE = character(1)))
  ref[is.na(ref)] <- 0L; alt[is.na(alt)] <- 0L
  r <- data.frame(ref = ref, alt = alt, stringsAsFactors = FALSE) # r$ref = CS allele; r$alt = P allele
  tibble(gene = gene_id_vec, ref = r$ref, alt = r$alt) %>%
    group_by(gene) %>%
    summarise(ref = sum(ref), alt = sum(alt), .groups = "drop") %>%
    mutate(sample = sname)
}

gene_counts_by_sample <- bind_rows(lapply(hyb_cols, get_gene_counts_for_sample))
gene_counts_by_sample <- gene_counts_by_sample[!(gene_counts_by_sample$gene %in% rownames(ase_poi_genes2)),] ## removing genes that have different effects depending on cross direction
gene_counts_by_sample <- gene_counts_by_sample %>%
  group_by(gene, sample) %>%
  summarise(ref = sum(ref, na.rm = TRUE),alt = sum(alt, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from  = sample,values_from = c(ref, alt),names_glue  = "{sample}_{.value}",values_fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames("gene")

## note, already tried to test if differ by direction of cross, but no significant hits
sample_info.edger <- factor(c( rep("ref", 5), rep("alt", 5)))
edgeR.DGElist <- DGEList(counts = gene_counts_by_sample, group = sample_info.edger)
keep <- rowSums( cpm(edgeR.DGElist) >= 2) >= 4
edgeR.DGElist <- edgeR.DGElist[keep,]
edgeR.DGElist$samples$lib.size <- colSums(edgeR.DGElist$counts)
edgeR.DGElist <- calcNormFactors(edgeR.DGElist, method = "TMM")
edgeR.DGElist$samples[,3] <- rep(calcNormFactors(DGEList(counts = gene_counts_by_sample[1:5]+gene_counts_by_sample[6:10]), method = "TMM")$samples[,3],2)
edgeR.DGElist$samples[,2] <- rep(calcNormFactors(DGEList(counts = gene_counts_by_sample[1:5]+gene_counts_by_sample[6:10]), method = "TMM")$samples[,2],2)
mm <- model.matrix(~0+edgeR.DGElist$samples$group, data = edgeR.DGElist$samples)
colnames(mm) <- levels(edgeR.DGElist$samples$group)
y <- voom(edgeR.DGElist, mm, plot = T)
fit <- lmFit(y, mm)

cpm_log <- cpm(edgeR.DGElist, log = TRUE)
cpm_nolog <- cpm(edgeR.DGElist, log = FALSE)
cpm_nolog_relative <- cpm_nolog/rowMeans(cpm_nolog)

Group <- edgeR.DGElist$samples[1]
Group <- as.factor(unlist(Group))
cpm_log_forpca <- cpm_log
pca <- prcomp(t(cpm_log_forpca), scale. = TRUE) ## do pca
## plot of pca with groups in ellipses
asepca <- fviz_pca_ind(pca,
                       col.ind = Group,
                       legend.title = "Allele",
                       repel = TRUE,
                       pointshape = 16,
                       pointsize  = 3,
                       mean.point = FALSE, 
                       title = ""
) + guides(color = guide_legend(override.aes = list(shape = 16, size = 3)))
asepca

pdf("ASE_PCA.pdf",height=3.5,width=4.5)
asepca
dev.off()

csscree <- fviz_screeplot(pca, ncp=10,title = "")
csscree

asetest <- eBayes(contrasts.fit(fit, contrast = c(1, -1))) ## alt upregulated, ref downregulated
top.table <- topTable(asetest, sort.by = "P", n = Inf) ## sort by most significantly DE genes
asetest_pvals <- topTable(asetest, sort.by = "none", n = Inf,p.value=1,lfc=0) ## get all genes with logFC and pvalues
length(which(top.table$adj.P.Val < 0.05)) ## how many significantly DE genes
write.csv(top.table,file="Ref_vs_Alt.csv")

DGEgenes <- rownames(subset(top.table, top.table$adj.P.Val < 0.05))
mat_DGEgenes <- cpm_nolog_relative[DGEgenes, ]

pdf("ASE_heatmap.pdf",height=3.5,width=4)
Heatmap(mat_DGEgenes, name = "Scaled CPM", show_row_names = FALSE, use_raster = F)
dev.off()

fit$design

asetest_pvals_sub <- asetest_pvals[c(1,5)]
asetest_pvals_sub$Sig <- F
asetest_pvals_sub[(abs(asetest_pvals_sub$logFC)> 0.58) & (asetest_pvals_sub$adj.P.Val < 0.05), ]$Sig  <- T
asetest_pvals_sub <- asetest_pvals_sub[c(-2)]
colnames(asetest_pvals_sub) <- c("H_FC", "H")
asetest_pvals_sub$gene <- rownames(asetest_pvals_sub)

all.CSvP <- read.csv(file = "CSvP all genes.csv",row.names = 1)
all.CSvP <- all.CSvP[c(1,5)]
all.CSvP$Sig <- F
all.CSvP[(abs(all.CSvP$logFC)> 0.58) & (all.CSvP$adj.P.Val < 0.05), ]$Sig  <- T
all.CSvP <- all.CSvP[c(-2)]
colnames(all.CSvP) <- c("P_FC", "P")
all.CSvP$gene <- rownames(all.CSvP)

all.CS_PvCSxP <- read.csv(file = "CS_PvCSxP all genes.csv",row.names = 1)
all.CS_PvCSxP <- all.CS_PvCSxP[c(1,5)]
all.CS_PvCSxP$Sig <- F
all.CS_PvCSxP[(abs(all.CS_PvCSxP$logFC)> 0.58) & (all.CS_PvCSxP$adj.P.Val < 0.05), ]$Sig  <- T
all.CS_PvCSxP <- all.CS_PvCSxP[c(-2)]
colnames(all.CS_PvCSxP) <- c("T_FC", "T")
all.CS_PvCSxP$gene <- rownames(all.CS_PvCSxP)

all_genes <- inner_join(asetest_pvals_sub,all.CSvP,by="gene")
all_genes <- inner_join(all_genes,all.CS_PvCSxP,by="gene")

counts_HPT <- rbind(
  H = table(all_genes$H),
  P = table(all_genes$P),
  T = table(all_genes$T)
)
counts_HPT <- as.data.frame(as.table(counts_HPT))

pdf("Overall_patterns_limma.pdf",height=3.5,width=5)
ggplot(counts_HPT, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col(position = "stack") +
  geom_text(aes(label = paste0(round(100 * Freq / ave(Freq, Var1, FUN = sum), 1), "%")),position = position_stack(vjust = 0.5),size = 3) +
  labs(x = NULL, y = "Number of Genes", fill = "") +
  scale_x_discrete(breaks = c("H", "P", "T"), labels = c(H = "ASE", P = "DE between\nparents", T = "DE in parents\n& hybrids")) +
  theme_minimal(base_size = 12)
dev.off()

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

classified_limma_A <- classified_limma %>% filter(grepl("^TraesCS[0-9]+A", gene))

counts_limma_A <- as.data.frame(table(classified_limma_A$category))
names(counts_limma_A) <- c("category", "count")
props_limma_A <- prop.table(table(classified_limma_A$category))
df_limma_A <- counts_limma_A %>%
  mutate(prop = as.numeric(props_limma_A[as.character(category)]), category = factor(category, levels = category))

pdf(file="limma_classification_A_subgenome.pdf",height=3.5,width=4)
ggplot(df_limma_A, aes(x = category, y = prop, fill = category)) +
  geom_col() +
  geom_text(aes(label = percent(prop, accuracy = 0.1)),vjust = -0.3, size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494"), guide = "none") +
  scale_x_discrete(labels = function(x) {
    x <- gsub("\\b[Cc]is\\b", "<i>cis</i>", x, perl = TRUE)
    x <- gsub("\\b[Tt]rans\\b", "<i>trans</i>", x, perl = TRUE)
    x
  }) +
  labs(x = NULL, y = "Proportion of genes", title = paste("A:",nrow(classified_limma_A))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +
  coord_cartesian(clip = "off")
dev.off()

classified_limma_B <- classified_limma %>% filter(grepl("^TraesCS[0-9]+B", gene))

counts_limma_B <- as.data.frame(table(classified_limma_B$category))
names(counts_limma_B) <- c("category", "count")
props_limma_B <- prop.table(table(classified_limma_B$category))
df_limma_B <- counts_limma_B %>%
  mutate(prop = as.numeric(props_limma_B[as.character(category)]), category = factor(category, levels = category))

pdf(file="limma_classification_B_subgenome.pdf",height=3.5,width=4)
ggplot(df_limma_B, aes(x = category, y = prop, fill = category)) +
  geom_col() +
  geom_text(aes(label = percent(prop, accuracy = 0.1)),vjust = -0.3, size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494"), guide = "none") +
  scale_x_discrete(labels = function(x) {
    x <- gsub("\\b[Cc]is\\b", "<i>cis</i>", x, perl = TRUE)
    x <- gsub("\\b[Tt]rans\\b", "<i>trans</i>", x, perl = TRUE)
    x
  }) +
  labs(x = NULL, y = "Proportion of genes", title = paste("B:",nrow(classified_limma_B))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +
  coord_cartesian(clip = "off")
dev.off()

classified_limma_D <- classified_limma %>% filter(grepl("^TraesCS[0-9]+D", gene))

counts_limma_D <- as.data.frame(table(classified_limma_D$category))
names(counts_limma_D) <- c("category", "count")
props_limma_D <- prop.table(table(classified_limma_D$category))
df_limma_D <- counts_limma_D %>%
  mutate(prop = as.numeric(props_limma_D[as.character(category)]), category = factor(category, levels = category))

pdf(file="limma_classification_D_subgenome.pdf",height=3.5,width=4)
ggplot(df_limma_D, aes(x = category, y = prop, fill = category)) +
  geom_col() +
  geom_text(aes(label = percent(prop, accuracy = 0.1)),vjust = -0.3, size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494"), guide = "none") +
  scale_x_discrete(labels = function(x) {
    x <- gsub("\\b[Cc]is\\b", "<i>cis</i>", x, perl = TRUE)
    x <- gsub("\\b[Tt]rans\\b", "<i>trans</i>", x, perl = TRUE)
    x
  }) +
  labs(x = NULL, y = "Proportion of genes", title = paste("D:",nrow(classified_limma_D))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +
  coord_cartesian(clip = "off")
dev.off()

classified_all <- inner_join(classified_McManus,classified_limma,by="gene")
dim(classified_all)

df_tab  <- as.data.frame(table(classified_all$category.x, classified_all$category.y))
names(df_tab) <- c("x_cat", "y_cat", "n")
ital_words <- function(x) {
  x <- gsub("\\b[Cc]is\\b", "<i>cis</i>", x, perl = TRUE)
  x <- gsub("\\b[Tt]rans\\b", "<i>trans</i>", x, perl = TRUE)
  x
}
p_tab <- ggplot(df_tab, aes(x = y_cat, y = x_cat, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n)) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  scale_x_discrete(labels = ital_words) +
  scale_y_discrete(labels = ital_words) +
  labs(y = "Beta-binomial classification", x = "Linear model classification", fill = "Count", title = "") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_markdown(angle = 30, hjust = 1),axis.text.y = element_markdown())

ggsave("category_crosstab_heatmap.pdf", p_tab, width = 7, height = 4, device = cairo_pdf)

classified_all <- inner_join(classified_McManus,classified_limma,by="gene")
dim(classified_all)

classified_all <- classified_all %>%
  mutate(category = ifelse(category.x == category.y, category.x, "Ambiguous"))
classified_all <- classified_all[c(1,4)]
classified_all_A <- classified_all %>% filter(grepl("^TraesCS[0-9]+A", gene))

counts_all_A <- as.data.frame(table(classified_all_A$category))
names(counts_all_A) <- c("category", "count")
props_all_A <- prop.table(table(classified_all_A$category))
df_all_A <- counts_all_A %>%
  mutate(prop = as.numeric(props_all_A[as.character(category)]), category = factor(category, levels = category))

pdf(file="all_classification_A_subgenome.pdf",height=3.5,width=4)
ggplot(df_all_A, aes(x = category, y = prop, fill = category)) +
  geom_col() +
  geom_text(aes(label = percent(prop, accuracy = 0.1)),vjust = -0.3, size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = c("#66C2A5","#E78AC3","#A6D854","#FFD92F","#E5C494"), guide = "none") +
  scale_x_discrete(labels = function(x) {
    x <- gsub("\\b[Cc]is\\b", "<i>cis</i>", x, perl = TRUE)
    x <- gsub("\\b[Tt]rans\\b", "<i>trans</i>", x, perl = TRUE)
    x
  }) +
  labs(x = NULL, y = "Proportion of genes", title = paste("A:",nrow(classified_all_A))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +
  coord_cartesian(clip = "off")
dev.off()

classified_all_B <- classified_all %>% filter(grepl("^TraesCS[0-9]+B", gene))

counts_all_B <- as.data.frame(table(classified_all_B$category))
names(counts_all_B) <- c("category", "count")
props_all_B <- prop.table(table(classified_all_B$category))
df_all_B <- counts_all_B %>%
  mutate(prop = as.numeric(props_all_B[as.character(category)]), category = factor(category, levels = category))

pdf(file="all_classification_B_subgenome.pdf",height=3.5,width=4)
ggplot(df_all_B, aes(x = category, y = prop, fill = category)) +
  geom_col() +
  geom_text(aes(label = percent(prop, accuracy = 0.1)),vjust = -0.3, size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = c("#66C2A5","#E78AC3","#A6D854","#FFD92F","#E5C494"), guide = "none") +
  scale_x_discrete(labels = function(x) {
    x <- gsub("\\b[Cc]is\\b", "<i>cis</i>", x, perl = TRUE)
    x <- gsub("\\b[Tt]rans\\b", "<i>trans</i>", x, perl = TRUE)
    x
  }) +
  labs(x = NULL, y = "Proportion of genes", title = paste("B:",nrow(classified_all_B))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +
  coord_cartesian(clip = "off")
dev.off()

classified_all_D <- classified_all %>% filter(grepl("^TraesCS[0-9]+D", gene))

counts_all_D <- as.data.frame(table(classified_all_D$category))
names(counts_all_D) <- c("category", "count")
props_all_D <- prop.table(table(classified_all_D$category))
df_all_D <- counts_all_D %>%
  mutate(prop = as.numeric(props_all_D[as.character(category)]), category = factor(category, levels = category))

pdf(file="all_classification_D_subgenome.pdf",height=3.5,width=4)
ggplot(df_all_D, aes(x = category, y = prop, fill = category)) +
  geom_col() +
  geom_text(aes(label = percent(prop, accuracy = 0.1)),vjust = -0.3, size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = c("#66C2A5","#E78AC3","#A6D854","#FFD92F","#E5C494"), guide = "none") +
  scale_x_discrete(labels = function(x) {
    x <- gsub("\\b[Cc]is\\b", "<i>cis</i>", x, perl = TRUE)
    x <- gsub("\\b[Tt]rans\\b", "<i>trans</i>", x, perl = TRUE)
    x
  }) +
  labs(x = NULL, y = "Proportion of genes", title = paste("D:",nrow(classified_all_D))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +
  coord_cartesian(clip = "off")
dev.off()

cowplot()
## gene ontology analysis

all_go_subset <- subset(all_go, Gene %in% classified_McManus$gene)
degenes <- classified_McManus$gene
names(degenes) <- degenes
degenes <- setNames(rep(0L, length(degenes)), names(degenes))
degenes[names(degenes) %in% classified_McManus[classified_McManus$category %in% "Conserved",]$gene] <- 1
gene_length_subsest <- gene_length[gene_length$V5 %in% names(degenes),]
len_lookup <- setNames(as.numeric(gene_length$V6), gene_length$V5)
degenes_len <- len_lookup[names(degenes)]
pwf = nullp(degenes, bias.data = degenes_len, plot.fit = TRUE)
GO.wall_conserved = goseq(pwf, gene2cat = all_go_subset[1:2])
GO.wall_conserved$over_rep_padj=p.adjust(GO.wall_conserved$over_represented_pvalue, method="BH")
GO.wall_conserved$under_rep_padj=p.adjust(GO.wall_conserved$under_represented_pvalue, method="BH")
GO.wall_conserved_over <- GO.wall_conserved[GO.wall_conserved$over_rep_padj < 0.05,]
GO.wall_conserved_under <- GO.wall_conserved[GO.wall_conserved$under_rep_padj < 0.05,]

# Over-represented in Conserved
go_over_conserved <- GO.wall_conserved_over %>%
  transmute(term_label = paste0(term, " (", category, ")"), ontology, numDEInCat, numInCat, padj = over_rep_padj) %>%
  filter(!is.na(padj)) %>%
  mutate(log10p = -log10(pmax(padj, .Machine$double.xmin)), gene_ratio = numDEInCat / numInCat) %>%
  arrange(desc(log10p)) %>%
  slice_head(n = 12)

p_over_conserved <- ggplot(go_over_conserved,aes(x = log10p, y = fct_reorder(term_label, log10p), size = numDEInCat, color = ontology)) +
  geom_point() +
  labs(title = "", x = expression(-log[10]("adjusted P-value")), y = NULL, size = "# DE genes", color = "Ontology") +
  theme_minimal(base_size = 12)

# Under-represented in Conserved
go_under_conserved <- GO.wall_conserved_under %>%
  transmute(term_label = paste0(term, " (", category, ")"), ontology, numDEInCat, numInCat, padj = under_rep_padj) %>%
  filter(!is.na(padj)) %>%
  mutate(log10p = -log10(pmax(padj, .Machine$double.xmin)),  gene_ratio = numDEInCat / numInCat) %>%
  arrange(desc(log10p)) %>%
  slice_head(n = 12)

p_under_conserved <- ggplot(go_under_conserved, aes(x = log10p, y = fct_reorder(term_label, log10p), size = numDEInCat, color = ontology)) +
  geom_point() +
  labs(title = "Under-represented in Conserved", x = expression(-log[10]("adjusted p-value")), y = NULL, size = "# DE genes", color = "Ontology") +
  theme_minimal(base_size = 12)

pdf("ASE_GO.pdf",height=3,width=6.5)
p_over_conserved
dev.off()

### do ASE analysis comparing with 50:50 expectation
ase_files <- c("CSxP1.wasp.ase.tsv", "CSxP2.wasp.ase.tsv", "CSxP3.wasp.ase.tsv", "PxCS1.wasp.ase.tsv", "PxCS2.wasp.ase.tsv")
read_one_ase <- function(fp) {
  smp <- sub(".wasp.ase.tsv", "", basename(fp), ignore.case = TRUE)
  tb <- readr::read_tsv(fp, col_types = cols(contig = col_character(),position = col_integer(),refCount = col_integer(),altCount = col_integer()),progress = FALSE) %>%
    mutate(locus = paste0(contig, "_", position),!!smp := paste0(pmax(refCount, 0L), ",",pmax(altCount, 0L))) %>%
    select(locus, !!smp)
  tb
}
ase_list <- lapply(ase_files, read_one_ase)
datum2c <- purrr::reduce(ase_list, full_join, by = "locus") %>%
  arrange(locus)
# Replace NAs with "0,0"
for (j in setdiff(colnames(datum2c), "locus")) {
  datum2c[[j]][is.na(datum2c[[j]])] <- "0,0"
}
datum2c <- as.data.frame(datum2c,check.names = FALSE)
rownames(datum2c) <- datum2c$locus
datum2c$locus <- NULL
datum2c <- datum2c %>%
  rownames_to_column("site") %>%
  anti_join(toremove %>% transmute(site = paste0(contig, "_", as.integer(position))) %>% distinct(), by = "site") %>%
  column_to_rownames("site")
## convert part names in variant calls back to whole for gtf compatability
part_sizes <- read.table("part_sizes.txt", header = FALSE, col.names = c("part", "size"))
part_sizes <- part_sizes[-c(nrow(part_sizes)),]
ps2 <- part_sizes %>%
  filter(grepl("_part[12]$", part)) %>%
  mutate(base = sub("_part[12]$", "", part), part_num = ifelse(grepl("_part2$", part), 2L, 1L))
p1 <- ps2 %>% filter(part_num == 1L) %>% select(base, size1 = size) %>% distinct()
out <- ps2 %>% left_join(p1, by = "base") %>%
  mutate(offset = ifelse(part_num == 1L, 0L, size1))
offset_lookup <- setNames(out$offset, out$part)
df    <- datum2c %>% tibble::rownames_to_column("locus")
m <- stringr::str_match(df$locus, "^(Chr[^_]+)_part([12])_(\\d+)$")
hit <- !is.na(m[,1])
seqname <- sub("_(\\d+)$", "", df$locus)
pos <- as.integer(sub("^.*_(\\d+)$", "\\1", df$locus))
seqname[hit] <- m[hit, 2]
part_tag <- paste0(m[hit, 2], "_part", m[hit, 3])
off  <- unname(offset_lookup[part_tag]); off[is.na(off)] <- 0L
pos[hit] <- pos[hit] + off
coord <- tibble(seqname, pos)
rowRanges <- GRanges(seqnames = coord$seqname, ranges   = IRanges(coord$pos, width = 1), aseID    = df$locus, allele1  = rep("REF", nrow(df)), allele2  = rep("ALT", nrow(df)))
names(rowRanges) <- paste0(coord$seqname, ":", coord$pos)
#  Map SNPs to genes
genes_gr <- rtracklayer::import("genes_refseqv2_HC.gff3")
genes_gr <- genes_gr[genes_gr$type == "gene"]
id_field <- intersect(c("ID","gene_id","Name"), colnames(mcols(genes_gr)))[1]
mcols(genes_gr)$gid <- as.character(mcols(genes_gr)[[id_field]])
hits <- findOverlaps(rowRanges, genes_gr, select = "first")
mcols(rowRanges)$aseID <- mcols(genes_gr)$gid[hits]
keep_overlap <- !is.na(mcols(rowRanges)$aseID)
rowRanges <- rowRanges[keep_overlap]
df <- df[keep_overlap, , drop = FALSE]
gene_ids  <- mcols(rowRanges)$aseID
feature_ids <- names(rowRanges)
## MBASED with 50:50
run_one_sample <- function(sample_col, serial = TRUE) {
  x <- df[[sample_col]]
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "0,0"
  x <- gsub("\\s+", "", x)
  parts <- strsplit(x, ",", fixed = TRUE)
  ref <- as.integer(vapply(parts, function(v) if (length(v) >= 1) v[1] else "0", FUN.VALUE = character(1)))
  alt <- as.integer(vapply(parts, function(v) if (length(v) >= 2) v[2] else "0", FUN.VALUE = character(1)))
  ref[is.na(ref)] <- 0L; alt[is.na(alt)] <- 0L
  rna <- data.frame(ref = ref, alt = alt, stringsAsFactors = FALSE)
  tot <- rna$ref + rna$alt
  ok <- is.finite(tot) & tot > 0
  rr   <- rowRanges[ok]
  rref <- rna$ref[ok]
  ralt <- rna$alt[ok]
  feats <- names(rr)
  gids  <- mcols(rr)$aseID
  # Fixed 50:50 expectation (no DNA baseline)
  #pr <- rep(meds[meds$sample %in% sample_col,]$median_prop, length(rr))
  pr <- rep(0.5, length(rr))
  
  SE <- SummarizedExperiment(
    assays = list(
      lociAllele1Counts = matrix(rref, ncol = 1, dimnames = list(feats, sample_col)),
      lociAllele2Counts = matrix(ralt, ncol = 1, dimnames = list(feats, sample_col)),
      lociAllele1CountsNoASEProbs = matrix(pr,   ncol = 1, dimnames = list(feats, sample_col))
    ),
    rowRanges = rr
  )
  
  set.seed(1)
  bp <- if (serial) BiocParallel::SerialParam() else BiocParallel::bpparam()
  # numSim=0 for speed; increase for full inference
  mb <- runMBASED(SE, isPhased = FALSE, numSim = 10^6, BPPARAM = bp)
  res <- tibble::tibble(gene = rownames(assays(mb)$pValueASE), sample = sample_col, pval = assays(mb)$pValueASE[, 1], adjpval = p.adjust(assays(mb)$pValueASE[, 1], method = "BH"), MAF = assays(mb)$majorAlleleFrequency[, 1])
  dir_df <- tibble::tibble( snv = feats, gene = gids, rref = rref, ralt = ralt) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarize(delta = log2(((sum(rref)+0.0000001) / (sum(ralt)+0.0000001) )), rref = sum(rref) , ralt = sum(ralt), .groups = "drop") %>%
    dplyr::mutate(sample = sample_col)
  out <- dplyr::left_join(res, dplyr::select(dir_df, gene, sample, delta, rref, ralt),by = c("gene", "sample"))
  return(out)
}

hyb_cols <-  c("CSxP1","CSxP2","CSxP3","PxCS1","PxCS2")
#all_res <- lapply(hyb_cols, run_one_sample) %>%
#  dplyr::bind_rows()
all_res <- as_tibble(read.csv(file="all_res_CS_0.5.csv"))
length(unique(all_res$gene)) ## number of testable genes

## Some genes vary depending on direction of the cross

delta_mat <- all_res %>%
  select(gene, sample, delta) %>%
  group_by(gene, sample) %>%
  summarise(delta = dplyr::first(delta), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = delta) %>% 
  column_to_rownames("gene") %>% 
  as.matrix()

testdat2 <- delta_mat[complete.cases(delta_mat),]
testdat2 <- testdat2[apply(testdat2,1,sd)>0,]
ase_poi_genes2 <- as.data.frame(testdat2)
ase_poi_genes2 <- rbind(ase_poi_genes2[ase_poi_genes2$PxCS1 < -20 & ase_poi_genes2$PxCS2 < -20 & ase_poi_genes2$CSxP1 > -3 & ase_poi_genes2$CSxP2 > -3 & ase_poi_genes2$CSxP3 > -3,],
                        ase_poi_genes2[ase_poi_genes2$PxCS1 > 20 & ase_poi_genes2$PxCS2 > 20 & ase_poi_genes2$CSxP1 < 3 & ase_poi_genes2$CSxP2 < 3 & ase_poi_genes2$CSxP3 < 3,],
                        ase_poi_genes2[ase_poi_genes2$PxCS1 > -3 & ase_poi_genes2$PxCS2 > -3 & ase_poi_genes2$CSxP1 < -20 & ase_poi_genes2$CSxP2 < -20 & ase_poi_genes2$CSxP3 < -20,],
                        ase_poi_genes2[ase_poi_genes2$PxCS1 < 3 & ase_poi_genes2$PxCS2 < 3 & ase_poi_genes2$CSxP1 > 20 & ase_poi_genes2$CSxP2 > 20 & ase_poi_genes2$CSxP3 > 20,])
all_res <- all_res[!(all_res$gene %in% rownames(ase_poi_genes2)),] ## remove those genes from further analyses
csxp <- c("CSxP1","CSxP2","CSxP3")
pxcs <- c("PxCS1","PxCS2")
poi_ase_genes <- all_res %>%
  arrange(gene, sample) %>%
  distinct(gene, sample, .keep_all = TRUE) %>%
  filter(sample %in% c(csxp, pxcs)) %>%
  group_by(gene) %>%
  summarise(
    csxp_n    = sum(sample %in% csxp & !is.na(delta) & !is.na(adjpval)),
    pxcs_n    = sum(sample %in% pxcs & !is.na(delta) & !is.na(adjpval)),
    csxp_sig  = all(adjpval[sample %in% csxp] < 0.05, na.rm = TRUE),
    pxcs_sig  = all(adjpval[sample %in% pxcs] < 0.05, na.rm = TRUE),
    csxp_pos  = all(delta[sample %in% csxp] >  0, na.rm = TRUE),
    csxp_neg  = all(delta[sample %in% csxp] < 0, na.rm = TRUE),
    pxcs_pos  = all(delta[sample %in% pxcs] >  0, na.rm = TRUE),
    pxcs_neg  = all(delta[sample %in% pxcs] < 0, na.rm = TRUE),
    .groups = "drop") %>%
  filter(pxcs_n == 2,csxp_n >= 2,csxp_sig, pxcs_sig,
         (csxp_pos & pxcs_neg) | (csxp_neg & pxcs_pos)) %>%
  mutate(pattern = dplyr::case_when(csxp_pos & pxcs_neg ~ "CSxP +  /  PxCS -", csxp_neg & pxcs_pos ~ "CSxP -  /  PxCS +")) %>%
  select(gene)
ase_poi_genes2 <- rbind(ase_poi_genes2,testdat2[rownames(testdat2) %in% poi_ase_genes$gene,])
dim(ase_poi_genes2)
all_res <- all_res[!(all_res$gene %in% rownames(ase_poi_genes2)),] ## remove those genes from further analyses

aseout <- all_res %>%
  reframe(n_pos = sum(adjpval < 0.05 & delta >=  0.58, na.rm = TRUE), n_neg = sum(adjpval < 0.05 & delta <= -0.58, na.rm = TRUE), .by = gene) %>%
  mutate(hit = (n_pos >= 4 & n_neg == 0) | (n_neg >= 4 & n_pos == 0), direction = case_when(hit & n_pos >= 4 ~ "positive", hit & n_neg >= 4 ~ "negative", TRUE ~ NA_character_)) %>%
  select(gene, hit, direction)
colnames(aseout)[2:3] <- c("H","H_FC")
all.CSvP <- read.csv(file = "CSvP all genes.csv",row.names = 1)
all.CSvP <- all.CSvP[c(1,5)]
all.CSvP$Sig <- F
all.CSvP[(abs(all.CSvP$logFC)> 0.58) & (all.CSvP$adj.P.Val < 0.05), ]$Sig  <- T
all.CSvP <- all.CSvP[c(-2)]
colnames(all.CSvP) <- c("P_FC", "P")
all.CSvP$gene <- rownames(all.CSvP)
all.CS_PvCSxP <- read.csv(file = "CS_PvCSxP all genes.csv",row.names = 1)
all.CS_PvCSxP <- all.CS_PvCSxP[c(1,5)]
all.CS_PvCSxP$Sig <- F
all.CS_PvCSxP[(abs(all.CS_PvCSxP$logFC)> 0.58) & (all.CS_PvCSxP$adj.P.Val < 0.05), ]$Sig  <- T
all.CS_PvCSxP <- all.CS_PvCSxP[c(-2)]
colnames(all.CS_PvCSxP) <- c("T_FC", "T")
all.CS_PvCSxP$gene <- rownames(all.CS_PvCSxP)
all_genes <- inner_join(aseout,all.CSvP,by="gene")
all_genes <- inner_join(all_genes,all.CS_PvCSxP,by="gene")
classified <- all_genes %>%
  mutate(
    # convert to comparable signs
    p_sign = case_when(is.na(P_FC) ~ NA_integer_, P_FC > 0 ~  1L, P_FC < 0  ~ -1L, TRUE   ~  0L), 
    h_sign = case_when(H_FC == "positive" ~  1L, H_FC == "negative" ~ -1L, TRUE ~ NA_integer_)) %>%
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
    TRUE ~ "Ambiguous"
  )
  ) %>%
  select(gene, category)
classified_McManus_p0.5 <- classified
classified_McManus_p0.5 <- inner_join(classified_McManus,classified_McManus_p0.5,by="gene")

df_tab2  <- as.data.frame(table(classified_McManus_p0.5$category.x, classified_McManus_p0.5$category.y))
names(df_tab2) <- c("x_cat", "y_cat", "n")
ital_words <- function(x) {
  x <- gsub("\\b[Cc]is\\b", "<i>cis</i>", x, perl = TRUE)
  x <- gsub("\\b[Tt]rans\\b", "<i>trans</i>", x, perl = TRUE)
  x
}
p_tab2 <- ggplot(df_tab2, aes(x = y_cat, y = x_cat, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n)) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  scale_x_discrete(labels = ital_words) +
  scale_y_discrete(labels = ital_words) +
  labs(y = "Median Reference Allele Proportion", x = "50:50 expectation", fill = "Count", title = "") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_markdown(angle = 30, hjust = 1),axis.text.y = element_markdown())

ggsave("ase_stats_compare_heatmap.pdf", p_tab2, width = 7, height = 4, device = cairo_pdf)

## Now doing ASE analyses with PAR reference

ase_files <- c("CSxP1.par.ase.tsv", "CSxP2.par.ase.tsv", "CSxP3.par.ase.tsv", "PxCS1.par.ase.tsv", "PxCS2.par.ase.tsv")

read_one_ase <- function(fp) {
  smp <- sub(".par.ase.tsv", "", basename(fp), ignore.case = TRUE)
  tb <- readr::read_tsv(fp, col_types = cols(contig = col_character(),position = col_integer(),refCount = col_integer(),altCount = col_integer()),progress = FALSE) %>%
    mutate(
      locus = paste0(contig, "_", position),!!smp := paste0(pmax(refCount, 0L), ",",pmax(altCount, 0L))) %>%
    select(locus, !!smp)
  tb
}
ase_list <- lapply(ase_files, read_one_ase)

datum2c <- purrr::reduce(ase_list, full_join, by = "locus") %>%
  arrange(locus)

# Replace NAs with "0,0"
for (j in setdiff(colnames(datum2c), "locus")) {
  datum2c[[j]][is.na(datum2c[[j]])] <- "0,0"
}

datum2c <- as.data.frame(datum2c,check.names = FALSE)
rownames(datum2c) <- datum2c$locus
datum2c$locus <- NULL

datum2c <- datum2c %>%
  rownames_to_column("site") %>%
  anti_join(toremove_par %>% transmute(site = paste0(contig, "_", as.integer(position))) %>% distinct(), by = "site") %>%
  column_to_rownames("site")

## convert part names in variant calls back to whole for gtf compatability
part_sizes <- read.table("Paragon_part_chr_sizes.txt", header = FALSE, col.names = c("part", "size"))
part_sizes <- part_sizes[-c(nrow(part_sizes)),]
ps2 <- part_sizes %>%
  filter(grepl("_part[12]$", part)) %>%
  mutate(base = sub("_part[12]$", "", part), part_num = ifelse(grepl("_part2$", part), 2L, 1L))
p1 <- ps2 %>% filter(part_num == 1L) %>% select(base, size1 = size) %>% distinct()
out <- ps2 %>% left_join(p1, by = "base") %>%
  mutate(offset = ifelse(part_num == 1L, 0L, size1))
offset_lookup <- setNames(out$offset, out$part)
df <- datum2c %>% tibble::rownames_to_column("locus")
m <- stringr::str_match(df$locus, "^([^_]+)_part([12])_(\\d+)$")
hit <- !is.na(m[,1])
seqname <- sub("_(\\d+)$", "", df$locus)
pos <- as.integer(sub("^.*_(\\d+)$", "\\1", df$locus))
seqname[hit] <- m[hit, 2]
part_tag <- paste0(m[hit, 2], "_part", m[hit, 3])
off <- unname(offset_lookup[part_tag]); off[is.na(off)] <- 0L
pos[hit] <- pos[hit] + off
coord <- tibble(seqname, pos)
rowRanges <- GRanges(seqnames = coord$seqname, ranges = IRanges(coord$pos, width = 1), aseID = df$locus, allele1  = rep("REF", nrow(df)), allele2 = rep("ALT", nrow(df)))
names(rowRanges) <- paste0(coord$seqname, ":", coord$pos)

#  Map SNPs to genes
genes_gr <- rtracklayer::import("genes_par.gff3")
genes_gr <- genes_gr[genes_gr$type == "gene"]
id_field <- intersect(c("ID","gene_id","Name"), colnames(mcols(genes_gr)))[1]
mcols(genes_gr)$gid <- as.character(mcols(genes_gr)[[id_field]])
hits <- findOverlaps(rowRanges, genes_gr, select = "first")
mcols(rowRanges)$aseID <- mcols(genes_gr)$gid[hits]
keep_overlap <- !is.na(mcols(rowRanges)$aseID)
rowRanges <- rowRanges[keep_overlap]
df <- df[keep_overlap, , drop = FALSE]
gene_ids  <- mcols(rowRanges)$aseID
feature_ids <- names(rowRanges)

## MBASED
run_one_sample <- function(sample_col, serial = TRUE) {
  x <- df[[sample_col]]
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "0,0"
  x <- gsub("\\s+", "", x)
  parts <- strsplit(x, ",", fixed = TRUE)
  ref <- as.integer(vapply(parts, function(v) if (length(v) >= 1) v[1] else "0", FUN.VALUE = character(1)))
  alt <- as.integer(vapply(parts, function(v) if (length(v) >= 2) v[2] else "0", FUN.VALUE = character(1)))
  ref[is.na(ref)] <- 0L; alt[is.na(alt)] <- 0L
  rna <- data.frame(ref = ref, alt = alt, stringsAsFactors = FALSE)
  tot <- rna$ref + rna$alt
  ok <- is.finite(tot) & tot > 0
  rr   <- rowRanges[ok]
  rref <- rna$ref[ok]
  ralt <- rna$alt[ok]
  feats <- names(rr)
  gids  <- mcols(rr)$aseID
  
  # comparing against median values
  pr <- rep(meds_par[meds_par$sample %in% sample_col,]$median_prop, length(rr))
  #pr <- rep(0.5, length(rr))
  
  SE <- SummarizedExperiment(
    assays = list(
      lociAllele1Counts = matrix(rref, ncol = 1, dimnames = list(feats, sample_col)),
      lociAllele2Counts = matrix(ralt, ncol = 1, dimnames = list(feats, sample_col)),
      lociAllele1CountsNoASEProbs = matrix(pr,   ncol = 1, dimnames = list(feats, sample_col))
    ),
    rowRanges = rr
  )
  
  set.seed(1)
  bp <- if (serial) BiocParallel::SerialParam() else BiocParallel::bpparam()
  # numSim=0 for speed; increase for full inference
  mb <- runMBASED(SE, isPhased = FALSE, numSim = 0, BPPARAM = bp)
  res <- tibble::tibble(gene = rownames(assays(mb)$pValueASE), sample = sample_col, pval = assays(mb)$pValueASE[, 1], adjpval = p.adjust(assays(mb)$pValueASE[, 1], method = "BH"), MAF = assays(mb)$majorAlleleFrequency[, 1])
  dir_df <- tibble::tibble( snv = feats, gene = gids, rref = rref, ralt = ralt) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarize(delta = log2(((sum(rref)+0.0000001) / (sum(ralt)+0.0000001) )), rref = sum(rref) , ralt = sum(ralt), .groups = "drop") %>%
    dplyr::mutate(sample = sample_col)
  out <- dplyr::left_join(res, dplyr::select(dir_df, gene, sample, delta, rref, ralt),by = c("gene", "sample"))
  return(out)
}

hyb_cols <-  c("CSxP1","CSxP2","CSxP3","PxCS1","PxCS2")
#all_res <- lapply(hyb_cols, run_one_sample) %>%
#  dplyr::bind_rows()
#all_res$gene <- gsub("gene:","", all_res$gene)

all_res <- as_tibble(read.csv(file="all_res_PAR.csv"))
length(unique(all_res$gene)) ## number of testable genes

pval_mat <- all_res %>%
  select(gene, sample, adjpval) %>%
  group_by(gene, sample) %>%
  summarise(pval = dplyr::first(adjpval), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = pval) %>% 
  column_to_rownames("gene") %>% 
  as.matrix()

counts_mat <- all_res %>% 
  select(gene, sample, c(rref, ralt)) %>% 
  pivot_longer(cols = c(rref, ralt), names_to = "allele", values_to = "count") %>% 
  pivot_wider(id_cols = gene, names_from = c(sample, allele), values_from = count, names_glue = "{sample}_{allele}")%>% 
  column_to_rownames("gene") %>% 
  as.matrix()

testdat <- pval_mat[complete.cases(pval_mat),]
testdat <- testdat[apply(testdat,1,sd)>0,]
pdf("ase_persample_pval_par.pdf",height=3.5,width=4)
Heatmap(testdat, name = "P-value", show_row_names = FALSE, use_raster = F)
dev.off()

## Some genes vary depending on direction of the cross
delta_mat <- all_res %>%
  select(gene, sample, delta) %>%
  dplyr::group_by(gene, sample) %>%
  dplyr::summarise(delta = dplyr::first(delta), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = delta) %>% 
  column_to_rownames("gene") %>% 
  as.matrix()

testdat2 <- delta_mat[complete.cases(delta_mat),]
testdat2 <- testdat2[(apply(testdat2,1,sd)>0),]

ase_poi_genes2 <- as.data.frame(testdat2)
ase_poi_genes2 <- rbind(ase_poi_genes2[ase_poi_genes2$PxCS1 < -20 & ase_poi_genes2$PxCS2 < -20 & ase_poi_genes2$CSxP1 > -3 & ase_poi_genes2$CSxP2 > -3 & ase_poi_genes2$CSxP3 > -3,],
ase_poi_genes2[ase_poi_genes2$PxCS1 > 20 & ase_poi_genes2$PxCS2 > 20 & ase_poi_genes2$CSxP1 < 3 & ase_poi_genes2$CSxP2 < 3 & ase_poi_genes2$CSxP3 < 3,],
ase_poi_genes2[ase_poi_genes2$PxCS1 > -3 & ase_poi_genes2$PxCS2 > -3 & ase_poi_genes2$CSxP1 < -20 & ase_poi_genes2$CSxP2 < -20 & ase_poi_genes2$CSxP3 < -20,],
ase_poi_genes2[ase_poi_genes2$PxCS1 < 3 & ase_poi_genes2$PxCS2 < 3 & ase_poi_genes2$CSxP1 > 20 & ase_poi_genes2$CSxP2 > 20 & ase_poi_genes2$CSxP3 > 20,])

all_res <- all_res[!(all_res$gene %in% rownames(ase_poi_genes2)),] ## remove those genes from further analyses
csxp <- c("CSxP1","CSxP2","CSxP3")
pxcs <- c("PxCS1","PxCS2")
poi_ase_genes <- all_res %>%
  arrange(gene, sample) %>%
  distinct(gene, sample, .keep_all = TRUE) %>%
  filter(sample %in% c(csxp, pxcs)) %>%
  group_by(gene) %>%
  summarise(
    csxp_n    = sum(sample %in% csxp & !is.na(delta) & !is.na(adjpval)),
    pxcs_n    = sum(sample %in% pxcs & !is.na(delta) & !is.na(adjpval)),
    csxp_sig  = all(adjpval[sample %in% csxp] < 0.05, na.rm = TRUE),
    pxcs_sig  = all(adjpval[sample %in% pxcs] < 0.05, na.rm = TRUE),
    csxp_pos  = all(delta[sample %in% csxp] >  0, na.rm = TRUE),
    csxp_neg  = all(delta[sample %in% csxp] < 0, na.rm = TRUE),
    pxcs_pos  = all(delta[sample %in% pxcs] >  0, na.rm = TRUE),
    pxcs_neg  = all(delta[sample %in% pxcs] < 0, na.rm = TRUE),
    .groups = "drop") %>%
  filter(pxcs_n == 2,csxp_n >= 2,csxp_sig, pxcs_sig,
         (csxp_pos & pxcs_neg) | (csxp_neg & pxcs_pos)) %>%
  mutate(pattern = dplyr::case_when(csxp_pos & pxcs_neg ~ "CSxP +  /  PxCS -", csxp_neg & pxcs_pos ~ "CSxP -  /  PxCS +")) %>%
  select(gene)
ase_poi_genes2 <- rbind(ase_poi_genes2,delta_mat[rownames(delta_mat) %in% poi_ase_genes$gene,])
all_res <- all_res[!(all_res$gene %in% rownames(ase_poi_genes2)),] ## remove those genes from further analyses

write.csv(ase_poi_genes2,file="ase_poi_genes2_par.csv")

pdf("ase_cross_direction_par.pdf",height=12,width=5.5)
Heatmap(log2(pmax(counts_mat[rownames(counts_mat) %in% rownames(ase_poi_genes2), , drop = FALSE][complete.cases(counts_mat[rownames(counts_mat) %in% rownames(ase_poi_genes2), ]), , drop = FALSE], 1)),cluster_columns = F,cluster_rows = F, name = "Log2 (Allele counts)")
dev.off()

##  delta plot

testdat2 <- testdat2[!rownames(testdat2) %in% rownames(ase_poi_genes2),]

aseout <- all_res %>%
  reframe(n_pos = sum(adjpval < 0.05 & delta >=  0.58, na.rm = TRUE), n_neg = sum(adjpval < 0.05 & delta <= -0.58, na.rm = TRUE), .by = gene) %>%
  mutate(hit = (n_pos >= 4 & n_neg == 0) | (n_neg >= 4 & n_pos == 0), direction = case_when(hit & n_pos >= 4 ~ "positive", hit & n_neg >= 4 ~ "negative", TRUE ~ NA_character_)) %>%
  select(gene, hit, direction)

colnames(aseout)[2:3] <- c("H","H_FC")

all.CSvP <- read.csv(file = "CSvP all genes_PAR.csv",row.names = 1)
all.CSvP <- all.CSvP[c(1,5)]
all.CSvP$Sig <- F
all.CSvP[(abs(all.CSvP$logFC)> 0.58) & (all.CSvP$adj.P.Val < 0.05), ]$Sig  <- T
all.CSvP <- all.CSvP[c(-2)]
colnames(all.CSvP) <- c("P_FC", "P")
all.CSvP$gene <- rownames(all.CSvP)

all.CS_PvCSxP <- read.csv(file = "CS_PvCSxP all genes_PAR.csv",row.names = 1)
all.CS_PvCSxP <- all.CS_PvCSxP[c(1,5)]
all.CS_PvCSxP$Sig <- F
all.CS_PvCSxP[(abs(all.CS_PvCSxP$logFC)> 0.58) & (all.CS_PvCSxP$adj.P.Val < 0.05), ]$Sig  <- T
all.CS_PvCSxP <- all.CS_PvCSxP[c(-2)]
colnames(all.CS_PvCSxP) <- c("T_FC", "T")
all.CS_PvCSxP$gene <- rownames(all.CS_PvCSxP)

all_genes <- inner_join(aseout,all.CSvP,by="gene")
all_genes <- inner_join(all_genes,all.CS_PvCSxP,by="gene")

counts_HPT <- rbind(
  H = table(all_genes$H),
  P = table(all_genes$P),
  T = table(all_genes$T)
)
counts_HPT <- as.data.frame(as.table(counts_HPT))

pdf("Overall_patterns_PAR.pdf",height=3.5,width=5)
ggplot(counts_HPT, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col(position = "stack") +
  geom_text(aes(label = paste0(round(100 * Freq / ave(Freq, Var1, FUN = sum), 1), "%")),position = position_stack(vjust = 0.5),size = 3) +
  labs(x = NULL, y = "Number of Genes", fill = "") +
  scale_x_discrete(breaks = c("H", "P", "T"), labels = c(H = "ASE", P = "DE between\nparents", T = "DE in parents\n& hybrids")) +
  theme_minimal(base_size = 12)
dev.off()

classified <- all_genes %>%
  mutate(
    # convert to comparable signs
    p_sign = case_when(is.na(P_FC) ~ NA_integer_, P_FC > 0 ~  1L, P_FC < 0  ~ -1L, TRUE   ~  0L), 
    h_sign = case_when(H_FC == "positive" ~  1L, H_FC == "negative" ~ -1L, TRUE ~ NA_integer_)) %>%
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

classified_McManus <- classified
write.csv(classified_McManus,file="classified_McManus_par.csv",row.names = F)
counts_McManus <- as.data.frame(table(classified$category))
names(counts_McManus) <- c("category", "count")

props_McManus <- prop.table(table(classified$category))

df_McManus <- counts_McManus %>%
  mutate(prop = as.numeric(props_McManus[as.character(category)]), category = factor(category, levels = category))

pdf(file="McManus_classification_par.pdf",height=3.5,width=4)
ggplot(df_McManus, aes(x = category, y = prop, fill = category)) +
  geom_col() +
  geom_text(aes(label = percent(prop, accuracy = 0.1)),vjust = -0.3, size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494"), guide = "none") +
  scale_x_discrete(labels = function(x) {
    x <- gsub("\\b[Cc]is\\b", "<i>cis</i>", x, perl = TRUE)
    x <- gsub("\\b[Tt]rans\\b", "<i>trans</i>", x, perl = TRUE)
    x
  }) +
  labs(x = NULL, y = "Proportion of genes", title = paste("n=",nrow(classified_McManus))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +
  coord_cartesian(clip = "off")
dev.off()

## now using limma for detecting ASE, comparing ref and alt counts

gene_id_vec <- mcols(rowRanges)$aseID

get_gene_counts_for_sample <- function(sname) {
  x <- df[[sname]]
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "0,0"
  x <- gsub("\\s+", "", x)
  parts <- strsplit(x, ",", fixed = TRUE)
  ref <- as.integer(vapply(parts, function(v) if (length(v) >= 1) v[1] else "0", FUN.VALUE = character(1)))
  alt <- as.integer(vapply(parts, function(v) if (length(v) >= 2) v[2] else "0", FUN.VALUE = character(1)))
  ref[is.na(ref)] <- 0L; alt[is.na(alt)] <- 0L
  r <- data.frame(ref = ref, alt = alt, stringsAsFactors = FALSE) # r$ref = CS allele; r$alt = P allele
  
  tibble(gene = gene_id_vec, ref = r$ref, alt = r$alt) %>%
    group_by(gene) %>%
    summarise(ref = sum(ref), alt = sum(alt), .groups = "drop") %>%
    mutate(sample = sname)
}

gene_counts_by_sample <- bind_rows(lapply(hyb_cols, get_gene_counts_for_sample))
gene_counts_by_sample$gene <- gsub("gene:","", gene_counts_by_sample$gene)
gene_counts_by_sample <- gene_counts_by_sample[!(gene_counts_by_sample$gene %in% rownames(ase_poi_genes2)),] ## removing genes that have different effects depending on cross direction
gene_counts_by_sample <- gene_counts_by_sample %>%
  group_by(gene, sample) %>%
  summarise(ref = sum(ref, na.rm = TRUE),alt = sum(alt, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from  = sample,values_from = c(ref, alt),names_glue  = "{sample}_{.value}",values_fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames("gene")

sample_info.edger <- factor(c( rep("ref", 5), rep("alt", 5)))
edgeR.DGElist <- DGEList(counts = gene_counts_by_sample, group = sample_info.edger)
keep <- rowSums( cpm(edgeR.DGElist) >= 2) >= 4
edgeR.DGElist <- edgeR.DGElist[keep,]
edgeR.DGElist$samples$lib.size <- colSums(edgeR.DGElist$counts)
edgeR.DGElist <- calcNormFactors(edgeR.DGElist, method = "TMM")
edgeR.DGElist$samples[,3] <- rep(calcNormFactors(DGEList(counts = gene_counts_by_sample[1:5]+gene_counts_by_sample[6:10]), method = "TMM")$samples[,3],2)
edgeR.DGElist$samples[,2] <- rep(calcNormFactors(DGEList(counts = gene_counts_by_sample[1:5]+gene_counts_by_sample[6:10]), method = "TMM")$samples[,2],2)
mm <- model.matrix(~0+edgeR.DGElist$samples$group, data = edgeR.DGElist$samples)
colnames(mm) <- levels(edgeR.DGElist$samples$group)
y <- voom(edgeR.DGElist, mm, plot = T)
fit <- lmFit(y, mm)

cpm_log <- cpm(edgeR.DGElist, log = TRUE)
cpm_nolog <- cpm(edgeR.DGElist, log = FALSE)
cpm_nolog_relative <- cpm_nolog/rowMeans(cpm_nolog)

Group <- edgeR.DGElist$samples[1]
Group <- as.factor(unlist(Group))
cpm_log_forpca <- cpm_log
pca <- prcomp(t(cpm_log_forpca), scale. = TRUE) ## do pca
## plot of pca with groups in ellipses
asepca <- fviz_pca_ind(pca,
                       col.ind = Group,
                       legend.title = "Allele",
                       repel = TRUE,
                       pointshape = 16,
                       pointsize  = 3,
                       mean.point = FALSE, 
                       title = ""
) + guides(color = guide_legend(override.aes = list(shape = 16, size = 3)))
asepca

pdf("ASE_PCA_par.pdf",height=3.5,width=4.5)
asepca
dev.off()

csscree <- fviz_screeplot(pca, ncp=10,title = "")
csscree

asetest <- eBayes(contrasts.fit(fit, contrast = c(1, -1))) ## alt upregulated, ref downregulated
top.table <- topTable(asetest, sort.by = "P", n = Inf) ## sort by most significantly DE genes
asetest_pvals <- topTable(asetest, sort.by = "none", n = Inf,p.value=1,lfc=0) ## get all genes with logFC and pvalues
length(which(top.table$adj.P.Val < 0.05)) ## how many significantly DE genes
write.csv(top.table,file="Ref_vs_Alt_par.csv")

DGEgenes <- rownames(subset(top.table, top.table$adj.P.Val < 0.05))
mat_DGEgenes <- cpm_nolog_relative[DGEgenes, ]

pdf("ASE_heatmap_par.pdf",height=3.5,width=4)
Heatmap(mat_DGEgenes, name = "Scaled CPM", show_row_names = FALSE, use_raster = F)
dev.off()

fit$design

asetest_pvals_sub <- asetest_pvals[c(1,5)]
asetest_pvals_sub$Sig <- F
asetest_pvals_sub[(abs(asetest_pvals_sub$logFC)> 0.58) & (asetest_pvals_sub$adj.P.Val < 0.05), ]$Sig  <- T
asetest_pvals_sub <- asetest_pvals_sub[c(-2)]
colnames(asetest_pvals_sub) <- c("H_FC", "H")
asetest_pvals_sub$gene <- rownames(asetest_pvals_sub)

all.CSvP <- read.csv(file = "CSvP all genes_PAR.csv",row.names = 1)
all.CSvP <- all.CSvP[c(1,5)]
all.CSvP$Sig <- F
all.CSvP[(abs(all.CSvP$logFC)> 0.58) & (all.CSvP$adj.P.Val < 0.05), ]$Sig  <- T
all.CSvP <- all.CSvP[c(-2)]
colnames(all.CSvP) <- c("P_FC", "P")
all.CSvP$gene <- rownames(all.CSvP)

all.CS_PvCSxP <- read.csv(file = "CS_PvCSxP all genes_PAR.csv",row.names = 1)
all.CS_PvCSxP <- all.CS_PvCSxP[c(1,5)]
all.CS_PvCSxP$Sig <- F
all.CS_PvCSxP[(abs(all.CS_PvCSxP$logFC)> 0.58) & (all.CS_PvCSxP$adj.P.Val < 0.05), ]$Sig  <- T
all.CS_PvCSxP <- all.CS_PvCSxP[c(-2)]
colnames(all.CS_PvCSxP) <- c("T_FC", "T")
all.CS_PvCSxP$gene <- rownames(all.CS_PvCSxP)

all_genes <- inner_join(asetest_pvals_sub,all.CSvP,by="gene")
all_genes <- inner_join(all_genes,all.CS_PvCSxP,by="gene")

counts_HPT <- rbind(
  H = table(all_genes$H),
  P = table(all_genes$P),
  T = table(all_genes$T)
)
counts_HPT <- as.data.frame(as.table(counts_HPT))

pdf("Overall_patterns_limma_par.pdf",height=3.5,width=5)
ggplot(counts_HPT, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col(position = "stack") +
  geom_text(aes(label = paste0(round(100 * Freq / ave(Freq, Var1, FUN = sum), 1), "%")),position = position_stack(vjust = 0.5),size = 3) +
  labs(x = NULL, y = "Number of Genes", fill = "") +
  scale_x_discrete(breaks = c("H", "P", "T"), labels = c(H = "ASE", P = "DE between\nparents", T = "DE in parents\n& hybrids")) +
  theme_minimal(base_size = 12)
dev.off()

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
write.csv(classified_limma,file="classified_limma_par.csv",row.names = F)
counts_limma <- as.data.frame(table(classified$category))
names(counts_limma) <- c("category", "count")
props_limma <- prop.table(table(classified$category))
df_limma <- counts_limma %>%
  mutate(prop = as.numeric(props_limma[as.character(category)]), category = factor(category, levels = category))

pdf(file="limma_classification_par.pdf",height=3.5,width=4)
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

classified_all <- inner_join(classified_McManus,classified_limma,by="gene")
df_tab  <- as.data.frame(table(classified_all$category.x, classified_all$category.y))
names(df_tab) <- c("x_cat", "y_cat", "n")
ital_words <- function(x) {
  x <- gsub("\\b[Cc]is\\b", "<i>cis</i>", x, perl = TRUE)
  x <- gsub("\\b[Tt]rans\\b", "<i>trans</i>", x, perl = TRUE)
  x
}
p_tab <- ggplot(df_tab, aes(x = y_cat, y = x_cat, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n)) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  scale_x_discrete(labels = ital_words) +
  scale_y_discrete(labels = ital_words) +
  labs(y = "McManus et al. classification", x = "limma classification", fill = "Count", title = "") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_markdown(angle = 30, hjust = 1),axis.text.y = element_markdown())

ggsave("category_crosstab_heatmap_par.pdf", p_tab, width = 7, height = 4, device = cairo_pdf)

classified_all <- inner_join(classified_McManus,classified_limma,by="gene")
dim(classified_all)
classified_all <- classified_all %>%
  mutate(category = ifelse(category.x == category.y, category.x, "Ambiguous"))
classified_all <- classified_all[c(1,4)]
counts_all <- as.data.frame(table(classified_all$category))
names(counts_all) <- c("category", "count")
props_all <- prop.table(table(classified_all$category))
df_all <- counts_all %>%
  mutate(prop = as.numeric(props_all[as.character(category)]), category = factor(category, levels = category))
write.csv(classified_all,file="classified_all_par.csv",row.names = F)

pdf(file="McManus_limma_classification_par.pdf",height=3.5,width=4)
ggplot(df_all, aes(x = category, y = prop, fill = category)) +
  geom_col() +
  geom_text(aes(label = percent(prop, accuracy = 0.1)),vjust = -0.3, size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = c("#66C2A5","#E78AC3","#A6D854","#FFD92F","#E5C494"), guide = "none") +
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

## files for 1:1 mapping
SingleCopyOrthologues <- read.table(file="SingleCopyOrthologues_matrix.tsv")
SingleCopyOrthologues$V3 <- gsub("\\..*","",SingleCopyOrthologues$V3 )
SingleCopyOrthologues$V2 <- sub("^([^.]*\\.[^.]*)\\..*$", "\\1", SingleCopyOrthologues$V2)
SingleCopyOrthologues_unique <- SingleCopyOrthologues %>%
  add_count(V2, name = "nV2") %>%
  add_count(V3, name = "nV3") %>%
  filter(nV2 == 1, nV3 == 1) %>%
  select(-nV2, -nV3)

classified_all <- read.csv(file="classified_all.csv")
classified_all <- classified_all[classified_all$gene %in% SingleCopyOrthologues_unique$V3,]
classified_all_par <- read.csv(file="classified_all_par.csv")
classified_all_par <- classified_all_par[classified_all_par$gene %in% SingleCopyOrthologues_unique$V2,]
colnames(classified_all)[1] <- "V3"
colnames(classified_all_par)[1] <- "V2"
classified_all <- left_join(classified_all,SingleCopyOrthologues_unique[2:3],by="V3")
classified_all <- inner_join(classified_all,classified_all_par,by="V2")
table(classified_all$category.x,classified_all$category.y)

df_tab  <- as.data.frame(table(classified_all$category.x, classified_all$category.y))
names(df_tab) <- c("x_cat", "y_cat", "n")
ital_words <- function(x) {
  x <- gsub("\\b[Cc]is\\b", "<i>cis</i>", x, perl = TRUE)
  x <- gsub("\\b[Tt]rans\\b", "<i>trans</i>", x, perl = TRUE)
  x
}
p_tab <- ggplot(df_tab, aes(x = y_cat, y = x_cat, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n)) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  scale_x_discrete(labels = ital_words) +
  scale_y_discrete(labels = ital_words) +
  labs(y = "CS-based classification", x = "Paragon-based classification", fill = "Count", title = "") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_markdown(angle = 30, hjust = 1),axis.text.y = element_markdown())

ggsave("category_crosstab_heatmap_both_ref.pdf", p_tab, width = 7, height = 4, device = cairo_pdf)

classified_all <- classified_all[classified_all$category.x == classified_all$category.y,]
table(classified_all$category.x,classified_all$category.y)

counts_all <- as.data.frame(table(classified_all$category.x))
names(counts_all) <- c("category", "count")
props_all <- prop.table(table(classified_all$category.x))
df_all <- counts_all %>%
  mutate(prop = as.numeric(props_all[as.character(category)]), category = factor(category, levels = category))

pdf(file="Classification_both_ref.pdf",height=3.5,width=4)
ggplot(df_all, aes(x = category, y = prop, fill = category)) +
  geom_col() +
  geom_text(aes(label = percent(prop, accuracy = 0.1)),vjust = -0.3, size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = c("#66C2A5","#E78AC3","#A6D854","#FFD92F","#E5C494"), guide = "none") +
  scale_x_discrete(labels = function(x) {
    x <- gsub("\\b[Cc]is\\b", "<i>cis</i>", x, perl = TRUE)
    x <- gsub("\\b[Tt]rans\\b", "<i>trans</i>", x, perl = TRUE)
    x
  }) +
  labs(x = NULL, y = "Proportion of genes", title = paste("n=",nrow(classified_all))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +
  coord_cartesian(clip = "off")
dev.off()

classified_all <- classified_all[c(1,3,4)]
colnames(classified_all) <- c("CS_reference","Paragon_reference","Category")
write.table(classified_all,file="Regulatory_Classification.txt",row.names = F)

all_res <- read.csv(file="all_res_CS.csv")
all_res <- all_res[all_res$gene %in% SingleCopyOrthologues_unique$V3,]
all_res <- all_res %>%
  group_by(sample) %>%
  mutate(adjpval = p.adjust(pval, method = "BH")) %>%
  ungroup()
colnames(all_res) <- gsub("$","_CS_ref", colnames(all_res))
all_res <- inner_join(all_res,SingleCopyOrthologues_unique,by = c("gene_CS_ref" = "V3"))
all_res2 <- read.csv(file="ASE_test_par.csv")
all_res2 <- all_res2[all_res2$gene %in% SingleCopyOrthologues_unique$V2,]
all_res2 <- all_res2 %>%
  group_by(sample) %>%
  mutate(adjpval = p.adjust(pval, method = "BH")) %>%
  ungroup()
colnames(all_res2) <- gsub("$","_PAR_ref", colnames(all_res2))
all_res2 <- inner_join(all_res2,SingleCopyOrthologues_unique,by = c("gene_PAR_ref" = "V2"))
all_res_both <- inner_join(all_res,all_res2,by = c(c("gene_CS_ref" = "V3"),c("sample_CS_ref" = "sample_PAR_ref")) )
dim(all_res_both[all_res_both$adjpval_CS_ref < 0.05 & all_res_both$adjpval_PAR_ref < 0.05,])
dim(all_res_both)
write.csv(all_res_both,file="ASE_test_one_one_mapping.csv", row.names = F)

p_ase_compare1 <- ggscatter(data=all_res_both,x="rref_CS_ref",y="ralt_PAR_ref", add = "reg.line", conf.int = F, cor.coef = TRUE, cor.method = "pearson",shape = 1)+
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(x="CS allele in CS reference",y="CS allele in Paragon reference")
p_ase_compare2 <-ggscatter(data=all_res_both,x="rref_PAR_ref",y="ralt_CS_ref", add = "reg.line", conf.int = F, cor.coef = TRUE, cor.method = "pearson",shape = 1)+
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, linetype = 2)+
  labs(x="Paragon allele in Paragon reference",y="Paragon allele in CS reference")
p_ase_compare3 <-ggscatter(data=all_res_both[all_res_both$delta_CS_ref < 20 & all_res_both$delta_PAR_ref > -20,],x="delta_CS_ref",y="delta_PAR_ref", add = "reg.line", conf.int = F, cor.coef = TRUE, cor.method = "pearson",shape = 1) +
  geom_abline(slope = -1, intercept = 0, linetype = 2)
p_ase_compare4 <-ggplot(all_res_both %>% filter(pval_CS_ref > 0, pval_PAR_ref > 0) %>% transmute(x = -log10(pval_CS_ref), y = -log10(pval_PAR_ref)), aes(x, y)) +
  geom_point(shape = 1) +
  stat_cor(method = "pearson",label.x.npc = "right", label.y.npc = "top",hjust = 1, vjust = 1, size = 4) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(x = "-log10(CS_ref p-value)", y = "-log10(PAR_ref p-value)")

pdf("ASE_compare_McManus.pdf",height=7,width=7)
plot_grid(p_ase_compare1,p_ase_compare2,p_ase_compare3,p_ase_compare4,ncol=2)
dev.off()

Ref_vs_Alt <- read.csv(file="Ref_vs_Alt.csv",row.names = 1)
Ref_vs_Alt_par <- read.csv(file="Ref_vs_Alt_par.csv",row.names = 1)
Ref_vs_Alt <- Ref_vs_Alt %>% rownames_to_column(var = "CS_id")
Ref_vs_Alt <- inner_join( Ref_vs_Alt, SingleCopyOrthologues_unique,by = c("CS_id" = "V3"))
Ref_vs_Alt_par <- Ref_vs_Alt_par %>% rownames_to_column(var = "PAR_id")
Ref_vs_Alt_par <- inner_join(Ref_vs_Alt_par, SingleCopyOrthologues_unique,by = c("PAR_id" = "V2"))
Ref_vs_Alt <- inner_join(Ref_vs_Alt,Ref_vs_Alt_par,by = c("CS_id" = "V3"))
dim(Ref_vs_Alt)
dim(Ref_vs_Alt[Ref_vs_Alt$adj.P.Val.x < 0.05 & Ref_vs_Alt$adj.P.Val.y < 0.05,])

p_ase_limma_compare1 <- ggscatter(data=Ref_vs_Alt,x="AveExpr.x",y="AveExpr.y", add = "reg.line", conf.int = F, cor.coef = TRUE, cor.method = "pearson",shape = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(x = "Mean expression using CS reference", y = "Mean expression using\nParagon reference")
p_ase_limma_compare2 <- ggscatter(data=Ref_vs_Alt,x="logFC.x",y="logFC.y", add = "reg.line", conf.int = F, cor.coef = TRUE, cor.method = "pearson",shape = 1) +
  geom_abline(slope = -1, intercept = 0, linetype = 2) +
  labs(x = "log2FC using CS reference", y = "log2FC using Paragon reference")
p_ase_limma_compare3 <-ggplot(Ref_vs_Alt %>% filter(P.Value.x > 0, P.Value.y > 0) %>% transmute(x = -log10(P.Value.x), y = -log10(P.Value.y)), aes(x, y)) +
  geom_point(shape = 1) +
  stat_cor(method = "pearson",label.x.npc = "right", label.y.npc = "top",hjust = 1, vjust = 1, size = 4) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(x = "-log10(CS_ref p-value)", y = "-log10(PAR_ref p-value)")
pdf("ASE_compare_limma.pdf",height=3.5,width=10)
plot_grid(p_ase_limma_compare1,p_ase_limma_compare2,p_ase_limma_compare3,ncol=3,labels="AUTO")
dev.off()

## get cpm estimates for ASE genes, first with CS ref
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
cpm_log <- cpm(edgeR.DGElist, log = TRUE)
cpm_nolog <- cpm(edgeR.DGElist, log = FALSE)
colnames(cpm_log) <- sub("_.*","",colnames(cpm_log))
colnames(cpm_log) <- sub("_.*","",colnames(cpm_log))
cpm_nolog_relative <- cpm_nolog/rowMeans(cpm_nolog)
colnames(cpm_nolog_relative) <- sub("_.*","",colnames(cpm_nolog_relative))
colnames(cpm_nolog_relative) <- sub("_.*","",colnames(cpm_nolog_relative))
colnames(cpm_nolog) <- sub("_.*","",colnames(cpm_nolog))
colnames(cpm_nolog) <- sub("_.*","",colnames(cpm_nolog))
colnames(cpm_nolog_relative) <- gsub("PXCS2","PxCS2",colnames(cpm_nolog_relative))
pdf(file="Trans_cpm.pdf",width=4,height=2.5)
Heatmap(cpm_nolog_relative[rownames(cpm_nolog_relative) %in% classified_all[classified_all$category.x == "Trans only",]$V3,],
        show_row_names = F, cluster_columns = T, name = "Scaled CPM",  use_raster = F,row_title = "Trans only Genes")
dev.off()
pdf(file="Cis_cpm.pdf",width=4,height=2.5)
Heatmap(cpm_nolog_relative[rownames(cpm_nolog_relative) %in% classified_all[classified_all$category.x == "Cis only",]$V3,],
        show_row_names = F, cluster_columns = T, name = "Scaled CPM",  use_raster = F,row_title = "Cis only Genes")
dev.off()

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
cpm_log <- cpm(edgeR.DGElist, log = TRUE)
cpm_nolog <- cpm(edgeR.DGElist, log = FALSE)
colnames(cpm_log) <- sub("_.*","",colnames(cpm_log))
colnames(cpm_log) <- sub("_.*","",colnames(cpm_log))
cpm_nolog_relative <- cpm_nolog/rowMeans(cpm_nolog)
colnames(cpm_nolog_relative) <- sub("_.*","",colnames(cpm_nolog_relative))
colnames(cpm_nolog_relative) <- sub("_.*","",colnames(cpm_nolog_relative))
colnames(cpm_nolog) <- sub("_.*","",colnames(cpm_nolog))
colnames(cpm_nolog) <- sub("_.*","",colnames(cpm_nolog))
colnames(cpm_nolog_relative) <- gsub("PXCS2","PxCS2",colnames(cpm_nolog_relative))
pdf(file="Trans_cpm_par.pdf",width=4,height=2.5)
Heatmap(cpm_nolog_relative[rownames(cpm_nolog_relative) %in% classified_all[classified_all$category.y == "Trans only",]$V2,],
        show_row_names = F, cluster_columns = T, name = "Scaled CPM",  use_raster = F,row_title = "Trans only Genes")
dev.off()
pdf(file="Cis_cpm_par.pdf",width=4,height=2.5)
Heatmap(cpm_nolog_relative[rownames(cpm_nolog_relative) %in% classified_all[classified_all$category.y == "Cis only",]$V2,],
        show_row_names = F, cluster_columns = T, name = "Scaled CPM",  use_raster = F,row_title = "Cis only Genes")
dev.off()

##  harper dataset for cis only and trans only genes

read.counts <- read.table("cs_csxp_count.tsv", header = TRUE) 
read.counts <- read.counts[!grepl("LC$", rownames(read.counts)), , drop = FALSE]
read.counts <- read.counts[-c(12)]
CSxPvPxCS_sig_genes <- read.csv(file = "CSxPvPxCS sig genes.csv")
read.counts <- read.counts[!rownames(read.counts) %in% CSxPvPxCS_sig_genes$X,]
sample_info.edger <- factor(c( rep("CS", 3), rep("CSxP", 3), rep("P", 3), rep("CSxP", 2),rep("CS", 1),rep("F2", 46),rep("P",1))) ### treatment as grouping variables
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
colnames(cpm_nolog_relative) <- gsub("PXCS2","PxCS2",colnames(cpm_nolog_relative))

D <- cpm_nolog[rownames(cpm_nolog) %in% classified_all[classified_all$category.x == "Cis only",]$V3,]
D2 <- cbind(log2(abs((rowMeans(D[,c(4:6,10:11)]) - rowMeans(D[,1:3]) ) / (rowMeans(D[,c(4:6,10:11)]) - rowMeans(D[,7:9])))),"Cis only")
D <- cpm_nolog[rownames(cpm_nolog) %in% classified_all[classified_all$category.x == "Trans only",]$V3,]
D2 <- rbind(D2,cbind(log2(abs((rowMeans(D[,c(4:6,10:11)]) - rowMeans(D[,1:3]) ) / (rowMeans(D[,c(4:6,10:11)]) - rowMeans(D[,7:9])))),"Trans only"))
D2 <- as.data.frame(D2)
D2 <- D2 %>% mutate(V1 = as.numeric(V1))
counts <- D2 %>%
  group_by(V2) %>%
  summarise(n = n(), .groups = "drop")
pdf(file="cis_trans_spread_F1.pdf",height=2.5,width=3)
ggplot(D2, aes(x = V2, y = V1)) +
  geom_violin() +
  scale_x_discrete(labels = c("Cis only"   = expression(italic(Cis)~"only"),"Trans only" = expression(italic(Trans)~"only"))) +
  geom_text(data = counts,aes(x = V2, y = 11, label = paste0("n=", n)),fontface = "bold") +
  labs(x="Genes",y = expression(log[2]~frac(plain("Difference from ")*italic(CS),plain("Difference from ")*italic(Paragon)) ))
dev.off()

D <- cpm_nolog[rownames(cpm_nolog) %in% classified_all[classified_all$category.x == "Cis only",]$V3,]
D <- D[,-c(1:11)]
D2 <- cbind(log2(abs((rowMeans(D[,c(2:47)]) - D[,1]) / (rowMeans(D[,c(2:47)]) - D[,48]))),"Cis only")
D <- cpm_nolog[rownames(cpm_nolog) %in% classified_all[classified_all$category.x == "Trans only",]$V3,]
D2 <- rbind(D2,cbind(log2(abs((rowMeans(D[,c(2:47)]) - D[,1]) / (rowMeans(D[,c(2:47)]) - D[,48]))),"Trans only"))
D2 <- as.data.frame(D2)
D2 <- D2 %>% mutate(V1 = as.numeric(V1))
counts <- D2 %>%
  group_by(V2) %>%
  summarise(n = n(), .groups = "drop")
pdf(file="cis_trans_spread_F2.pdf",height=2.5,width=3)
ggplot(D2, aes(x = V2, y = V1)) +
  geom_violin() +
  scale_x_discrete(labels = c("Cis only"   = expression(italic(Cis)~"only"),"Trans only" = expression(italic(Trans)~"only"))) +
  geom_text(data = counts,aes(x = V2, y = 12, label = paste0("n=", n)),fontface = "bold") +
  labs(x="Genes",y = expression(log[2]~frac(plain("Difference from ")*italic(CS),plain("Difference from ")*italic(Paragon)) ))
dev.off()
