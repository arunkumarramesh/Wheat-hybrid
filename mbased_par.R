library(readr)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(rtracklayer)
library(GenomicRanges)
library(IRanges)
library(S4Vectors)
library(SummarizedExperiment)
library(BiocParallel)
library(MBASED)


files <- c("P1.par.ase.tsv", "P2.par.ase.tsv", "P3.par.ase.tsv", "CS1.par.ase.tsv", "CS2.par.ase.tsv", "CS3.par.ase.tsv")

all_ase_parents <- map_dfr(files, function(f) {
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(sample = sub(".par.ase.tsv", "", basename(f), perl = TRUE))
})

toremove_par <- all_ase_parents[(all_ase_parents$sample %in% c("CS1","CS2","CS3") ) & (all_ase_parents$refCount/all_ase_parents$totalCount > 0.6),]


files <- c("CSxP1.par.ase.tsv", "CSxP2.par.ase.tsv", "CSxP3.par.ase.tsv", "PxCS1.par.ase.tsv", "PxCS2.par.ase.tsv")

all_ase <- map_dfr(files, function(f) {
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(sample = sub(".par.ase.tsv", "", basename(f), perl = TRUE))
})

all_ase_filtered <- all_ase %>%
  anti_join(toremove_par %>% distinct(contig, position), by = c("contig", "position"))

all_ase_filtered %>% group_by(sample) %>% summarise(median_prop = median(refCount/totalCount))
meds_par <- all_ase_filtered %>% group_by(sample) %>% summarise(median_prop = median(refCount/totalCount))


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
      lociAllele1Counts           = matrix(rref, ncol = 1, dimnames = list(feats, sample_col)),
      lociAllele2Counts           = matrix(ralt, ncol = 1, dimnames = list(feats, sample_col)),
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
all_res <- lapply(hyb_cols, run_one_sample) %>%
  dplyr::bind_rows()
all_res$gene <- gsub("gene:","", all_res$gene)

write.csv(all_res,file = "all_res_PAR.csv",row.names = F)
