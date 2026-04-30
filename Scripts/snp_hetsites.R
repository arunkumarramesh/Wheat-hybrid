
library(vcfR)
library(factoextra)
library(pheatmap)
library(Polychrome)

datum <- read.vcfR(file="wheat_sites_merged_cs.vcf")
datum <- as.data.frame(extract.gt(datum))
colnames(datum) <- gsub("PXCS2","PxCS2", gsub(".sort.cs.bam","", gsub("_.*","", colnames(datum))))
chrids_datum <- rownames(datum)
datum[datum == "0/0"] <- 0
datum[datum == "0|0"] <- 0
datum[datum == "0/1"] <- 1
datum[datum == "0|1"] <- 1
datum[datum == "1|0"] <- 1
datum[datum == "0/2"] <- 9
datum[datum == "0|2"] <- 9
datum[datum == "0/3"] <- 9
datum[datum == "0|3"] <- 9
datum[datum == "1/2"] <- 9
datum[datum == "1|2"] <- 9
datum[datum == "1/3"] <- 9
datum[datum == "1|3"] <- 9
datum[datum == "2/3"] <- 9
datum[datum == "2|3"] <- 9
datum[datum == "1/1"] <- 2
datum[datum == "1|1"] <- 2
datum[datum == "2/2"] <- 9
datum[datum == "2|2"] <- 9
datum <- as.data.frame(apply(datum,2,as.numeric))
rownames(datum) <- chrids_datum
datum <- datum[complete.cases(datum),]
datum <- datum %>%
  mutate(across(everything(), ~ replace(.x, !(.x %in% 0:2), NA_real_)))
datum <- datum[complete.cases(datum),]
datum <- datum[apply(datum,1,sd) > 0,]
pca <- prcomp(t(datum), scale. = TRUE) ## do pca
## plot of pca with groups in ellipses
CSpca <- fviz_pca_ind(pca,
                      col.ind = c(rep("CS",3),rep("CSxP",3),rep("P",3),rep("CSxP",3)), # color by groups
                      palette = c("#0072B2", "#E69F00", "#CC79A7"),
                      legend.title = "Genotypes",
                      repel = TRUE,
                      pointshape = 16,                                  # filled circles
                      pointsize  = 3,
                      mean.point = FALSE, 
                      title = ""
) + guides(color = guide_legend(override.aes = list(shape = 16, size = 3)))

pdf("SNPs_PCA_first_pass.pdf",height=3.5,width=4.5)
CSpca
dev.off()

datum <- read.vcfR(file="wheat_ase_het_snps_filtered.vcf")
datum <- as.data.frame(extract.gt(datum))
colnames(datum) <- gsub("PXCS2","PxCS2", gsub(".sort.cs.bam","", gsub("_.*","", colnames(datum))))
chrids_datum <- rownames(datum)
datum[datum == "0/0"] <- 0
datum[datum == "0|0"] <- 0
datum[datum == "0/1"] <- 1
datum[datum == "0|1"] <- 1
datum[datum == "1|0"] <- 1
datum[datum == "0/2"] <- 9
datum[datum == "0|2"] <- 9
datum[datum == "0/3"] <- 9
datum[datum == "0|3"] <- 9
datum[datum == "1/2"] <- 9
datum[datum == "1|2"] <- 9
datum[datum == "1/3"] <- 9
datum[datum == "1|3"] <- 9
datum[datum == "2/3"] <- 9
datum[datum == "2|3"] <- 9
datum[datum == "1/1"] <- 2
datum[datum == "1|1"] <- 2
datum[datum == "2/2"] <- 9
datum[datum == "2|2"] <- 9
datum <- as.data.frame(apply(datum,2,as.numeric))
rownames(datum) <- chrids_datum
datum <- datum[complete.cases(datum),]
filtered_list <- datum[rowSums(datum[1:3]) < 1,]
filtered_list <- filtered_list[rowSums(filtered_list[7:9]) == 6,]
A <- intersect(c("CS1","CS2","CS3","P1","P2","P3"), names(datum))
B <- intersect(c("CSxP1","CSxP2","CSxP3","PxCS1","PxCS2"), names(datum))
filtered_list <- filtered_list %>%
  dplyr::filter(
    !if_any(all_of(A), ~ .x %in% c(1, 9)),     # no 1 or 9 in A
    !if_any(all_of(B), ~ .x == 9),             # no 9 in B
  )
filtered_list <- as.data.frame(rownames(filtered_list))
filtered_list$`rownames(filtered_list)` <- gsub("ChrUnknown","Chr_Unknown",filtered_list$`rownames(filtered_list)`)
filtered_list <- filtered_list %>%
  tidyr::extract(`rownames(filtered_list)`,
                 into = c("chr", "position"),
                 regex = '^([^_]+_[^_]+)_(.+)$',
                 remove = FALSE, convert = TRUE)
filtered_list$chr <- gsub("Chr_Unknown","ChrUnknown",filtered_list$chr)
write.table(filtered_list[2:3],"filtered_set_CS.txt",row.names = F, col.names = F, quote = F)

datum <- datum %>%
  mutate(across(everything(), ~ replace(.x, !(.x %in% 0:2), NA_real_)))
datum <- datum[complete.cases(datum),]
datum <- datum[apply(datum,1,sd) > 0,]
pca <- prcomp(t(datum), scale. = TRUE) ## do pca
## plot of pca with groups in ellipses
CSpca <- fviz_pca_ind(pca,
                      col.ind = c(rep("CS",3),rep("CSxP",3),rep("P",3),rep("CSxP",2)), # color by groups
                      palette = c("#0072B2", "#E69F00", "#CC79A7"),
                      legend.title = "Genotypes",
                      repel = TRUE,
                      pointshape = 16,                                  # filled circles
                      pointsize  = 3,
                      mean.point = FALSE, 
                      title = paste("n=",nrow(datum)," sites",sep="")
) + guides(color = guide_legend(override.aes = list(shape = 16, size = 3)))

pdf("SNPs_PCA.pdf",height=3.5,width=4.5)
CSpca
dev.off()

csscree <- fviz_screeplot(pca, ncp=10,title = "")
csscree

## now with PAR reference

datum <- read.vcfR(file="par_ase_het_snps_filtered.vcf")
datum <- as.data.frame(extract.gt(datum))
colnames(datum) <- gsub("PXCS2","PxCS2", gsub(".sort.cs.bam","", gsub("_.*","", colnames(datum))))
chrids_datum <- rownames(datum)

datum[datum == "0/0"] <- 0
datum[datum == "0|0"] <- 0
datum[datum == "0/1"] <- 1
datum[datum == "0|1"] <- 1
datum[datum == "1|0"] <- 1
datum[datum == "0/2"] <- 9
datum[datum == "0|2"] <- 9
datum[datum == "0/3"] <- 9
datum[datum == "0|3"] <- 9
datum[datum == "1/2"] <- 9
datum[datum == "1|2"] <- 9
datum[datum == "1/3"] <- 9
datum[datum == "1|3"] <- 9
datum[datum == "2/3"] <- 9
datum[datum == "2|3"] <- 9
datum[datum == "1/1"] <- 2
datum[datum == "1|1"] <- 2
datum[datum == "2/2"] <- 9
datum[datum == "2|2"] <- 9
datum <- as.data.frame(apply(datum,2,as.numeric))
rownames(datum) <- chrids_datum

datum <- datum[complete.cases(datum),]
filtered_list <- datum[rowSums(datum[1:3]) == 6,]
filtered_list <- filtered_list[rowSums(filtered_list[7:9]) < 1,]
A <- intersect(c("CS1","CS2","CS3","P1","P2","P3"), names(datum))
B <- intersect(c("CSxP1","CSxP2","CSxP3","PxCS1","PxCS2"), names(datum))
filtered_list <- filtered_list %>%
  dplyr::filter(
    !if_any(all_of(A), ~ .x %in% c(1, 9)),     # no 1 or 9 in A
    !if_any(all_of(B), ~ .x == 9),             # no 9 in B
  )

filtered_list <- as.data.frame(rownames(filtered_list))
filtered_list <- filtered_list %>%
  tidyr::extract(`rownames(filtered_list)`,
                 into = c("chr", "position"),
                 regex = '^([^_]+_[^_]+)_(.+)$',   # grab 1st_two_fields and the rest
                 remove = FALSE, convert = TRUE)

write.table(filtered_list[2:3],"filtered_set_PAR.txt",row.names = F, col.names = F, quote = F)

datum <- datum %>%
  mutate(across(everything(), ~ replace(.x, !(.x %in% 0:2), NA_real_)))
datum <- datum[complete.cases(datum),]
datum <- datum[apply(datum,1,sd) > 0,]
pca <- prcomp(t(datum), scale. = TRUE) ## do pca
## plot of pca with groups in ellipses
PARpca <- fviz_pca_ind(pca,
                      col.ind = c(rep("CS",3),rep("CSxP",3),rep("P",3),rep("CSxP",2)), # color by groups
                      palette = c("#0072B2", "#E69F00", "#CC79A7"),
                      legend.title = "Genotypes",
                      repel = TRUE,
                      pointshape = 16,                                  # filled circles
                      pointsize  = 3,
                      mean.point = FALSE, 
                      title = paste("n=",nrow(datum)," sites",sep="")
) + guides(color = guide_legend(override.aes = list(shape = 16, size = 3)))

pdf("SNPs_PCA_PAR.pdf",height=3.5,width=4.5)
PARpca
dev.off()

csscree <- fviz_screeplot(pca, ncp=10,title = "")
csscree

