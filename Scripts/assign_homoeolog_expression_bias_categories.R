# Aim is to classify triads as balanced, dominant or suppressed like Ricardo did

library(tidyverse)
library(fields)

tpms <- read.table(file="cs_tpm.tsv",sep = "\t", header=T)
homologies <- read.csv(file="homoeologs_1_1_1_synt_and_non_synt.csv")
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

## drop  homologies lacking a unique match for all three
homologies <- homologies %>%
  mutate(row_id = row_number()) %>%
  # long join once
  select(row_id, A, B, D) %>%
  pivot_longer(c(A, B, D), names_to = "subgenome", values_to = "v11") %>%
  left_join(map_unique, by = "v11") %>%
  group_by(row_id) %>%
  filter(all(!is.na(v21))) %>%               # keep only rows where A,B,D each have a unique value
  ungroup() %>%
  select(row_id, subgenome, v21) %>%
  pivot_wider(names_from = subgenome, values_from = v21) %>%
  inner_join(homologies %>% mutate(row_id = row_number()), by = "row_id") %>%
  select(-A.y, -B.y, -D.y) %>%
  rename(A = A.x, B = B.x, D = D.x) %>%
  select(-row_id)

# now I want to only keep genes which are expressed >0.5 tpm in at least 3 samples

tpms <- tpms[rowSums(tpms > 0.5) >= 3, ]

# now only keep HC genes
HC_expr <- tpms[!grepl("LC",row.names(tpms)),]
dim(HC_expr)

hc_genes_to_use <- data.frame(gene=rownames(HC_expr))
nrow(hc_genes_to_use)

# which triads have at least 1 homoeolog expressed >0.5 tpm?
head(homologies)

# make into long format and add column whether each gene is "hc_genes_to_use"
long_homoeologs <- homologies %>%
  gather(homoeolog, gene, A:D) %>%
  mutate("gene_to_use" = gene %in% hc_genes_to_use$gene ) 
head(long_homoeologs)
dim(long_homoeologs)
nrow(long_homoeologs[long_homoeologs$gene_to_use == "TRUE",])

# does each triad have an expressed gene?
group_expr <- long_homoeologs %>%
  group_by(group_id) %>%
  summarise (n_true = sum(gene_to_use == "TRUE"),
             n_false = sum(gene_to_use == "FALSE"))
head(group_expr)

# select only groups which have at least 1 homoeolog expressed

expressed_groups <- group_expr[group_expr$n_true >0 ,]
head(expressed_groups)
dim(expressed_groups)

# now I want to combine the expression data with this list of expressed groups
head(expressed_groups)

head(long_homoeologs)
dim(long_homoeologs)

long_homoeologs_to_use <- long_homoeologs[long_homoeologs$group_id %in% expressed_groups$group_id,]
head(long_homoeologs_to_use)
dim(long_homoeologs_to_use)

# do a test to make a matrix for tpm expression for one sample:
head(tpms)

tpm_long_homoeologs_to_use <- merge(long_homoeologs_to_use, tpms, by.x="gene", by.y=0) # add the tpm values for the long_homoeologs_to_use
head(tpm_long_homoeologs_to_use)
dim(tpm_long_homoeologs_to_use)

# now select just 1 sample and calculate the relative expression of A, B, D for each triad:
test_sample <-
  tpm_long_homoeologs_to_use %>%
  select(group_id, homoeolog, P3_RNA_MKRN250026356.1A_22VTNMLT4_L4_CS) %>%
  spread(homoeolog,  P3_RNA_MKRN250026356.1A_22VTNMLT4_L4_CS)
head(test_sample)

# now calculate relative ABD
test_sample$total <- test_sample$A + test_sample$B + test_sample$D
head(test_sample)

test_sample$A_rel <- test_sample$A/test_sample$total
test_sample$B_rel <- test_sample$B/test_sample$total
test_sample$D_rel <- test_sample$D/test_sample$total

head(test_sample)
# only keep triads with a sum >0.5 tpm
test_sample <- test_sample[test_sample$total >0.5,]
dim(test_sample)
head(test_sample)


## how did Ricardo actually assign the homoeologs to bias categories?

### NB to get this to work I will need to make a "test_mat" for each sample then I can just loop through!

# make the "ideal" categories
centers<-t(matrix(c(0.33,0.33,0.33,1,0,0,0,1,0,0,0,1,0,0.5,0.5,0.5,0,0.5,0.5,0.5,0), nrow=3))
colnames(centers)<-c("A","B","D")
rownames(centers)<-c("Central","A.dominant","B.dominant","D.dominant","A.suppressed","B.suppressed","D.suppressed")
head(centers)

test_mat <- as.matrix(test_sample[,c("A_rel","B_rel","D_rel")])
head(test_mat)
is.matrix(test_mat)
rownames(test_mat)<- test_sample$group_id
colnames(test_mat) <- c("A","B","D")
head(test_mat)

expectation_distance <- rdist(test_mat,centers) # this calculate the euclidian distance for each triad to each category
head(expectation_distance)

colnames(expectation_distance)<-c("Central",
                                  "A.dominant",  "B.dominant",  "D.dominant",
                                  "A.suppressed","B.suppressed","D.suppressed")

head(expectation_distance)
rownames(expectation_distance)<-rownames(test_mat) # add back in the triad names
head(expectation_distance)
X <- as.matrix(expectation_distance)
storage.mode(X) <- "numeric"        # ensure numeric
mins <- max.col(-X, ties.method = "first")
head(mins)
clust_desc<-colnames(expectation_distance) 
head(clust_desc)
name_mins<-clust_desc[mins]
head(name_mins) # give the categories names

general_desc<-c("Central","Dominant",  "Dominant",  "Dominant",
                "Suppressed","Suppressed","Suppressed")

general_name_mins<-general_desc[mins] # add the general category names too

head(general_name_mins)

output_df <- cbind(test_mat,name_mins,general_name_mins) # add together this information about dominance for each triad
head(output_df)


### now do this for each sample:
head(tpm_long_homoeologs_to_use)

list_of_samples <- c(colnames(tpm_long_homoeologs_to_use[,10:ncol(tpm_long_homoeologs_to_use)]))
list_of_samples

# make output dataframe:
head(output_df)
output_df_all_samples <- data.frame(A= numeric(), B= numeric(), D= numeric(), name_mins=character(), general_name_mins=character(),
                                    group_id = numeric(), sample=character())     
output_df_all_samples

for(sample in list_of_samples) {
  print(sample)
  # now select just 1 sample and calculate the relative expression of A, B, D for each triad:
  test_sample <-
    tpm_long_homoeologs_to_use %>%
    select(group_id, homoeolog, sample) %>%
    spread(homoeolog,  sample)
  head(test_sample)
  
  # now calculate relative ABD
  test_sample$total <- test_sample$A + test_sample$B + test_sample$D
  head(test_sample)
  
  test_sample$A_rel <- test_sample$A/test_sample$total
  test_sample$B_rel <- test_sample$B/test_sample$total
  test_sample$D_rel <- test_sample$D/test_sample$total
  
  head(test_sample)
  dim(test_sample)
  # only keep triads with a sum >0.5 tpm
  test_sample <- test_sample[test_sample$total >0.5,] 
  
  
  ## how did Ricardo actually assign the homoeologs to bias categories?
  
  ### NB to get this to work I will need to make a "test_mat" for each sample then I can just loop through!
  
  # make the "ideal" categories
  centers<-t(matrix(c(0.33,0.33,0.33,1,0,0,0,1,0,0,0,1,0,0.5,0.5,0.5,0,0.5,0.5,0.5,0), nrow=3))
  colnames(centers)<-c("A","B","D")
  rownames(centers)<-c("Central","A.dominant","B.dominant","D.dominant","A.suppressed","B.suppressed","D.suppressed")
  #head(centers)
  
  test_mat <- as.matrix(test_sample[,c("A_rel","B_rel","D_rel")])
  #head(test_mat)
  rownames(test_mat)<- test_sample$group_id
  colnames(test_mat) <- c("A","B","D")
  #head(test_mat)
  
  expectation_distance <- rdist(test_mat,centers) # this calculate the euclidian distance for each triad to each category
  head(expectation_distance)
  
  colnames(expectation_distance)<-c("Central",
                                    "A.dominant",  "B.dominant",  "D.dominant",
                                    "A.suppressed","B.suppressed","D.suppressed")
  
  head(expectation_distance)
  rownames(expectation_distance)<-rownames(test_mat) # add back in the triad names
  head(expectation_distance)
  X <- as.matrix(expectation_distance)
  storage.mode(X) <- "numeric"        # ensure numeric
  mins <- max.col(-X, ties.method = "first")
  #head(mins)
  clust_desc<-colnames(expectation_distance) 
  #head(clust_desc)
  name_mins<-clust_desc[mins]
  #head(name_mins) # give the categories names
  
  general_desc<-c("Central","Dominant",  "Dominant",  "Dominant",
                  "Suppressed","Suppressed","Suppressed")
  
  general_name_mins<-general_desc[mins] # add the general category names too
  
  #head(general_name_mins)
  
  output_mat <- cbind(test_mat,name_mins,general_name_mins)# add together this information about dominance for each triad
  output_df <- as.data.frame(output_mat)
  head(output_df)
  output_df$group_id <- rownames(output_df) # make group_id into a column
  output_df$sample <- sample # add which sample this is as a column
  
  head(output_df)
  
  output_df_all_samples <- rbind(output_df_all_samples, output_df) # puts all samples into a big table
  
}

head(output_df_all_samples)
dim(output_df_all_samples)

output_df_all_samples <- output_df_all_samples %>%
  mutate(group_id = trimws(as.character(group_id))) %>%
  filter(
    !is.na(group_id),
    group_id != "",
    !grepl("NA", group_id, ignore.case = TRUE)   # remove partial NA matches
  ) %>%
  add_count(group_id, name = "gcount") %>%
  filter(gcount == 12) %>%
  select(-gcount)

write.csv(file="bias_category_all_samples.csv", output_df_all_samples, row.names = F)


### now do this for each sample AND keep original values:
head(tpm_long_homoeologs_to_use)

list_of_samples <- c(colnames(tpm_long_homoeologs_to_use[,10:ncol(tpm_long_homoeologs_to_use)]))
list_of_samples

# make output dataframe:
head(output_df)
output_df_all_samples <- data.frame(A_tpm= numeric(), B_tpm= numeric(), D_tpm= numeric(), A= numeric(), B= numeric(), D= numeric(), name_mins=character(), general_name_mins=character(),
                                    group_id = numeric(), sample=character())     
output_df_all_samples

for(sample in list_of_samples) {
  print(sample)
  # now select just 1 sample and calculate the relative expression of A, B, D for each triad:
  test_sample <-
    tpm_long_homoeologs_to_use %>%
    select(group_id, homoeolog, sample) %>%
    spread(homoeolog,  sample)
  head(test_sample)
  
  # now calculate relative ABD
  test_sample$total <- test_sample$A + test_sample$B + test_sample$D
  head(test_sample)
  
  test_sample$A_rel <- test_sample$A/test_sample$total
  test_sample$B_rel <- test_sample$B/test_sample$total
  test_sample$D_rel <- test_sample$D/test_sample$total
  
  head(test_sample)
  dim(test_sample)
  # only keep triads with a sum >0.5 tpm
  test_sample <- test_sample[test_sample$total >0.5,] 
  
  
  ## how did Ricardo actually assign the homoeologs to bias categories?
  
  ### NB to get this to work I will need to make a "test_mat" for each sample then I can just loop through!
  
  # make the "ideal" categories
  centers<-t(matrix(c(0.33,0.33,0.33,1,0,0,0,1,0,0,0,1,0,0.5,0.5,0.5,0,0.5,0.5,0.5,0), nrow=3))
  colnames(centers)<-c("A","B","D")
  rownames(centers)<-c("Central","A.dominant","B.dominant","D.dominant","A.suppressed","B.suppressed","D.suppressed")
  #head(centers)
  
  test_mat <- as.matrix(test_sample[,c("A_rel","B_rel","D_rel")])
  #head(test_mat)
  rownames(test_mat)<- test_sample$group_id
  colnames(test_mat) <- c("A","B","D")
  #head(test_mat)
  
  expectation_distance <- rdist(test_mat,centers) # this calculate the euclidian distance for each triad to each category
  head(expectation_distance)
  
  colnames(expectation_distance)<-c("Central",
                                    "A.dominant",  "B.dominant",  "D.dominant",
                                    "A.suppressed","B.suppressed","D.suppressed")
  
  head(expectation_distance)
  rownames(expectation_distance)<-rownames(test_mat) # add back in the triad names
  head(expectation_distance)
  X <- as.matrix(expectation_distance)
  storage.mode(X) <- "numeric"        # ensure numeric
  mins <- max.col(-X, ties.method = "first")
  #head(mins)
  clust_desc<-colnames(expectation_distance) 
  #head(clust_desc)
  name_mins<-clust_desc[mins]
  #head(name_mins) # give the categories names
  
  general_desc<-c("Central","Dominant",  "Dominant",  "Dominant",
                  "Suppressed","Suppressed","Suppressed")
  
  general_name_mins<-general_desc[mins] # add the general category names too
  
  #head(general_name_mins)
  
  output_mat <- cbind("A_tpm" =test_sample$A, "B_tpm" =test_sample$B, "D_tpm"=test_sample$D, test_mat,name_mins,general_name_mins)# add together this information about dominance for each triad
  output_df <- as.data.frame(output_mat)
  head(output_df)
  output_df$group_id <- rownames(output_df) # make group_id into a column
  output_df$sample <- sample # add which sample this is as a column
  
  head(output_df)
  
  output_df_all_samples <- rbind(output_df_all_samples, output_df) # puts all samples into a big table
  
}

head(output_df_all_samples)
dim(output_df_all_samples)

output_df_all_samples <- output_df_all_samples %>%
  mutate(group_id = trimws(as.character(group_id))) %>%
  filter(
    !is.na(group_id),
    group_id != "",
    !grepl("NA", group_id, ignore.case = TRUE)   # remove partial NA matches
  ) %>%
  add_count(group_id, name = "gcount") %>%
  filter(gcount == 12) %>%
  select(-gcount)

write.csv(file="bias_category_all_samples_inc_orig_expr.csv", output_df_all_samples, row.names = F)
write.csv(homologies,file="homologies.csv",row.names = F)


library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(scales)
library(cowplot)
library(purrr)
library(rlang)
library(edgeR)
library(agricolae)
library(ggtext)

homologies <- read.csv(file="homologies.csv")

bias_categories <- read.csv(file="bias_category_all_samples_inc_orig_expr.csv")
bias_categories$sample <- gsub("_.*","",bias_categories$sample)
bias_categories <- bias_categories[!bias_categories$sample %in% "PxCS3",]
bias_categories$sample <- gsub("PXCS2","PxCS2",bias_categories$sample)
bias_categories$CV <- apply(bias_categories[4:6],1,sd)/apply(bias_categories[4:6],1,mean)
bias_categories$triad_tpm <- bias_categories$A_tpm+bias_categories$B_tpm+bias_categories$D_tpm
bias_categories <- bias_categories %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample), genotype = ifelse(genotype == "PxCS", "CSxP", genotype)) %>% mutate(sample = factor(sample, levels = c("CS1","CS2","CS3","CSxP1","CSxP2","CSxP3","PxCS1","PxCS2","P1","P2","P3")))
bias_categories$group_id <- as.numeric(gsub("X","",bias_categories$group_id))

## get cpm estimates for ASE genes with CS ref

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
colnames(cpm_nolog) <- gsub("PXCS2","PxCS2",colnames(cpm_nolog))
colnames(cpm_nolog_relative) <- gsub("PXCS2","PxCS2",colnames(cpm_nolog_relative))

ids <- rownames(cpm_nolog)
homologies_kept <- homologies %>%
  mutate(across(c(A, B, D), as.character)) %>%
  filter(A %in% ids, B %in% ids, D %in% ids)

bias_categories_cpm <- bias_categories[bias_categories$group_id %in% homologies_kept$group_id,]

attach_cpm <- function(bias_categories_cpm, cpm_nolog, homologies_kept) {
  triads <- homologies_kept[, c("group_id","A","B","D")]
  names(triads) <- c("group_id","A_gene","B_gene","D_gene")
  bias2 <- merge(bias_categories_cpm, triads, by = "group_id", all.x = TRUE, sort = FALSE)
  
  rn <- rownames(cpm_nolog)
  cn <- colnames(cpm_nolog)
  ci <- match(as.character(bias2$sample), cn)
  
  ai <- match(as.character(bias2$A_gene), rn)
  bi <- match(as.character(bias2$B_gene), rn)
  di <- match(as.character(bias2$D_gene), rn)
  
  bias2$A_cpm <- cpm_nolog[cbind(ai, ci)]
  bias2$B_cpm <- cpm_nolog[cbind(bi, ci)]
  bias2$D_cpm <- cpm_nolog[cbind(di, ci)]
  
  bias2
}

bias_categories_cpm <- attach_cpm(bias_categories_cpm, cpm_nolog, homologies_kept)
bias_categories_cpm$triad_cpm <- bias_categories_cpm$A_cpm+bias_categories_cpm$B_cpm+bias_categories_cpm$D_cpm

hebbias_plot_a <- ggplot(data=bias_categories, aes(x = sample, y = CV, fill = genotype)) +
  geom_boxplot() +
  geom_text(data = bias_categories %>% group_by(sample, genotype) %>% summarise(med = median(CV, na.rm = TRUE), .groups = "drop"),
    aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  geom_text(data = HSD.test(aov(CV ~ sample, data = bias_categories), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 0.8, label = groups),hjust = -0.3) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  coord_cartesian(clip = "off") +
  labs(title=paste("n=",nrow(bias_categories[bias_categories$sample %in% "CS1",])," triads",sep="")) +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

pdf(file="hebbias_cv_plot.pdf",height=2.5,width=7)
hebbias_plot_a
dev.off()

hebbias_plot_b <- ggplot(data=bias_categories_cpm,aes(x = sample, y = log10(triad_cpm), fill = genotype)) +
  geom_boxplot() +
  geom_text(data = bias_categories_cpm %>% 
              group_by(sample, genotype) %>% summarise(med = median(log10(triad_cpm), na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  #geom_text(data = HSD.test(aov(triad_cpm ~ sample, data = bias_categories_cpm), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
  #          aes(sample, y = 2, label = groups),hjust = -0.3) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

biasdf <- as.data.frame.matrix( prop.table(table(bias_categories$sample,bias_categories$general_name_mins), margin = 1)) %>%
  tibble::rownames_to_column("sample") %>%
  tidyr::pivot_longer(-sample, names_to = "bias_class", values_to = "prop") %>%
  mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample), genotype = ifelse(genotype == "PxCS", "CSxP", genotype)) %>%
  mutate(sample = factor(sample, levels = c("CS1","CS2","CS3","CSxP1","CSxP2","CSxP3","PxCS1","PxCS2","P1","P2","P3")))

hebbias_plot_c <- ggplot(biasdf, aes(x = sample, y = prop, fill = genotype)) +
  geom_col(width = 0.85) +
  facet_wrap(~ bias_class, nrow = 1,scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype") +
  labs(x = "Sample", y = "Proportion of triads") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.placement = "outside",
        strip.background = element_blank())

biasdf_narrow <- as.data.frame.matrix( prop.table(table(bias_categories$sample,bias_categories$name_mins), margin = 1)) %>%
  tibble::rownames_to_column("sample") %>%
  tidyr::pivot_longer(-sample, names_to = "bias_class", values_to = "prop") %>%
  mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample), genotype = ifelse(genotype == "PxCS", "CSxP", genotype)) %>%
  mutate(sample = factor(sample, levels = c("CS1","CS2","CS3","CSxP1","CSxP2","CSxP3","PxCS1","PxCS2","P1","P2","P3")))  %>%
  mutate(bias_class = factor(bias_class, levels = c("Central","A.dominant","B.dominant","D.dominant","A.suppressed","B.suppressed","D.suppressed")))

hebbias_plot_d <- ggplot(biasdf_narrow, aes(x = sample, y = prop, fill = genotype)) +
  geom_col(width = 0.85) +
  facet_wrap(~ bias_class, ncol = 4,scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype") +
  labs(x = "Sample", y = "Proportion of triads") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.placement = "outside",
        strip.background = element_blank())

pdf(file="hebbias_plot.pdf",height=9,width=8)
plot_grid(hebbias_plot_b,hebbias_plot_c,hebbias_plot_d,ncol=1,labels="AUTO")
dev.off()

data <- read.csv(file="CSvP all genes.csv",row.names = 1)
data_sig <- data[data$adj.P.Val < 0.05 & data$logFC > 0.58,]
data_sig2 <- data[data$adj.P.Val < 0.05 & data$logFC < -0.58,]
ids <- rownames(data)
ids_sig <- rownames(data_sig)
ids_sig2 <- rownames(data_sig2)

homologies_kept <- homologies %>%
  mutate(across(c(A, B, D), as.character)) %>%
  filter(A %in% ids, B %in% ids, D %in% ids)

hits <- homologies_kept %>%
  mutate(hit_A = A %in% ids_sig,hit_B = B %in% ids_sig,hit_D = D %in% ids_sig) %>%
  mutate(from = paste0(ifelse(hit_A, "A", ""),ifelse(hit_B, "B", ""),ifelse(hit_D, "D", "")),from = ifelse(from == "", "none", from),n_hits = (hit_A + hit_B + hit_D))

hits2 <- homologies_kept %>%
  mutate(hit_A = A %in% ids_sig2,hit_B = B %in% ids_sig2,hit_D = D %in% ids_sig2) %>%
  mutate(from = paste0(ifelse(hit_A, "A", ""),ifelse(hit_B, "B", ""),ifelse(hit_D, "D", "")),from = ifelse(from == "", "none", from),n_hits = (hit_A + hit_B + hit_D))

hits$direction <- "Up"
hits2$direction <- "Down"
hits_both <- rbind(hits,hits2)
write.csv(hits_both,file="triads_CS_P.csv",row.names = F)

combo_levels <- c("A","B","D","AB","AD","BD","ABD","none")
combo_counts <- hits %>%
  dplyr::count(from, name = "n_rows") %>%
  mutate(from = factor(from, levels = combo_levels)) %>%
  arrange(from)
colnames(combo_counts) <- c("Subgenome","Number_DE_homoeologs")
combo_counts_mod_1 <- combo_counts %>%
  mutate(Subgenome = factor(Subgenome, levels = c("A","B","D","AB","AD","BD","ABD","none"))) %>%
  filter(Subgenome != "none")
combo_counts_mod_1$Direction <- "CS"

## pay attention to what direction means here. Quite important. Direction means change in individual homoeolog direction. So for a triad in the AB categories, say if A is overexpressed and B is underexpressed, it will show up in different categories as single changes. Not sure if this actually happens.

combo_levels <- c("A","B","D","AB","AD","BD","ABD","none")
combo_counts <- hits2 %>%
  dplyr::count(from, name = "n_rows") %>%
  mutate(from = factor(from, levels = combo_levels)) %>%
  arrange(from)
colnames(combo_counts) <- c("Subgenome","Number_DE_homoeologs")
combo_counts_mod_2 <- combo_counts %>%
  mutate(Subgenome = factor(Subgenome, levels = c("A","B","D","AB","AD","BD","ABD","none"))) %>%
  filter(Subgenome != "none")
combo_counts_mod_2$Direction <- "Paragon"
combo_counts_mod <- rbind(combo_counts_mod_1,combo_counts_mod_2)

total_n    <- sum(combo_counts$Number_DE_homoeologs, na.rm = TRUE)
non_none_n <- sum(combo_counts_mod$Number_DE_homoeologs, na.rm = TRUE)
pct_non    <- 100 * non_none_n / total_n

pdf("DE_subgenome_cs_p.pdf",height=3,width=5)
ggplot(combo_counts_mod,
       aes(x = Subgenome, y = Number_DE_homoeologs, fill = Direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = Number_DE_homoeologs),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  scale_fill_manual(values = c(CS = "#0072B2", Paragon =  "#CC79A7")) +
  labs(
    x = "Subgenome of DE homoeologs",
    y = "Number of DE homoeologs",
    title = sprintf("%.1f%% of %s triads have DE homoeologs", pct_non, total_n),
    fill = "Direction"
  ) +
  theme_minimal(base_size = 12)
dev.off()

buckets <- homologies_kept %>%
  mutate(n_hits = rowSums(across(c(A, B, D), ~ .x %in% ids_sig))) %>%
  group_split(n_hits)  # list of tibbles for 0,1,2,3 hits (if present)
buckets2 <- homologies_kept %>%
  mutate(n_hits = rowSums(across(c(A, B, D), ~ .x %in% ids_sig2))) %>%
  group_split(n_hits)  # list of tibbles for 0,1,2,3 hits (if present)

hebbias_cv_plot_a <- ggplot(data=bias_categories[bias_categories$group_id %in% buckets[[1]]$group_id,], aes(x = sample, y = CV, fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(CV ~ sample, data = bias_categories[bias_categories$group_id %in% buckets[[1]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 0.7, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories[bias_categories$group_id %in% buckets[[1]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(CV, na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_cv_plot_b <- ggplot(data=bias_categories[bias_categories$group_id %in% buckets[[2]]$group_id,], aes(x = sample, y = CV, fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(CV ~ sample, data = bias_categories[bias_categories$group_id %in% buckets[[2]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 1, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories[bias_categories$group_id %in% buckets[[2]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(CV, na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_cv_plot_c <- ggplot(data=bias_categories[bias_categories$group_id %in% buckets[[3]]$group_id,], aes(x = sample, y = CV, fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(CV ~ sample, data = bias_categories[bias_categories$group_id %in% buckets[[3]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 1.2, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories[bias_categories$group_id %in% buckets[[3]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(CV, na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_cv_plot_d <- ggplot(data=bias_categories[bias_categories$group_id %in% buckets[[4]]$group_id,], aes(x = sample, y = CV, fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(CV ~ sample, data = bias_categories[bias_categories$group_id %in% buckets[[4]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 0.8, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories[bias_categories$group_id %in% buckets[[4]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(CV, na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_cv_plot_e <- ggplot(data=bias_categories[bias_categories$group_id %in% buckets2[[2]]$group_id,], aes(x = sample, y = CV, fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(CV ~ sample, data = bias_categories[bias_categories$group_id %in% buckets2[[2]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 0.8, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories[bias_categories$group_id %in% buckets2[[2]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(CV, na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_cv_plot_f <- ggplot(data=bias_categories[bias_categories$group_id %in% buckets2[[3]]$group_id,], aes(x = sample, y = CV, fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(CV ~ sample, data = bias_categories[bias_categories$group_id %in% buckets2[[3]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 0.8, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories[bias_categories$group_id %in% buckets2[[3]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(CV, na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_cv_plot_g <- ggplot(data=bias_categories[bias_categories$group_id %in% buckets2[[4]]$group_id,], aes(x = sample, y = CV, fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(CV ~ sample, data = bias_categories[bias_categories$group_id %in% buckets2[[4]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 0.8, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories[bias_categories$group_id %in% buckets2[[4]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(CV, na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

pdf(file="hebbias_cv_cs_p.pdf",height=16,width=7)
plot_grid(hebbias_cv_plot_a,hebbias_cv_plot_b,hebbias_cv_plot_c,hebbias_cv_plot_d,hebbias_cv_plot_e,hebbias_cv_plot_f,hebbias_cv_plot_g,ncol=1,labels="AUTO")
dev.off()


hebbias_avecpm_plot_a <- ggplot(data=bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[1]]$group_id,], aes(x = sample, y = log10(triad_cpm), fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(triad_cpm ~ sample, data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[1]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 2, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[1]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(log10(triad_cpm), na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_avecpm_plot_b <- ggplot(data=bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[2]]$group_id,], aes(x = sample, y = log10(triad_cpm), fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(triad_cpm ~ sample, data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[2]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 2, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[2]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(log10(triad_cpm), na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_avecpm_plot_c <- ggplot(data=bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[3]]$group_id,], aes(x = sample, y = log10(triad_cpm), fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(triad_cpm ~ sample, data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[3]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 2, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[3]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(log10(triad_cpm), na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_avecpm_plot_d <- ggplot(data=bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[4]]$group_id,], aes(x = sample, y = log10(triad_cpm), fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(triad_cpm ~ sample, data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[4]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 2, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[4]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(log10(triad_cpm), na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_avecpm_plot_e <- ggplot(data=bias_categories_cpm[bias_categories_cpm$group_id %in% buckets2[[2]]$group_id,], aes(x = sample, y = log10(triad_cpm), fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(triad_cpm ~ sample, data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets2[[2]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 2, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets2[[2]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(log10(triad_cpm), na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_avecpm_plot_f <- ggplot(data=bias_categories_cpm[bias_categories_cpm$group_id %in% buckets2[[3]]$group_id,], aes(x = sample, y = log10(triad_cpm), fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(triad_cpm ~ sample, data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets2[[3]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 2, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets2[[3]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(log10(triad_cpm), na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_avecpm_plot_g <- ggplot(data=bias_categories_cpm[bias_categories_cpm$group_id %in% buckets2[[4]]$group_id,], aes(x = sample, y = log10(triad_cpm), fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(triad_cpm ~ sample, data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets2[[4]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 2, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets2[[4]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(log10(triad_cpm), na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

pdf(file="hebbias_avecpm_cs_p.pdf",height=14,width=7)
plot_grid(hebbias_avecpm_plot_a,hebbias_avecpm_plot_b,hebbias_avecpm_plot_c,hebbias_avecpm_plot_d,hebbias_avecpm_plot_e,hebbias_avecpm_plot_f,hebbias_avecpm_plot_g,ncol=1,labels="AUTO")
dev.off()


data <- read.csv(file="CS_PvCSxP all genes.csv",row.names = 1)
data_sig <- data[data$adj.P.Val < 0.05 & data$logFC > 0.58,]
data_sig2 <- data[data$adj.P.Val < 0.05 & data$logFC < -0.58,]
ids <- rownames(data)
ids_sig <- rownames(data_sig)
ids_sig2 <- rownames(data_sig2)

homologies_kept <- homologies %>%
  mutate(across(c(A, B, D), as.character)) %>%
  filter(A %in% ids, B %in% ids, D %in% ids)

hits <- homologies_kept %>%
  mutate(hit_A = A %in% ids_sig,hit_B = B %in% ids_sig,hit_D = D %in% ids_sig) %>%
  mutate(from = paste0(ifelse(hit_A, "A", ""),ifelse(hit_B, "B", ""),ifelse(hit_D, "D", "")),from = ifelse(from == "", "none", from),n_hits = (hit_A + hit_B + hit_D))

hits2 <- homologies_kept %>%
  mutate(hit_A = A %in% ids_sig2,hit_B = B %in% ids_sig2,hit_D = D %in% ids_sig2) %>%
  mutate(from = paste0(ifelse(hit_A, "A", ""),ifelse(hit_B, "B", ""),ifelse(hit_D, "D", "")),from = ifelse(from == "", "none", from),n_hits = (hit_A + hit_B + hit_D))

hits$direction <- "Up"
hits2$direction <- "Down"
hits_both <- rbind(hits,hits2)
write.csv(hits_both,file="triads_hybrids_parents.csv",row.names = F)

combo_levels <- c("A","B","D","AB","AD","BD","ABD","none")
combo_counts <- hits %>%
  dplyr::count(from, name = "n_rows") %>%
  mutate(from = factor(from, levels = combo_levels)) %>%
  arrange(from)
colnames(combo_counts) <- c("Subgenome","Number_DE_homoeologs")
combo_counts_mod_1 <- combo_counts %>%
  mutate(Subgenome = factor(Subgenome, levels = c("A","B","D","AB","AD","BD","ABD","none"))) %>%
  filter(Subgenome != "none")
combo_counts_mod_1$Direction <- "Hybrids"

combo_levels <- c("A","B","D","AB","AD","BD","ABD","none")
combo_counts <- hits2 %>%
  dplyr::count(from, name = "n_rows") %>%
  mutate(from = factor(from, levels = combo_levels)) %>%
  arrange(from)
colnames(combo_counts) <- c("Subgenome","Number_DE_homoeologs")
combo_counts_mod_2 <- combo_counts %>%
  mutate(Subgenome = factor(Subgenome, levels = c("A","B","D","AB","AD","BD","ABD","none"))) %>%
  filter(Subgenome != "none")
combo_counts_mod_2$Direction <- "Parents"
combo_counts_mod <- rbind(combo_counts_mod_1,combo_counts_mod_2)

total_n    <- sum(combo_counts$Number_DE_homoeologs, na.rm = TRUE)
non_none_n <- sum(combo_counts_mod$Number_DE_homoeologs, na.rm = TRUE)
pct_non    <- 100 * non_none_n / total_n

pdf("DE_subgenome_CS_PvCSxP.pdf",height=3,width=6)
ggplot(combo_counts_mod,aes(x = Subgenome, y = Number_DE_homoeologs, fill = Direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = Number_DE_homoeologs),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  scale_fill_manual(values = c(Hybrids = "#E69F00", Parents = "#0072B2")) +
  labs(
    x = "Subgenome of DE homoeologs",
    y = "Number of DE homoeologs",
    title = sprintf("%.1f%% of %s triads have DE homoeologs", pct_non, total_n),
    fill = "Direction"
  ) +
  theme_minimal(base_size = 12)
dev.off()


# Check if number of triads in DE homoeologs in hybrids in unusual
# n_rows: number of rows, p: proportion of cells to sample, trials: number of independent repeats

sample_full_rows <- function(n_rows, p = NULL, n_cols = 3, trials = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  N <- n_rows * n_cols
  k <- round(p * N)
  
  coords <- expand.grid(row = seq_len(n_rows), col = seq_len(n_cols))
  counts <- integer(trials)
  
  for (t in seq_len(trials)) {
    idx <- sample.int(N, size = k, replace = FALSE)
    picks <- coords[idx, , drop = FALSE]
    row_counts <- tabulate(picks$row, nbins = n_rows)
    counts[t] <- sum(row_counts == n_cols)
  }
  
  list(
    n_rows = n_rows,
    n_cols = n_cols,
    total_cells = N,
    sampled_cells = k,
    full_rows_count = counts,
    mean_full_rows = mean(counts)
  )
}

# 10000 trials to estimate the distribution
res2 <- sample_full_rows(n_rows = 7892, p = 0.42, n_cols = 3, trials = 10000, seed = 123)
mean(res2$full_rows_count)

pdf(file="Number_triads_expected.pdf",height=2.5,width=4)
ggplot(data.frame(full_rows_count = res2$full_rows_count), aes(x = full_rows_count)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 1790, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Expected number of triads\nwith three DE homoeologs",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14)
dev.off()


library(dplyr)
library(stringr)

# Euclidean distance in 3D
e3 <- function(a1,b1,d1, a2,b2,d2) sqrt((a1 - a2)^2 + (b1 - b2)^2 + (d1 - d2)^2)
is_cs_rep <- function(s) s %in% c("CS1","CS2","CS3")
is_p_rep  <- function(s) s %in% c("P1","P2","P3")
cs_centroid <- bias_categories %>%
  filter(is_cs_rep(sample)) %>%
  group_by(group_id) %>%
  summarise(A_cs = mean(A, na.rm = TRUE),B_cs = mean(B, na.rm = TRUE),D_cs = mean(D, na.rm = TRUE),.groups = "drop")
p_centroid <- bias_categories %>%
  filter(is_p_rep(sample)) %>%
  group_by(group_id) %>%
  summarise(A_p = mean(A, na.rm = TRUE),B_p = mean(B, na.rm = TRUE),D_p = mean(D, na.rm = TRUE),.groups = "drop")
bc2 <- bias_categories %>%
  left_join(cs_centroid, by = "group_id") %>%
  left_join(p_centroid,  by = "group_id")
bc2 <- bc2 %>%
  mutate(dist_CS_P = ifelse(!is.na(A_cs) & !is.na(A_p),e3(A_cs, B_cs, D_cs, A_p, B_p, D_p),NA_real_))
bc2 <- bc2 %>%
  mutate(dist_to_CS = ifelse(!is_cs_rep(sample) & !is_p_rep(sample) & !is.na(A_cs),e3(A, B, D, A_cs, B_cs, D_cs),NA_real_))
bc2 <- bc2 %>%
  mutate(dist_to_P = ifelse(!is_cs_rep(sample) & !is_p_rep(sample) & !is.na(A_p),e3(A, B, D, A_p, B_p, D_p),NA_real_))
bc2 <- left_join(bc2,hits[c(4,14)],by="group_id")

library(LSD)

pdf(file="bias_distance.pdf",height=6,width=6)
par(mfrow = c(2,2))

with(bc2 %>% filter(!is.na(n_hits), n_hits > 0),  heatscatter(dist_to_CS, dist_to_P,xlim=c(0,1), ylim=c(0,1),xlab="Bias difference from CS", ylab="Bias difference from Paragon",main="With DE homoeologs"))
abline(0,1)

with(bc2 %>% filter(!is.na(n_hits), n_hits < 1)  , heatscatter(dist_to_CS, dist_to_P,xlim=c(0,1), ylim=c(0,1),xlab="Bias difference from CS", ylab="Bias difference from Paragon",main="Without DE homoeologs"))
abline(0,1)

with(rbind(bc2 %>% filter(!is.na(n_hits), n_hits > 0) %>% select(dist_CS_P, value = dist_to_CS),bc2 %>% filter(!is.na(n_hits), n_hits > 0) %>% select(dist_CS_P, value = dist_to_P)), heatscatter(dist_CS_P, value,xlim = c(0,1), ylim = c(0,1),xlab = "Parental bias difference",ylab = "Bias difference from parents",main = "With DE homoeologs"))
abline(0,1)

with(rbind(bc2 %>% filter(!is.na(n_hits), n_hits < 1) %>% select(dist_CS_P, value = dist_to_CS),bc2 %>% filter(!is.na(n_hits), n_hits < 1) %>% select(dist_CS_P, value = dist_to_P)), heatscatter(dist_CS_P, value,xlim = c(0,1), ylim = c(0,1),xlab = "Parental bias difference",ylab = "Bias difference from parents",main = "Without DE homoeologs"))
abline(0,1)

dev.off()

## check if bias distance between parents differs between DE genes vs non DE genes, boxplot

pdf(file="bias_distance_parents_de.pdf",height=3.5,width=3)
ggplot(data = (bc2 %>% filter(!is.na(n_hits)) %>% mutate(n_hits_f = factor(n_hits))), aes(x = n_hits_f, y = dist_CS_P)) +
  geom_boxplot(fill = "#69b3a2", color = "black", outlier.shape = 16, outlier.size = 1.5) +
  stat_summary(fun = median, geom = "text",aes(label = sprintf("%.2f", after_stat(y))),vjust = -0.5, size = 3.3) +
  geom_text(data = HSD.test(aov(dist_CS_P ~ n_hits_f, data = (bc2 %>% filter(!is.na(n_hits)) %>% mutate(n_hits_f = factor(n_hits)))), "n_hits_f", group = TRUE)$groups %>% rownames_to_column("n_hits_f") %>% as_tibble(), aes(n_hits_f, y = 1.15, label = groups),vjust = -0.2, size = 4, fontface = "bold") +
  labs(x = "Number of DE homoeologs",y = "Parental bias difference") +
  theme_minimal(base_size = 13) +
  ylim(0,1.25)
dev.off()


buckets <- homologies_kept %>%
  mutate(n_hits = rowSums(across(c(A, B, D), ~ .x %in% ids_sig))) %>%
  group_split(n_hits)  # list of tibbles for 0,1,2,3 hits (if present)
buckets2 <- homologies_kept %>%
  mutate(n_hits = rowSums(across(c(A, B, D), ~ .x %in% ids_sig2))) %>%
  group_split(n_hits)  # list of tibbles for 0,1,2,3 hits (if present)

hebbias_cv_plot_a <- ggplot(data=bias_categories[bias_categories$group_id %in% buckets[[1]]$group_id,], aes(x = sample, y = CV, fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(CV ~ sample, data = bias_categories[bias_categories$group_id %in% buckets[[1]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 0.7, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories[bias_categories$group_id %in% buckets[[1]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(CV, na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_cv_plot_b <- ggplot(data=bias_categories[bias_categories$group_id %in% buckets[[2]]$group_id,], aes(x = sample, y = CV, fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(CV ~ sample, data = bias_categories[bias_categories$group_id %in% buckets[[2]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 1, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories[bias_categories$group_id %in% buckets[[2]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(CV, na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_cv_plot_c <- ggplot(data=bias_categories[bias_categories$group_id %in% buckets[[3]]$group_id,], aes(x = sample, y = CV, fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(CV ~ sample, data = bias_categories[bias_categories$group_id %in% buckets[[3]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 1.2, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories[bias_categories$group_id %in% buckets[[3]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(CV, na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_cv_plot_d <- ggplot(data=bias_categories[bias_categories$group_id %in% buckets[[4]]$group_id,], aes(x = sample, y = CV, fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(CV ~ sample, data = bias_categories[bias_categories$group_id %in% buckets[[4]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 0.8, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories[bias_categories$group_id %in% buckets[[4]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(CV, na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  xlab("") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_cv_plot_e <- ggplot(data=bias_categories[bias_categories$group_id %in% buckets2[[2]]$group_id,], aes(x = sample, y = CV, fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(CV ~ sample, data = bias_categories[bias_categories$group_id %in% buckets2[[2]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 0.8, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories[bias_categories$group_id %in% buckets2[[2]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(CV, na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_cv_plot_f <- ggplot(data=bias_categories[bias_categories$group_id %in% buckets2[[3]]$group_id,], aes(x = sample, y = CV, fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(CV ~ sample, data = bias_categories[bias_categories$group_id %in% buckets2[[3]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 0.8, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories[bias_categories$group_id %in% buckets2[[3]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(CV, na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_cv_plot_g <- ggplot(data=bias_categories[bias_categories$group_id %in% buckets2[[4]]$group_id,], aes(x = sample, y = CV, fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(CV ~ sample, data = bias_categories[bias_categories$group_id %in% buckets2[[4]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 0.8, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories[bias_categories$group_id %in% buckets2[[4]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(CV, na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  xlab("") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

pdf(file="hebbias_cv_CS_PvCSxP.pdf",height=14,width=7)
plot_grid(hebbias_cv_plot_a,hebbias_cv_plot_b,hebbias_cv_plot_c,hebbias_cv_plot_d,hebbias_cv_plot_e,hebbias_cv_plot_f,hebbias_cv_plot_g,ncol=1,labels="AUTO")
dev.off()

hebbias_avecpm_plot_a <- ggplot(data=bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[1]]$group_id,], aes(x = sample, y = log10(triad_cpm), fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(triad_cpm ~ sample, data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[1]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 2, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[1]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(log10(triad_cpm), na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_avecpm_plot_b <- ggplot(data=bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[2]]$group_id,], aes(x = sample, y = log10(triad_cpm), fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(triad_cpm ~ sample, data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[2]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 2, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[2]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(log10(triad_cpm), na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_avecpm_plot_c <- ggplot(data=bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[3]]$group_id,], aes(x = sample, y = log10(triad_cpm), fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(triad_cpm ~ sample, data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[3]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 2, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[3]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(log10(triad_cpm), na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_avecpm_plot_d <- ggplot(data=bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[4]]$group_id,], aes(x = sample, y = log10(triad_cpm), fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(triad_cpm ~ sample, data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[4]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 2.5, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets[[4]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(log10(triad_cpm), na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  coord_cartesian(clip = "off") +
  xlab("") +
  ylab("log10 (Triad CPM)") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_avecpm_plot_e <- ggplot(data=bias_categories_cpm[bias_categories_cpm$group_id %in% buckets2[[2]]$group_id,], aes(x = sample, y = log10(triad_cpm), fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(triad_cpm ~ sample, data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets2[[2]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 2, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets2[[2]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(log10(triad_cpm), na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_avecpm_plot_f <- ggplot(data=bias_categories_cpm[bias_categories_cpm$group_id %in% buckets2[[3]]$group_id,], aes(x = sample, y = log10(triad_cpm), fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(triad_cpm ~ sample, data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets2[[3]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 2, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets2[[3]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(log10(triad_cpm), na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_avecpm_plot_g <- ggplot(data=bias_categories_cpm[bias_categories_cpm$group_id %in% buckets2[[4]]$group_id,], aes(x = sample, y = log10(triad_cpm), fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(triad_cpm ~ sample, data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets2[[4]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 2.5, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories_cpm[bias_categories_cpm$group_id %in% buckets2[[4]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(log10(triad_cpm), na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # headroom for labels
  coord_cartesian(clip = "off") +
  xlab("") +
  ylab("log10 (Triad CPM)") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

pdf(file="hebbias_avecpm_CS_PvCSxP.pdf",height=14,width=7)
plot_grid(hebbias_avecpm_plot_a,hebbias_avecpm_plot_b,hebbias_avecpm_plot_c,hebbias_avecpm_plot_d,hebbias_avecpm_plot_e,hebbias_avecpm_plot_f,hebbias_avecpm_plot_g,ncol=1,labels="AUTO")
dev.off()

pdf(file="ABD_CS_PvCSxP.pdf",height=5,width=13)
plot_grid(hebbias_avecpm_plot_d,hebbias_cv_plot_d,hebbias_avecpm_plot_g,hebbias_cv_plot_g,ncol=2, labels="AUTO")
dev.off()

## cannot do similar anayses with ASE genes as too few are in triads

data <- read.csv(file="all_res_CS.csv")
ids <- data$gene
homologies_kept <- homologies %>%
  mutate(across(c(A, B, D), as.character)) %>%
  filter(A %in% ids, B %in% ids, D %in% ids)
dim(homologies_kept)

data <- read.csv(file="Ref_vs_Alt.csv",row.names = 1)
ids <- rownames(data)
homologies_kept <- homologies %>%
  mutate(across(c(A, B, D), as.character)) %>%
  filter(A %in% ids, B %in% ids, D %in% ids)
dim(homologies_kept)

## but get bias estimates for the ASE categories

classified <- read.csv(file="classified_all.csv")
ids <- classified$gene
homologies_kept <- homologies %>%
  mutate(across(c(A, B, D), as.character)) %>%
  filter(A %in% ids, B %in% ids, D %in% ids)
dim(homologies_kept)
cat_map <- setNames(classified$category, classified$gene)
homologies_kept <- homologies_kept %>%
  mutate(A_cat = cat_map[A],B_cat = cat_map[B],D_cat = cat_map[D])
homologies_kept <- homologies_kept[-c(5:9)]
colnames(homologies_kept)[1:3] <- c("A_id","B_id","C_id")

lev <- c("Ambiguous","Cis + trans","Cis × trans","Cis only","Compensatory","Conserved","Trans only")
ital <- function(x){
  x <- gsub("\\b[Cc]is\\b","<i>cis</i>",x,perl=TRUE)
  x <- gsub("\\b[Tt]rans\\b","<i>trans</i>",x,perl=TRUE)
  x
}

mk_df <- function(x, y, pair){
  as.data.frame(table(x, y)) %>%
    dplyr::rename(x_cat = x, y_cat = y, n = Freq) %>%
    mutate(pair = pair)
}
df_compare_subgenomes <- bind_rows(mk_df(homologies_kept$A_cat, homologies_kept$B_cat, "A vs B"),mk_df(homologies_kept$A_cat, homologies_kept$D_cat, "A vs D"),mk_df(homologies_kept$B_cat, homologies_kept$D_cat, "B vs D")) %>%
  mutate(x_cat = factor(x_cat, levels = lev),y_cat = factor(y_cat, levels = lev))

pdf(file="traid_ase_bias.pdf",height=4,width=7)
ggplot(df_compare_subgenomes, aes(x = y_cat, y = x_cat, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n), size = 3) +
  facet_wrap(~ pair, nrow = 1) +
  scale_fill_gradient(low = "white", high = "steelblue", name = "Count") +
  labs(x = "", y = "",
       title = "") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_markdown(angle = 90, hjust = 1),
    axis.text.y = element_markdown()
  ) +
  scale_x_discrete(labels = ital) +
  scale_y_discrete(labels = ital) +
  coord_equal()
dev.off()

bias_categories_ase_sub <- bias_categories[bias_categories$group_id %in% homologies_kept$group_id,]
bias_categories_ase_sub <- left_join(bias_categories_ase_sub,homologies_kept,by="group_id")
bias_categories_ase_sub <- bias_categories_ase_sub[!(bias_categories_ase_sub$A_cat %in% "Ambiguous"),]
bias_categories_ase_sub <- bias_categories_ase_sub[!(bias_categories_ase_sub$B_cat %in% "Ambiguous"),]
bias_categories_ase_sub <- bias_categories_ase_sub[!(bias_categories_ase_sub$D_cat %in% "Ambiguous"),]
bias_categories_ase_sub <- bias_categories_ase_sub %>%
  rowwise() %>%
  mutate(cat_canonical = paste(as.character(sort(factor(c(A_cat, B_cat, D_cat), levels = c("Cis only", "Cis + trans", "Cis × trans", "Trans only", "Compensatory", "Conserved")))),collapse = ", ")) %>%
  ungroup()

label_cis_trans_md <- function(x) {
  str_replace_all(x,regex("\\b(cis|trans)\\b", ignore_case = TRUE),~ paste0("<i>", .x, "</i>"))
}


cat_counts <- bias_categories_ase_sub %>%
  dplyr::count(cat_canonical, name = "n")

cv_labels <- HSD.test(aov(CV ~ cat_canonical, data = bias_categories_ase_sub),"cat_canonical",group = TRUE)$groups %>%
  rownames_to_column("cat_canonical") %>%
  left_join(cat_counts, by = "cat_canonical") %>%
  mutate(label = paste0(groups, "\n", "n=", n))

homologies_McManus_plot_a <- ggplot(data = bias_categories_ase_sub,aes(x = cat_canonical, y = CV)) +
  geom_boxplot() +
  xlab("") +
  ylab("CV") +
  geom_text(data = cv_labels,aes(cat_canonical, y = 1.3, label = label),vjust = 0.5) +
  coord_flip() +
  scale_x_discrete(labels = label_cis_trans_md) +
  theme(axis.text.y = element_markdown(),strip.placement = "outside",strip.background = element_blank())

tpm_labels <- HSD.test(aov(triad_tpm ~ cat_canonical, data = bias_categories_ase_sub),"cat_canonical",group = TRUE)$groups %>%
  rownames_to_column("cat_canonical") %>%
  left_join(cat_counts, by = "cat_canonical") %>%
  mutate(label = paste0(groups, "\n", "n=", n))

homologies_McManus_plot_b <- ggplot(data = bias_categories_ase_sub,aes(x = cat_canonical, y = triad_tpm)) +
  geom_boxplot() +
  xlab("") +
  ylab("Triad TPM") +
  ylim(0, 180) +
  geom_text(data = tpm_labels,aes(cat_canonical, y = 120, label = label), vjust = 0.5) +
  coord_flip() +
  scale_x_discrete(labels = label_cis_trans_md) +
  theme(axis.text.y = element_markdown(),strip.placement = "outside",strip.background = element_blank())

pdf(file="triad_ase_cv.pdf",height=2,width=5)
homologies_McManus_plot_a
dev.off()

pdf(file="triad_ase_tpm.pdf",height=2,width=5)
homologies_McManus_plot_b
dev.off()