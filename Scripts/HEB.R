library(dplyr)
library(ggplot2)
library(tibble)
library(scales)
library(cowplot)
library(edgeR)
library(agricolae)

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
total_genes <- read.csv(file = "CSvP all genes.csv")
read.counts <- read.counts[rownames(read.counts) %in% total_genes$X,]
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
res2 <- sample_full_rows(n_rows = 6591 , p = 0.436, n_cols = 3, trials = 10000, seed = 123)
mean(res2$full_rows_count)

pdf(file="Number_triads_expected.pdf",height=2.5,width=4)
ggplot(data.frame(full_rows_count = res2$full_rows_count), aes(x = full_rows_count)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 1500, color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Expected number of triads\nwith three DE homoeologs", y = "Frequency") +
  theme_minimal(base_size = 14)
dev.off()

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

hebbias_cv_plot_d <- ggplot(data=bias_categories[bias_categories$group_id %in% buckets[[4]]$group_id,], aes(x = sample, y = CV, fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(CV ~ sample, data = bias_categories[bias_categories$group_id %in% buckets[[4]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 0.8, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories[bias_categories$group_id %in% buckets[[4]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(CV, na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  
  coord_cartesian(clip = "off") +
  ylab("HEB") +
  xlab("") +
  scale_fill_manual(values = c(CS = "#0072B2", CSxP = "#E69F00", P = "#CC79A7"), name = "Genotype")

hebbias_cv_plot_g <- ggplot(data=bias_categories[bias_categories$group_id %in% buckets2[[4]]$group_id,], aes(x = sample, y = CV, fill = genotype)) +
  geom_boxplot() +
  geom_text(data = HSD.test(aov(CV ~ sample, data = bias_categories[bias_categories$group_id %in% buckets2[[4]]$group_id,]), "sample", group = TRUE)$groups %>% rownames_to_column("sample") %>% as_tibble() %>% mutate(genotype = sub("^(CSxP|PxCS|CS|P).*", "\\1", sample),genotype = ifelse(genotype == "PxCS", "CSxP", genotype)), 
            aes(sample, y = 0.8, label = groups),hjust = -0.3) +
  geom_text(data = bias_categories[bias_categories$group_id %in% buckets2[[4]]$group_id, ] %>% group_by(sample, genotype) %>% summarise(med = median(CV, na.rm = TRUE), .groups = "drop"),
            aes(x = sample, y = med, label = number(med, accuracy = 0.01)),vjust = -0.4, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  coord_cartesian(clip = "off") +
  xlab("") +
  ylab("HEB") +
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

pdf(file="ABD_CS_PvCSxP.pdf",height=7.5,width=7)
plot_grid(hebbias_avecpm_plot_d,hebbias_avecpm_plot_g,hebbias_cv_plot_g,ncol=1, labels="AUTO")
dev.off()

pdf(file="ABD_HEB.pdf",height=2.5,width=7)
hebbias_cv_plot_d
dev.off()
