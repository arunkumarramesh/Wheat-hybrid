# Aim is to run combine samples to gene expression level from transcript level

# Philippa Borrill
#06-03-2018 # updated 14-03-2019, updated 22-02-2021

#Steps will be:

#1: Summarise counts per gene (rather than transcript) using tximport for the studies which are to be included in the manuscript

#2: Summarise tpm per gene (rather than transcript) using tximport for the studies which are to be included in the manuscript

##### #1: Summarise counts per gene ########## 

#BiocManager::install("tximportData")
library(tximportData)
library(readr)
library(tximport)
library(rhdf5)

# read in pre-constructed tx2gene table (transcript to gene table)
tx2gene <- read.csv("transcript_to_gene_refseqv2.1.csv", header=T)
head(tx2gene)

# make vector pointing to the kallisto results files   ########
samples <- read.table("cs_kallisto_samplenames.txt", header=F)
samples

files <- file.path(samples$V1, "abundance.tsv", fsep ="/")
files
names(files) <- paste0(samples$V1)
head(files)
all(file.exists(files))

# read in the files and sum per gene
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
names(txi)

head(txi$abundance)

# move into directory where I will save this analysis

# to see counts summarised per gene
head(txi$counts)
colnames(txi$counts)

# save counts summarised per gene
write.table(txi$counts, file="cs_count.tsv",sep = "\t")

# save tpm summarised per gene
write.table(txi$abundance, file="cs_tpm.tsv",sep = "\t")

# calculate average gene length across all samples
gene_lengths <- as.data.frame(rowMeans(txi$length))
head(gene_lengths)
colnames(gene_lengths) <- c("length")
head(gene_lengths)
#save length per gene
write.csv(gene_lengths, file="cs_gene_lengths.csv")


# read in pre-constructed tx2gene table (transcript to gene table)
tx2gene <- read.csv("transcript_to_gene_paragon.GCA949126075v1.csv", header=T)
head(tx2gene)

# make vector pointing to the kallisto results files   ########
samples <- read.table("par_kallisto_samplenames.txt", header=F)
samples

files <- file.path(samples$V1, "abundance.tsv", fsep ="/")
files
names(files) <- paste0(samples$V1)
head(files)
all(file.exists(files))

# read in the files and sum per gene
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
names(txi)

head(txi$abundance)

# move into directory where I will save this analysis

# to see counts summarised per gene
head(txi$counts)
colnames(txi$counts)

# save counts summarised per gene
write.table(txi$counts, file="par_count.tsv",sep = "\t")

# save tpm summarised per gene
write.table(txi$abundance, file="par_tpm.tsv",sep = "\t")

# calculate average gene length across all samples
gene_lengths <- as.data.frame(rowMeans(txi$length))
head(gene_lengths)
colnames(gene_lengths) <- c("length")
head(gene_lengths)
#save length per gene
write.csv(gene_lengths, file="par_gene_lengths.csv")


## Repeat for Harper et al data

# Aim is to run combine samples to gene expression level from transcript level

# Philippa Borrill
#06-03-2018 # updated 14-03-2019, updated 22-02-2021

#Steps will be:

#1: Summarise counts per gene (rather than transcript) using tximport for the studies which are to be included in the manuscript

#2: Summarise tpm per gene (rather than transcript) using tximport for the studies which are to be included in the manuscript

##### #1: Summarise counts per gene ########## 

#BiocManager::install("tximportData")
library(tximportData)
library(readr)
library(tximport)
library(rhdf5)

# read in pre-constructed tx2gene table (transcript to gene table)
tx2gene <- read.csv("transcript_to_gene_refseqv2.1.csv", header=T)
head(tx2gene)

# make vector pointing to the kallisto results files   ########
samples <- read.table("cs_kallisto_samplenames.txt", header=F)
samples

samples2 <- read.table("samples_trimmed_CS", header=F)
samples2

samples <- rbind(samples,samples2)

files <- file.path(samples$V1, "abundance.tsv", fsep ="/")
files
names(files) <- paste0(samples$V1)
head(files)
all(file.exists(files))

# read in the files and sum per gene
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
names(txi)

head(txi$abundance)

# move into directory where I will save this analysis

# to see counts summarised per gene
head(txi$counts)
colnames(txi$counts)

# save counts summarised per gene
write.table(txi$counts, file="cs_csxp_count.tsv",sep = "\t")

# save tpm summarised per gene
write.table(txi$abundance, file="cs_csxp_tpm.tsv",sep = "\t")

# calculate average gene length across all samples
gene_lengths <- as.data.frame(rowMeans(txi$length))
head(gene_lengths)
colnames(gene_lengths) <- c("length")
head(gene_lengths)
#save length per gene
write.csv(gene_lengths, file="cs_csxp_gene_lengths.csv")


# read in pre-constructed tx2gene table (transcript to gene table)
tx2gene <- read.csv("transcript_to_gene_paragon.GCA949126075v1.csv", header=T)
head(tx2gene)

# make vector pointing to the kallisto results files   ########
samples <- read.table("par_kallisto_samplenames.txt", header=F)
samples

samples2 <- read.table("samples_trimmed_PAR", header=F)
samples2

samples <- rbind(samples,samples2)

files <- file.path(samples$V1, "abundance.tsv", fsep ="/")
files
names(files) <- paste0(samples$V1)
head(files)
all(file.exists(files))

# read in the files and sum per gene
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
names(txi)

head(txi$abundance)

# move into directory where I will save this analysis

# to see counts summarised per gene
head(txi$counts)
colnames(txi$counts)

# save counts summarised per gene
write.table(txi$counts, file="par_csxp_count.tsv",sep = "\t")

# save tpm summarised per gene
write.table(txi$abundance, file="par_csxp_tpm.tsv",sep = "\t")

# calculate average gene length across all samples
gene_lengths <- as.data.frame(rowMeans(txi$length))
head(gene_lengths)
colnames(gene_lengths) <- c("length")
head(gene_lengths)
#save length per gene
write.csv(gene_lengths, file="par_csxp_gene_lengths.csv")

