
1. Trim RNA-seq reads
```
for file in *_1.fq.gz; do java -jar /proj/popgen/a.ramesh/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 20 $file ${file/_1.fq.gz/_2.fq.gz} ${file/_1.fq.gz/_1.paired.fq.gz} ${file/_1.fq.gz/_1.unpaired.fq.gz} ${file/_1.fq.gz/_2.paired.fq.gz} ${file/_1.fq.gz/_2.unpaired.fq.gz} ILLUMINACLIP:adaptors_novogene.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36; done
```
2. Obtained Chinese Spring and Paragon Reference Transcriptomes

```
wget https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v2.1/iwgsc_refseqv2.1_gene_annotation_200916.zip
wget  https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/fasta/triticum_aestivum_paragon/cdna/Triticum_aestivum_paragon.GCA949126075v1.cdna.all.fa.gz

cat iwgsc_refseqv2.1_annotation_200916_LC_mrna.fasta iwgsc_refseqv2.1_annotation_200916_HC_mrna.fasta > ../iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna.fasta
cat iwgsc_refseqv2.1_annotation_200916_LC_pep.fasta iwgsc_refseqv2.1_annotation_200916_HC_pep.fasta > ../iwgsc_refseqv2.1_annotation_200916_HC_LC_pep.fasta

```
3. Map reads to transcriptome references
```
for file in *_trimmed.fq.gz; do /proj/popgen/a.ramesh/software/kallisto/build/src/kallisto quant -i iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna_index  -o ${file/.fq.gz/_CS} --single -l 200 -s 20 -t 20  $file ; done
for file in *_trimmed.fq.gz; do /proj/popgen/a.ramesh/software/kallisto/build/src/kallisto quant -i Triticum_aestivum_paragon.GCA949126075v1.cdna.all_index  -o ${file/.fq.gz/_PAR} --single -l 200 -s 20 -t 20  $file ; done

grep '>' iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna.fasta | cut -d ' ' -f 1 | sed 's/>//' >transnames
sed 's/\..*//g' transnames | paste -d ',' transnames - >transcript_to_gene_refseqv2.1.csv
ls -d *CS/ | sed 's/\///' >cs_kallisto_samplenames.txt

grep '>' Triticum_aestivum_paragon.GCA949126075v1.cdna.all.fa | cut -d ' ' -f 1 | sed 's/>//' >transnames2
sed 's/\.[0-9]\+$//' transnames2 | paste -d ',' transnames2 - >transcript_to_gene_paragon.GCA949126075v1.csv
ls -d *PAR/ | sed 's/\///' >par_kallisto_samplenames.txt

for file in *_1.paired.fq.gz; do /proj/popgen/a.ramesh/software/kallisto/build/src/kallisto quant -i iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna_index -o ${file/_1.paired.fq.gz/_CS} -t 20  $file ${file/_1.paired.fq.gz/_2.paired.fq.gz} ; done
for file in *_1.paired.fq.gz; do /proj/popgen/a.ramesh/software/kallisto/build/src/kallisto quant -i Triticum_aestivum_paragon.GCA949126075v1.cdna.all_index -o ${file/_1.paired.fq.gz/_PAR} -t 20  $file ${file/_1.paired.fq.gz/_2.paired.fq.gz} ; done
```
4. Obtain 1-1 orthologs for Chinese Spring and Paragon Reference Transcripts
```
mkdir orthofinder_input_2_A
mkdir orthofinder_input_2_B
mkdir orthofinder_input_2_D
cp iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC_mrna.fasta .

awk '/^>TraesCS[1-7]A/{p=1; print; next} /^>/{p=0} p' iwgsc_refseqv2.1_annotation_200916_HC_mrna.fasta >iwgsc_refseqv2.1_annotation_200916_HC_mrna_A.fasta
awk '/^>TraesCS[1-7]A/{p=1; print; next} /^>/{p=0} p' iwgsc_refseqv2.1_annotation_200916_HC_mrna.fasta >iwgsc_refseqv2.1_annotation_200916_HC_mrna_B.fasta
awk '/^>TraesCS[1-7]A/{p=1; print; next} /^>/{p=0} p' iwgsc_refseqv2.1_annotation_200916_HC_mrna.fasta >iwgsc_refseqv2.1_annotation_200916_HC_mrna_D.fasta

awk '/^>.*:[1-7]A:/{p=1; print; next} /^>/{p=0} p' Triticum_aestivum_paragon.GCA949126075v1.cdna.all.fa > Triticum_aestivum_paragon.GCA949126075v1.cdna.all_A.fa
awk '/^>.*:[1-7]B:/{p=1; print; next} /^>/{p=0} p' Triticum_aestivum_paragon.GCA949126075v1.cdna.all.fa > Triticum_aestivum_paragon.GCA949126075v1.cdna.all_B.fa
awk '/^>.*:[1-7]D:/{p=1; print; next} /^>/{p=0} p' Triticum_aestivum_paragon.GCA949126075v1.cdna.all.fa > Triticum_aestivum_paragon.GCA949126075v1.cdna.all_D.fa

mv iwgsc_refseqv2.1_annotation_200916_HC_mrna_A.fasta orthofinder_input_2_A/
mv iwgsc_refseqv2.1_annotation_200916_HC_mrna_B.fasta orthofinder_input_2_B/
mv iwgsc_refseqv2.1_annotation_200916_HC_mrna_D.fasta orthofinder_input_2_D/
mv Triticum_aestivum_paragon.GCA949126075v1.cdna.all_A.fa orthofinder_input_2_A/
mv Triticum_aestivum_paragon.GCA949126075v1.cdna.all_B.fa orthofinder_input_2_B/
mv Triticum_aestivum_paragon.GCA949126075v1.cdna.all_D.fa orthofinder_input_2_D

export PATH="/proj/popgen/a.ramesh/software/ncbi-blast-2.16.0+/bin:$PATH"

python3 /proj/popgen/a.ramesh/software/OrthoFinder-2.5.5/orthofinder.py -f orthofinder_input_2_A/ -t 20 -d
python3 /proj/popgen/a.ramesh/software/OrthoFinder-2.5.5/orthofinder.py -f orthofinder_input_2_B/ -t 20 -d
python3 /proj/popgen/a.ramesh/software/OrthoFinder-2.5.5/orthofinder.py -f orthofinder_input_2_D/ -t 20 -d

cd /proj/popgen/a.ramesh/projects/wheat/orthofinder_input_2_A/OrthoFinder/Results_Jul05_1/Orthogroups
awk 'NR==FNR{keep[$1]; next} ($1 in keep)' Orthogroups_SingleCopyOrthologues.txt Orthogroups.tsv > SingleCopyOrthologues_matrix_A.tsv

cd /proj/popgen/a.ramesh/projects/wheat/orthofinder_input_2_B/OrthoFinder/Results_Jul05_1/Orthogroups
awk 'NR==FNR{keep[$1]; next} ($1 in keep)' Orthogroups_SingleCopyOrthologues.txt Orthogroups.tsv > SingleCopyOrthologues_matrix_B.tsv

cd /proj/popgen/a.ramesh/projects/wheat/orthofinder_input_2_D/OrthoFinder/Results_Jul05_1/Orthogroups
awk 'NR==FNR{keep[$1]; next} ($1 in keep)' Orthogroups_SingleCopyOrthologues.txt Orthogroups.tsv > SingleCopyOrthologues_matrix_D.tsv

cp /proj/popgen/a.ramesh/projects/wheat/orthofinder_input_2_A/OrthoFinder/Results_Jul05_1/Orthogroups/SingleCopyOrthologues_matrix_A.tsv .
cp /proj/popgen/a.ramesh/projects/wheat/orthofinder_input_2_B/OrthoFinder/Results_Jul05_1/Orthogroups/SingleCopyOrthologues_matrix_B.tsv .
cp /proj/popgen/a.ramesh/projects/wheat/orthofinder_input_2_D/OrthoFinder/Results_Jul05_1/Orthogroups/SingleCopyOrthologues_matrix_D.tsv .
cat SingleCopyOrthologues_matrix_A.tsv SingleCopyOrthologues_matrix_B.tsv SingleCopyOrthologues_matrix_D.tsv > SingleCopyOrthologues_matrix.tsv

```
5. Summerise counts per gene
```
# Aim is to run combine samples to gene expression level from transcript level
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

# Aim is to run combine samples to gene expression level from transcript level

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
```


