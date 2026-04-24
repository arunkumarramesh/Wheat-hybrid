
1. Trim RNA-seq reads
```
for file in *_1.fq.gz; do java -jar /software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 20 $file ${file/_1.fq.gz/_2.fq.gz} ${file/_1.fq.gz/_1.paired.fq.gz} ${file/_1.fq.gz/_1.unpaired.fq.gz} ${file/_1.fq.gz/_2.paired.fq.gz} ${file/_1.fq.gz/_2.unpaired.fq.gz} ILLUMINACLIP:adaptors_novogene.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36; done
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
for file in *_trimmed.fq.gz; do /software/kallisto/build/src/kallisto quant -i iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna_index  -o ${file/.fq.gz/_CS} --single -l 200 -s 20 -t 20  $file ; done
for file in *_trimmed.fq.gz; do /software/kallisto/build/src/kallisto quant -i Triticum_aestivum_paragon.GCA949126075v1.cdna.all_index  -o ${file/.fq.gz/_PAR} --single -l 200 -s 20 -t 20  $file ; done

grep '>' iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna.fasta | cut -d ' ' -f 1 | sed 's/>//' >transnames
sed 's/\..*//g' transnames | paste -d ',' transnames - >transcript_to_gene_refseqv2.1.csv
ls -d *CS/ | sed 's/\///' >cs_kallisto_samplenames.txt

grep '>' Triticum_aestivum_paragon.GCA949126075v1.cdna.all.fa | cut -d ' ' -f 1 | sed 's/>//' >transnames2
sed 's/\.[0-9]\+$//' transnames2 | paste -d ',' transnames2 - >transcript_to_gene_paragon.GCA949126075v1.csv
ls -d *PAR/ | sed 's/\///' >par_kallisto_samplenames.txt

for file in *_1.paired.fq.gz; do /software/kallisto/build/src/kallisto quant -i iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna_index -o ${file/_1.paired.fq.gz/_CS} -t 20  $file ${file/_1.paired.fq.gz/_2.paired.fq.gz} ; done
for file in *_1.paired.fq.gz; do /software/kallisto/build/src/kallisto quant -i Triticum_aestivum_paragon.GCA949126075v1.cdna.all_index -o ${file/_1.paired.fq.gz/_PAR} -t 20  $file ${file/_1.paired.fq.gz/_2.paired.fq.gz} ; done
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

export PATH="/software/ncbi-blast-2.16.0+/bin:$PATH"

python3 /software/OrthoFinder-2.5.5/orthofinder.py -f orthofinder_input_2_A/ -t 20 -d
python3 /software/OrthoFinder-2.5.5/orthofinder.py -f orthofinder_input_2_B/ -t 20 -d
python3 /software/OrthoFinder-2.5.5/orthofinder.py -f orthofinder_input_2_D/ -t 20 -d

cd /projects/wheat/orthofinder_input_2_A/OrthoFinder/Results_Jul05_1/Orthogroups
awk 'NR==FNR{keep[$1]; next} ($1 in keep)' Orthogroups_SingleCopyOrthologues.txt Orthogroups.tsv > SingleCopyOrthologues_matrix_A.tsv

cd /projects/wheat/orthofinder_input_2_B/OrthoFinder/Results_Jul05_1/Orthogroups
awk 'NR==FNR{keep[$1]; next} ($1 in keep)' Orthogroups_SingleCopyOrthologues.txt Orthogroups.tsv > SingleCopyOrthologues_matrix_B.tsv

cd /projects/wheat/orthofinder_input_2_D/OrthoFinder/Results_Jul05_1/Orthogroups
awk 'NR==FNR{keep[$1]; next} ($1 in keep)' Orthogroups_SingleCopyOrthologues.txt Orthogroups.tsv > SingleCopyOrthologues_matrix_D.tsv

cp /projects/wheat/orthofinder_input_2_A/OrthoFinder/Results_Jul05_1/Orthogroups/SingleCopyOrthologues_matrix_A.tsv .
cp /projects/wheat/orthofinder_input_2_B/OrthoFinder/Results_Jul05_1/Orthogroups/SingleCopyOrthologues_matrix_B.tsv .
cp /projects/wheat/orthofinder_input_2_D/OrthoFinder/Results_Jul05_1/Orthogroups/SingleCopyOrthologues_matrix_D.tsv .
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
6. Obtain genome references
```
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/fasta/triticum_aestivum_paragon/dna_index/Triticum_aestivum_paragon.GCA949126075v1.dna.toplevel.fa.gz
wget https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Assemblies/v2.1/iwgsc_refseqv2.1_assembly.fa.zip
unzip iwgsc_refseqv2.1_assembly.fa.zip
```
7. Split genome references in half
```
TH=450000000

cut -f1,2 Triticum_aestivum_paragon.GCA949126075v1.dna.toplevel.fa.gz.fai | while read CTG LEN; do
  if [ "$LEN" -gt "$TH" ]; then
    MID=$(( LEN/2 ))
    samtools faidx Triticum_aestivum_paragon.GCA949126075v1.dna.toplevel.fa.gz ${CTG}:1-${MID} > ${CTG}_part1.fa
    samtools faidx Triticum_aestivum_paragon.GCA949126075v1.dna.toplevel.fa.gz ${CTG}:$((MID+1))-${LEN}  > ${CTG}_part2.fa
    sed -i "1s/>.*/>${CTG}_part1/" ${CTG}_part1.fa
    sed -i "1s/>.*/>${CTG}_part2/" ${CTG}_part2.fa
  fi
done

samtools faidx Triticum_aestivum_paragon.GCA949126075v1.dna.toplevel.fa.gz  $(awk '$1 ~ /^scaffold_/ { print $1 }' Triticum_aestivum_paragon.GCA949126075v1.dna.toplevel.fa.gz.fai)  > scaffolds_only.fa
cat 1*.fa 2*.fa 3*.fa 4*.fa 5*.fa 6*.fa 7*.fa scaffolds_only.fa > Paragon_part.fa
samtools faidx Paragon_part.fa

TH=450000000

cut -f1,2 iwgsc_refseqv2.1_assembly.fa.fai | while read CTG LEN; do
  if [ "$LEN" -gt "$TH" ]; then
    MID=$(( LEN/2 ))
    samtools faidx iwgsc_refseqv2.1_assembly.fa ${CTG}:1-${MID}   > ${CTG}_part1.fa
    samtools faidx iwgsc_refseqv2.1_assembly.fa ${CTG}:$((MID+1))-${LEN}  > ${CTG}_part2.fa
    sed -i "1s/>.*/>${CTG}_part1/" ${CTG}_part1.fa
    sed -i "1s/>.*/>${CTG}_part2/" ${CTG}_part2.fa
  fi
done

samtools faidx iwgsc_refseqv2.1_assembly.fa  $(awk '$1 ~ /^ChrUnknown/ { print $1 }' iwgsc_refseqv2.1_assembly.fa.fai)  > scaffolds_only.fa
cat Chr*fa scaffolds_only.fa > iwgsc_refseqv2.1_part.fa
samtools faidx iwgsc_refseqv2.1_part.fa

```

8. Map reads to genome reference
```

/software/hisat2-2.2.1/hisat2-build -p 20  iwgsc_refseqv2.1_part.fa iwgsc_refseqv2.1_part
/software/hisat2-2.2.1/hisat2-build -p 20  Paragon_part.fa Paragon_part

java -jar /software/picard.jar CreateSequenceDictionary -R  iwgsc_refseqv2.1_part.fa -O iwgsc_refseqv2.1_part.dict
java -jar /software/picard.jar CreateSequenceDictionary -R   Paragon_part.fa -O Paragon_part.dict

for file in *_1.paired.fq.gz; do /software/hisat2-2.2.1/hisat2  -p 20 -x iwgsc_refseqv2.1_part -1 $file -2 ${file/_1.paired.fq.gz/_2.paired.fq.gz} -S ${file/_1.paired.fq.gz/.cs.sam} ; done

for file in *.cs.sam; do samtools sort -n -@ 4 -O bam -o ${file/.cs.sam/.sortname.cs.bam} $file; done
for file in *.sortname.cs.bam; do samtools fixmate -m $file ${file/.sortname.cs.bam/.fixmate.cs.bam}; done
for file in *.fixmate.cs.bam; do samtools sort -@ 4 -O bam -o ${file/.fixmate.cs.bam/.sort.cs.bam} $file; done
for file in *.sort.cs.bam ; do java -jar /software/picard.jar  AddOrReplaceReadGroups -I $file -O ${file/.sort.cs.bam/.readgroup.cs.bam} -LB species -PL illumina -PU 1 -SM $file; done
for file in *.readgroup.cs.bam ; do java -jar /software/picard.jar  MarkDuplicates -I $file -O ${file/.readgroup.cs.bam/_marked.cs.bam} -M ${file/.readgroup.cs.bam/_metrics.cs.txt}; done
for file in *_marked.cs.bam ; do java -jar /software/picard.jar BuildBamIndex -I $file; done
for file in *_marked.cs.bam ; do /software/gatk-4.3.0.0/gatk SplitNCigarReads -R iwgsc_refseqv2.1_part.fa -OBI F -I $file -O  ${file/_marked.cs.bam/_split.cs.bam}  ; done
for file in *_split.cs.bam ; do java -jar /software/picard.jar BuildBamIndex -I $file; done
# interval list is  the chr names
for file in *_split.cs.bam ; do /software/gatk-4.3.0.0/gatk HaplotypeCaller -R iwgsc_refseqv2.1_part.fa -I $file -O ${file/_split.cs.bam/.g.vcf.gz} -ERC GVCF -L interval.list ; done

ls *.g.vcf.gz | tail -n +3 >samples
/software/gatk-4.3.0.0/gatk --java-options "-Xmx45g -Xms1g" GenomicsDBImport -V CS1_RNA_MKRN250026357-1A_22VTNMLT4_L3.g.vcf.gz -V CS2_RNA_MKRN250026358-1A_22VTNMLT4_L4.g.vcf.gz --genomicsdb-workspace-path genomicsdb --tmp-dir /projects/wheat/tmp -L interval.list
cat samples | while read line; do /software/gatk-4.3.0.0/gatk --java-options "-Xmx80g" GenomicsDBImport --genomicsdb-update-workspace-path genomicsdb --tmp-dir /projects/wheat/tmp -V $line; done

/software/gatk-4.3.0.0/gatk  --java-options "-Xmx45g" GenotypeGVCFs -R iwgsc_refseqv2.1_part.fa -V gendb://genomicsdb -G StandardAnnotation -O wheat.ase.output.vcf.gz -L interval.list
/software/gatk-4.3.0.0/gatk SelectVariants -V  wheat.ase.output.vcf.gz -select-type SNP -O  wheat.ase.snps.vcf.gz

for file in *_1.paired.fq.gz; do /software/hisat2-2.2.1/hisat2 -p 20 -x Paragon_part -1 $file -2 ${file/_1.paired.fq.gz/_2.paired.fq.gz} | samtools sort -@ 8 -o ${file/_1.paired.fq.gz/.par.sorted.bam} - ; done

for file in *.par.sorted.bam; do samtools sort -n -@ 8 -O bam -o ${file/.par.sorted.bam/.sortname.par.bam} $file; done
for file in *.sortname.par.bam; do samtools fixmate -m $file ${file/.sortname.par.bam/.fixmate.par.bam}; done
for file in *.fixmate.par.bam; do samtools sort -@ 8 -O bam -o ${file/.fixmate.par.bam/.sort.par.bam} $file; done
for file in *.sort.par.bam ; do java -jar /software/picard.jar  AddOrReplaceReadGroups -I $file -O ${file/.sort.par.bam/.readgroup.par.bam} -LB species -PL illumina -PU 1 -SM $file; done
for file in *.readgroup.par.bam ; do java -jar /software/picard.jar  MarkDuplicates -I $file -O ${file/.readgroup.par.bam/_marked.par.bam} -M ${file/.readgroup.par.bam/_metrics.par.txt}; done
for file in *_marked.par.bam ; do java -jar /software/picard.jar BuildBamIndex -I $file; done
for file in *_marked.par.bam ; do /software/gatk-4.3.0.0/gatk SplitNCigarReads -R Paragon_part.fa -OBI F -I $file -O  ${file/_marked.par.bam/_split.par.bam}  ; done
for file in *_split.par.bam ; do java -jar /software/picard.jar BuildBamIndex -I $file; done
# interval list is  the chr names
for file in *_split.par.bam ; do /software/gatk-4.3.0.0/gatk HaplotypeCaller -R Paragon_part.fa -I $file -O ${file/_split.par.bam/.par.g.vcf.gz} -ERC GVCF -L interval_par.list ; done
/software/gatk-4.3.0.0/gatk --java-options "-Xmx45g -Xms1g" GenomicsDBImport -V CS1_RNA_MKRN250026357-1A_22VTNMLT4_L3.par.g.vcf.gz  -V CS2_RNA_MKRN250026358-1A_22VTNMLT4_L4.par.g.vcf.gz  --genomicsdb-workspace-path genomicsdb_par --tmp-dir /projects/wheat/tmp -L interval_par.list
cat samples_par | while read line; do /software/gatk-4.3.0.0/gatk --java-options "-Xmx45g" GenomicsDBImport --genomicsdb-update-workspace-path genomicsdb_par --tmp-dir /projects/wheat/tmp -V $line; done
/software/gatk-4.3.0.0/gatk  --java-options "-Xmx45g" GenotypeGVCFs -R Paragon_part.fa -V gendb://genomicsdb_par -G StandardAnnotation -O par.ase.output.vcf.gz -L interval_par.list
/software/gatk-4.3.0.0/gatk SelectVariants -V  par.ase.output.vcf.gz -select-type SNP -O  par.ase.snps.vcf.gz

```
