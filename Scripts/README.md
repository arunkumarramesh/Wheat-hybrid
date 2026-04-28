
1. Trim RNA-seq reads
```
for file in *_1.fq.gz; do java -jar /software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 20 $file ${file/_1.fq.gz/_2.fq.gz} ${file/_1.fq.gz/_1.paired.fq.gz} ${file/_1.fq.gz/_1.unpaired.fq.gz} ${file/_1.fq.gz/_2.paired.fq.gz} ${file/_1.fq.gz/_2.unpaired.fq.gz} ILLUMINACLIP:adaptors_novogene.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36; done
```
2. Obtained Chinese Spring and Paragon Reference Transcriptomes

```
wget https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v2.1/iwgsc_refseqv2.1_gene_annotation_200916.zip
wget  https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/fasta/triticum_aestivum_paragon/cdna/Triticum_aestivum_paragon.GCA949126075v1.cdna.all.fa.gz

```
3. Map reads to transcriptome references
```
cat iwgsc_refseqv2.1_annotation_200916_LC_mrna.fasta iwgsc_refseqv2.1_annotation_200916_HC_mrna.fasta > ../iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna.fasta
/software/kallisto/build/src/kallisto index -i iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna_index -t 20 iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna.fasta
/software/kallisto/build/src/kallisto index -i Triticum_aestivum_paragon.GCA949126075v1.cdna.all_index -t 20 Triticum_aestivum_paragon.GCA949126075v1.cdna.all.fa

grep '>' iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna.fasta | cut -d ' ' -f 1 | sed 's/>//' >transnames
sed 's/\..*//g' transnames | paste -d ',' transnames - >transcript_to_gene_refseqv2.1.csv

grep '>' Triticum_aestivum_paragon.GCA949126075v1.cdna.all.fa | cut -d ' ' -f 1 | sed 's/>//' >transnames2
sed 's/\.[0-9]\+$//' transnames2 | paste -d ',' transnames2 - >transcript_to_gene_paragon.GCA949126075v1.csv


for file in *_1.paired.fq.gz; do /software/kallisto/build/src/kallisto quant -i iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna_index -o ${file/_1.paired.fq.gz/_CS} -t 20  $file ${file/_1.paired.fq.gz/_2.paired.fq.gz} ; done
for file in *_1.paired.fq.gz; do /software/kallisto/build/src/kallisto quant -i Triticum_aestivum_paragon.GCA949126075v1.cdna.all_index -o ${file/_1.paired.fq.gz/_PAR} -t 20  $file ${file/_1.paired.fq.gz/_2.paired.fq.gz} ; done

ls -d *CS/ | sed 's/\///' >cs_kallisto_samplenames.txt
ls -d *PAR/ | sed 's/\///' >par_kallisto_samplenames.txt
```
4. Obtain 1-1 orthologs for Chinese Spring and Paragon Reference Transcripts
```
mkdir orthofinder_input_2_A
mkdir orthofinder_input_2_B
mkdir orthofinder_input_2_D
cp iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC_mrna.fasta .

awk '/^>TraesCS[1-7]A/{p=1; print; next} /^>/{p=0} p' iwgsc_refseqv2.1_annotation_200916_HC_mrna.fasta >iwgsc_refseqv2.1_annotation_200916_HC_mrna_A.fasta
awk '/^>TraesCS[1-7]B/{p=1; print; next} /^>/{p=0} p' iwgsc_refseqv2.1_annotation_200916_HC_mrna.fasta >iwgsc_refseqv2.1_annotation_200916_HC_mrna_B.fasta
awk '/^>TraesCS[1-7]D/{p=1; print; next} /^>/{p=0} p' iwgsc_refseqv2.1_annotation_200916_HC_mrna.fasta >iwgsc_refseqv2.1_annotation_200916_HC_mrna_D.fasta

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

cd /projects/wheat/orthofinder_input_2_B/OrthoFinder/Results_Apr24/Orthogroups
awk 'NR==FNR{keep[$1]; next} ($1 in keep)' Orthogroups_SingleCopyOrthologues.txt Orthogroups.tsv > SingleCopyOrthologues_matrix_B.tsv

cd /projects/wheat/orthofinder_input_2_D/OrthoFinder/Results_Apr24/Orthogroups
awk 'NR==FNR{keep[$1]; next} ($1 in keep)' Orthogroups_SingleCopyOrthologues.txt Orthogroups.tsv > SingleCopyOrthologues_matrix_D.tsv

cp /projects/wheat/orthofinder_input_2_A/OrthoFinder/Results_Jul05_1/Orthogroups/SingleCopyOrthologues_matrix_A.tsv .
cp /projects/wheat/orthofinder_input_2_B/OrthoFinder/Results_Apr24/Orthogroups/SingleCopyOrthologues_matrix_B.tsv .
cp /projects/wheat/orthofinder_input_2_D/OrthoFinder/Results_Apr24/Orthogroups/SingleCopyOrthologues_matrix_D.tsv .
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
#  samplenames with trailing "_trimmed_CS"
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

#  samplenames with trailing "_trimmed_PAR"
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
ls *.par.g.vcf.gz | tail -n +3 > samples_par
cat samples_par | while read line; do /software/gatk-4.3.0.0/gatk --java-options "-Xmx45g" GenomicsDBImport --genomicsdb-update-workspace-path genomicsdb_par --tmp-dir /projects/wheat/tmp -V $line; done
/software/gatk-4.3.0.0/gatk  --java-options "-Xmx45g" GenotypeGVCFs -R Paragon_part.fa -V gendb://genomicsdb_par -G StandardAnnotation -O par.ase.output.vcf.gz -L interval_par.list
/software/gatk-4.3.0.0/gatk SelectVariants -V  par.ase.output.vcf.gz -select-type SNP -O  par.ase.snps.vcf.gz

```

9. Also split annotation files by part
```
awk '$3 == "exon" {print $1, $4, $5}'  iwgsc_refseqv2.1_annotation_200916_HC.gff3 > iwgsc_refseqv2.1_annotation_200916_HC_exon.bed

awk '$3 == "gene"' iwgsc_refseqv2.1_annotation_200916_HC.gff3 | cut -f 1,3-5,9 | sed -e 's/;.*//' -e 's/ID=//' > gene.gff3
cut -f 1-2 iwgsc_refseqv2.1_part.fa.fai > iwgsc_refseqv2.1_part_chr_sizes.txt
sed -i '/ChrUnknown/d' iwgsc_refseqv2.1_part_chr_sizes.txt
python3 split_bed.py
grep 'ChrUnknown' iwgsc_refseqv2.1_annotation_200916_HC_exon.bed | sed 's/ /\t/g'| cat iwgsc_refseqv2.1_annotation_200916_HC_exon_part.bed - > iwgsc_refseqv2.1_annotation_200916_HC_exon_unknown_part.bed
gffread iwgsc_refseqv2.1_annotation_200916_HC.gff3 -T -o iwgsc_refseqv2.1_annotation_200916_HC.gtf 
python3 split_gff.py
grep '^ChrUnknown' iwgsc_refseqv2.1_annotation_200916_HC.gtf | cat iwgsc_refseqv2.1_annotation_200916_HC_part.gtf - >iwgsc_refseqv2.1_annotation_200916_HC_unknown_part.gtf
awk '$3 == "transcript"' iwgsc_refseqv2.1_annotation_200916_HC_unknown_part.gtf > transcripts_part.gtf
cut -f 3 SingleCopyOrthologues_matrix.tsv | tail -n +2 | sed 's/\..*//'  | grep -F -f - transcripts_part.gtf | cut -f 1,4,5  | sort -u >one_one_orthologs.bed

grep '>' Triticum_aestivum_paragon.GCA949126075v1.cdna.all.fa | awk 'match($0,/primary_assembly:[^:]+:([^:]+:[0-9]+:[0-9]+):/,m){print m[1]}'  | sed 's/:/\t/g' > paragon_cdna.bed
awk '$3 == "exon" {print $1, $4, $5}'  Triticum_aestivum_paragon.GCA949126075v1.62.gff3 >  Triticum_aestivum_paragon.GCA949126075v1.62_exon.bed
grep 'scaffold' Triticum_aestivum_paragon.GCA949126075v1.62_exon.bed | sed 's/ /\t/g'| cat Triticum_aestivum_paragon.GCA949126075v1.62_exon_part.bed - > Triticum_aestivum_paragon.GCA949126075v1.62_scaf_exon_part.bed
gffread Triticum_aestivum_paragon.GCA949126075v1.62.gff3 -T -o Triticum_aestivum_paragon.GCA949126075v1.62.gtf 
cut -f 1-2 Paragon_part.fa.fai >Paragon_part_chr_sizes.txt
sed -i '/scaffold/d' Paragon_part_chr_sizes.txt
python3 split_gff_par.py
python3 split_bed_par.py
grep '^scaffold' Triticum_aestivum_paragon.GCA949126075v1.62.gtf | cat Triticum_aestivum_paragon.GCA949126075v1.62_part.gtf - >Triticum_aestivum_paragon.GCA949126075v1.62_scaf_part.gtf
awk '$3 == "exon" {print $1, $4, $5}'  Triticum_aestivum_paragon.GCA949126075v1.62.gff3 >  Triticum_aestivum_paragon.GCA949126075v1.62_exon.bed
awk '$3 == "transcript"' Triticum_aestivum_paragon.GCA949126075v1.62_scaf_part.gtf > transcripts_par_part.gtf
awk -F'\t' 'BEGIN{OFS="\t"} /^#/ || $3=="gene" {print}' Triticum_aestivum_paragon.GCA949126075v1.62.gff3 > genes_par.gff3
cut -f 2 SingleCopyOrthologues_matrix.tsv | tail -n +2 | sed 's/\.[^.]*$//' | grep -F -f - transcripts_par_part.gtf | cut -f 1,4,5 | sort -u  >one_one_orthologs_par.bed
awk -F'\t' 'BEGIN{OFS="\t"} /^#/ || $3=="gene" {print}' iwgsc_refseqv2.1_annotation_200916_HC.gff3 > genes_refseqv2_HC.gff3
```
split_bed.py
```
chrom_sizes = {}
with open("iwgsc_refseqv2.1_part_chr_sizes.txt") as f:
    for line in f:
        chrom, size = line.strip().split("\t")
        base_chrom = chrom.replace("_part1", "").replace("_part2", "")
        chrom_sizes.setdefault(base_chrom, [0, 0])
        if "part1" in chrom:
            chrom_sizes[base_chrom][0] = int(size)
        else:
            chrom_sizes[base_chrom][1] = int(size)

with open("iwgsc_refseqv2.1_annotation_200916_HC_exon.bed") as f, open("iwgsc_refseqv2.1_annotation_200916_HC_exon_part.bed", "w") as out:
    for line in f:
        chrom, start, end = line.strip().split()
        start = int(start)
        end = int(end)
        part1_size = chrom_sizes[chrom][0]

        if start < part1_size:
            out.write(f"{chrom}_part1\t{start}\t{end}\n")
        else:
            out.write(f"{chrom}_part2\t{start - part1_size}\t{end - part1_size}\n")

```

split_gff.py
```
sizes_file = "../iwgsc_refseqv2.1_part_chr_sizes.txt"
in_gtf = "iwgsc_refseqv2.1_annotation_200916_HC.gtf"
out_gtf = "iwgsc_refseqv2.1_annotation_200916_HC_part.gtf"

chrom_sizes = {}
with open(sizes_file) as f:
    for line in f:
        chrom, size = line.strip().split("\t")[:2]
        base = chrom.replace("_part1", "").replace("_part2", "")
        chrom_sizes.setdefault(base, [None, None])
        if chrom.endswith("_part1"):
            chrom_sizes[base][0] = int(size)
        elif chrom.endswith("_part2"):
            chrom_sizes[base][1] = int(size)

def with_split_attr(attr_str, tag):
    s = attr_str.rstrip()
    if not s.endswith(";"):
        s += ";"
    return s + f' split "{tag}";'

with open(in_gtf) as fin, open(out_gtf, "w") as fout:
    for line in fin:
        if line.startswith("#") or not line.strip():
            fout.write(line)
            continue
        cols = line.rstrip("\n").split("\t")
        chrom = cols[0]
        start = int(cols[3])
        end = int(cols[4])
        part1_size = chrom_sizes[chrom][0]

        if end <= part1_size:
            cols[0] = f"{chrom}_part1"
            fout.write("\t".join(cols) + "\n")
        elif start > part1_size:
            cols[0] = f"{chrom}_part2"
            cols[3] = str(start - part1_size)
            cols[4] = str(end - part1_size)
            fout.write("\t".join(cols) + "\n")
        else:
            left = cols.copy()
            left[0] = f"{chrom}_part1"
            left[4] = str(part1_size)
            left[8] = with_split_attr(left[8], "left")
            fout.write("\t".join(left) + "\n")

            right = cols.copy()
            right[0] = f"{chrom}_part2"
            right[3] = "1"
            right[4] = str(end - part1_size)
            right[8] = with_split_attr(right[8], "right")
            fout.write("\t".join(right) + "\n")

```
split_gff_par.py
```
sizes_file = "Paragon_part_chr_sizes.txt"
in_gtf = "Triticum_aestivum_paragon.GCA949126075v1.62.gtf"
out_gtf = "Triticum_aestivum_paragon.GCA949126075v1.62_part.gtf"

chrom_sizes = {}
with open(sizes_file) as f:
    for line in f:
        chrom, size = line.strip().split("\t")[:2]
        base = chrom.replace("_part1", "").replace("_part2", "")
        chrom_sizes.setdefault(base, [None, None])
        if chrom.endswith("_part1"):
            chrom_sizes[base][0] = int(size)
        elif chrom.endswith("_part2"):
            chrom_sizes[base][1] = int(size)

def with_split_attr(attr_str, tag):
    s = attr_str.rstrip()
    if not s.endswith(";"):
        s += ";"
    return s + f' split "{tag}";'

with open(in_gtf) as fin, open(out_gtf, "w") as fout:
    for line in fin:
        if line.startswith("#") or not line.strip():
            fout.write(line)
            continue
        cols = line.rstrip("\n").split("\t")
        chrom = cols[0]
        start = int(cols[3])
        end = int(cols[4])
        part1_size = chrom_sizes[chrom][0]

        if end <= part1_size:
            cols[0] = f"{chrom}_part1"
            fout.write("\t".join(cols) + "\n")
        elif start > part1_size:
            cols[0] = f"{chrom}_part2"
            cols[3] = str(start - part1_size)
            cols[4] = str(end - part1_size)
            fout.write("\t".join(cols) + "\n")
        else:
            left = cols.copy()
            left[0] = f"{chrom}_part1"
            left[4] = str(part1_size)
            left[8] = with_split_attr(left[8], "left")
            fout.write("\t".join(left) + "\n")

            right = cols.copy()
            right[0] = f"{chrom}_part2"
            right[3] = "1"
            right[4] = str(end - part1_size)
            right[8] = with_split_attr(right[8], "right")
            fout.write("\t".join(right) + "\n")

```
split_bed_par.py
```
chrom_sizes={}
with open("Paragon_part_chr_sizes.txt") as f:
    for line in f:
        chrom,size=line.strip().split("\t")
        base=chrom.replace("_part1","").replace("_part2","")
        chrom_sizes.setdefault(base,[0,0])
        if "part1" in chrom:
            chrom_sizes[base][0]=int(size)
        else:
            chrom_sizes[base][1]=int(size)

with open("Triticum_aestivum_paragon.GCA949126075v1.62_exon.bed") as f,open("Triticum_aestivum_paragon.GCA949126075v1.62_exon_part.bed","w") as out:
    for line in f:
        chrom,start,end=line.strip().split()
        start=int(start);end=int(end)
        part1_size=chrom_sizes[chrom][0]
        if start<part1_size:
            out.write(f"{chrom}_part1\t{start}\t{end}\n")
        else:
            out.write(f"{chrom}_part2\t{start-part1_size}\t{end-part1_size}\n")
```
10. Identify sites for ASE analyses
```
/proj/popgen/a.ramesh/software/bcftools-1.16/bcftools view  -R iwgsc_refseqv2.1_annotation_200916_HC_exon_unknown_part.bed -i 'QUAL>=20 && N_ALT>=1 && COUNT(GT!="mis" && FMT/DP>=20 && FMT/GQ>=20)>0' -Oz -o wheat_ase_het_snps_filtered.vcf.gz wheat.ase.snps.vcf.gz
/software/bcftools-1.16/bcftools view  -i 'QUAL>=10 && N_ALT>=1 && COUNT(GT!="mis" && FMT/DP>=10 )>0' -Oz -o wheat_het_snps_filtered.vcf.gz wheat.ase.snps.vcf.gz
/software/htslib-1.16/tabix -p vcf wheat_het_snps_filtered.vcf.gz
gunzip wheat_het_snps_filtered.vcf.gz

/proj/popgen/a.ramesh/software/bcftools-1.16/bcftools view  -R Triticum_aestivum_paragon.GCA949126075v1.62_scaf_exon_part.bed -i 'QUAL>=20 && N_ALT>=1 && COUNT(GT!="mis" && FMT/DP>=20 && FMT/GQ>=20)>0' -Oz -o par_ase_het_snps_filtered.vcf.gz par.ase.snps.vcf.gz
/software/bcftools-1.16/bcftools view  -i 'QUAL>=10 && N_ALT>=1 && COUNT(GT!="mis" && FMT/DP>=10 )>0' -Oz -o par_het_snps_filtered.vcf.gz par.ase.snps.vcf.gz
/software/htslib-1.16/tabix -p vcf par_het_snps_filtered.vcf.gz
gunzip par_het_snps_filtered.vcf.gz

```

11. Include a merged file with PxCS3 to show no discrepency when genotype calls are included. only at gene expression level
```
/proj/popgen/a.ramesh/software/gatk-4.3.0.0/gatk GenotypeGVCFs  -R iwgsc_refseqv2.1_part.fa  -V PxCS3_RNA_MKRN250026362-1A_22VTNMLT4_L3.g.vcf.gz  -L wheat_ase_het_snps_filtered.vcf.gz  --include-non-variant-sites true  -O PxCS3_on_wheat_sites.vcf.gz
/proj/popgen/a.ramesh/software/bcftools-1.16/bcftools merge -m all -Oz -o wheat_sites_merged_cs.vcf.gz wheat_ase_het_snps_filtered.vcf.gz PxCS3_on_wheat_sites.vcf.gz
gunzip wheat_sites_merged_cs.vcf.gz
```

12. Only keep sites that are not heterozygous in parents and are biallelic. Also do PCA for SNPs.
```
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
```

13. Map reads using WASP for ASE
```
/software/vcftools-vcftools-581c231/bin/vcftools --vcf wheat_ase_het_snps_filtered.vcf --positions filtered_set_CS.txt --recode --recode-INFO-all --out wheat_ase_snps_het
/software/htslib-1.16/bgzip wheat_ase_snps_het.recode.vcf
/software/htslib-1.16/tabix wheat_ase_snps_het.recode.vcf.gz

/software/vcftools-vcftools-581c231/bin/vcftools --vcf par_ase_het_snps_filtered.vcf --positions filtered_set_PAR.txt --recode --recode-INFO-all --out par_ase_snps_het
/software/htslib-1.16/bgzip par_ase_snps_het.recode.vcf
/software/htslib-1.16/tabix par_ase_snps_het.recode.vcf.gz

/software/STAR-2.7.10b/bin/Linux_x86_64_static/STAR --runThreadN 16 --runMode genomeGenerate  --genomeDir star_index --genomeFastaFiles iwgsc_refseqv2.1_part.fa  --sjdbGTFfile iwgsc_refseqv2.1_annotation_200916_HC_unknown_part.gtf --sjdbOverhang 100 --limitGenomeGenerateRAM 48889586954

for file in *_1.paired.fq.gz; do /software/STAR-2.7.10b/bin/Linux_x86_64_static/STAR --runThreadN 16 --genomeDir star_index --readFilesIn $file ${file/_1.paired.fq.gz/_2.paired.fq.gz} --readFilesCommand zcat  --varVCFfile wheat_het_snps_filtered.vcf  --waspOutputMode SAMtag  --outSAMtype BAM SortedByCoordinate  --outFilterMultimapNmax 1  --outSAMattrRGline ID:$file SM:$file PL:ILLUMINA LB:lib1 PU:unit1 --outSAMattributes NH HI AS nM NM MD jM jI rB MC vA vG vW  --outFileNamePrefix ${file/_1.paired.fq.gz/} ; done
for file in *Aligned.sortedByCoord.out.bam; do samtools view -@ 10 -h $file | awk 'BEGIN{OFS="\t"} /^@/{print;next} {vw=""; for(i=12;i<=NF;i++) if($i~/^vW:i:/){split($i,a,":"); vw=a[3]; break} if(vw=="" || vw==1) print}' | samtools sort -@ 8 -o ${file/Aligned.sortedByCoord.out.bam/.wasp.bam} ; done
for file in *.wasp.bam ; do java -jar /software/picard.jar  MarkDuplicates -I $file -O ${file/.wasp.bam/.ase.bam} -M ${file/.wasp.bam/_metrics.ase.txt}; done
for file in *.ase.bam ; do java -jar /software/picard.jar BuildBamIndex -I $file; done

/software/STAR-2.7.10b/bin/Linux_x86_64_static/STAR --runThreadN 16 --runMode genomeGenerate  --genomeDir star_index_par --genomeFastaFiles Paragon_part.fa  --sjdbGTFfile Triticum_aestivum_paragon.GCA949126075v1.62_scaf_part.gtf --sjdbOverhang 100 --limitGenomeGenerateRAM 48889586954

for file in *_1.paired.fq.gz; do /software/STAR-2.7.10b/bin/Linux_x86_64_static/STAR --runThreadN 16 --genomeDir star_index_par --readFilesIn $file ${file/_1.paired.fq.gz/_2.paired.fq.gz} --readFilesCommand zcat  --varVCFfile par_het_snps_filtered.vcf  --waspOutputMode SAMtag  --outSAMtype BAM SortedByCoordinate  --outFilterMultimapNmax 1  --outSAMattrRGline ID:$file SM:$file PL:ILLUMINA LB:lib1 PU:unit1 --outSAMattributes NH HI AS nM NM MD jM jI rB MC vA vG vW  --outFileNamePrefix ${file/_1.paired.fq.gz/} ; done
for file in *Aligned.sortedByCoord.out.bam; do samtools view -@ 10 -h $file | awk 'BEGIN{OFS="\t"} /^@/{print;next} {vw=""; for(i=12;i<=NF;i++) if($i~/^vW:i:/){split($i,a,":"); vw=a[3]; break} if(vw=="" || vw==1) print}' | samtools sort -@ 8 -o ${file/Aligned.sortedByCoord.out.bam/.wasp.par.bam} ; done

for file in *.wasp.par.bam ; do java -jar /software/picard.jar  MarkDuplicates -I $file -O ${file/.wasp.par.bam/.ase.par.bam} -M ${file/.wasp.par.bam/_metrics.ase.par.txt}; done
for file in *.ase.par.bam ; do java -jar /software/picard.jar BuildBamIndex -I $file; done
```

14. Profile ASE using GATK
```
/software/gatk-4.3.0.0/gatk ASEReadCounter -R iwgsc_refseqv2.1_part.fa -I CSxP1_MKRN250026363-1A_22VTNMLT4_L3.ase.bam -V wheat_ase_snps_het.recode.vcf.gz -O CSxP1.wasp.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10
/software/gatk-4.3.0.0/gatk ASEReadCounter -R iwgsc_refseqv2.1_part.fa -I CSxP2_MKRN250026364-1A_22VTNMLT4_L3.ase.bam -V wheat_ase_snps_het.recode.vcf.gz -O CSxP2.wasp.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10
/software/gatk-4.3.0.0/gatk ASEReadCounter -R iwgsc_refseqv2.1_part.fa -I CSxP3_RNA_MKRN250033262-1A_22VTNMLT4_L3.ase.bam -V wheat_ase_snps_het.recode.vcf.gz -O CSxP3.wasp.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10
/software/gatk-4.3.0.0/gatk ASEReadCounter -R iwgsc_refseqv2.1_part.fa -I PxCS1_MKRN250026360-1A_22VTNMLT4_L3.ase.bam -V wheat_ase_snps_het.recode.vcf.gz -O PxCS1.wasp.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10
/software/gatk-4.3.0.0/gatk ASEReadCounter -R iwgsc_refseqv2.1_part.fa -I PXCS2_RNA_MKRN250033261-1A_22VTNMLT4_L4.ase.bam -V wheat_ase_snps_het.recode.vcf.gz -O PxCS2.wasp.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10

/software/gatk-4.3.0.0/gatk ASEReadCounter -R iwgsc_refseqv2.1_part.fa -I CS1_RNA_MKRN250026357-1A_22VTNMLT4_L3.ase.bam -V wheat_ase_snps_het.recode.vcf.gz -O CS1.wasp.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10
/software/gatk-4.3.0.0/gatk ASEReadCounter -R iwgsc_refseqv2.1_part.fa -I CS2_RNA_MKRN250026358-1A_22VTNMLT4_L4.ase.bam -V wheat_ase_snps_het.recode.vcf.gz -O CS2.wasp.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10
/software/gatk-4.3.0.0/gatk ASEReadCounter -R iwgsc_refseqv2.1_part.fa -I CS3_RNA_MKRN250026359-1A_22VTNMLT4_L3.ase.bam -V wheat_ase_snps_het.recode.vcf.gz -O CS3.wasp.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10
/software/gatk-4.3.0.0/gatk ASEReadCounter -R iwgsc_refseqv2.1_part.fa -I P1_RNA_MKRN250026354-1A_22VTNMLT4_L3.ase.bam -V wheat_ase_snps_het.recode.vcf.gz -O P1.wasp.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10
/software/gatk-4.3.0.0/gatk ASEReadCounter -R iwgsc_refseqv2.1_part.fa -I P2_RNA_MKRN250026355-1A_22VTNMLT4_L3.ase.bam -V wheat_ase_snps_het.recode.vcf.gz -O P2.wasp.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10
/software/gatk-4.3.0.0/gatk ASEReadCounter -R iwgsc_refseqv2.1_part.fa -I P3_RNA_MKRN250026356-1A_22VTNMLT4_L4.ase.bam -V wheat_ase_snps_het.recode.vcf.gz -O P3.wasp.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10


/software/gatk-4.3.0.0/gatk ASEReadCounter -R Paragon_part.fa -I CSxP1_MKRN250026363-1A_22VTNMLT4_L3.ase.par.bam -V par_ase_snps_het.recode.vcf.gz -O CSxP1.par.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10
/software/gatk-4.3.0.0/gatk ASEReadCounter -R Paragon_part.fa -I CSxP2_MKRN250026364-1A_22VTNMLT4_L3.ase.par.bam -V par_ase_snps_het.recode.vcf.gz -O CSxP2.par.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10
/software/gatk-4.3.0.0/gatk ASEReadCounter -R Paragon_part.fa -I CSxP3_RNA_MKRN250033262-1A_22VTNMLT4_L3.ase.par.bam -V par_ase_snps_het.recode.vcf.gz -O CSxP3.par.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10
/software/gatk-4.3.0.0/gatk ASEReadCounter -R Paragon_part.fa -I PxCS1_MKRN250026360-1A_22VTNMLT4_L3.ase.par.bam -V par_ase_snps_het.recode.vcf.gz -O PxCS1.par.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10
/software/gatk-4.3.0.0/gatk ASEReadCounter -R Paragon_part.fa -I PXCS2_RNA_MKRN250033261-1A_22VTNMLT4_L4.ase.par.bam -V par_ase_snps_het.recode.vcf.gz -O PxCS2.par.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10

/software/gatk-4.3.0.0/gatk ASEReadCounter -R Paragon_part.fa -I CS1_RNA_MKRN250026357-1A_22VTNMLT4_L3.ase.par.bam -V par_ase_snps_het.recode.vcf.gz -O CS1.par.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10
/software/gatk-4.3.0.0/gatk ASEReadCounter -R Paragon_part.fa -I CS2_RNA_MKRN250026358-1A_22VTNMLT4_L4.ase.par.bam -V par_ase_snps_het.recode.vcf.gz -O CS2.par.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10
/software/gatk-4.3.0.0/gatk ASEReadCounter -R Paragon_part.fa -I CS3_RNA_MKRN250026359-1A_22VTNMLT4_L3.ase.par.bam -V par_ase_snps_het.recode.vcf.gz -O CS3.par.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10
/software/gatk-4.3.0.0/gatk ASEReadCounter -R Paragon_part.fa -I P1_RNA_MKRN250026354-1A_22VTNMLT4_L3.ase.par.bam -V par_ase_snps_het.recode.vcf.gz -O P1.par.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10
/software/gatk-4.3.0.0/gatk ASEReadCounter -R Paragon_part.fa -I P2_RNA_MKRN250026355-1A_22VTNMLT4_L3.ase.par.bam -V par_ase_snps_het.recode.vcf.gz -O P2.par.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10
/software/gatk-4.3.0.0/gatk ASEReadCounter -R Paragon_part.fa -I P3_RNA_MKRN250026356-1A_22VTNMLT4_L4.ase.par.bam -V par_ase_snps_het.recode.vcf.gz -O P3.par.ase.tsv  --min-mapping-quality 20 --min-base-quality 20   --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE --min-depth 10

```

15. Run differential expression tests using de_wheat.R

16. Compare homoeolog expression bias for differentially expressed genes using assign_homoeolog_expression_bias_categories.R

17. Compare allele specific expression using ase_test.R

18. Process Harper et al. data
```

/proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/prefetch --option-file harper_acc -O harper_data
for file in *.sra; do /proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/fastq-dump --gzip --split-3  $file; done

mv SRR3031950.fastq.gz CS.fastq.gz
mv SRR3031953.fastq.gz P.fastq.gz
mv SRR2983146.fastq.gz CSxP04.fastq.gz
mv SRR2983147.fastq.gz CSxP05.fastq.gz
mv SRR2983148.fastq.gz CSxP17.fastq.gz
mv SRR2983149.fastq.gz CSxP19.fastq.gz
mv SRR2983150.fastq.gz CSxP20.fastq.gz
mv SRR2983151.fastq.gz CSxP22.fastq.gz
mv SRR2983152.fastq.gz CSxP24.fastq.gz
mv SRR2983153.fastq.gz CSxP25.fastq.gz
mv SRR2983154.fastq.gz CSxP26.fastq.gz
mv SRR2983155.fastq.gz CSxP29.fastq.gz
mv SRR2983156.fastq.gz CSxP30.fastq.gz
mv SRR2983157.fastq.gz CSxP34.fastq.gz
mv SRR2983158.fastq.gz CSxP06.fastq.gz
mv SRR2983159.fastq.gz CSxP35.fastq.gz
mv SRR2983160.fastq.gz CSxP37.fastq.gz
mv SRR2983161.fastq.gz CSxP38.fastq.gz
mv SRR2983172.fastq.gz CSxP41.fastq.gz
mv SRR2983174.fastq.gz CSxP42.fastq.gz
mv SRR2983175.fastq.gz CSxP46.fastq.gz
mv SRR2983176.fastq.gz CSxP47.fastq.gz
mv SRR2983177.fastq.gz CSxP49.fastq.gz
mv SRR2983179.fastq.gz CSxP53.fastq.gz
mv SRR2983182.fastq.gz CSxP54.fastq.gz
mv SRR2983183.fastq.gz CSxP07.fastq.gz
mv SRR2983185.fastq.gz CSxP56.fastq.gz
mv SRR2983186.fastq.gz CSxP58.fastq.gz
mv SRR2983187.fastq.gz CSxP59.fastq.gz
mv SRR2983189.fastq.gz CSxP61.fastq.gz
mv SRR2983190.fastq.gz CSxP62.fastq.gz
mv SRR2983192.fastq.gz CSxP65.fastq.gz
mv SRR2983194.fastq.gz CSxP66.fastq.gz
mv SRR2983195.fastq.gz CSxP67.fastq.gz
mv SRR2983196.fastq.gz CSxP73.fastq.gz
mv SRR2983197.fastq.gz CSxP76.fastq.gz
mv SRR2983198.fastq.gz CSxP11.fastq.gz
mv SRR2983199.fastq.gz CSxP78.fastq.gz
mv SRR2983200.fastq.gz CSxP81.fastq.gz
mv SRR2983201.fastq.gz CSxP84.fastq.gz
mv SRR2983202.fastq.gz CSxP87.fastq.gz
mv SRR2983315.fastq.gz CSxP89.fastq.gz
mv SRR2983317.fastq.gz CSxP91.fastq.gz
mv SRR2983318.fastq.gz CSxP93.fastq.gz
mv SRR2983319.fastq.gz CSxP12.fastq.gz
mv SRR2983320.fastq.gz CSxP13.fastq.gz
mv SRR2983321.fastq.gz CSxP14.fastq.gz
mv SRR2983322.fastq.gz CSxP15.fastq.gz
mv SRR2983323.fastq.gz CSxP16.fastq.gz

for file in *.fastq.gz; do java -jar /proj/popgen/a.ramesh/software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 -threads 20 $file ${file/.fastq.gz/_trimmed.fq.gz} ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36; done

for file in *_trimmed.fq.gz; do /proj/popgen/a.ramesh/software/kallisto/build/src/kallisto quant -i iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna_index  -o ${file/.fq.gz/_CS} --single -l 200 -s 20 -t 20  $file ; done
for file in *_trimmed.fq.gz; do /proj/popgen/a.ramesh/software/kallisto/build/src/kallisto quant -i Triticum_aestivum_paragon.GCA949126075v1.cdna.all_index  -o ${file/.fq.gz/_PAR} --single -l 200 -s 20 -t 20  $file ; done

```
