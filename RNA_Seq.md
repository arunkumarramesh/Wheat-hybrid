
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
#/proj/popgen/a.ramesh/software/kallisto/build/src/kallisto index -i iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna_index -t 20 iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna.fasta
#/proj/popgen/a.ramesh/software/kallisto/build/src/kallisto index -i Triticum_aestivum_paragon.GCA949126075v1.cdna.all_index -t 20 Triticum_aestivum_paragon.GCA949126075v1.cdna.all.fa

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
ls *.par.g.vcf.gz | tail -n +3 > samples_par
cat samples_par | while read line; do /software/gatk-4.3.0.0/gatk --java-options "-Xmx45g" GenomicsDBImport --genomicsdb-update-workspace-path genomicsdb_par --tmp-dir /projects/wheat/tmp -V $line; done
/software/gatk-4.3.0.0/gatk  --java-options "-Xmx45g" GenotypeGVCFs -R Paragon_part.fa -V gendb://genomicsdb_par -G StandardAnnotation -O par.ase.output.vcf.gz -L interval_par.list
/software/gatk-4.3.0.0/gatk SelectVariants -V  par.ase.output.vcf.gz -select-type SNP -O  par.ase.snps.vcf.gz

```

9. Also split annotation files by part
```
awk '$3 == "exon" {print $1, $4, $5}'  iwgsc_refseqv2.1_annotation_200916_HC.gff3 > iwgsc_refseqv2.1_annotation_200916_HC_exon.bed
awk '$3 == "transcript"' iwgsc_refseqv2.1_annotation_200916_HC_unknown_part.gtf > transcripts_part.gtf
cut -f 3 SingleCopyOrthologues_matrix.tsv | tail -n +2 | sed 's/\..*//'  | grep -F -f - transcripts_part.gtf | cut -f 1,4,5  | sort -u >one_one_orthologs.bed
awk '$3 == "gene"' iwgsc_refseqv2.1_annotation_200916_HC.gff3 | cut -f 1,3-5,9 | sed -e 's/;.*//' -e 's/ID=//' > gene.gff3
cut -f 1-2 iwgsc_refseqv2.1_part.fa.fai > iwgsc_refseqv2.1_part_chr_sizes.txt
sed -i '/ChrUnknown/d' iwgsc_refseqv2.1_part_chr_sizes.txt
python3 split_bed.py
grep 'ChrUnknown' iwgsc_refseqv2.1_annotation_200916_HC_exon.bed | sed 's/ /\t/g'| cat iwgsc_refseqv2.1_annotation_200916_HC_exon_part.bed - > iwgsc_refseqv2.1_annotation_200916_HC_exon_unknown_part.bed
gffread iwgsc_refseqv2.1_annotation_200916_HC.gff3 -T -o iwgsc_refseqv2.1_annotation_200916_HC.gtf 
python3 split_gff.py
grep '^ChrUnknown' iwgsc_refseqv2.1_annotation_200916_HC.gtf | cat iwgsc_refseqv2.1_annotation_200916_HC_part.gtf - >iwgsc_refseqv2.1_annotation_200916_HC_unknown_part.gtf


grep '>' Triticum_aestivum_paragon.GCA949126075v1.cdna.all.fa | awk 'match($0,/primary_assembly:[^:]+:([^:]+:[0-9]+:[0-9]+):/,m){print m[1]}'  | sed 's/:/\t/g' > paragon_cdna.bed
awk '$3 == "exon" {print $1, $4, $5}'  Triticum_aestivum_paragon.GCA949126075v1.62.gff3 >  Triticum_aestivum_paragon.GCA949126075v1.62_exon.bed
grep 'scaffold' Triticum_aestivum_paragon.GCA949126075v1.62_exon.bed | sed 's/ /\t/g'| cat Triticum_aestivum_paragon.GCA949126075v1.62_exon_part.bed - > Triticum_aestivum_paragon.GCA949126075v1.62_scaf_exon_part.bed
gffread Triticum_aestivum_paragon.GCA949126075v1.62.gff3 -T -o Triticum_aestivum_paragon.GCA949126075v1.62.gtf 
cut -f 1-2 Paragon_part.fa.fai >Paragon_part_chr_sizes.txt
sed -i '/scaffold/d' Paragon_part_chr_sizes.txt
python3 split_gff_par.py
grep '^scaffold' Triticum_aestivum_paragon.GCA949126075v1.62.gtf | cat Triticum_aestivum_paragon.GCA949126075v1.62_part.gtf - >Triticum_aestivum_paragon.GCA949126075v1.62_scaf_part.gtf
awk '$3 == "exon" {print $1, $4, $5}'  Triticum_aestivum_paragon.GCA949126075v1.62.gff3 >  Triticum_aestivum_paragon.GCA949126075v1.62_exon.bed
awk '$3 == "transcript"' Triticum_aestivum_paragon.GCA949126075v1.62_scaf_part.gtf > transcripts_par_part.gtf
awk -F'\t' 'BEGIN{OFS="\t"} /^#/ || $3=="gene" {print}' Triticum_aestivum_paragon.GCA949126075v1.62.gff3 > genes_par.gff3
cut -f 2 SingleCopyOrthologues_matrix.tsv | tail -n +2 | sed 's/\.[^.]*$//' | grep -F -f - transcripts_par_part.gtf | cut -f 1,4,5 | sort -u  >one_one_orthologs_par.bed
python3 split_bed_par.py
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
10. Identify heterzygous sites
```
/software/bcftools-1.16/bcftools view  -R iwgsc_refseqv2.1_annotation_200916_HC_exon_unknown_part.bed -i 'QUAL>=20 && N_ALT>=1 && COUNT(GT!="mis" && FMT/DP>=20 && FMT/GQ>=20)>0' -Oz -o wheat_ase_het_snps_filtered.vcf.gz wheat.ase.snps.vcf.gz
/software/bcftools-1.16/bcftools view  -i 'QUAL>=10 && N_ALT>=1 && COUNT(GT!="mis" && FMT/DP>=10 )>0' -Oz -o wheat_het_snps_filtered.vcf.gz wheat.ase.snps.vcf.gz
/software/htslib-1.16/tabix -p vcf wheat_het_snps_filtered.vcf.gz
gunzip wheat_het_snps_filtered.vcf.gz

/software/vcftools-vcftools-581c231/bin/vcftools --vcf wheat_ase_het_snps_filtered.vcf --positions filtered_set_CS.txt --recode --recode-INFO-all --out wheat_ase_snps_het
/software/htslib-1.16/bgzip wheat_ase_snps_het.recode.vcf
/software/htslib-1.16/tabix wheat_ase_snps_het.recode.vcf.gz

/software/bcftools-1.16/bcftools view  -R Triticum_aestivum_paragon.GCA949126075v1.62_scaf_exon_part.bed -i 'QUAL>=20 && N_ALT>=1 && COUNT(GT!="mis" && FMT/DP>=20 && FMT/GQ>=20)>0' -Oz -o par_ase_het_snps_filtered.vcf.gz par.ase.snps.vcf.gz
/software/bcftools-1.16/bcftools view  -i 'QUAL>=10 && N_ALT>=1 && COUNT(GT!="mis" && FMT/DP>=10 )>0' -Oz -o par_het_snps_filtered.vcf.gz par.ase.snps.vcf.gz
/software/htslib-1.16/tabix -p vcf par_het_snps_filtered.vcf.gz
gunzip par_het_snps_filtered.vcf.gz

/software/vcftools-vcftools-581c231/bin/vcftools --vcf par_ase_het_snps_filtered.vcf --positions filtered_set_PAR.txt --recode --recode-INFO-all --out par_ase_snps_het
/software/htslib-1.16/bgzip par_ase_snps_het.recode.vcf
/software/htslib-1.16/tabix par_ase_snps_het.recode.vcf.gz
```

11. Profile ASE
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
