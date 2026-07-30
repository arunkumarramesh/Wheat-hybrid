
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

grep '>' iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna.fasta | cut -d ' ' -f 1 | sed 's/>//' >transnames
sed 's/\..*//g' transnames | paste -d ',' transnames - >transcript_to_gene_refseqv2.1.csv

for file in *_1.paired.fq.gz; do /software/kallisto/build/src/kallisto quant -i iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna_index -o ${file/_1.paired.fq.gz/_CS} -t 20  $file ${file/_1.paired.fq.gz/_2.paired.fq.gz} ; done
ls -d *CS/ | sed 's/\///' >cs_kallisto_samplenames.txt
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

cd /projects/wheat/orthofinder_input_2_A/OrthoFinder/Results_Jul05_1/Orthologues/Orthologues_iwgsc_refseqv2.1_annotation_200916_HC_mrna_A
awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) {if ($i ~ /iwgsc_refseqv2.1_annotation_200916_HC_mrna_A/) iwgsc=i; if ($i ~ /Triticum_aestivum_paragon.GCA949126075v1.cdna.all_A/) paragon=i} next} $iwgsc != "" && $paragon != "" {n=split($iwgsc, genes, /, */); for (j=1; j<=n; j++) {gsub(/\r/,"",genes[j]); gsub(/[[:space:]]+$/,"",genes[j]); sub(/\.[0-9]+$/,"",genes[j]); print genes[j]}}'  iwgsc_refseqv2.1_annotation_200916_HC_mrna_A__v__Triticum_aestivum_paragon.GCA949126075v1.cdna.all_A.tsv >A_CS_Par_match.txt
cp  A_CS_Par_match.txt ../../../../../

cd /projects/wheat/orthofinder_input_2_B/OrthoFinder/Results_Apr24/Orthologues/Orthologues_iwgsc_refseqv2.1_annotation_200916_HC_mrna_B
awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) {if ($i ~ /iwgsc_refseqv2.1_annotation_200916_HC_mrna_B/) iwgsc=i; if ($i ~ /Triticum_aestivum_paragon.GCA949126075v1.cdna.all_B/) paragon=i} next} $iwgsc != "" && $paragon != "" {n=split($iwgsc, genes, /, */); for (j=1; j<=n; j++) {gsub(/\r/,"",genes[j]); gsub(/[[:space:]]+$/,"",genes[j]); sub(/\.[0-9]+$/,"",genes[j]); print genes[j]}}'  iwgsc_refseqv2.1_annotation_200916_HC_mrna_B__v__Triticum_aestivum_paragon.GCA949126075v1.cdna.all_B.tsv  >B_CS_Par_match.txt
cp  B_CS_Par_match.txt ../../../../../

cd /projects/wheat/orthofinder_input_2_D/OrthoFinder/Results_Apr24/Orthologues/Orthologues_iwgsc_refseqv2.1_annotation_200916_HC_mrna_D
awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) {if ($i ~ /iwgsc_refseqv2.1_annotation_200916_HC_mrna_D/) iwgsc=i; if ($i ~ /Triticum_aestivum_paragon.GCA949126075v1.cdna.all_D/) paragon=i} next} $iwgsc != "" && $paragon != "" {n=split($iwgsc, genes, /, */); for (j=1; j<=n; j++) {gsub(/\r/,"",genes[j]); gsub(/[[:space:]]+$/,"",genes[j]); sub(/\.[0-9]+$/,"",genes[j]); print genes[j]}}'  iwgsc_refseqv2.1_annotation_200916_HC_mrna_D__v__Triticum_aestivum_paragon.GCA949126075v1.cdna.all_D.tsv  >D_CS_Par_match.txt
cp  D_CS_Par_match.txt ../../../../../

cd ../../../../../
cat A_CS_Par_match.txt B_CS_Par_match.txt D_CS_Par_match.txt >CS_Par_match.txt
```
5. Summarise counts per gene using [`combine_sample_tpm.R`](./combine_sample_tpm.R). Script from Philippa Borrill. 

6. Run differential expression tests using [`de_wheat.R`](./de_wheat.R). Based on preliminary analyses done by Cris https://github.com/crisforgiarini/Data-Gene-expression-and-methylation-in-intraspecific-hybrids-of-hexaploidy-wheat

7. Identify genes showing allele-specific expression using combined transcriptome reference. Then do ASE analyses using [`ase_new.R`](./ase_new.R)
```
python3 longest_transcript_ref.py iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna.fasta  iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna_longest.fasta
python3 longest_transcript_ref.py Triticum_aestivum_paragon.GCA949126075v1.cdna.all.fa Triticum_aestivum_paragon.GCA949126075v1.cdna.longest.fa
cat iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna_longest.fasta Triticum_aestivum_paragon.GCA949126075v1.cdna.longest.fa >iwgsc2.1_paragon.fasta
samtools faidx iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna_longest.fasta
samtools faidx Triticum_aestivum_paragon.GCA949126075v1.cdna.longest.fa
hisat2-build -p 20 iwgsc2.1_paragon.fasta iwgsc2.1_paragon
for file in *_1.paired.fq.gz; do hisat2  -p 20  --no-spliced-alignment --no-softclip --score-min L,0,0  --no-mixed  --no-discordant -k 100 --secondary -x iwgsc2.1_paragon -1 $file -2 ${file/_1.paired.fq.gz/_2.paired.fq.gz} -S ${file/_1.paired.fq.gz/.cs.sam} ; done
for file in *.cs.sam; do samtools sort -n -@ 10 -O bam -o ${file/.cs.sam/.sortname.cs.bam} $file; done
for file in *.sortname.cs.bam; do samtools view $file| python3 count_unique_reads.py > ${file/sortname.cs.bam/tsv} ; done
```

8. Compare homoeolog expression bias for differentially expressed genes using [`assign_homoeolog_expression_bias_categories.R`](./assign_homoeolog_expression_bias_categories.R) and [`HEB.R`](./HEB.R)

9. Trim bisulfite reads
```
for file in P-1_R1.fq.gz; do java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 -threads 20 $file ${file/_R1.fq.gz/_R2.fq.gz} ${file/_R1.fq.gz/_1.paired.fq.gz} ${file/_R1.fq.gz/_1.unpaired.fq.gz} ${file/_R1.fq.gz/_2.paired.fq.gz} ${file/_R1.fq.gz/_2.unpaired.fq.gz} ILLUMINACLIP:TruSeq3-PE_sailgene.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36; done
```

10. Index genome, map bisulfite reads and deduplicate
```
# genome folder has 161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta
bismark_genome_preparation --hisat2 --verbose --parallel 5 genome

for file in *.paired.fq.gz ; do bismark --multicore 4 --hisat2 --genome_folder genome -1 $file -2 ${file/_1.paired.fq.gz/_2.paired.fq.gz}  ; done
for file in *_bismark_hisat2_pe.bam; do deduplicate_bismark -p --bam $file ; done
```

11. Extract methylation counts
```
for file in *.deduplicated.bam; do bismark_methylation_extractor --multicore 4 --gzip --bedGraph --buffer_size 280G --CX --genome_folder genome $file; done
for file in *.deduplicated.bam; do coverage2cytosine --gzip --genome_folder genome --coverage_threshold 1 --CX -o ${file/.paired_bismark_hisat2_pe.deduplicated.bismark.cov.gz/} $file ; done
```

12. Merge methylation counts from the three replicate libraries
```
./merge_cx_reports.sh P-1_1.CX_report.txt.gz P-2_1.CX_report.txt.gz P-3_1.CX_report.txt.gz P_combined.CX_report.txt.gz
./merge_cx_reports.sh CS-1_1.CX_report.txt.gz CS-2_1.CX_report.txt.gz CS-3_1.CX_report.txt.gz CS_combined.CX_report.txt.gz
./merge_cx_reports.sh CSxP-1_1.CX_report.txt.gz CSxP-2_1.CX_report.txt.gz CSxP-3_1.CX_report.txt.gz CSxP_combined.CX_report.txt.gz
```

13. From merged methylation count files, split into files containing seperate cytosine contexts using split_cx_report.sh
```
./split_cx_report.sh P_combined.CX_report.txt.gz
./split_cx_report.sh CS_combined.CX_report.txt.gz
./split_cx_report.sh CSxP_combined.CX_report.txt.gz
```
14. For each CG methylation site pair (consecutive sites), it sums the methylated and unmethylated counts across both strands
```
./collapse_cg_symmetric.sh P_combined.CX_report.CG_symmetric.txt.gz P_combined.CG_symmetric_collapsed.txt.gz
./collapse_cg_symmetric.sh CS_combined.CX_report.CG_symmetric.txt.gz CS_combined.CG_symmetric_collapsed.txt.gz
./collapse_cg_symmetric.sh CSxP_combined.CX_report.CG_symmetric.txt.gz CSxP_combined.CG_symmetric_collapsed.txt.gz
```

15. For each CHG methylation site pair (two sites apart), it sums the methylated and unmethylated counts across both strands
```
./collapse_chg_symmetric.sh P_combined.CX_report.CHG_symmetric.txt.gz P_combined.CHG_symmetric_collapsed.txt.gz
./collapse_chg_symmetric.sh CS_combined.CX_report.CHG_symmetric.txt.gz CS_combined.CHG_symmetric_collapsed.txt.gz
./collapse_chg_symmetric.sh CSxP_combined.CX_report.CHG_symmetric.txt.gz CSxP_combined.CHG_symmetric_collapsed.txt.gz
```

16. Merge methylation counts from all three samples into a single file
```
./methylation_merge.sh CS_combined.CX_report.CHH.txt.gz CSxP_combined.CX_report.CHH.txt.gz P_combined.CX_report.CHH.txt.gz merged_CHH_sites.txt.gz
./methylation_merge.sh CS_combined.CX_report.CHG_other.txt.gz CSxP_combined.CX_report.CHG_other.txt.gz P_combined.CX_report.CHG_other.txt.gz merged_CHG_other_sites.txt.gz
./methylation_merge.sh CS_combined.CX_report.CG_other.txt.gz CSxP_combined.CX_report.CG_other.txt.gz P_combined.CX_report.CG_other.txt.gz merged_CG_other_sites.txt.gz

./methylation_merge_sym.sh CS_combined.CG_symmetric_collapsed.txt.gz CSxP_combined.CG_symmetric_collapsed.txt.gz P_combined.CG_symmetric_collapsed.txt.gz merged_CG_symmetric.txt.gz
./methylation_merge_sym.sh CS_combined.CHG_symmetric_collapsed.txt.gz CSxP_combined.CHG_symmetric_collapsed.txt.gz P_combined.CHG_symmetric_collapsed.txt.gz merged_CHG_symmetric.txt.gz

```

17. Convert chromosome part coordinates into full genome coordinates for methylation sites
```
awk 'BEGIN{FS=OFS="\t"} NR==FNR{c[$1]=$4; o[$1]=$5; next} FNR==1{print; next} {$2=$2+o[$1]; $1=c[$1]; print}' 161010_Chinese_Spring_v1.0_pseudomolecules_parts_to_chr.bed <(zcat merged_CG_symmetric.txt.gz) | gzip > merged_CG_symmetric_fullchr.txt.gz
awk 'BEGIN{FS=OFS="\t"} NR==FNR{c[$1]=$4; o[$1]=$5; next} FNR==1{print; next} {$2=$2+o[$1]; $1=c[$1]; print}' 161010_Chinese_Spring_v1.0_pseudomolecules_parts_to_chr.bed <(zcat merged_CHG_symmetric.txt.gz) | gzip > merged_CHG_symmetric_fullchr.txt.gz
awk 'BEGIN{FS=OFS="\t"} NR==FNR{c[$1]=$4; o[$1]=$5; next} FNR==1{print; next} {$2=$2+o[$1]; $1=c[$1]; print}' 161010_Chinese_Spring_v1.0_pseudomolecules_parts_to_chr.bed <(zcat merged_CHH_sites.txt.gz) | gzip > merged_CHH_fullchr.txt.gz

```

18. Identify SNPs differences between CS and Paragon reference genomes. Whole genome alignments from Ensembl Plants
```
cd taes_iwgsc.v.tapa_gca949126075v1.lastz_net/
for file in taes*.maf; do python3 maf_snps_cs_paragon.py $file >${file/.maf/_snps.txt}; done
cat *_snps.txt | sed '/cs_src/d' > ../cs_par_snps.txt
awk 'BEGIN{FS=OFS="\t"} $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $3==$2+1{chr=$1; sub(/^triticum_aestivum\./,"chr",chr); print chr,$2,$3}' cs_par_snps.txt | LC_ALL=C sort -k1,1 -k2,2n -k3,3n -u > cs_par_all_snps.bed

```

19. Remove methylation sites with SNPs between CS and Paragon reference genomes. The classify sites based on inheritance categories. Classification scheme based on scripts developed by Asena: https://github.com/AsenaArdaman/Hybrid_inheritance_models.
```
(zcat merged_CG_symmetric_fullchr.txt.gz | head -n 1; zcat merged_CG_symmetric_fullchr.txt.gz | awk 'BEGIN{FS=OFS="\t"} FNR>1{print $1,$2-1,$2+1,$0}' | bedtools intersect -sorted -a - -b cs_par_all_snps.bed -v | cut -f4-) | gzip -c > merged_CG_symmetric_all.txt.gz

(zcat merged_CHG_symmetric_fullchr.txt.gz | head -n 1; zcat merged_CHG_symmetric_fullchr.txt.gz | awk 'BEGIN{FS=OFS="\t"} FNR>1{print $1,$2-1,$2+2,$0}' | bedtools intersect -sorted -a - -b cs_par_all_snps.bed -v | cut -f4-) | gzip -c > merged_CHG_symmetric_all.txt.gz

python3 remove_CHH_snps.py

Rscript boman_classification_cg.R
Rscript boman_classification_chg.R
Rscript boman_classification_chh.R

```

20. Convert IWGSC v1.1 gene annotation into BED files for CDS for the longest transcript and 1 kb promoter regions, while replacing v1.1 gene IDs with their high-confidence v2.1 gene IDs using [`bed_intervals.sh`](./bed_intervals.sh)

21. Subset CDS regions from methylation sites. Remove any duplicate positions. Data available on https://doi.org/10.6084/m9.figshare.32144041.
```
(printf "chr\tpos\tpct_CS\tcov_CS\tpct_CSxP\tcov_CSxP\tpct_P\tcov_P\tgene_id\n"; zcat merged_CG_symmetric_all.txt.gz | awk 'BEGIN{FS=OFS="\t"} NR>1{print $1,$2-1,$2,$0}' | bedtools intersect -a stdin -b CDS.bed  -wa -wb | awk 'BEGIN{FS=OFS="\t"}{print $4,$5,$6,$7,$8,$9,$10,$11,$15}') | gzip > merged_CG_symmetric_CDS.txt.gz
zcat merged_CG_symmetric_CDS.txt.gz | awk 'NR==1 || !seen[$1 FS $2]++' | gzip > tmp && mv tmp merged_CG_symmetric_CDS.txt.gz

(printf "chr\tpos\tpct_CS\tcov_CS\tpct_CSxP\tcov_CSxP\tpct_P\tcov_P\tgene_id\n"; zcat merged_CHG_symmetric_all.txt.gz | awk 'BEGIN{FS=OFS="\t"} NR>1{print $1,$2-1,$2,$0}' | bedtools intersect -a stdin -b CDS.bed  -wa -wb | awk 'BEGIN{FS=OFS="\t"}{print $4,$5,$6,$7,$8,$9,$10,$11,$15}') | gzip > merged_CHG_symmetric_CDS.txt.gz
zcat merged_CHG_symmetric_CDS.txt.gz | awk 'NR==1 || !seen[$1 FS $2]++' | gzip > tmp && mv tmp merged_CHG_symmetric_CDS.txt.gz

(printf "chr\tpos\tstrand\tpct_CS\tcov_CS\tpct_CSxP\tcov_CSxP\tpct_P\tcov_P\tgene_id\n"; zcat merged_CHH_all.txt.gz | awk 'BEGIN{FS=OFS="\t"} NR>1{print $1,$2-1,$2,$0}' | bedtools intersect -a stdin -b CDS.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{print $4,$5,$6,$7,$8,$9,$10,$11,$12,$16}') | gzip > merged_CHH_all_CDS.txt.gz
zcat merged_CHH_all_CDS.txt.gz | awk 'NR==1 || !seen[$1 FS $2]++' | gzip > tmp && mv tmp merged_CHH_all_CDS.txt.gz

#old
#python3 subset_chh_by_cds.py CDS.bed merged_CHH_all.txt.gz merged_CHH_all_CDS.txt.gz
```

22. Plot methylation results using [`boman_classification_gene.R`](./boman_classification_gene.R) and [`gbM_wheat.R`](./gbM_wheat.R)

