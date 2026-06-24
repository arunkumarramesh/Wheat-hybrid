## File Descriptions

### Reference and annotation files

| File | Description |
|---|---|
| [`adaptors_novogene.fa`](adaptors_novogene.fa) | Adapter sequences used for RNA-seq read trimming. |
| [`TruSeq3-PE_sailgene.fa`](TruSeq3-PE_sailgene.fa) | Adapter sequences used for methyl-seq read trimming. |
| [`SingleCopyOrthologues_matrix.tsv`](SingleCopyOrthologues_matrix.tsv) | One-to-one orthologues identified between Chinese Spring IWGSC RefSeq v2.1 and Paragon GCA949126075v1. |
| [`BP.csv`](BP.csv), [`MF.csv`](MF.csv), [`CC.csv`](CC.csv) | Gene Ontology annotation files for IWGSC RefSeq v2.1 genes. |
| [`homoeologs_1_1_1_synt_and_non_synt.csv`](homoeologs_1_1_1_synt_and_non_synt.csv) | Triad classification from Ramírez-González et al., *Science*, DOI: 10.1126/science.aar6089. |
| [`bias_category_all_samples_inc_orig_expr.csv`](bias_category_all_samples_inc_orig_expr.csv) | Triad identification and homoeolog expression-bias categories. |
| [`iwgsc_refseq_all_correspondances.csv`](iwgsc_refseq_all_correspondances.csv) | v1.1 to v2.1 IWGSC RefSeq gene IDs from urgi.versailles.inrae.fr. |
| [`161010_Chinese_Spring_v1.0_pseudomolecules_parts_to_chr.bed`](161010_Chinese_Spring_v1.0_pseudomolecules_parts_to_chr.bed)  | Chromosome part sizes for IWGSC v1.0. |
| [`CDS.bed`](CDS.bed)  | CDS interval for longest transcript |
| [`promoter1kb.bed`](promoter1kb.bed)  | 1Kb Promoter interval |
| [`cs_par_snps_50kb_counts.tsv`](cs_par_snps_50kb_counts.tsv)  | Number of SNVs per 50kb between IWGSC v1.0 and GCA949126075v1 reference.|
| [`chromatin_states.txt.zip`](chromatin_states.txt.zip)  | File containing IWGSC v1.0 coordinates for chromatin states 1-4 and 13 defined in Li et al 2019 https://doi.org/10.1186/s13059-019-1746-8 |
| [`choulet_gene_tissue.tsv`](choulet_gene_tissue.tsv) | Gene IDs and expression breadth classification |
| [`core_genes.txt`](core_genes.txt) | Core wheat genes identified in https://doi.org/10.1038/s41467-025-64046-1. |


### Expression quantification files

| File | Description |
|---|---|
| [`cs_count.tsv`](cs_count.tsv) | Gene read counts estimated using the Chinese Spring IWGSC RefSeq v2.1 reference. |
| [`cs_tpm.tsv`](cs_tpm.tsv) | Transcript-per-million estimates generated using the Chinese Spring IWGSC RefSeq v2.1 reference. |
| [`cs_gene_lengths.csv`](cs_gene_lengths.csv) | Gene lengths for the Chinese Spring IWGSC RefSeq v2.1 reference. |


### Differential expression results

| File | Description |
|---|---|
| [`CSvP all genes.csv`](CSvP%20all%20genes.csv) | Differential expression results for Chinese Spring vs Paragon using IWGSC RefSeq v2.1. |
| [`CS_PvCSxP all genes.csv`](CS_PvCSxP%20all%20genes.csv) | Differential expression results for parents vs hybrids using IWGSC RefSeq v2.1. |
| [`triads_CS_P.csv`](triads_CS_P.csv) | Triads with one or more homoeologs differentially expressed between Chinese Spring and Paragon using IWGSC RefSeq v2.1. “Up” indicates higher expression in hybrids; “Down” indicates higher expression in parents. |
| [`triads_hybrids_parents.csv`](triads_hybrids_parents.csv) | Triads with one or more homoeologs differentially expressed between mid-parental estimates and hybrids using IWGSC RefSeq v2.1. “Up” indicates higher expression in CS; “Down” indicates higher expression in Paragon. |

### Allele-specific expression results

| File or folder | Description |
|---|---|
| [`ASE data/`](ASE%20data/) | Folder containing Chinese Spring and Paragon read counts for hybrid genotypes after mapping to a combined parental transcriptome reference. |
| [`Ref_vs_Alt.csv`](Ref_vs_Alt.csv) | Linear-model test for allele-specific expression. |
| [`classified_limma.csv`](classified_limma.csv) | Regulatory classifications based on linear-model ASE tests. |


### Methylation files

| File | Description |
|---|---|
| [`gene_level_TE_CG.tsv`](gene_level_TE_CG.tsv) | Mean CG methylation of transposable elements near bread wheat genes. |
| [`gene_level_TE_CHG.tsv`](gene_level_TE_CHG.tsv) | Mean CHG methylation of transposable elements near bread wheat genes. |
| [`gene_level_TE_CHH.tsv`](gene_level_TE_CHH.tsv) | Mean CHH methylation of transposable elements near bread wheat genes. |
| [`CDS_meth_pct.txt`](CDS_meth_pct.txt) | Mean CG, CHG and CHH methylation across CDS of bread wheat genes. |
| [`promoter_meth_pct.txt`](promoter_meth_pct.txt) | Mean CG, CHG and CHH methylation across promoters of bread wheat genes. |
| [`gene_body_methylation.tsv`](gene_body_methylation.tsv) | Classification of genes as gbM or non-gbM in Chinese Spring, Paragon, and their hybrid. |
