## File descriptions

### Reference and annotation files

| File                                                                                                                                 | Description                                                                                                                               |
| ------------------------------------------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------- |
| [`161010_Chinese_Spring_v1.0_pseudomolecules_parts_to_chr.bed`](161010_Chinese_Spring_v1.0_pseudomolecules_parts_to_chr.bed)         | Coordinates and sizes used to convert Chinese Spring IWGSC RefSeq v1.0 pseudomolecule parts to complete chromosomes.                      |
| [`CDS.bed`](CDS.bed)                                                                                                                 | Coding-sequence intervals for the longest transcript of each Chinese Spring gene.                                                         |
| [`SingleCopyOrthologues_matrix.tsv`](SingleCopyOrthologues_matrix.tsv)                                                               | One-to-one orthologues between Chinese Spring IWGSC RefSeq v2.1 and Paragon GCA949126075v1.                                               |
| [`Triticum_aestivum_paragon.GCA949126075v1.cdna.longest.fa.fai`](Triticum_aestivum_paragon.GCA949126075v1.cdna.longest.fa.fai)       | FASTA index for the longest Paragon cDNA sequence assigned to each gene.                                                                  |
| [`iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna_longest.fasta.fai`](iwgsc_refseqv2.1_annotation_200916_HC_LC_mrna_longest.fasta.fai) | FASTA index for the longest Chinese Spring IWGSC RefSeq v2.1 mRNA sequence assigned to each gene.                                         |
| [`homoeologs_1_1_1_synt_and_non_synt.csv`](homoeologs_1_1_1_synt_and_non_synt.csv)                                                   | Syntenic and non-syntenic 1:1:1 wheat homoeologous triads classified by Ramírez-González et al., *Science*, DOI: 10.1126/science.aar6089. |
| [`bias_category_all_samples_inc_orig_expr.csv`](bias_category_all_samples_inc_orig_expr.csv)                                         | Homoeolog expression values and expression-bias classifications for triads across the parental and hybrid samples.                        |
| [`iwgsc_refseq_all_correspondances.csv.zip`](iwgsc_refseq_all_correspondances.csv.zip)                                               | Table with IWGSC RefSeq v1.1 and v2.1 gene identifiers, obtained from URGI.                                                |
| [`core_genes.txt`](core_genes.txt)                                                                                                   | List of core wheat genes identified in DOI: 10.1038/s41467-025-64046-1.                                                                   |
| [`TruSeq3-PE_sailgene.fa`](TruSeq3-PE_sailgene.fa)                                                                                   | Adapter sequences used for methyl-seq read trimming.                                                                                      |
| [`adaptors_novogene.fa`](adaptors_novogene.fa)                                                                                       | Adapter sequences used for RNA-seq read trimming.                                                                                         |

### Expression quantification files

| File                                         | Description                                                                                               |
| -------------------------------------------- | --------------------------------------------------------------------------------------------------------- |
| [`cs_count.tsv`](cs_count.tsv)               | Gene-level RNA-seq read counts estimated using the Chinese Spring IWGSC RefSeq v2.1 reference.            |
| [`cs_tpm.tsv`](cs_tpm.tsv)                   | Gene-level transcripts-per-million values estimated using the Chinese Spring IWGSC RefSeq v2.1 reference. |
| [`cs_gene_lengths.csv`](cs_gene_lengths.csv) | Gene lengths used for expression quantification with the Chinese Spring IWGSC RefSeq v2.1 reference.      |

### Differential expression results

| File                                                       | Description                                                                                                                                                                                                                              |
| ---------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [`CSvP all genes.csv`](CSvP%20all%20genes.csv)             | Differential expression results comparing Chinese Spring and Paragon using Chinese Spring IWGSC RefSeq v2.1 gene identifiers.                                                                                                            |
| [`CS_PvCSxP all genes.csv`](CS_PvCSxP%20all%20genes.csv)   | Differential expression results comparing the mid-parental expression estimate with the Chinese Spring × Paragon hybrids using Chinese Spring IWGSC RefSeq v2.1 gene identifiers.                                                        |
| [`triads_CS_P.csv`](triads_CS_P.csv)                       | Homoeologous triads containing at least one gene differentially expressed between Chinese Spring and Paragon. “Up” indicates higher expression in Chinese Spring, whereas “Down” indicates higher expression in Paragon.                 |
| [`triads_hybrids_parents.csv`](triads_hybrids_parents.csv) | Homoeologous triads containing at least one gene differentially expressed between the mid-parental estimate and the hybrids. “Up” indicates higher expression in the hybrids, whereas “Down” indicates higher expression in the parents. |
| [`additive_genes.txt`](additive_genes.txt)                 | Genes classified as showing additive expression in the hybrids relative to parental expression.                                                                                                                                          |
| [`dominant_genes.txt`](dominant_genes.txt)                 | Genes classified as showing expression dominance, with hybrid expression resembling one of the parents.                                                                                                                                  |
| [`transgressive_genes.txt`](transgressive_genes.txt)       | Genes for which hybrid expression was significantly higher or lower than expression in both parents.                                                                                                                                     |
| [`nonDE_genes.txt`](nonDE_genes.txt)                       | Genes without a significant expression difference between the hybrids and the parental comparison.                                                                                                                                       |

### Allele-specific expression results

| File or folder                     | Description                                                                                                                                                  |
| ---------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| [`ASE data/`](ASE%20data/)         | Folder containing Chinese Spring- and Paragon-specific read counts in the hybrid samples after mapping reads to a combined parental transcriptome reference. |
| [`Ref_vs_Alt.csv`](Ref_vs_Alt.csv) | Results of linear-model tests comparing Chinese Spring and Paragon allele-specific expression in the hybrids.                                                |

### Homoeolog expression and triad summaries

| File                                                       | Description                                                                                                 |
| ---------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------- |
| [`triads_CS_P.csv`](triads_CS_P.csv)                       | Triad-level summary of homoeolog expression differences between Chinese Spring and Paragon.                 |
| [`triads_hybrids_parents.csv`](triads_hybrids_parents.csv) | Triad-level summary of homoeolog expression differences between the hybrids and the mid-parental estimates. |

### DNA methylation files

| File                                                     | Description                                                                                                               |
| -------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| [`gene_body_methylation.tsv`](gene_body_methylation.tsv) | Classification of genes as gene-body methylated or non-gene-body methylated in Chinese Spring, Paragon, and their hybrid. |
