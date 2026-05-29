#!/usr/bin/env bash

INPUT="$1"
OUTPUT="$2"
TE_BED="${3:-TEs.bed}"

(
printf "chr\tstart\tend\tte_id\tte_strand\tte_class\tte_family\tte_consensus\tte_consensus_pct\tte_status\tte_copie\tte_compo\tte_post\tte_length\tsubgenome\tnearest_gene_id\tdistance_to_gene\tcentromere_status\tcentromere_interval\tn_sites\tpct_CS\tcov_CS\tpct_CSxP\tcov_CSxP\tpct_P\tcov_P\n"

zcat "$INPUT" |
awk 'BEGIN{FS=OFS="\t"} NR>1{print $1,$2-1,$2,$0}' |
bedtools intersect -a stdin -b "$TE_BED" -wa -wb |
awk 'BEGIN{FS=OFS="\t"}
{
    key=$12 OFS $13 OFS $14 OFS $15 OFS $16 OFS $17 OFS $18 OFS $19 OFS $20 OFS $21 OFS $22 OFS $23 OFS $24 OFS $25 OFS $26 OFS $27 OFS $28 OFS $29 OFS $30
    if(!(key in seen)){
        seen[key]=1
        order[++n]=key
    }
    nsites[key]++
    meth_CS[key]+=$6*$7/100
    cov_CS[key]+=$7
    meth_CSxP[key]+=$8*$9/100
    cov_CSxP[key]+=$9
    meth_P[key]+=$10*$11/100
    cov_P[key]+=$11
}
END{
    for(i=1;i<=n;i++){
        key=order[i]
        pct_CS=(cov_CS[key]>0 ? 100*meth_CS[key]/cov_CS[key] : "NA")
        pct_CSxP=(cov_CSxP[key]>0 ? 100*meth_CSxP[key]/cov_CSxP[key] : "NA")
        pct_P=(cov_P[key]>0 ? 100*meth_P[key]/cov_P[key] : "NA")
        print key,nsites[key],pct_CS,cov_CS[key],pct_CSxP,cov_CSxP[key],pct_P,cov_P[key]
    }
}'
) | gzip > "$OUTPUT"
