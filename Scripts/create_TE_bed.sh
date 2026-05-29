#!/usr/bin/env bash

TE_GFF="iwgsc_refseqv1.0_TransposableElements_2017Mar13.gff3"
GENE_GFF="IWGSC_v1.1_HC_20170706.gff3"
CENTROMERES="centromeres.bed"
CORRESPONDENCE="../iwgsc_refseq_all_correspondances.csv"

awk -v OFS="\t" '
BEGIN{FS="\t"}
function get_attr(attrs,key,n,i,a,kv){
    n=split(attrs,a,";");
    for(i=1;i<=n;i++){
        gsub(/^ +| +$/,"",a[i]);
        split(a[i],kv,"=");
        if(kv[1]==key)return kv[2]
    }
    return "."
}
$0!~/^#/ && ($3=="match" || $3=="repeat_region") && $1~/^[Cc]hr[1-7][ABD]$/ {
    chr=$1; sub(/^chr/,"Chr",chr);
    id=get_attr($9,"ID"); compo=get_attr($9,"compo"); copie=get_attr($9,"copie");
    post=get_attr($9,"post"); status=get_attr($9,"status");
    split(compo,c,/ +/);
    te_consensus=c[1]; te_consensus_pct=c[2];
    te_family=te_consensus; sub(/\.[0-9]+$/,"",te_family);
    split(te_family,a,"_"); te_class=a[1];
    gsub(/[ \t]+/,",",compo); gsub(/[ \t]+/,",",post);
    subgenome=substr(chr,length(chr),1);
    te_length=$5-$4+1;
    print chr,$4-1,$5,id,".",$7,te_class,te_family,te_consensus,te_consensus_pct,status,copie,compo,post,te_length,subgenome
}' "$TE_GFF" | sort -k1,1 -k2,2n > TEs.metadata.bed

awk -v OFS="\t" '
BEGIN{FS="\t"}
function get_attr(attrs,key,n,i,a,kv){
    n=split(attrs,a,";");
    for(i=1;i<=n;i++){
        gsub(/^ +| +$/,"",a[i]);
        split(a[i],kv,"=");
        if(kv[1]==key)return kv[2]
    }
    return "."
}
$0!~/^#/ && $3=="gene" && $1~/^[Cc]hr[1-7][ABD]$/ {
    chr=$1; sub(/^chr/,"Chr",chr);
    id=get_attr($9,"ID");
    print chr,$4-1,$5,id,".",$7
}' "$GENE_GFF" | sort -k1,1 -k2,2n > genes.bed

bedtools closest -a TEs.metadata.bed -b genes.bed -d -t first > TEs.with_nearest_gene.raw.bed

awk -v OFS="\t" '
FNR==NR{
    chr=$1; n[chr]++;
    cen_start[chr,n[chr]]=$2;
    cen_end[chr,n[chr]]=$3;
    cen_id[chr,n[chr]]=$4;
    next
}
{
    chr=$1; mid=int(($2+$3)/2);
    centromere_status="non_centromeric";
    centromere_interval=".";
    for(i=1;i<=n[chr];i++){
        if(mid>=cen_start[chr,i] && mid<=cen_end[chr,i]){
            centromere_status="centromeric";
            centromere_interval=cen_id[chr,i];
            break
        }
    }
    print $1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$20,$23,centromere_status,centromere_interval
}' "$CENTROMERES" TEs.with_nearest_gene.raw.bed > TEs.with_nearest_gene.centromere.bed

awk -v OFS="\t" '
FNR==NR{
    n=split($0,a,/[[:space:]]+/);
    if(a[1]=="v1.0" && a[2]=="v1.1" && a[3]=="v2.1")next;
    if(a[2]!="-" && a[3]!="-" && a[2]!~/LC$/ && a[3]!~/LC$/)map[a[2]]=a[3];
    next
}
{
    gene_v11=$16;
    gene_v21=(gene_v11 in map ? map[gene_v11] : ".");
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,gene_v11,gene_v21,$17,$18,$19
}' "$CORRESPONDENCE" TEs.with_nearest_gene.centromere.bed > TEs.bed
