#!/usr/bin/env bash
set -euo pipefail

gff="IWGSC_v1.1_HC_20170706.gff3"
map="iwgsc_refseq_all_correspondances.csv"
core="shell_v1.txt"
out_cds="CDS_shell.bed"
out_prom="promoter1kb_shell.bed"

awk -v OFS="\t" -v map="$map" -v out_cds="$out_cds" -v out_prom="$out_prom" '
BEGIN { FS = "[[:space:]]+" }
function get_attr(attrs, key,   n, i, a, kv) {
    n = split(attrs, a, ";")
    for (i = 1; i <= n; i++) {
        split(a[i], kv, "=")
        if (kv[1] == key) return kv[2]
    }
    return ""
}
FNR == NR {
    core_gene[$1] = 1
    next
}
FILENAME == map {
    if ($1 == "v1.0" && $2 == "v1.1" && $3 == "v2.1") next
    if (!($2 in core_gene) || $2 ~ /LC$/ || $3 ~ /LC$/ || $2 == "-" || $3 == "-") next
    map_v11_to_v21[$2] = $3
    next
}
/^##sequence-region/ {
    chr_min[$2] = $3
    chr_max[$2] = $4
    next
}
/^#/ { next }
$3 == "gene" {
    gene_id_v11 = get_attr($9, "ID")
    gene_id_v21 = map_v11_to_v21[gene_id_v11]
    if (gene_id_v21 == "") next
    gene_v11_to_v21[gene_id_v11] = gene_id_v21

    if ($7 == "+") {
        prom_start_1 = $4 - 1000
        prom_end_1 = $4 - 1
    } else if ($7 == "-") {
        prom_start_1 = $5 + 1
        prom_end_1 = $5 + 1000
    } else next

    if ($1 in chr_min) {
        if (prom_start_1 < chr_min[$1]) prom_start_1 = chr_min[$1]
        if (prom_end_1 > chr_max[$1]) prom_end_1 = chr_max[$1]
    }
    if (prom_start_1 <= prom_end_1) print $1, prom_start_1 - 1, prom_end_1, gene_id_v21, ".", $7 > out_prom
    next
}
$3 == "mRNA" || $3 == "transcript" {
    tx_id = get_attr($9, "ID")
    parent_gene_v11 = get_attr($9, "Parent")
    if (tx_id != "" && parent_gene_v11 != "") tx_to_gene_v11[tx_id] = parent_gene_v11
    next
}
$3 == "CDS" {
    parent_tx = get_attr($9, "Parent")
    gene_id_v11 = tx_to_gene_v11[parent_tx]
    gene_id_v21 = gene_v11_to_v21[gene_id_v11]
    if (parent_tx == "" || gene_id_v11 == "" || gene_id_v21 == "") next
    cds_len[parent_tx] += $5 - $4 + 1
    tx_gene_v21[parent_tx] = gene_id_v21
    n_cds[parent_tx]++
    cds_line[parent_tx, n_cds[parent_tx]] = $1 OFS ($4 - 1) OFS $5 OFS gene_id_v21 OFS "." OFS $7
    next
}
END {
    for (tx in cds_len) {
        gene = tx_gene_v21[tx]
        if (!(gene in best_tx) || cds_len[tx] > best_len[gene]) {
            best_tx[gene] = tx
            best_len[gene] = cds_len[tx]
        }
    }
    for (gene in best_tx) {
        tx = best_tx[gene]
        for (i = 1; i <= n_cds[tx]; i++) print cds_line[tx, i] > out_cds
    }
}
' "$core" "$map" "$gff"


gff="IWGSC_v1.1_HC_20170706.gff3"
map="iwgsc_refseq_all_correspondances.csv"
core="cloud_v1.txt"
out_cds="CDS_cloud.bed"
out_prom="promoter1kb_cloud.bed"

awk -v OFS="\t" -v map="$map" -v out_cds="$out_cds" -v out_prom="$out_prom" '
BEGIN { FS = "[[:space:]]+" }
function get_attr(attrs, key,   n, i, a, kv) {
    n = split(attrs, a, ";")
    for (i = 1; i <= n; i++) {
        split(a[i], kv, "=")
        if (kv[1] == key) return kv[2]
    }
    return ""
}
FNR == NR {
    core_gene[$1] = 1
    next
}
FILENAME == map {
    if ($1 == "v1.0" && $2 == "v1.1" && $3 == "v2.1") next
    if (!($2 in core_gene) || $2 ~ /LC$/ || $3 ~ /LC$/ || $2 == "-" || $3 == "-") next
    map_v11_to_v21[$2] = $3
    next
}
/^##sequence-region/ {
    chr_min[$2] = $3
    chr_max[$2] = $4
    next
}
/^#/ { next }
$3 == "gene" {
    gene_id_v11 = get_attr($9, "ID")
    gene_id_v21 = map_v11_to_v21[gene_id_v11]
    if (gene_id_v21 == "") next
    gene_v11_to_v21[gene_id_v11] = gene_id_v21

    if ($7 == "+") {
        prom_start_1 = $4 - 1000
        prom_end_1 = $4 - 1
    } else if ($7 == "-") {
        prom_start_1 = $5 + 1
        prom_end_1 = $5 + 1000
    } else next

    if ($1 in chr_min) {
        if (prom_start_1 < chr_min[$1]) prom_start_1 = chr_min[$1]
        if (prom_end_1 > chr_max[$1]) prom_end_1 = chr_max[$1]
    }
    if (prom_start_1 <= prom_end_1) print $1, prom_start_1 - 1, prom_end_1, gene_id_v21, ".", $7 > out_prom
    next
}
$3 == "mRNA" || $3 == "transcript" {
    tx_id = get_attr($9, "ID")
    parent_gene_v11 = get_attr($9, "Parent")
    if (tx_id != "" && parent_gene_v11 != "") tx_to_gene_v11[tx_id] = parent_gene_v11
    next
}
$3 == "CDS" {
    parent_tx = get_attr($9, "Parent")
    gene_id_v11 = tx_to_gene_v11[parent_tx]
    gene_id_v21 = gene_v11_to_v21[gene_id_v11]
    if (parent_tx == "" || gene_id_v11 == "" || gene_id_v21 == "") next
    cds_len[parent_tx] += $5 - $4 + 1
    tx_gene_v21[parent_tx] = gene_id_v21
    n_cds[parent_tx]++
    cds_line[parent_tx, n_cds[parent_tx]] = $1 OFS ($4 - 1) OFS $5 OFS gene_id_v21 OFS "." OFS $7
    next
}
END {
    for (tx in cds_len) {
        gene = tx_gene_v21[tx]
        if (!(gene in best_tx) || cds_len[tx] > best_len[gene]) {
            best_tx[gene] = tx
            best_len[gene] = cds_len[tx]
        }
    }
    for (gene in best_tx) {
        tx = best_tx[gene]
        for (i = 1; i <= n_cds[tx]; i++) print cds_line[tx, i] > out_cds
    }
}
' "$core" "$map" "$gff"


gff="IWGSC_v1.1_HC_20170706.gff3"
map="iwgsc_refseq_all_correspondances.csv"
core="seedling_specific_v1.txt"
out_cds="CDS_seedling.bed"
out_prom="promoter1kb_seedling.bed"

awk -v OFS="\t" -v map="$map" -v out_cds="$out_cds" -v out_prom="$out_prom" '
BEGIN { FS = "[[:space:]]+" }
function get_attr(attrs, key,   n, i, a, kv) {
    n = split(attrs, a, ";")
    for (i = 1; i <= n; i++) {
        split(a[i], kv, "=")
        if (kv[1] == key) return kv[2]
    }
    return ""
}
FNR == NR {
    core_gene[$1] = 1
    next
}
FILENAME == map {
    if ($1 == "v1.0" && $2 == "v1.1" && $3 == "v2.1") next
    if (!($2 in core_gene) || $2 ~ /LC$/ || $3 ~ /LC$/ || $2 == "-" || $3 == "-") next
    map_v11_to_v21[$2] = $3
    next
}
/^##sequence-region/ {
    chr_min[$2] = $3
    chr_max[$2] = $4
    next
}
/^#/ { next }
$3 == "gene" {
    gene_id_v11 = get_attr($9, "ID")
    gene_id_v21 = map_v11_to_v21[gene_id_v11]
    if (gene_id_v21 == "") next
    gene_v11_to_v21[gene_id_v11] = gene_id_v21

    if ($7 == "+") {
        prom_start_1 = $4 - 1000
        prom_end_1 = $4 - 1
    } else if ($7 == "-") {
        prom_start_1 = $5 + 1
        prom_end_1 = $5 + 1000
    } else next

    if ($1 in chr_min) {
        if (prom_start_1 < chr_min[$1]) prom_start_1 = chr_min[$1]
        if (prom_end_1 > chr_max[$1]) prom_end_1 = chr_max[$1]
    }
    if (prom_start_1 <= prom_end_1) print $1, prom_start_1 - 1, prom_end_1, gene_id_v21, ".", $7 > out_prom
    next
}
$3 == "mRNA" || $3 == "transcript" {
    tx_id = get_attr($9, "ID")
    parent_gene_v11 = get_attr($9, "Parent")
    if (tx_id != "" && parent_gene_v11 != "") tx_to_gene_v11[tx_id] = parent_gene_v11
    next
}
$3 == "CDS" {
    parent_tx = get_attr($9, "Parent")
    gene_id_v11 = tx_to_gene_v11[parent_tx]
    gene_id_v21 = gene_v11_to_v21[gene_id_v11]
    if (parent_tx == "" || gene_id_v11 == "" || gene_id_v21 == "") next
    cds_len[parent_tx] += $5 - $4 + 1
    tx_gene_v21[parent_tx] = gene_id_v21
    n_cds[parent_tx]++
    cds_line[parent_tx, n_cds[parent_tx]] = $1 OFS ($4 - 1) OFS $5 OFS gene_id_v21 OFS "." OFS $7
    next
}
END {
    for (tx in cds_len) {
        gene = tx_gene_v21[tx]
        if (!(gene in best_tx) || cds_len[tx] > best_len[gene]) {
            best_tx[gene] = tx
            best_len[gene] = cds_len[tx]
        }
    }
    for (gene in best_tx) {
        tx = best_tx[gene]
        for (i = 1; i <= n_cds[tx]; i++) print cds_line[tx, i] > out_cds
    }
}
' "$core" "$map" "$gff"



gff="IWGSC_v1.1_HC_20170706.gff3"
map="iwgsc_refseq_all_correspondances.csv"
core="broad_v1.txt"
out_cds="CDS_broad.bed"
out_prom="promoter1kb_broad.bed"

awk -v OFS="\t" -v map="$map" -v out_cds="$out_cds" -v out_prom="$out_prom" '
BEGIN { FS = "[[:space:]]+" }
function get_attr(attrs, key,   n, i, a, kv) {
    n = split(attrs, a, ";")
    for (i = 1; i <= n; i++) {
        split(a[i], kv, "=")
        if (kv[1] == key) return kv[2]
    }
    return ""
}
FNR == NR {
    core_gene[$1] = 1
    next
}
FILENAME == map {
    if ($1 == "v1.0" && $2 == "v1.1" && $3 == "v2.1") next
    if (!($2 in core_gene) || $2 ~ /LC$/ || $3 ~ /LC$/ || $2 == "-" || $3 == "-") next
    map_v11_to_v21[$2] = $3
    next
}
/^##sequence-region/ {
    chr_min[$2] = $3
    chr_max[$2] = $4
    next
}
/^#/ { next }
$3 == "gene" {
    gene_id_v11 = get_attr($9, "ID")
    gene_id_v21 = map_v11_to_v21[gene_id_v11]
    if (gene_id_v21 == "") next
    gene_v11_to_v21[gene_id_v11] = gene_id_v21

    if ($7 == "+") {
        prom_start_1 = $4 - 1000
        prom_end_1 = $4 - 1
    } else if ($7 == "-") {
        prom_start_1 = $5 + 1
        prom_end_1 = $5 + 1000
    } else next

    if ($1 in chr_min) {
        if (prom_start_1 < chr_min[$1]) prom_start_1 = chr_min[$1]
        if (prom_end_1 > chr_max[$1]) prom_end_1 = chr_max[$1]
    }
    if (prom_start_1 <= prom_end_1) print $1, prom_start_1 - 1, prom_end_1, gene_id_v21, ".", $7 > out_prom
    next
}
$3 == "mRNA" || $3 == "transcript" {
    tx_id = get_attr($9, "ID")
    parent_gene_v11 = get_attr($9, "Parent")
    if (tx_id != "" && parent_gene_v11 != "") tx_to_gene_v11[tx_id] = parent_gene_v11
    next
}
$3 == "CDS" {
    parent_tx = get_attr($9, "Parent")
    gene_id_v11 = tx_to_gene_v11[parent_tx]
    gene_id_v21 = gene_v11_to_v21[gene_id_v11]
    if (parent_tx == "" || gene_id_v11 == "" || gene_id_v21 == "") next
    cds_len[parent_tx] += $5 - $4 + 1
    tx_gene_v21[parent_tx] = gene_id_v21
    n_cds[parent_tx]++
    cds_line[parent_tx, n_cds[parent_tx]] = $1 OFS ($4 - 1) OFS $5 OFS gene_id_v21 OFS "." OFS $7
    next
}
END {
    for (tx in cds_len) {
        gene = tx_gene_v21[tx]
        if (!(gene in best_tx) || cds_len[tx] > best_len[gene]) {
            best_tx[gene] = tx
            best_len[gene] = cds_len[tx]
        }
    }
    for (gene in best_tx) {
        tx = best_tx[gene]
        for (i = 1; i <= n_cds[tx]; i++) print cds_line[tx, i] > out_cds
    }
}
' "$core" "$map" "$gff"


