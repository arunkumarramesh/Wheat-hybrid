#!/usr/bin/env bash
set -euo pipefail

infile="$1"
prefix="${2:-${infile%.txt.gz}}"

gzip -dc "$infile" | awk \
    -v chh="${prefix}.CHH.txt.gz" \
    -v chg_sym="${prefix}.CHG_symmetric.txt.gz" \
    -v chg_other="${prefix}.CHG_other.txt.gz" \
    -v cg_sym="${prefix}.CG_symmetric.txt.gz" \
    -v cg_other="${prefix}.CG_other.txt.gz" '
BEGIN {
    OFS = "\t"
    chh_cmd = "gzip > " chh
    chg_sym_cmd = "gzip > " chg_sym
    chg_other_cmd = "gzip > " chg_other
    cg_sym_cmd = "gzip > " cg_sym
    cg_other_cmd = "gzip > " cg_other
}

function print5(cmd, c1, c2, c3, c4, c5) {
    print c1, c2, c3, c4, c5 | cmd
}

function print_cov(cmd, c1, c2, c3, c4, c5) {
    if (c4 + c5 > 2) print5(cmd, c1, c2, c3, c4, c5)
}

$6 == "CHH" {
    print_cov(chh_cmd, $1, $2, $3, $4, $5)
    next
}

$6 == "CHG" {
    if (have_chg && p_chg_chr == $1 && p_chg_strand == "+" && $3 == "-" && ($2 - p_chg_pos) == 2) {
        print5(chg_sym_cmd, p_chg_chr, p_chg_pos, p_chg_strand, p_chg_m, p_chg_u)
        print5(chg_sym_cmd, $1, $2, $3, $4, $5)
        have_chg = 0
    } else {
        if (have_chg) print_cov(chg_other_cmd, p_chg_chr, p_chg_pos, p_chg_strand, p_chg_m, p_chg_u)
        p_chg_chr = $1; p_chg_pos = $2; p_chg_strand = $3; p_chg_m = $4; p_chg_u = $5
        have_chg = 1
    }
    next
}

$6 == "CG" {
    if (have_cg && p_cg_chr == $1 && p_cg_strand == "+" && $3 == "-" && ($2 - p_cg_pos) == 1) {
        print5(cg_sym_cmd, p_cg_chr, p_cg_pos, p_cg_strand, p_cg_m, p_cg_u)
        print5(cg_sym_cmd, $1, $2, $3, $4, $5)
        have_cg = 0
    } else {
        if (have_cg) print_cov(cg_other_cmd, p_cg_chr, p_cg_pos, p_cg_strand, p_cg_m, p_cg_u)
        p_cg_chr = $1; p_cg_pos = $2; p_cg_strand = $3; p_cg_m = $4; p_cg_u = $5
        have_cg = 1
    }
}

END {
    if (have_chg) print_cov(chg_other_cmd, p_chg_chr, p_chg_pos, p_chg_strand, p_chg_m, p_chg_u)
    if (have_cg) print_cov(cg_other_cmd, p_cg_chr, p_cg_pos, p_cg_strand, p_cg_m, p_cg_u)
}
'
