#!/usr/bin/env bash
set -euo pipefail

f1="$1"
f2="$2"
f3="$3"
out="$4"

SEP=$'\034'

LC_ALL=C join -t $'\t' \
    <(
        gzip -dc "$f1" | awk -v OFS='\t' -v sep="$SEP" '
        {
            cov = $4 + $5
            pct = (cov > 0 ? ($4 / cov) * 100 : 0)
            printf "%s%s%s%s%s\t%.6f\t%d\n", $1, sep, $2, sep, $3, pct, cov
        }'
    ) \
    <(
        gzip -dc "$f2" | awk -v OFS='\t' -v sep="$SEP" '
        {
            cov = $4 + $5
            pct = (cov > 0 ? ($4 / cov) * 100 : 0)
            printf "%s%s%s%s%s\t%.6f\t%d\n", $1, sep, $2, sep, $3, pct, cov
        }'
    ) |
LC_ALL=C join -t $'\t' - \
    <(
        gzip -dc "$f3" | awk -v OFS='\t' -v sep="$SEP" '
        {
            cov = $4 + $5
            pct = (cov > 0 ? ($4 / cov) * 100 : 0)
            printf "%s%s%s%s%s\t%.6f\t%d\n", $1, sep, $2, sep, $3, pct, cov
        }'
    ) |
awk -v OFS='\t' -v sep="$SEP" '
BEGIN {
    print "chr","pos","strand","pct_CS","cov_CS","pct_CSxP","cov_CSxP","pct_P","cov_P"
}
{
    split($1, a, sep)
    print a[1], a[2], a[3], $2, $3, $4, $5, $6, $7
}' | gzip > "$out"
