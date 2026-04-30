#!/usr/bin/env bash

f1="$1"
f2="$2"
f3="$3"
out="$4"

SEP=$'\034'

LC_ALL=C join -t $'\t' \
    <(
        gzip -dc "$f1" | awk -v OFS='\t' -v sep="$SEP" '
        {
            cov = $3 + $4
            if (cov > 2)
                printf "%s%s%s\t%.6f\t%d\n", $1, sep, $2, ($3/cov)*100, cov
        }'
    ) \
    <(
        gzip -dc "$f2" | awk -v OFS='\t' -v sep="$SEP" '
        {
            cov = $3 + $4
            if (cov > 2)
                printf "%s%s%s\t%.6f\t%d\n", $1, sep, $2, ($3/cov)*100, cov
        }'
    ) |
LC_ALL=C join -t $'\t' - \
    <(
        gzip -dc "$f3" | awk -v OFS='\t' -v sep="$SEP" '
        {
            cov = $3 + $4
            if (cov > 2)
                printf "%s%s%s\t%.6f\t%d\n", $1, sep, $2, ($3/cov)*100, cov
        }'
    ) |
awk -v OFS='\t' -v sep="$SEP" '
BEGIN {
    print "chr","pos","pct_CS","cov_CS","pct_CSxP","cov_CSxP","pct_P","cov_P"
}
{
    split($1, a, sep)
    print a[1], a[2], $2, $3, $4, $5, $6, $7
}' | gzip > "$out"
