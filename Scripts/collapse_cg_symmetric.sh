#!/usr/bin/env bash

infile="$1"
outfile="${2:-${infile%.txt.gz}.collapsed.txt.gz}"

gzip -dc "$infile" | awk '
BEGIN { OFS="\t" }
NR % 2 == 1 {
    chr=$1; pos=$2; m=$4; u=$5
    next
}
NR % 2 == 0 {
    print chr, pos, m + $4, u + $5
}
' | gzip > "$outfile"
