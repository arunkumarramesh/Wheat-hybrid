#!/usr/bin/env bash

f1="$1"
f2="$2"
f3="$3"
out="$4"

LC_ALL=C sort -m -k1,1 -k2,2n -k3,3 -k6,6 -k7,7 \
    <(gzip -dc "$f1") \
    <(gzip -dc "$f2") \
    <(gzip -dc "$f3") |
awk '
BEGIN { OFS = "\t" }

NR == 1 {
    chr=$1; pos=$2; strand=$3
    meth=$4; unmeth=$5
    ctx=$6; tri=$7
    next
}

{
    if ($1==chr && $2==pos && $3==strand && $6==ctx && $7==tri) {
        meth   += $4
        unmeth += $5
    } else {
        print chr, pos, strand, meth, unmeth, ctx, tri
        chr=$1; pos=$2; strand=$3
        meth=$4; unmeth=$5
        ctx=$6; tri=$7
    }
}

END {
    if (NR > 0)
        print chr, pos, strand, meth, unmeth, ctx, tri
}
' | gzip > "$out"
