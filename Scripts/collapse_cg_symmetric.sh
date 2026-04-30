#!/usr/bin/env bash

infile="$1"
outfile="${2:-${infile%.txt.gz}.collapsed.txt.gz}"

gzip -dc "$infile" | awk '
BEGIN { OFS="\t" }
!h {
    chr=$1; pos=$2; strand=$3; m=$4; u=$5; h=1
    next
}
chr==$1 && strand=="+" && $3=="-" && ($2-pos)==1 {
    print chr, pos, m+$4, u+$5
    h=0
    next
}
{
    print "Error: unexpected non-paired rows:", chr, pos, strand, m, u, "/", $1, $2, $3, $4, $5 > "/dev/stderr"
    exit 1
}
END {
    if (h) {
        print "Error: unpaired final row:", chr, pos, strand, m, u > "/dev/stderr"
        exit 1
    }
}
' | gzip > "$outfile"
