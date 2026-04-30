#!/usr/bin/env bash
set -euo pipefail

infile="$1"
outfile="${2:-${infile%.txt.gz}.collapsed.txt.gz}"

gzip -dc "$infile" | awk '
BEGIN { OFS = "\t" }

{
    if (!have_prev) {
        p_chr = $1; p_pos = $2; p_m = $4; p_u = $5
        have_prev = 1
    } else {
        print p_chr, p_pos, p_m + $4, p_u + $5
        have_prev = 0
    }
}
' | gzip > "$outfile"
