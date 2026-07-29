#!/usr/bin/env python3

import gzip
import sys

chh_file = "merged_CHH_fullchr.txt.gz"
snp_file = "cs_par_all_snps.bed"
output_file = "merged_CHH_all.txt.gz"

def read_next_snp(handle):
    for line in handle:
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        return fields[0], int(fields[1]), int(fields[2])
    return None

with open(snp_file, "r") as snps, gzip.open(chh_file, "rt") as chh, gzip.open(output_file, "wt") as out:
    out.write(next(chh))

    current_snp = read_next_snp(snps)

    for line_number, line in enumerate(chh, start=2):
        fields = line.rstrip("\n").split("\t")

        if len(fields) < 3:
            raise ValueError(f"Invalid CHH line {line_number}: {line.rstrip()}")

        chrom = fields[0]
        pos = int(fields[1])
        strand = fields[2]

        if strand == "+":
            site_start = pos - 1
            site_end = pos + 2
        elif strand == "-":
            site_start = max(0, pos - 3)
            site_end = pos
        else:
            raise ValueError(f"Invalid strand at line {line_number}: {strand}")

        while current_snp is not None and (
            current_snp[0] < chrom
            or (current_snp[0] == chrom and current_snp[2] <= site_start)
        ):
            current_snp = read_next_snp(snps)

        overlaps_snp = (
            current_snp is not None
            and current_snp[0] == chrom
            and current_snp[1] < site_end
            and current_snp[2] > site_start
        )

        if not overlaps_snp:
            out.write(line)
