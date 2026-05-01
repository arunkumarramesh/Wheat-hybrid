#!/usr/bin/env python3

import sys
import gzip
from collections import defaultdict

bed_file = sys.argv[1]
input_file = sys.argv[2]
output_file = sys.argv[3]

intervals = defaultdict(list)

with open(bed_file, "rt") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        gene_id = fields[3]
        intervals[chrom].append((start, end, gene_id))

for chrom in intervals:
    intervals[chrom].sort(key=lambda x: (x[0], x[1]))

pointers = {chrom: 0 for chrom in intervals}
active = {chrom: [] for chrom in intervals}

with gzip.open(input_file, "rt") as fin, gzip.open(output_file, "wt") as fout:
    header = fin.readline().rstrip("\n")
    fout.write(header + "\tgene_id\n")

    for line in fin:
        line = line.rstrip("\n")
        if not line:
            continue

        fields = line.split()
        chrom = fields[0]
        pos = int(fields[1])

        if chrom not in intervals:
            continue

        chrom_intervals = intervals[chrom]
        i = pointers[chrom]
        chrom_active = active[chrom]

        while i < len(chrom_intervals) and chrom_intervals[i][0] < pos:
            chrom_active.append(chrom_intervals[i])
            i += 1
        pointers[chrom] = i

        chrom_active = [iv for iv in chrom_active if pos <= iv[1]]
        active[chrom] = chrom_active

        genes = []
        for start, end, gene_id in chrom_active:
            if start < pos <= end:
                genes.append(gene_id)

        if genes:
            seen = set()
            uniq_genes = []
            for g in genes:
                if g not in seen:
                    uniq_genes.append(g)
                    seen.add(g)

            fout.write(line + "\t" + ",".join(uniq_genes) + "\n")
