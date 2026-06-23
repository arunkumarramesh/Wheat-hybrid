#!/usr/bin/env python3

import sys
import re

in_fa = sys.argv[1]
out_fa = sys.argv[2]

best = {}
header = None
seq = []

def save_record():
    tx = header[1:].split()[0]
    gene = re.sub(r"\.[0-9]+$", "", tx)
    s = "".join(seq)

    if gene not in best or len(s) > len(best[gene][1]):
        best[gene] = (header, s)

with open(in_fa) as f:
    for line in f:
        line = line.rstrip()

        if line.startswith(">"):
            if header is not None:
                save_record()
            header = line
            seq = []
        else:
            seq.append(line)

save_record()

with open(out_fa, "w") as out:
    for gene in sorted(best):
        header, seq = best[gene]
        out.write(header + "\n")
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60] + "\n")
