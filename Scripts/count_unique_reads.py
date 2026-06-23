#!/usr/bin/env python3

import sys
import re
from collections import defaultdict

counts = defaultdict(int)
qname = None
hits = defaultdict(set)

def good(f):
    flag = int(f[1])
    return (
        flag & 2 and
        not flag & 4 and
        not flag & 8 and
        not flag & 2048 and
        "NM:i:0" in f[11:] and
        not re.search(r"[SHIDN]", f[5])
    )

def finish():
    valid = [tx for tx, mates in hits.items() if mates == {1, 2}]
    if len(valid) == 1:
        counts[valid[0]] += 1

for line in sys.stdin:
    f = line.rstrip().split("\t")

    if qname is not None and f[0] != qname:
        finish()
        hits = defaultdict(set)

    qname = f[0]

    if good(f):
        if int(f[1]) & 64:
            hits[f[2]].add(1)
        elif int(f[1]) & 128:
            hits[f[2]].add(2)

finish()

for tx in sorted(counts):
    print(tx, counts[tx], sep="\t")
