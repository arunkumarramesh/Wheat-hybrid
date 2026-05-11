#!/usr/bin/env python3

import sys
import gzip
from collections import defaultdict

bed_file = sys.argv[1]
input_file = sys.argv[2]
output_file = sys.argv[3]

te_cols = ["te_id","te_strand","te_class","te_family","te_consensus","te_consensus_pct","te_status","te_copie","te_compo","te_post","te_length","subgenome","nearest_gene_id_v1.1","nearest_gene_id_v2.1","distance_to_gene","centromere_status","centromere_interval"]

intervals = defaultdict(list)

with open(bed_file, "rt") as f:
    for line in f:
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        meta = fields[3:20]
        if len(meta) != len(te_cols):
            sys.exit(f"ERROR: expected 20 columns in TE BED, but found {len(fields)} columns:\n{line}")
        key = (chrom,start,end,tuple(meta))
        intervals[chrom].append((start,end,key,meta))

for chrom in intervals:
    intervals[chrom].sort(key=lambda x: (x[0],x[1]))

pointers = {chrom: 0 for chrom in intervals}
active = {chrom: [] for chrom in intervals}

summary = {}
order = []

with gzip.open(input_file, "rt") as fin:
    header = fin.readline().rstrip("\n").split()
    idx = {name: i for i,name in enumerate(header)}
    required = ["chr","pos","pct_CS","cov_CS","pct_CSxP","cov_CSxP","pct_P","cov_P"]
    missing = [x for x in required if x not in idx]
    if missing:
        sys.exit("ERROR: missing required columns in methylation file header: " + ",".join(missing))

    for line in fin:
        line = line.rstrip("\n")
        if not line:
            continue
        fields = line.split()
        chrom = fields[idx["chr"]]
        pos = int(fields[idx["pos"]])
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
        for start,end,key,meta in chrom_active:
            if start < pos <= end:
                if key not in summary:
                    summary[key] = {"n_sites": 0,"meth_CS": 0.0,"cov_CS": 0.0,"meth_CSxP": 0.0,"cov_CSxP": 0.0,"meth_P": 0.0,"cov_P": 0.0}
                    order.append(key)
                pct_CS = float(fields[idx["pct_CS"]])
                cov_CS = float(fields[idx["cov_CS"]])
                pct_CSxP = float(fields[idx["pct_CSxP"]])
                cov_CSxP = float(fields[idx["cov_CSxP"]])
                pct_P = float(fields[idx["pct_P"]])
                cov_P = float(fields[idx["cov_P"]])
                summary[key]["n_sites"] += 1
                summary[key]["meth_CS"] += pct_CS * cov_CS / 100.0
                summary[key]["cov_CS"] += cov_CS
                summary[key]["meth_CSxP"] += pct_CSxP * cov_CSxP / 100.0
                summary[key]["cov_CSxP"] += cov_CSxP
                summary[key]["meth_P"] += pct_P * cov_P / 100.0
                summary[key]["cov_P"] += cov_P

with gzip.open(output_file, "wt") as fout:
    fout.write("chr\tstart\tend\t" + "\t".join(te_cols) + "\tn_sites\tpct_CS\tcov_CS\tpct_CSxP\tcov_CSxP\tpct_P\tcov_P\n")
    for key in order:
        chrom,start,end,meta_tuple = key
        meta = list(meta_tuple)
        s = summary[key]
        pct_CS = 100.0 * s["meth_CS"] / s["cov_CS"] if s["cov_CS"] > 0 else "NA"
        pct_CSxP = 100.0 * s["meth_CSxP"] / s["cov_CSxP"] if s["cov_CSxP"] > 0 else "NA"
        pct_P = 100.0 * s["meth_P"] / s["cov_P"] if s["cov_P"] > 0 else "NA"
        fout.write("\t".join(map(str,[chrom,start,end] + meta + [s["n_sites"],pct_CS,s["cov_CS"],pct_CSxP,s["cov_CSxP"],pct_P,s["cov_P"]])) + "\n")