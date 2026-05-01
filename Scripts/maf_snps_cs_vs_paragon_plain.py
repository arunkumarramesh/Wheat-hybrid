#!/usr/bin/env python3
import sys
import gzip

CS_PREFIXES = ["triticum_aestivum.", "taes_iwgsc."]
PAR_PREFIXES = ["triticum_aestivum_paragon.", "tapa_gca949126075v1."]
VALID = set("ACGT")

def open_maybe_gzip(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

def starts_with_any(s, prefixes):
    return any(s.startswith(p) for p in prefixes)

def parse_s_line(line):
    f = line.rstrip("\n").split()
    return {"src": f[1], "start": int(f[2]), "size": int(f[3]), "strand": f[4], "src_size": int(f[5]), "text": f[6].upper()}

def block_to_snps(cs, par, out):
    cs_text, par_text = cs["text"], par["text"]
    cs_pos = cs["start"]
    par_pos = par["start"] if par["strand"] == "+" else par["src_size"] - (par["start"] + par["size"])

    for a, b in zip(cs_text, par_text):
        cs_consumes, par_consumes = a != "-", b != "-"

        if cs_consumes and par_consumes and a in VALID and b in VALID and a != b:
            out.write(f"{cs['src']}\t{cs_pos}\t{cs_pos + 1}\t{a}\t{b}\t{par['src']}\t{par_pos}\t{par['strand']}\n")

        if cs_consumes: cs_pos += 1
        if par_consumes: par_pos += 1

def process_file(path, out):
    current_block = []

    def flush_block(lines):
        s_rows = [parse_s_line(line) for line in lines if line.startswith("s ")]
        if len(s_rows) < 2: return

        cs = par = None
        for rec in s_rows:
            if cs is None and starts_with_any(rec["src"], CS_PREFIXES): cs = rec
            elif par is None and starts_with_any(rec["src"], PAR_PREFIXES): par = rec

        if cs is None or par is None:
            if len(s_rows) == 2: cs, par = s_rows
            else: return

        block_to_snps(cs, par, out)

    with open_maybe_gzip(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                flush_block(current_block); current_block = []; continue
            if line.startswith("#"): continue
            if line.startswith("a"):
                flush_block(current_block); current_block = [line]; continue
            current_block.append(line)

        flush_block(current_block)

def main():
    sys.stdout.write("cs_src\tcs_pos0\tcs_pos1\tcs_base\tpar_base\tpar_src\tpar_pos0\tpar_strand\n")
    process_file(sys.argv[1], sys.stdout)

if __name__ == "__main__":
    main()
