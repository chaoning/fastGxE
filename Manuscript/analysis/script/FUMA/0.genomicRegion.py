import numpy as np
import pandas as pd
import sys

in_file = sys.argv[1]
out_file = sys.argv[2]

with open(in_file) as fin, open(out_file, "w") as fout:
    head_line = fin.readline()
    head_lst = head_line.split()
    fout.write(",".join(["GenomicLocus", "chrom", "start",	"end",	"nSNPs", "nGWASSNPs", "IndSigSNPs", "LeadSNPs"]) + "\n")
    for line in fin:
        lst = line.split()
        chrom = lst[head_lst.index("chr")]
        start = lst[head_lst.index("start")]
        end = lst[head_lst.index("end")]
        nSNPs = lst[head_lst.index("nSNPs")]
        nGWASSNPs = lst[head_lst.index("nGWASSNPs")]
        IndSigSNPs = lst[head_lst.index("IndSigSNPs")]
        LeadSNPs = lst[head_lst.index("LeadSNPs")]
        GenomicLocus = f"{chrom}:{start}-{end}"
        fout.write(",".join([GenomicLocus, chrom, start, end, nSNPs, nGWASSNPs, IndSigSNPs, LeadSNPs]) + "\n")
