import pandas as pd
import numpy as np
import sys

in_file = sys.argv[1]
bed_file = sys.argv[2]
out_file = sys.argv[3]

df = pd.read_csv(in_file)
dfBed = pd.read_csv(bed_file, sep="\s+", header=None)

dfBed.columns = ["chrom_hg38", "start_hg38", "end_hg38", "GenomicLocus"]
dfBed = dfBed.iloc[:, 1:].copy()

dfR = pd.merge(df, dfBed, on="GenomicLocus", how="left")
GenomicLocus_hg38_lst = []
for chrom, start_hg38, end_hg38 in zip(dfR["chrom"], dfR["start_hg38"], dfR["end_hg38"]):
    GenomicLocus_hg38_lst.append(f"{chrom}:{start_hg38}-{end_hg38}")
dfR["GenomicLocus_hg38"] = GenomicLocus_hg38_lst
dfR.to_csv(out_file, index=False)
