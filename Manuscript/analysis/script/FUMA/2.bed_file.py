import pandas as pd
import numpy as np
import sys

in_file = sys.argv[1]
out_file = sys.argv[2]

df = pd.read_csv(in_file)
dfR = df.iloc[:, [1, 2, 3, 0]].copy()

for i in range(dfR.shape[0]):
    chrom = dfR.iloc[i, 0]
    chrom = f"chr{chrom}"
    dfR.iloc[i, 0] = chrom


dfR.to_csv(out_file, sep="\t", index=False, header=False)
