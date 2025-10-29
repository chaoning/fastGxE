import pandas as pd
import numpy as np

dfLoci = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/GxE_signal_h2/FUMA/GenomicRiskLoci.txt",
                     sep=r"\s+")

df = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/GxE_signal_h2/2.1_real_data_5e-8.csv")

for i in range(df.shape[0]):
    chrom = df.loc[i, "chrom"]
    pos = df.loc[i, "base"]
    dfSub = dfLoci[(dfLoci["chr"] == chrom) & (dfLoci["start"] <= pos) & (dfLoci["end"] >= pos)]
    if dfSub.shape[0] > 0:
        df.loc[i, "locus"] = ":".join(list(np.array(dfSub["GenomicLocus"].values, dtype=str)))
    else:
        df.loc[i, "locus"] = "NA"


dfSub = df[df["h2gxe"] > 0.12/100]

print(dfSub.shape[0], "\n", dfSub["locus"].value_counts(), "\n", dfSub["locus"].value_counts().shape[0])

dfSub = df[(df["h2gxe"] <= 0.12/100) & (df["h2gxe"] > 0.04/100)]

print(dfSub.shape[0], "\n", dfSub["locus"].value_counts(), "\n", dfSub["locus"].value_counts().shape[0])

dfSub = df[(df["h2gxe"] <= 0.04/100) & (df["h2gxe"] > 0.02/100)]

print(dfSub.shape[0], "\n", dfSub["locus"].value_counts(), "\n", dfSub["locus"].value_counts().shape[0])


dfSub = df[(df["h2gxe"] <= 0.02/100)]

print(dfSub.shape[0], "\n", dfSub["locus"].value_counts(), "\n", dfSub["locus"].value_counts().shape[0])

