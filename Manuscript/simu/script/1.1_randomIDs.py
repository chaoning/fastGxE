import pandas as pd
import numpy as np

dfE = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/envi.csv", usecols=["eid"], dtype=str)

dfFam = pd.read_csv("/net/zootopia/disk1/chaon/data/UKB/imp/ukb_imp_info.qc.fam", sep=r"\s+",
                    header=None, names=["fid", "eid"], usecols=[0, 1], dtype=str)

dfm = pd.merge(dfFam, dfE, on="eid", how="inner")

for nID in [50_000, 100_000, 200_000, 300_000, 400_000]:
    dfs = dfm.sample(n=nID, replace=False, random_state=1)
    dfs.to_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/random_{nID}.fam", sep="\t", index=False, 
               header=False)
