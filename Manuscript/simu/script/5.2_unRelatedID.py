import pandas as pd
import numpy as np

dfU = pd.read_csv("/net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp_rm.agrm.id",
                  sep=r"\s+", header=None, names=["IID"])

for nID in [50000, 100000, 200000, 300000, 400000]:
    df = pd.read_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/random_{nID}.fam", sep=r"\s+", header=None,
                     names=["FID", "IID"])
    dfm = pd.merge(df, dfU, on="IID")
    print(f"Number of individuals in {nID} random sample: {len(dfm)}")
    dfm.to_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/random_{nID}_unRel.agrm.id", sep="\t",
                header=False, index=False)
