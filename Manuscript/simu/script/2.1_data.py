import sys
import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import sys

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno/")
dfBim = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink.100K.chip.bim", sep=r"\s+", header=None)
nSNP = dfBim.shape[0]
with open("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink.100K.chip.anno", "w") as fout:
    for i in range(nSNP):
        fout.write("1\n")

dfFam = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink.100K.chip.fam", sep=r"\s+", header=None,
                    usecols=[0, 1], names=["FID", "IID"])

# prefix = f"h2_add_30.h2_gxe_{gxe}.h2_nxe_15.sample_100000"
prefix = sys.argv[1]
dfE = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/envi.csv")
dfE.rename(columns={"eid": "FID"}, inplace=True)
dfE.insert(1, "IID", dfE["FID"])
for rep in tqdm(range(1, 101)):
    dfE2 = pd.merge(dfFam, dfE, on=["FID", "IID"], how="left")
    dfE2.to_csv(f"../pheno_GENIE/{prefix}.E.{rep}.txt", sep=" ", index=False)
    file = f"{prefix}.rep_{rep}.txt"
    df = pd.read_csv(file, sep=r"\s+")
    df.columns = ["FID", "IID", "Pheno"]
    df = pd.merge(dfFam, df, on=["FID", "IID"], how="left")
    df.to_csv(f"../pheno_GENIE/{prefix}.{rep}.txt", sep=" ", index=False)
