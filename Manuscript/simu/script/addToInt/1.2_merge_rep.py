import sys
import pandas as pd

prefix = sys.argv[1]
dfE = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/envi.csv")
dfP = pd.read_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_addToint/{prefix}.rep_1.txt",
                    sep=r"\s+")
dfP = pd.merge(dfP.iloc[:, [1]], dfE, on="eid", how="left")

for rep in range(1, 101):
    file = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_addToint/{prefix}.rep_{rep}.txt"
    dfP[f"trait{rep}"] = pd.read_csv(file, sep=r"\s+").iloc[:, 2].copy()

if dfP.isnull().values.any():
    print("There are null values in the merged DataFrame. Please check the input files.")
    sys.exit(1)

dfP.to_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_addToint_fastgxe/{prefix}.txt",
           index=False, sep=" ")
