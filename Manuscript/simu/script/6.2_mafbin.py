import pandas as pd
import numpy as np

df = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink.frq", sep=r"\s+")


maf_bin = [
    [0.01, 0.05],
    [0.05, 0.1],
    [0.1, 0.2],
    [0.2, 0.3],
    [0.3, 0.4],
    [0.4, 0.5]
]

with open("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/QQplot_maf/mafbin.txt", "w") as f:
    for i in range(len(maf_bin)):
        f.write(f"{maf_bin[i][0]} {maf_bin[i][1]}\n")


for i in range(len(maf_bin)):
    df_maf = df[(df["MAF"] >= maf_bin[i][0]) & (df["MAF"] < maf_bin[i][1])]
    df_maf.to_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/QQplot_maf/mafbin_{i+1}.frq", sep="\t", index=False)
    print(f"MAF bin {i+1}: {len(df_maf)} SNPs")
