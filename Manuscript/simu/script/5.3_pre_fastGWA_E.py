import pandas as pd

df = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/envi.csv")
df.insert(0, 'eid0', df["eid"])
df.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_fastGWA/cov.txt",
          sep=" ", index=False, header=False
)

dfE = df.loc[:, ["eid0", "eid"]]
for i in range(df.shape[1] - 2):
    dfE2 = dfE.copy()
    dfE2.insert(2, 'E', df.iloc[:, i + 2])
    dfE2.to_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_fastGWA/E.{i+1}.txt",
                sep=" ", index=False, header=False
    )
