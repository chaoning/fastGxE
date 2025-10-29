import numpy as np
import pandas as pd
import sys


df = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/PT_PIP.csv")
print(df.shape[0])
df["trait"] = np.array(df["trait"], dtype=str)

trait_lst = []
leading_lst = []
envi_lst = []
pip_lst = []
for i in range(df.shape[0]):
    dfSub = df.iloc[i, -42:].copy()
    if df.iloc[i, -43] < 0:
        dfSub = -dfSub
    envi = list(dfSub[(dfSub > 0.5) | (dfSub < -0.5)].index)
    pip_lst.extend(list(dfSub[(dfSub > 0.5) | (dfSub < -0.5)]))
    envi_lst.extend(envi)
    num_envi = len(envi)
    trait_lst.extend([df.iloc[i, 13]] * num_envi)
    leading_lst.extend([df.iloc[i, 14]] * num_envi)

dfRes = pd.DataFrame({
    "trait": trait_lst,
    "leadingSNP": leading_lst,
    "envi": envi_lst,
    "pip": pip_lst
})

dfRes.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/PT.PIP.sum.loci.envi.csv", index=False)

df = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/BB_PIP.csv")
print(df.shape[0])

trait_lst = []
leading_lst = []
envi_lst = []
pip_lst = []
for i in range(df.shape[0]):
    dfSub = df.iloc[i, -42:].copy()
    if df.iloc[i, -43] < 0:
        dfSub = -dfSub
    envi = list(dfSub[(dfSub > 0.5) | (dfSub < -0.5)].index)
    pip_lst.extend(list(dfSub[(dfSub > 0.5) | (dfSub < -0.5)]))
    envi_lst.extend(envi)
    num_envi = len(envi)
    trait_lst.extend([df.iloc[i, 13]] * num_envi)
    leading_lst.extend([df.iloc[i, 14]] * num_envi)

dfRes = pd.DataFrame({
    "trait": trait_lst,
    "leadingSNP": leading_lst,
    "envi": envi_lst,
    "pip": pip_lst
})

dfRes.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/BB.PIP.sum.loci.envi.csv", index=False)

