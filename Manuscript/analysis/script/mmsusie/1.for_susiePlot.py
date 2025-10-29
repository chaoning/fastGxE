import numpy as np
import pandas as pd
import sys


df = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/PT_PIP.csv")
print(df.shape[0])
# df["Alcohol"] = -df["Alcohol"]
df["trait"] = np.array(df["trait"], dtype=str)
index_df = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/trait32.txt",
                       sep="\t")
index_dct = dict(zip(index_df.iloc[:, 0], index_df.iloc[:, 2]))
traitName_lst = [index_dct[trait] for trait in df["trait"]]
df["TraitName"] = traitName_lst


locus_trait_lst = []
for locus, trait in zip(df["Genes"], df["TraitName"]):
    locus_trait_lst.append(f"{locus}:{trait}")

df = df.iloc[:, -43:].copy()
df.index = locus_trait_lst

df.loc[df.loc[:, "Main"] < 0, :] = -df.loc[df.loc[:, "Main"] < 0, :]

bool_arr = (df < 0)
df[bool_arr] = -df[bool_arr]
df[(df < 0.5)] = 0
df[(df >= 0.9)] = 2
df[(df >= 0.5) & (df < 0.9)] = 1
df[bool_arr] = -df[bool_arr]


df.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/PT_PIP_plot.csv")


df = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/BB_PIP.csv")
print(df.shape[0])
# df["Alcohol"] = -df["Alcohol"]
index_df = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/trait67.txt",
                       sep="\t")
index_dct = dict(zip(index_df.iloc[:, 0], index_df.iloc[:, 2]))
traitName_lst = [index_dct[trait] for trait in df["trait"]]
df["TraitName"] = traitName_lst


locus_trait_lst = []
for locus, trait, leadsnp in zip(df["Genes"], df["TraitName"], df["trait_leading_snp"]):
    locus_trait_lst.append(f"{locus}:{trait}:{leadsnp}")

df = df.iloc[:, -43:].copy()
df.index = locus_trait_lst

df.loc[df.loc[:, "Main"] < 0, :] = -df.loc[df.loc[:, "Main"] < 0, :]

bool_arr = (df < 0)
df[bool_arr] = -df[bool_arr]
df[(df < 0.5)] = 0
df[(df >= 0.9)] = 2
df[(df >= 0.5) & (df < 0.9)] = 1
df[bool_arr] = -df[bool_arr]


df.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/BB_PIP_plot.csv")
