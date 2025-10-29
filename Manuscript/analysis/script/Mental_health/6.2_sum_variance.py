import numpy as np
import pandas as pd
from scipy.stats import norm
import re
from scipy.stats import chi2
from tqdm import tqdm



file = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/trait.txt"
df = pd.read_csv(file, sep="\t", header=0)

sample_size_lst = []
for trait in tqdm(df.iloc[:, 0]):
    file = f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/pheno/{trait}.txt"
    dfT = pd.read_csv(file, header=0)
    sample_size_lst.append(dfT.shape[0])

res_lst = []
for trait in tqdm(df.iloc[:, 0]):
    file = f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/mom/{trait}.var"
    arr = np.loadtxt(file)
    h2_add = arr[0] / np.sum(arr)
    h2_gxe = arr[1] / np.sum(arr)
    h2_nxe = arr[2] / np.sum(arr)
    h2_e = arr[3] / np.sum(arr)
    res_lst.append([trait, h2_add, h2_gxe, h2_nxe, h2_e])

df_res1 = pd.DataFrame(res_lst, columns=["trait", "h2_add", "h2_gxe", "h2_nxe", "h2_e"])


res_lst = []
for trait in tqdm(df.iloc[:, 0]):
    file = f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/testGxE/{trait}.var"
    arr = np.loadtxt(file)
    arr = arr[:, 0]
    h2_add = arr[0] / np.sum(arr)
    h2_gxe = arr[1] / np.sum(arr)
    h2_nxe = arr[2] / np.sum(arr)
    h2_e = arr[3] / np.sum(arr)
    res_lst.append([trait, h2_add, h2_gxe, h2_nxe, h2_e])
df_res2 = pd.DataFrame(res_lst, columns=["trait", "h2_add_REML", "h2_gxe_REML", "h2_nxe_REML", "h2_e_REML"])


df_res = pd.merge(df_res1, df_res2, on="trait")

df_res["sample_size"] = sample_size_lst

df_res.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/results/mom/h2_mom.csv", index=False)

'''
# -----------------------------
# Summary (numeric cols only)
# -----------------------------
num_df = df_res[["h2_add", "h2_gxe", "h2_nxe", "h2_e"]]
print("Quantiles (25/50/75%):")
print(num_df.quantile([0.25, 0.50, 0.75]).round(4))

print("\nDescribe:")
print(num_df.describe(percentiles=[0.25, 0.5, 0.75]).round(4))

print("\nMin / Max for each component:")
mins = num_df.min().round(4)
maxs = num_df.max().round(4)
print("Min:\n", mins)
print("Max:\n", maxs)
'''
