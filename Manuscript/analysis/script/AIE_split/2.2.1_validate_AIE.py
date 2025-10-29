import logging
import sys

import numpy as np
import pandas as pd
import os

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/AIE_split/mmsusie/")
trait = sys.argv[1]
snp = sys.argv[2]


file = f"{trait}.{snp}.alpha.txt"
alpha = np.loadtxt(file)
alpha = alpha.reshape(-1, 42)

file = f"{trait}.{snp}.mu.txt"
mu = np.loadtxt(file)
mu = mu.reshape(-1, 42)

effect = np.sum(alpha * mu, axis=0)


dfP = pd.read_csv(f"../pheno/{trait}.part2.txt", sep=r"\s+")


score = np.dot(dfP.loc[:, "Age":"Confide"].values, effect)
print(dfP.loc[:, "Age":"Confide"].values.shape, effect.shape, score.shape)

dfP_out = dfP.copy()
dfP_out["score"] = score

try:
    # Try to divide the samples into 5 quantile-based groups by score
    dfP_out["score_group"] = pd.qcut(dfP_out["score"], q=5, labels=False)
    
    # If fewer than 5 groups are created (due to duplicate values), fall back to random assignment
    if dfP_out["score_group"].nunique() < 5:
        raise ValueError("Less than 5 unique groups formed by qcut")
except Exception as e:
    print(f"qcut failed due to: {e}. Falling back to random grouping.")
    
    # Randomly assign each sample to one of 5 groups (0–4)
    dfP_out["score_group"] = np.random.randint(0, 5, size=len(dfP_out))


for i in range(5):
    dfP_out[dfP_out["score_group"] == i].to_csv(f"../split5/{trait}.{snp}.score_group_{i+1}.txt", sep="\t", index=False)
