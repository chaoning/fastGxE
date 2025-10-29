import logging
import sys

import numpy as np
import pandas as pd
import os

wk_dir = sys.argv[1]
out_dir = sys.argv[2]
os.chdir(wk_dir)

prefix = "h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000"

logging.basicConfig(level=logging.INFO)
rep = int(sys.argv[3])
input_window_index = int(sys.argv[4])

file = f"{rep}.{input_window_index}.susie.alpha.txt"
alpha = np.loadtxt(file)
alpha = alpha.reshape(-1, 40)

file = f"{rep}.{input_window_index}.susie.mu.txt"
mu = np.loadtxt(file)
mu = mu.reshape(-1, 40)

effect = np.sum(alpha * mu, axis=0)


dfP = pd.read_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_fastgxe/{prefix}.txt", sep=r"\s+")


score = np.dot(dfP.loc[:, "age":"alcohol_frequency"].values, effect)
print(dfP.loc[:, "age":"alcohol_frequency"].values.shape, effect.shape, score.shape)

dfP_out = dfP.loc[:, :"alcohol_frequency"].copy()
dfP_out["trait"] = dfP[f"trait{rep}"]
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

os.chdir(out_dir)

for i in range(5):
    dfP_out[dfP_out["score_group"] == i].to_csv(f"{rep}.{input_window_index}.score_group_{i+1}.txt", sep="\t", index=False)
