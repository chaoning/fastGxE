import random
import pandas as pd
import numpy as np
import sys
import os
import re
import glob
from tqdm import tqdm
sys.path.append("/net/zootopia/disk1/chaon/WORK/structLMM/")
from saddle_point import pchisqsum
from acat import acat

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_fastgxe_testgxe_take_home/")
# prefix = "power_snp_h2_gxe_0.005.num_envi_10.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000"
prefix = sys.argv[1]
top_envi = int(sys.argv[2])  # 1, 5, 10
num_envi_power = int(sys.argv[3])  # 1, 5, 10
dfP = pd.read_csv(f"../pheno_take_home_fastgxe/{prefix}.txt", sep=r"\s+")
envi_corr = dfP.loc[:, "age":"alcohol_frequency"].corr().values
print(envi_corr)
q = envi_corr.shape[0]
beta_cols = [f"beta{i+1}" for i in range(q)]
se_cols = [f"se{i+1}" for i in range(q)]
p_cols = [f"p{i+1}" for i in range(q)]

p_combined_lst = []
for rep in tqdm(range(1, 101)):
    random.seed(rep)
    try:
        powerSNP = pd.read_csv(f"../pheno_take_home/{prefix}.rep_{rep}.power_snp_name.txt", sep=r"\s+", header=None).iloc[0, 0]
    except Exception as e:
        p_combined_lst.append(np.nan)
        continue
    try:
        power_envi_idx = pd.read_csv(f"../pheno_take_home/{prefix}.rep_{rep}.power_snp_envi.txt", sep=r"\s+").iloc[:, 0].values.tolist()
    except Exception as e:
        p_combined_lst.append(np.nan)
        continue
    pattern = re.compile(rf"{prefix}\.rep_{rep}\.[0-9]+_[0-9]+\.res")
    all_files = glob.glob(f"*")
    matched_files = [f for f in all_files if pattern.search(f)]
    if len(matched_files) != 1:
        p_combined_lst.append(np.nan)
        continue
    print(matched_files[0])
    dfRes = pd.read_csv(matched_files[0], sep=r"\s+")
    row = dfRes.loc[dfRes["SNP"] == powerSNP]
    try:
        row = row.iloc[0]
    except Exception as e:
        p_combined_lst.append(np.nan)
        continue
    beta_vals = row[beta_cols].to_numpy(dtype=float, copy=False)
    se_vals   = row[se_cols].to_numpy(dtype=float, copy=False)
    p_vals    = row[p_cols].to_numpy(dtype=float, copy=False)
    dfRes = pd.DataFrame({
        "order": range(q),
        "beta": beta_vals,
        "se": se_vals,
        "p": p_vals
    })
    if top_envi <= num_envi_power:
        idx = random.sample(power_envi_idx, top_envi)
    else:
        idx = power_envi_idx + random.sample(list(set(range(q)) - set(power_envi_idx)), top_envi - num_envi_power)
    idx = sorted(idx)
    dfRes_selected = dfRes.iloc[idx, :]
    z_scores = dfRes_selected["beta"] / dfRes_selected["se"]
    scores = np.sum(z_scores ** 2)
    idx = dfRes_selected["order"].values
    envi_corr_sub = envi_corr[np.ix_(idx, idx)]
    a = np.linalg.eigvalsh(envi_corr_sub)
    p_sad = pchisqsum(np.array([scores]), [1]*len(a), a, lower_tail=True, method="saddlepoint")[0]
    p_acat = acat(dfRes_selected["p"].values)[1]
    p_combined = acat([p_sad, p_acat])[1]
    p_combined_lst.append(p_combined)

dfR = pd.DataFrame({
    "rep": range(1, 101),
    "p_gxe": p_combined_lst
})
dfR.to_csv(f"../res_take_home_whichE/{prefix}.random_{top_envi}.csv", index=False)
