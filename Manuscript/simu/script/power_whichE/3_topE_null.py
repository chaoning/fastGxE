import pandas as pd
import numpy as np
import sys
sys.path.append("/net/zootopia/disk1/chaon/WORK/structLMM/")
from saddle_point import pchisqsum
from acat import acat
import os
from tqdm import tqdm

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_fastgxe_testgxe/")

prefix = "h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000"
rep = sys.argv[1]
topE = int(sys.argv[2])  # 1, 2, 5, 10, 20, 30, 40
in_file = f"{prefix}.rep_{rep}.res"
out_file = f"../res_fastgxe_testgxe_topE/{prefix}.rep_{rep}.top{topE}.csv"

beta_cols = [f"beta{i+1}" for i in range(40)]
se_cols = [f"se{i+1}" for i in range(40)]
p_cols = [f"p{i+1}" for i in range(40)]

dfP = pd.read_csv(f"../pheno_fastgxe/{prefix}.txt", sep=r"\s+")
envi_corr = dfP.loc[:, "age":"alcohol_frequency"].corr().values

df = pd.read_csv(in_file, sep=r"\s+")
p_combined_lst = []
for i in tqdm(range(df.shape[0])):
    row = df.iloc[i]
    df_row = pd.DataFrame({
        "order": range(40),
        "beta": row[beta_cols].values,
        "se": row[se_cols].values,
        "p": row[p_cols].values
    })
    df_row.sort_values(by="p", inplace=True)
    dfRes_selected = df_row.iloc[:topE, :]
    z_scores = dfRes_selected["beta"].values / dfRes_selected["se"].values
    scores = np.sum(z_scores ** 2)
    idx = dfRes_selected["order"].values
    envi_corr_sub = envi_corr[np.ix_(idx, idx)]
    a = np.linalg.eigvalsh(envi_corr_sub)
    p_sad = pchisqsum(np.array([scores]), [1]*len(a), a, lower_tail=True, method="saddlepoint")[0]
    p_acat = acat(dfRes_selected["p"].values)[1]
    p_combined = acat([p_sad, p_acat])[1]
    p_combined_lst.append(p_combined)

df_out = pd.DataFrame({
    "p_combined": p_combined_lst
})

df_out.to_csv(out_file, index=False)
