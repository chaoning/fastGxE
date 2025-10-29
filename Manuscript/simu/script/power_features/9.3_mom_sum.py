import pandas as pd
import numpy as np
import os

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_baselineLD_mom/")

res_lst = []
for i in range(100):
    file = f"h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000.rep_{i+1}.var"
    arr = np.loadtxt(file)
    res_lst.append(arr / np.sum(arr))

df_res = pd.DataFrame(res_lst, columns=["h2_add", "h2_gxe", "h2_nxe", "h2_e"])
df_res.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_baselineLD_mom_sum/mom_sum.csv", index=False)
