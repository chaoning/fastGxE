import numpy as np
import os
import re
import glob
import pandas as pd

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_rGE_fastgxe/")

power_snp_h2_add = 0.2
prefix="h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000"

dct = {}
for rGE in [0, 0.001, 0.01, 0.1, 0.5]:
    p_lst = []
    for rep in range(1, 101):
        pattern = re.compile(rf"power_snp_h2_add_{power_snp_h2_add}\.rGE_{rGE}\.{prefix}\.rep_{rep}\.[0-9]+_[0-9]+\.res")
        all_files = glob.glob(f"*")
        matched_files = [f for f in all_files if pattern.search(f)]
        if len(matched_files) != 1:
            p_lst.append(np.nan)
            continue
        file = matched_files[0]
        try:
            df = pd.read_csv(file, sep=r"\s+")
            p = df.iloc[0, -1]
            p_lst.append(p)
        except Exception as e:
            print(f"Error reading {file}: {e}")
            p_lst.append(np.nan)
            continue
    dct[f"rGE_{rGE}"] = p_lst

dfR = pd.DataFrame(dct)

dfR.to_csv(f"../res_rGE_sum/power_snp_h2_add_{power_snp_h2_add}.{prefix}.csv", index=False)
