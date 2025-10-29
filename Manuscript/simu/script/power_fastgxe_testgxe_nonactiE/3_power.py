import numpy as np
import pandas as pd
import sys
import glob
import re
import os

# prefix = "h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000"
num_snp = 1_000_000

file = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/QQ_plot/h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000_quantile.csv"
df = pd.read_csv(file)
p_cutoff_fastGxE = df.iloc[1, 1] / num_snp
p_cutoff_fastGxE_noNxE = df.iloc[1, 2] / num_snp
p_cutoff_StructLMM = df.iloc[1, 3] / num_snp
p_cutoff_fastGWA_GE = df.iloc[1, 4] / num_snp



h2_add = 30
h2_gxe = 5
h2_nxe = 15
num_envi_power = 30

reps = range(1, 101)

res_lst = []
os.chdir(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_fastgxe_testgxe_nonactiE/")
for nonE in [0, 10, 20, 30]:
    lst = []
    for rep in reps:
        pattern = re.compile(rf"nonE_{nonE}\.rep_{rep}\.[0-9]+_[0-9]+\.res")
        all_files = glob.glob(f"*")
        matched_files = [f for f in all_files if pattern.search(f)]
        if len(matched_files) != 1:
            continue
        file = matched_files[0]
        try:
            df = pd.read_csv(file, sep=r"\s+")
            p = df.iloc[:, -1].min()
            lst.append(p < p_cutoff_fastGxE)
        except Exception as e:
            print(f"Error reading {file}: {e}")
            continue
    print(f"nonE {nonE}: {len(lst)}")
    power_val = np.sum(lst) / len(lst)
    res_lst.append([nonE, power_val])

dfR = pd.DataFrame(res_lst, columns=["nonE", "power_val"])

dfR.to_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/nonE.csv", index=False)
print(dfR)
