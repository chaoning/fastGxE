import numpy as np
import pandas as pd
import sys
import glob
import re
import os
from tqdm import tqdm

# prefix = "h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000"
p_cutoff_fastGxE = 5.0e-8


h2_add = 30
h2_gxe = 5
h2_nxe = 15



os.chdir(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_take_home_extend/")

power_snp_h2_gxe_values = [0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14]
sample_sizes = [50000, 100000, 200000, 300000, 400000, 800000, 1200000, 1600000, 2000000]
num_envi_power_lst = [1, 2, 5, 10, 20, 30, 40]

res_lst = []
for sample_size in sample_sizes:
    for power_snp_h2_gxe in power_snp_h2_gxe_values:
        for num_envi_power in num_envi_power_lst:
            file = f"{sample_size}.h2_{power_snp_h2_gxe}.E_{num_envi_power}.csv"
            df_tmp = pd.read_csv(file)
            power_val = np.sum(df_tmp["p_gxe"] < p_cutoff_fastGxE) / len(df_tmp)
            res_lst.append([sample_size, power_snp_h2_gxe, num_envi_power, power_val])

dfR = pd.DataFrame(res_lst, columns=["sample_size", "power_snp_h2_gxe", "num_envi_power", "power_val"])

dfR.to_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/take_home_extend.csv", index=False)

