import pandas as pd
import numpy as np
import os
import sys
base_dir = "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power"
out_dir = f"{base_dir}/res_filter_sum"
os.chdir(out_dir)
num_envi_power = sys.argv[1]

df_lst = []
for rGI in [-0.9, -0.3, 0, 0.3, 0.9]:
	for add_gxe_ratio in [0.006, 0.02, 0.25, 1, 2, 4]:
		file = f"rGI_{rGI}.add_gxe_ratio_{add_gxe_ratio}.num_envi_{num_envi_power}.power_summary.csv"
		df = pd.read_csv(file)
		df_lst.append(df)

df = pd.concat(df_lst, ignore_index=True)

df.to_csv(f"num_envi_{num_envi_power}.power_summary_all.csv", index=False)
