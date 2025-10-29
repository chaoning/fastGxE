import numpy as np
import pandas as pd
import os
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_take_home_whichE_sum/")
res_lst = []
for power_snp_h2_gxe in [0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14]:
	for sample_size in [50000, 100000, 200000, 300000, 400000]:
		for num_envi_power in [1, 2, 5, 10, 20, 30, 40]:
			for topE in [1, 2, 5, 10, 20, 30, 40]:
				# Read the results for each combination
				df = pd.read_csv(f"../res_take_home_whichE/power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi_power}.h2_add_30.h2_gxe_5.h2_nxe_15.sample_{sample_size}.random_{topE}.csv")
				df.dropna(inplace=True)
				power_val = df[df["p_gxe"] < 5e-8].shape[0] / df.shape[0]
				res_lst.append((power_snp_h2_gxe, sample_size, num_envi_power, topE, power_val))

df_res = pd.DataFrame(res_lst, columns=["power_snp_h2_gxe", "sample_size", "num_envi_power", "topE", "power"])
df_res.to_csv("summary_results_random.csv", index=False)
