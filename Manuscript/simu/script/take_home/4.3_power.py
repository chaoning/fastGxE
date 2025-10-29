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



os.chdir(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_fastgxe_testgxe_nEtested/")

rndE_lst = [1, 2, 5, 10, 20, 30, 40]
reps = range(1, 101)
power_snp_h2_gxe = 0.08
sample_sizes = [50000, 100000, 200000, 300000, 400000]
num_envi_power_lst = [2, 30]

res_lst = []
res_lst_rm = []
for sample_size in sample_sizes:
    for num_envi_power in num_envi_power_lst:
        for rndE in rndE_lst:
            lst = []
            lst_rm = []
            for rep in reps:
                lead_snp_file = f"../pheno_take_home/power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi_power}.h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_{sample_size}.rep_{rep}.power_snp_name.txt"

                pattern = re.compile(rf"rndE{rndE}\.power_snp_h2_gxe_{power_snp_h2_gxe}\.num_envi_{num_envi_power}\.h2_add_{h2_add}\.h2_gxe_{h2_gxe}\.h2_nxe_{h2_nxe}\.sample_{sample_size}\.rep_{rep}\.[0-9]+_[0-9]+\.res")
                all_files = glob.glob(f"*")
                matched_files = [f for f in all_files if pattern.search(f)]
                if len(matched_files) != 1:
                    continue
                file = matched_files[0]
                try:
                    lead_snp = pd.read_csv(lead_snp_file, sep=r"\s+", header=None).iloc[0, 0]
                    df = pd.read_csv(file, sep=r"\s+")
                    p = df.iloc[:, -1].min()
                    lst.append(p < p_cutoff_fastGxE)

                    df = df[df["SNP"] != lead_snp]
                    p_rm = df.iloc[:, -1].min()
                    lst_rm.append(p_rm < p_cutoff_fastGxE)
                except Exception as e:
                    print(f"Error reading {file}: {e}")
                    continue
            print(f"sample size {sample_size}, power_snp_h2_gxe {power_snp_h2_gxe}, rndE {rndE}, num_envi_power {num_envi_power}: {len(lst)}")
            power_val = np.sum(lst) / len(lst)
            power_val_rm = np.sum(lst_rm) / len(lst_rm)
            res_lst.append([sample_size, num_envi_power, rndE, power_val])
            res_lst_rm.append([sample_size, num_envi_power, rndE, power_val_rm])

dfR = pd.DataFrame(res_lst, columns=["sample_size", "num_envi_power", "rndE", "power_val"])
dfR_rm = pd.DataFrame(res_lst_rm, columns=["sample_size", "num_envi_power", "rndE", "power_val_rm"])

dfR.to_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/take_home_rndE.csv", index=False)
dfR_rm.to_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/take_home_rndE_rm.csv", index=False)
