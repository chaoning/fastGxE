import numpy as np
import pandas as pd
import sys
import glob
import re
import os

# prefix = "h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000"
prefix = sys.argv[1]
num_snp = 1_000_000

file = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/QQ_plot/{prefix}_quantile.csv"
df = pd.read_csv(file)
p_cutoff_fastGxE = df.iloc[1, 1] / num_snp
p_cutoff_fastGxE_noNxE = df.iloc[1, 2] / num_snp
p_cutoff_StructLMM = df.iloc[1, 3] / num_snp
p_cutoff_fastGWA_GE = df.iloc[1, 4] / num_snp

num_envi = 30

power_lst = []
power_snp_h2_gxe_lst = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14]

# fastGxE
os.chdir(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_fastgxe_testgxe/")
for power_snp_h2_gxe in power_snp_h2_gxe_lst:
    lst = []
    for rep in range(1, 101):
        pattern = re.compile(rf"power_snp_h2_gxe_{power_snp_h2_gxe}\.num_envi_{num_envi}\.{prefix}\.rep_{rep}\.[0-9]+_[0-9]+\.res")
        all_files = glob.glob(f"*")
        matched_files = [f for f in all_files if pattern.search(f)]
        if len(matched_files) != 1:
            continue
        file = matched_files[0]
        try:
            df = pd.read_csv(file, sep=r"\s+")
            p = df.iloc[0, -1]
            lst.append(p < p_cutoff_fastGxE)
        except Exception as e:
            print(f"Error reading {file}: {e}")
            continue
    print(f"power_snp_h2_gxe_{power_snp_h2_gxe}: {len(lst)}")
    power_lst.append(["fastGxE", power_snp_h2_gxe, np.sum(lst) / len(lst)])


# fastGxE_noNxE
os.chdir(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_fastgxe_testgxe/")
for power_snp_h2_gxe in power_snp_h2_gxe_lst:
    lst = []
    for rep in range(1, 101):
        pattern = re.compile(rf"power_snp_h2_gxe_{power_snp_h2_gxe}\.num_envi_{num_envi}\.{prefix}\.rep_{rep}\.noNxE\.[0-9]+_[0-9]+\.res")
        all_files = glob.glob(f"*")
        matched_files = [f for f in all_files if pattern.search(f)]
        if len(matched_files) != 1:
            continue
        file = matched_files[0]
        try:
            df = pd.read_csv(file, sep=r"\s+")
            p = df.iloc[0, -1]
            lst.append(p < p_cutoff_fastGxE_noNxE)
        except Exception as e:
            print(f"Error reading {file}: {e}")
            continue
    print(f"power_snp_h2_gxe_{power_snp_h2_gxe}: {len(lst)}")
    power_lst.append(["fastGxE-noNxE", power_snp_h2_gxe, np.sum(lst) / len(lst)])


# StructLMM
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_structlmm/")
for power_snp_h2_gxe in power_snp_h2_gxe_lst:
    lst = []
    for rep in range(1, 101):
        file = f"power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{prefix}.{rep}.csv"
        try:
            df = pd.read_csv(file)
            p = df.iloc[0, -1]
            lst.append(p < p_cutoff_StructLMM)
        except Exception as e:
            print(f"Error reading {file}: {e}")
            continue
    print(f"power_snp_h2_gxe_{power_snp_h2_gxe}: {len(lst)}")
    power_lst.append(["StructLMM", power_snp_h2_gxe, np.sum(lst) / len(lst)])


# fastGWA_GE
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_fastGWA_h2gxe/")
for power_snp_h2_gxe in power_snp_h2_gxe_lst:
    lst = []
    for rep in range(1, 101):
        p_lst = []
        for envi in range(1, 41):
            file = f"power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{prefix}.{envi}.{rep}.fastGWA.res"
            try:
                df = pd.read_csv(file, sep=r"\s+")
                p = df.iloc[0, -1]
                p_lst.append(p)
            except Exception as e:
                print(f"Error reading {file}: {e}")
                continue
        if len(p_lst) == 40:
            lst.append(np.min(p_lst) * 40 < p_cutoff_fastGWA_GE)
    print(f"power_snp_h2_gxe_{power_snp_h2_gxe}: {len(lst)}")
    power_lst.append(["fastGWA-GE", power_snp_h2_gxe, np.sum(lst) / len(lst)])

dfR = pd.DataFrame(power_lst, columns=["method", "snp_h2_gxe", "power"])

dfR.to_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/snp_h2_gxe.{prefix}.csv", index=False)
