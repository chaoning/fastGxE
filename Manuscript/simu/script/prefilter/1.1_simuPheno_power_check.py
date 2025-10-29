import pandas as pd
import numpy as np
import os
import re

pattern = re.compile(
    r"(?:INFO\s+-\s+)?\[After normalization\]\s+Additive Var=([0-9.eE+-]+)\s+"
    r"GxE Var=([0-9.eE+-]+)\s+"
    r"r\(Add,\s*sum\(beta_GxE·G\)\)=([0-9.eE+-]+)"
)


os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/prefilter/Aerr/simu_power/")
k = 0
res_lst = []
for rGI in [-0.9, -0.3, 0, 0.3, 0.9]:
    for add_gxe_ratio in [0.005, 0.25, 1, 2, 4]:
        additive_sum = 0
        gxe_sum = 0
        correlation_sum = 0
        for rep in range(1, 101):
            k += 1
            file = f"e2.{k}.txt"
            with open(file, "r") as f:
                for line in f:
                    match = pattern.search(line)
                    if match:
                        additive = float(match.group(1))
                        gxe = float(match.group(2))
                        correlation = float(match.group(3))
                        additive_sum += additive
                        gxe_sum += gxe
                        correlation_sum += correlation
        additive_mean = additive_sum / 100
        gxe_mean = gxe_sum / 100
        correlation_mean = correlation_sum / 100
        res_lst.append([rGI, add_gxe_ratio, additive_mean, gxe_mean, correlation_mean])

df_res = pd.DataFrame(res_lst, columns=["rGI", "Add_GxE_Ratio", "Additive", "GxE", "Correlation"])
df_res.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/prefilter/1.1_simuPheno_power_check.e2.csv", index=False)
