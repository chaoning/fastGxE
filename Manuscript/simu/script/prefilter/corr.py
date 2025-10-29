import os
import pandas as pd
import sys

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/pheno_filter/")

rGI = 0.05
df_lst= []
for i in range(100):
    file = f"rGI_{rGI}.power_snp_h2_gxe_0.08.num_envi_30.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000.rep_{i+1}.power_snp_eff.txt"
    df = pd.read_csv(file, sep=r"\s+", header=None)
    df_lst.append(df)

df = pd.concat(df_lst)
