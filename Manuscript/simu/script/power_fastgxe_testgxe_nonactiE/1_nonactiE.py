import pandas as pd
import numpy as np
import os
import random


power_snp_h2_gxe=0.08
num_envi=10
prefix="h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000"


df = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/envi.csv", nrows=1)
envi40_set = set(df.loc[:, "age":"alcohol_frequency"].columns)
print(envi40_set)
print(len(envi40_set))

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/pheno/")
for rep in range(1, 101):
    envi_file = f"power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{prefix}.rep_{rep}.power_snp_envi.txt"
    envi_lst = pd.read_csv(envi_file, sep=r"\s+").iloc[:, 1].tolist()
    envi_nonactiE_lst = list(envi40_set - set(envi_lst))
    for n in [0, 10, 20, 30]:
        rndE_lst = set(random.sample(envi_nonactiE_lst, n)) | set(envi_lst)
        rndE_lst = list(rndE_lst)
        with open(f"../res_fastgxe_testgxe_nonactiE/power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{prefix}.rep_{rep}.rndE_{n}.txt", "w") as fout:
            for envi in rndE_lst:
                fout.write(f"{envi}\n")

