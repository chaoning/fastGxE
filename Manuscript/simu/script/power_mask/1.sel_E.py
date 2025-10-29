import pandas as pd
import numpy as np
import random
import os

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/pheno/")

dataE = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/envi.csv", 
                    nrows=1)
envi_lst = dataE.columns.tolist()[1:]
dfE = pd.DataFrame({
    0: range(len(envi_lst)),
    1: envi_lst
})



for rep in range(1, 101):
    file = f"power_snp_h2_gxe_0.08.num_envi_30.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000.rep_{rep}.power_snp_envi.txt"
    df = pd.read_csv(file, sep=r"\s+")
    index_lst = list(df.iloc[:, 0])
    for rmN in [1, 2, 5, 10, 20, 30]:
        rm_index_lst = random.sample(index_lst, rmN)
        dfE_filtered = dfE[~dfE[0].isin(rm_index_lst)].copy()
        out_file = f"../res_mask_fastGxE/envi.rm{rmN}.{rep}.txt"
        dfE_filtered.to_csv(out_file, sep=" ", header=False, index=False)

