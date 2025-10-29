import numpy as np
import os
import re
import glob
import pandas as pd

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_addToint_fastgxe/")

power_snp_h2_add = 0.2
prefix="h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000"

dct = {}
for unknown_h2SNP_h2E in [0.02, 0.08, 0.16]:
        for R2 in [0, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9]:
            for R2GE in [0]:
                p_lst = []
                for rep in range(1, 101):
                    file = f"h2_{unknown_h2SNP_h2E}.R2_{R2}.R2GE_{R2GE}.{prefix}.rep_{rep}.res"
                    try:
                        df = pd.read_csv(file, sep=r"\s+")
                        p = df.iloc[0, -1]
                        p_lst.append(p)
                    except Exception as e:
                        print(f"Error reading {file}: {e}")
                        p_lst.append(np.nan)
                        continue
                dct[f"h2_{unknown_h2SNP_h2E}.R2_{R2}.R2GE_{R2GE}"] = p_lst
        for R2 in [0.5]:
            for R2GE in [0.01, 0.1, 0.3]:
                p_lst = []
                for rep in range(1, 101):
                    file = f"h2_{unknown_h2SNP_h2E}.R2_{R2}.R2GE_{R2GE}.{prefix}.rep_{rep}.res"
                    try:
                        df = pd.read_csv(file, sep=r"\s+")
                        p = df.iloc[0, -1]
                        p_lst.append(p)
                    except Exception as e:
                        print(f"Error reading {file}: {e}")
                        p_lst.append(np.nan)
                        continue
                dct[f"h2_{unknown_h2SNP_h2E}.R2_{R2}.R2GE_{R2GE}"] = p_lst

dfR = pd.DataFrame(dct)

dfR.to_csv(f"../res_addToint_fastgxe_sum/sum.csv", index=False)
