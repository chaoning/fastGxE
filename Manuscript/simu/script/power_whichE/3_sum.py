import numpy as np
import pandas as pd
import os
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_fastgxe_testgxe_topE/")
prefix = "h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000"

quantile_lst = []
for topE in [1, 2, 5, 10, 20, 30, 40]:
    df_lst = []
    for rep in range(1, 101):
        try:
            in_file = f"{prefix}.rep_{rep}.top{topE}.csv"
            df = pd.read_csv(in_file)
            # Perform your analysis and summation here
            df_lst.append(df)
        except FileNotFoundError:
            print(f"File {in_file} not found, skipping.")
            continue
    df_combined = pd.concat(df_lst, axis=0)
    quantile_lst.append(df_combined.iloc[:, 0].quantile(0.05))

df_res = pd.DataFrame({
    "topE": [1, 2, 5, 10, 20, 30, 40],
    "quantile_0.05": quantile_lst
})

df_res.to_csv("../res_fastgxe_testgxe_topE_sum/quantile_results.csv", index=False)
