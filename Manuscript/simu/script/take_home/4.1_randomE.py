import pandas as pd
import numpy as np
import os
import random




df = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/envi.csv", nrows=1)
envi40_lst = df.loc[:, "age":"alcohol_frequency"].columns.to_list()

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_fastgxe_testgxe_nEtested/")
for rep in range(1, 101):
    for n in [1, 2, 5, 10, 20, 30, 40]:
        rndE_lst = random.sample(envi40_lst, n)
        with open(f"E.rep_{rep}.rnd_{n}.txt", "w") as fout:
            for envi in rndE_lst:
                fout.write(f"{envi}\n")
