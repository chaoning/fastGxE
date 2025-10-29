import pandas as pd
import numpy as np
import os
from scipy.stats import chi2
import glob
import re
import sys


if len(sys.argv) == 1:
    os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/validate_AIE_test_main/")
else:
    os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/validate_AIE_test_main_random/")


arr = np.array([[1, -1, 0, 0, 0],
                [0, 1, -1, 0, 0],
                [0, 0, 1, -1, 0],
                [0, 0, 0, 1, -1]])

p_lst = []
for rep in range(1, 11):
    for input_window_index in range(1, 139):
        beta_lst = []
        se_lst = []
        for group in range(1, 6):
            file = f"{rep}.{input_window_index}.score_group_{group}.[0-9]+_[0-9]+.res"
            # Glob to list all candidate files
            pattern = f"{rep}.{input_window_index}.score_group_{group}.*_*.res"
            file_list = glob.glob(pattern)
            regex = re.compile(rf"^{rep}\.{input_window_index}\.score_group_{group}\.\d+_\d+\.res$")
            for file in file_list:
                if regex.match(os.path.basename(file)):  # Ensure exact pattern match
                    print(file)
                    df = pd.read_csv(file, sep=r"\s+")
                    beta_lst.append(df.iloc[0, -3])
                    se_lst.append(df.iloc[0, -2])
        if len(beta_lst) == 5 and len(se_lst) == 5:
            beta = np.array(beta_lst)
            se = np.array(se_lst)
            chi2_stat = arr @ beta @ (np.linalg.inv(arr @ np.diag(se**2) @ arr.T) @ arr @ beta)
            p = chi2.sf(chi2_stat, 4)
            p_lst.append([rep, input_window_index, chi2_stat, p])

p_df = pd.DataFrame(p_lst, columns=["rep", "input_window_index", "chi2_stat", "p"])

if len(sys.argv) == 1:
    p_df.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/10.2_stratified_test.csv", index=False)
else:
    p_df.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/10.2_stratified_test_random.csv", index=False)
