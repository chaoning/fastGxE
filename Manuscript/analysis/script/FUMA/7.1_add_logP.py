import pandas as pd
import sys
import numpy as np
from tqdm import tqdm
import os

in_file = sys.argv[1]
test_dir = sys.argv[2]
os.chdir(test_dir)

dfGRL = pd.read_csv(f"{in_file}.csv", dtype=str)

p_gxe_log_lst = []
for trait, snp in tqdm(zip(dfGRL["trait"], dfGRL["trait_leading_snp"])):
    if trait == "WHRadjBMI" or trait == "78":
        p_gxe_log_lst.append(np.nan)
        continue
    df_lst = []
    for i in range(10):
        file = f"{trait}.10_{i+1}.res"
        df = pd.read_csv(file, sep=r"\s+", usecols=["SNP", "p_gxe"])
        df_lst.append(df)
    df = pd.concat(df_lst, ignore_index=True)
    p_gxe_log_lst.append(df[df["SNP"] == snp].iloc[0, -1])

dfGRL["p_gxe_log"] = p_gxe_log_lst
dfGRL.to_csv(f"{in_file}_logP.csv", index=False, na_rep="NA")
