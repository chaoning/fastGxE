import numpy as np
import pandas as pd
from scipy.stats import norm
import re
from scipy.stats import chi2
from tqdm import tqdm
import sys
import os

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/testGxE/")
trait = sys.argv[1]

df_lst = []

for i in range(10):
    file = f"{trait}.10_{i+1}.res"
    df = pd.read_csv(file, header=0, sep=r"\s+", usecols=["chrom", "SNP", "cm", "base", "p_main", "p_gxe"])
    df_lst.append(df)

df = pd.concat(df_lst, ignore_index=True)

df.to_csv(f"{trait}.res", index=False, sep="\t", na_rep="NA")

