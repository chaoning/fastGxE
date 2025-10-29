import pandas as pd
import numpy as np
import os
import glob
import re
from tqdm import tqdm
import sys
from concurrent.futures import ProcessPoolExecutor
import sys


os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/testGxE/")

df = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/trait.txt",
                 sep="\t")

trait_lst = list(df.iloc[:, 0])
print("The number of traits: {}".format(len(trait_lst)))


def process_trait(trait):
    df_lst = []
    for i in range(10):
        df = pd.read_csv(f"{trait}.10_{i+1}.res", sep=r"\s+")
        df_lst.append(df)
    df = pd.concat(df_lst, ignore_index=True)
    df["trait"] = trait
    print(trait, df.shape[0])
    df = df[df["p_gxe"] < 5e-8].copy()
    return df


def parallel_process_trait_list(trait_list):
    with ProcessPoolExecutor(max_workers=36) as executor:
        results = executor.map(process_trait, trait_list)
        df_lst = list(results)
    return df_lst


df_lst = parallel_process_trait_list(trait_lst)
df = pd.concat(df_lst)
df.to_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/results/signal/signal.txt",
          sep="\t", index=False, na_rep="NA")
