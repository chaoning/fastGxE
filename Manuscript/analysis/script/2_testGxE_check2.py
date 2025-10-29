import pandas as pd
import numpy as np
import os
import glob
import re
from tqdm import tqdm
import sys
from concurrent.futures import ProcessPoolExecutor
import sys
from pandas.api.types import CategoricalDtype



os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/testGxE/")

df = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/trait32.txt",
                 sep="\t")

trait_lst = list(df.iloc[:, 0])
print("The number of traits: {}".format(len(trait_lst)))


def process_trait(trait):
    lst = []
    for i in range(10):
        try:
            df = pd.read_csv(f"{trait}.10_{i+1}.res", sep=r"\s+")
            lst.append(df.shape[0])
        except Exception as e:
            print(f"Error processing {trait}.10_{i+1}.res: {e}")
            lst.append(0)
    res_df = pd.DataFrame({
        "trait": [trait] * len(lst),
        "replicate": list(range(1, 11)),
        "n": lst
    })
    return  res_df


def parallel_process_trait_list(trait_list):
    with ProcessPoolExecutor(max_workers=30) as executor:
        results = executor.map(process_trait, trait_list)
        res_lst = list(results)
    return res_lst


res_lst = parallel_process_trait_list(trait_lst)

df_res = pd.concat(res_lst, ignore_index=True)


trait_order = CategoricalDtype(categories=trait_lst, ordered=True)


df_res["trait"] = df_res["trait"].astype(trait_order)
df_res = df_res.sort_values(["trait", "replicate"]).reset_index(drop=True)

df_res.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/physical_trait_replicate_n.txt", sep="\t", index=False)
