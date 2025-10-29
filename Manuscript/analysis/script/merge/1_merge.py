import pandas as pd
import os
import sys
from tqdm import tqdm

trait = sys.argv[1]
in_dir = sys.argv[2]
out_dir = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/merge/"
os.chdir(in_dir)
df_lst = []
# col_lst = ["chrom", "SNP", "base", "allele1", "allele2", "af", "missing", "p_main", "p_gxe"]
for i in tqdm(range(10)):
    file = f"{trait}.10_{i+1}.res"
    df = pd.read_csv(file, sep=r"\s+")
    df_lst.append(df)

df_res = pd.concat(df_lst, axis=0, ignore_index=True)

os.chdir(out_dir)
df_res.to_csv(f"{trait}.gz", 
           sep="\t", index=False, na_rep="NA", compression="gzip")
