import pandas as pd
from tqdm import tqdm

file = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/trait32.txt"
df = pd.read_csv(file, sep="\t")
df = df[(df["FieldID"] != "78") & (df["FieldID"] != "WHRadjBMI")].copy()

df_lst = []
for trait in tqdm(df.iloc[:, 0], desc="Processing traits"):
    print(trait)
    file = f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/log_trans/vs/{trait}.filtered.csv"
    df_trait = pd.read_csv(file)
    if df_trait.shape[0] == 0:
        continue
    df_trait["trait"] = trait
    df_lst.append(df_trait)

df = pd.concat(df_lst, axis=0)

df.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/testGxE_log/5_sum_signal.txt", 
          index=False, sep="\t")
