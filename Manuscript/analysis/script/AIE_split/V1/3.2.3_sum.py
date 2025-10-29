import os
import pandas as pd
import numpy as np
from scipy.stats import chi2

def process_trait_specific(summary_csv, out_csv, trait_type):
    dfS = pd.read_csv(summary_csv)
    trait_lst = dfS.iloc[:, 14].to_list()
    snp_lst = dfS.iloc[:, 15].to_list()

    os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/AIE_split/test_main/")
    df_lst = []
    for trait, snp in zip(trait_lst, snp_lst):
        file = f"{trait}.{snp}.txt"
        df = pd.read_csv(file, sep=r"\s+")
        df_lst.append(df)

    df_res = pd.concat(df_lst, ignore_index=True)
    df_res = pd.concat(
        [dfS.loc[:, ["GenomicLocus_hg38", "Genes", "trait", "trait_leading_snp"]],
         df_res],
        axis=1
    )
    print(f"[{trait_type}] p < 5e-8:", (df_res["p"] < 5e-8).sum())
    print(f"[{trait_type}] p < 0.05 / N:", (df_res["p"] < 0.05 / df_res.shape[0]).sum())
    print(f"[{trait_type}] p < 0.05:", (df_res["p"] < 0.05).sum())

    p_col = [f"p{i}" for i in range(1, 6)]
    df_res["p_main_cutoff1"] = (df_res[p_col] < (5e-8 / 5)).sum(axis=1)
    df_res["p_main_cutoff2"] = (df_res[p_col] < (0.05 / df_res.shape[0] / 5)).sum(axis=1)
    df_res["p_main_cutoff3"] = (df_res[p_col] < 0.05/5).sum(axis=1)

    df_res.to_csv(out_csv, index=False)

# Process PT
process_trait_specific(
    "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/Physical_measures/result/SGEM/FUMA/IINT/GenomicRiskLoci.hg38.keygene.TraitSpecific.csv",
    "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/AIE_split/3.2.3_sum.PT_trait_specific.csv",
    "PT"
)

# Process BB
process_trait_specific(
    "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/Biological_samples/result/SGEM/FUMA/IINT/GenomicRiskLoci.hg38.keygene.TraitSpecific.Methods.addP.csv",
    "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/AIE_split/3.2.3_sum.BB_trait_specific.csv",
    "BB"
)
