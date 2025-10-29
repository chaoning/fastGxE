import numpy as np
import pandas as pd
import sys
import glob
import re
import os

prefix = "h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000"
num_snp = 1_000_000

vf_bins = ["exonic", "UTR5", "UTR3", "intronic", "upstream", "downstream", "intergenic"]

p_cutoff_fastGxE_dct = {}
p_cutoff_fastGxE_noNxE_dct = {}
p_cutoff_StructLMM_dct = {}
p_cutoff_fastGWA_GE_dct = {}

for vf in vf_bins:
    file = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/QQplot_vf/variant_function.{vf}.txt.{prefix}_quantile.csv"
    df = pd.read_csv(file)
    p_cutoff_fastGxE_dct[vf] = df.iloc[1, 1] / num_snp
    p_cutoff_fastGxE_noNxE_dct[vf] = df.iloc[1, 2] / num_snp
    p_cutoff_StructLMM_dct[vf] = df.iloc[1, 3] / num_snp
    p_cutoff_fastGWA_GE_dct[vf] = df.iloc[1, 4] / num_snp


power_snp_h2_gxe = 0.08
num_envi=30
dfR = pd.DataFrame({"vf_bin": vf_bins})

# fastGxE
os.chdir(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_fastgxe_testgxe_feature/")
power_lst = []
power_rm_lst = []
for vf_bin in vf_bins:
    lst = []
    lst_rm = []
    for rep in range(1, 101):
        lead_snp_file = f"../pheno_feature/variant_function.{vf_bin}.power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{prefix}.rep_{rep}.power_snp.txt"
        lead_snp = pd.read_csv(lead_snp_file, sep=r"\s+", header=None).iloc[0, 0]

        pattern = re.compile(rf"variant_function.{vf_bin}\.power_snp_h2_gxe_{power_snp_h2_gxe}\.num_envi_{num_envi}\.{prefix}\.rep_{rep}\.[0-9]+_[0-9]+\.res")
        all_files = glob.glob(f"*")
        matched_files = [f for f in all_files if pattern.search(f)]
        if len(matched_files) != 1:
            continue
        file = matched_files[0]
        try:
            df = pd.read_csv(file, sep=r"\s+")
            p = df.iloc[:, -1].min()
            lst.append(p < p_cutoff_fastGxE_dct[vf_bin])
            df = df[df["SNP"] != lead_snp]
            p_rm = df.iloc[:, -1].min()
            lst_rm.append(p_rm < p_cutoff_fastGxE_dct[vf_bin])
        except Exception as e:
            print(f"Error reading {file}: {e}")
            continue
    print(f"vf_bin_{vf_bin}: {len(lst)}")
    power_lst.append(np.sum(lst) / len(lst))
    power_rm_lst.append(np.sum(lst_rm) / len(lst_rm))

dfR["fastGxE"] = power_lst
dfR["fastGxE_rm"] = power_rm_lst

# StructLMM
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_structlmm_feature/")
power_lst = []
power_rm_lst = []
for vf_bin in vf_bins:
    lst = []
    lst_rm = []
    for rep in range(1, 101):
        lead_snp_file = f"../pheno_feature/variant_function.{vf_bin}.power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{prefix}.rep_{rep}.power_snp.txt"
        lead_snp = pd.read_csv(lead_snp_file, sep=r"\s+", header=None).iloc[0, 0]

        file = f"variant_function.{vf_bin}.power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{prefix}.{rep}.csv"
        try:
            df = pd.read_csv(file)
            p = df.iloc[:, -1].min()
            lst.append(p < p_cutoff_StructLMM_dct[vf_bin])
            df = df[df["SNP"] != lead_snp]
            p_rm = df.iloc[:, -1].min()
            lst_rm.append(p_rm < p_cutoff_StructLMM_dct[vf_bin])
        except Exception as e:
            print(f"Error reading {file}: {e}")
            continue
    print(f"vf_bin_{vf_bin}: {len(lst)}")
    power_lst.append(np.sum(lst) / len(lst))
    power_rm_lst.append(np.sum(lst_rm) / len(lst_rm))

dfR["StructLMM"] = power_lst
dfR["StructLMM_rm"] = power_rm_lst

# fastGWA_GE
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_fastGWA_feature/")
power_lst = []
power_rm_lst = []
for vf_bin in vf_bins:
    lst = []
    lst_rm = []
    for rep in range(1, 101):
        lead_snp_file = f"../pheno_feature/variant_function.{vf_bin}.power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{prefix}.rep_{rep}.power_snp.txt"
        lead_snp = pd.read_csv(lead_snp_file, sep=r"\s+", header=None).iloc[0, 0]

        p_lst = []
        p_lst_rm = []
        for envi in range(1, 41):
            file = f"variant_function.{vf_bin}.power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{prefix}.{envi}.{rep}.fastGWA.res"
            try:
                df = pd.read_csv(file, sep=r"\s+")
                p = df.iloc[:, -1].min()
                p_lst.append(p)
                df = df[df["SNP"] != lead_snp]
                p_rm = df.iloc[:, -1].min()
                p_lst_rm.append(p_rm)
            except Exception as e:
                print(f"Error reading {file}: {e}")
                continue
        if len(p_lst) == 40:
            lst.append(np.min(p_lst) * 40 < p_cutoff_fastGWA_GE_dct[vf_bin])
            lst_rm.append(np.min(p_lst_rm) * 40 < p_cutoff_fastGWA_GE_dct[vf_bin])
    print(f"vf_bin_{vf_bin}: {len(lst)}")
    power_lst.append(np.sum(lst) / len(lst))
    power_rm_lst.append(np.sum(lst_rm) / len(lst_rm))

dfR["fastGWA_GE"] = power_lst
dfR["fastGWA_GE_rm"] = power_rm_lst

dfR.to_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/vf.{prefix}.csv", index=False)
