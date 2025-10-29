import sys
sys.path.append('/net/zootopia/disk1/chaon/WORK/structLMM/')
from qqplot import qqplot, qqplot_pure, qqplot_frame

import pandas as pd
import numpy as np
from tqdm import tqdm
import os

prefix = sys.argv[1]

dfSNP = pd.read_csv(sys.argv[2], sep=r"\s+", usecols=["SNP"])

# fastGxE
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_fastgxe_testgxe/")
p_fastgxe_arr = np.array([])
p_fastgxe_noNxE_arr = np.array([])
for rep in tqdm(range(100)):
    try:
        file_path = f"{prefix}.rep_{rep + 1}.res"
        df = pd.read_csv(file_path, sep=r"\s+")
        df.dropna(inplace=True)
        df = pd.merge(dfSNP, df, on="SNP")
        p_fastgxe_arr = np.append(p_fastgxe_arr, df["p_gxe"])
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        continue

    try:
        file_path = f"{prefix}.rep_{rep + 1}.noNxE.res"
        df = pd.read_csv(file_path, sep=r"\s+")
        df.dropna(inplace=True)
        df = pd.merge(dfSNP, df, on="SNP")
        p_fastgxe_noNxE_arr = np.append(p_fastgxe_noNxE_arr, df["p_gxe"])
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        continue

# structlmm
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_structlmm/")
p_structlmm_arr = np.array([])
for rep in tqdm(range(100)):
    try:
        file_path = f"{prefix}.{rep + 1}.csv"
        df = pd.read_csv(file_path)
        df.dropna(inplace=True)
        df = pd.merge(dfSNP, df, on="SNP")
        p_structlmm_arr = np.append(p_structlmm_arr, df["p_sad"])
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        continue


#fastGWA-GE
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_fastGWA/")
p_fastGWA_arr = np.array([])


for rep in tqdm(range(100)):
    dfM = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22.bim", sep=r"\s+", header=None,
                    usecols=[1], names=["SNP"])
    dfM = dfSNP.copy()
    for envi in range(40):
        try:
            file_path = f"{prefix}.{envi+1}.{rep + 1}.fastGWA.res.gz"
            df = pd.read_csv(file_path, sep=r"\s+", usecols=["SNP", "P_G_by_E"])
            df.columns = ["SNP", f"P{envi + 1}"]
            dfM = pd.merge(dfM, df, on="SNP", how="left")
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            continue
    print(dfM.columns)
    dfPmin = dfM.iloc[:, 1:].min(axis=1) * 40
    dfPmin.dropna(inplace=True)
    p_fastGWA_arr = np.append(p_fastGWA_arr, dfPmin.values)

p_fastGWA_arr = np.sort(p_fastGWA_arr)[1:]

p_expected_arr = np.array([0.5, 0.05, 0.01, 0.001, 0.0001, 0.00001])
p_observed = np.zeros((len(p_expected_arr), 4))
for i in range(len(p_expected_arr)):
    p = p_expected_arr[i]
    p_observed[i, 0] = np.quantile(p_fastgxe_arr, p)
    p_observed[i, 1] = np.quantile(p_fastgxe_noNxE_arr, p)
    p_observed[i, 2] = np.quantile(p_structlmm_arr, p)
    p_observed[i, 3] = np.quantile(p_fastGWA_arr, p)

dfR = pd.DataFrame({
    "Expected": p_expected_arr,
    "fastGxE": np.log10(p_observed[:, 0]) / np.log10(p_expected_arr),
    "fastGxE-noNxE": np.log10(p_observed[:, 1]) / np.log10(p_expected_arr),
    "StructLMM": np.log10(p_observed[:, 2]) / np.log10(p_expected_arr),
    "fastGWA-GE": np.log10(p_observed[:, 3]) / np.log10(p_expected_arr)
})


p_arr_lst = [p_fastgxe_arr, p_structlmm_arr, p_fastGWA_arr]
color_lst = ["#1e466e", "#06948E",  "#e76254"]
label_lst = ["fastGxE", "structLMM", "fastGWA-GE"]
out_file = f"{sys.argv[2]}.{prefix}"
dfR.to_csv(out_file + ".csv", index=False)
qqplot(p_arr_lst, color_lst, label_lst, out_file, output_format="png")
qqplot_pure(p_arr_lst, color_lst, label_lst, out_file + "_pure", output_format="png")
qqplot_frame(p_arr_lst, color_lst, label_lst, out_file + "_frame", output_format="pdf")


dfR = pd.DataFrame({
    "Expected": p_expected_arr,
    "fastGxE": p_observed[:, 0],
    "fastGxE-noNxE": p_observed[:, 1],
    "StructLMM": p_observed[:, 2],
    "fastGWA-GE": p_observed[:, 3]
})

dfR.to_csv(out_file + "_quantile.csv", index=False)


