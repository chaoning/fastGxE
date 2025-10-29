import pandas as pd
import numpy as np
import os
import re
import glob

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/")
prefix = "h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000"
power_lst = []
num_envi = 30
for power_snp_h2_gxe in [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14]:
    power_mmSuSiE = 0
    for rep in range(1, 101):
        try:
            envi_file = f"./pheno/power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{prefix}.rep_{rep}.power_snp_envi.txt"
            envi_index_lst = pd.read_csv(envi_file, sep=r"\s+").iloc[:, 0].tolist()
            pattern = re.compile(rf"./res_fastgxe_testgxe/power_snp_h2_gxe_{power_snp_h2_gxe}\.num_envi_{num_envi}\.{prefix}\.rep_{rep}\.[0-9]+_[0-9]+\.res")
            all_files = glob.glob(f"./res_fastgxe_testgxe/*")
            matched_files = [f for f in all_files if pattern.search(f)]
            if len(matched_files) != 1:
                continue
            file = matched_files[0]
            df = pd.read_csv(file, sep=r"\s+")
            p = df.iloc[0, -1]
            if p < 5e-8:
                pip_file = f"./res_mmSuSiE/power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{rep}.mmsusie_out.pip.txt"
                df_pip = pd.read_csv(pip_file, sep=r"\s+", header=None)
                power_mmSuSiE += np.sum(df_pip.iloc[envi_index_lst, 0] > 0.5)
        except Exception as e:
            print(e)
            continue
    power_mmSuSiE /= (100 * num_envi)
    power_lst.append([power_snp_h2_gxe, num_envi, power_mmSuSiE, "mmSuSiE"])

    power_fastGWA = 0
    for rep in range(1, 101):
        try:
            envi_file = f"./pheno/power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{prefix}.rep_{rep}.power_snp_envi.txt"
            envi_index_lst = pd.read_csv(envi_file, sep=r"\s+").iloc[:, 0].tolist()
            for envi_index in envi_index_lst:
                file = f"./res_fastGWA_h2gxe/power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{prefix}.{envi_index+1}.{rep}.fastGWA.res"
                df = pd.read_csv(file, sep=r"\s+")
                p = df.iloc[0, -1]
                if p < 5e-8 / 40:
                    power_fastGWA += 1
        except Exception as e:
            print(e)
            continue
    power_fastGWA /= (100 * num_envi)
    print(f"{power_snp_h2_gxe}\t{num_envi}\t{power_mmSuSiE}\t{power_fastGWA}")
    power_lst.append([power_snp_h2_gxe, num_envi, power_fastGWA, "fastGWA-GE"])

power_snp_h2_gxe = 0.08
for num_envi in [1, 2, 5, 10, 20, 30, 40]:
    power_mmSuSiE = 0
    for rep in range(1, 101):
        try:
            envi_file = f"./pheno/power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{prefix}.rep_{rep}.power_snp_envi.txt"
            envi_index_lst = pd.read_csv(envi_file, sep=r"\s+").iloc[:, 0].tolist()
            pattern = re.compile(rf"./res_fastgxe_testgxe/power_snp_h2_gxe_{power_snp_h2_gxe}\.num_envi_{num_envi}\.{prefix}\.rep_{rep}\.[0-9]+_[0-9]+\.res")
            all_files = glob.glob(f"./res_fastgxe_testgxe/*")
            matched_files = [f for f in all_files if pattern.search(f)]
            if len(matched_files) != 1:
                continue
            file = matched_files[0]
            df = pd.read_csv(file, sep=r"\s+")
            p = df.iloc[0, -1]
            if p < 5e-8:
                pip_file = f"./res_mmSuSiE/power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{rep}.mmsusie_out.pip.txt"
                df_pip = pd.read_csv(pip_file, sep=r"\s+", header=None)
                power_mmSuSiE += np.sum(df_pip.iloc[envi_index_lst, 0] > 0.5)
        except Exception as e:
            print(e)
            continue
    power_mmSuSiE /= (100 * num_envi)
    power_lst.append([power_snp_h2_gxe, num_envi, power_mmSuSiE, "mmSuSiE"])

    power_fastGWA = 0
    for rep in range(1, 101):
        try:
            envi_file = f"./pheno/power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{prefix}.rep_{rep}.power_snp_envi.txt"
            envi_index_lst = pd.read_csv(envi_file, sep=r"\s+").iloc[:, 0].tolist()
            for envi_index in envi_index_lst:
                file = f"./res_fastGWA_nE/power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{prefix}.{envi_index+1}.{rep}.fastGWA.res"
                df = pd.read_csv(file, sep=r"\s+")
                p = df.iloc[0, -1]
                if p < 5e-8 / 40:
                    power_fastGWA += 1
        except Exception as e:
            print(e)
            continue
    power_fastGWA /= (100 * num_envi)
    print(f"{power_snp_h2_gxe}\t{num_envi}\t{power_mmSuSiE}\t{power_fastGWA}")
    power_lst.append([power_snp_h2_gxe, num_envi, power_fastGWA, "fastGWA-GE"])

df_power = pd.DataFrame(power_lst, columns=["power_snp_h2_gxe", "num_envi", "Power", "Method"])
df_power.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_mmSuSiE/power_summary.csv", index=False)

