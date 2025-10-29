import os

import numpy as np
import pandas as pd
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_GENIE/")

# Varying h2_gxe
h2_add = 30
h2_gxe = 5
h2_nxe = 15
sample_size = 100000
res_dct = {}
for h2_gxe in [1, 5, 10]:
    prefix = f"h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_{sample_size}"
    h2g_lst = []
    h2gxe_lst = []
    h2nxe_lst = []
    h2e_lst = []
    for rep in range(1, 101):
        file = f"{prefix}.{rep}.G_GxE_NxE.out2"
        try:
            with open(file) as f:
                h2g = 0
                h2gxe = 0
                h2nxe = 0
                h2e = 0
                total = 0
                for line in f:
                    if "Sigma^2_g[0]" in line:
                        arr = line.strip().split()
                        h2g += float(arr[2])
                    if "Sigma^2_gxe" in line:
                        arr = line.strip().split()
                        h2gxe += float(arr[2])
                    if "Sigma^2_nxe" in line:
                        arr = line.strip().split()
                        h2nxe += float(arr[2])
                    if "Sigma^2_e" in line:
                        arr = line.strip().split()
                        h2e += float(arr[2])
                    if "Sigma^2" in line:
                        arr = line.strip().split()
                        total += float(arr[2])
                h2g_lst.append(h2g / total)
                h2gxe_lst.append(h2gxe / total)
                h2nxe_lst.append(h2nxe / total)
                h2e_lst.append(h2e / total)
        except Exception as e:
            print(f"Error processing {file}: {e}")
            h2g_lst.append(np.nan)
            h2gxe_lst.append(np.nan)
            h2nxe_lst.append(np.nan)
            h2e_lst.append(np.nan)
    res_dct[f"h2_g_{h2_gxe}"] = h2g_lst
    res_dct[f"h2_gxe_{h2_gxe}"] = h2gxe_lst
    res_dct[f"h2_nxe_{h2_gxe}"] = h2nxe_lst
    res_dct[f"h2_e_{h2_gxe}"] = h2e_lst

dfR = pd.DataFrame(res_dct)
dfR.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_GENIE_sum/Varying_h2_gxe2.csv", index=False, na_rep="NA")

# Varying h2_nxe
h2_add = 30
h2_gxe = 5
h2_nxe = 15
sample_size = 100000
res_dct = {}
for h2_nxe in [0, 5, 15, 25]:
    prefix = f"h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_{sample_size}"
    h2g_lst = []
    h2gxe_lst = []
    h2nxe_lst = []
    h2e_lst = []
    for rep in range(1, 101):
        file = f"{prefix}.{rep}.G_GxE_NxE.out2"
        try:
            with open(file) as f:
                h2g = 0
                h2gxe = 0
                h2nxe = 0
                h2e = 0
                total = 0
                for line in f:
                    if "Sigma^2_g[0]" in line:
                        arr = line.strip().split()
                        h2g += float(arr[2])
                    if "Sigma^2_gxe" in line:
                        arr = line.strip().split()
                        h2gxe += float(arr[2])
                    if "Sigma^2_nxe" in line:
                        arr = line.strip().split()
                        h2nxe += float(arr[2])
                    if "Sigma^2_e" in line:
                        arr = line.strip().split()
                        h2e += float(arr[2])
                    if "Sigma^2" in line:
                        arr = line.strip().split()
                        total += float(arr[2])
                h2g_lst.append(h2g / total)
                h2gxe_lst.append(h2gxe / total)
                h2nxe_lst.append(h2nxe / total)
                h2e_lst.append(h2e / total)
        except Exception as e:
            print(f"Error processing {file}: {e}")
            h2g_lst.append(np.nan)
            h2gxe_lst.append(np.nan)
            h2nxe_lst.append(np.nan)
            h2e_lst.append(np.nan)
    res_dct[f"h2_g_{h2_nxe}"] = h2g_lst
    res_dct[f"h2_gxe_{h2_nxe}"] = h2gxe_lst
    res_dct[f"h2_nxe_{h2_nxe}"] = h2nxe_lst
    res_dct[f"h2_e_{h2_nxe}"] = h2e_lst

dfR = pd.DataFrame(res_dct)
dfR.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_GENIE_sum/Varying_h2_nxe2.csv", index=False, na_rep="NA")


# Varying h2_add
h2_add = 30
h2_gxe = 5
h2_nxe = 15
sample_size = 100000
res_dct = {}
for h2_add in [5, 30, 50]:
    prefix = f"h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_{sample_size}"
    h2g_lst = []
    h2gxe_lst = []
    h2nxe_lst = []
    h2e_lst = []
    for rep in range(1, 101):
        file = f"{prefix}.{rep}.G_GxE_NxE.out2"
        try:
            with open(file) as f:
                h2g = 0
                h2gxe = 0
                h2nxe = 0
                h2e = 0
                total = 0
                for line in f:
                    if "Sigma^2_g[0]" in line:
                        arr = line.strip().split()
                        h2g += float(arr[2])
                    if "Sigma^2_gxe" in line:
                        arr = line.strip().split()
                        h2gxe += float(arr[2])
                    if "Sigma^2_nxe" in line:
                        arr = line.strip().split()
                        h2nxe += float(arr[2])
                    if "Sigma^2_e" in line:
                        arr = line.strip().split()
                        h2e += float(arr[2])
                    if "Sigma^2" in line:
                        arr = line.strip().split()
                        total += float(arr[2])
                h2g_lst.append(h2g / total)
                h2gxe_lst.append(h2gxe / total)
                h2nxe_lst.append(h2nxe / total)
                h2e_lst.append(h2e / total)
        except Exception as e:
            print(f"Error processing {file}: {e}")
            h2g_lst.append(np.nan)
            h2gxe_lst.append(np.nan)
            h2nxe_lst.append(np.nan)
            h2e_lst.append(np.nan)
    res_dct[f"h2_g_{h2_add}"] = h2g_lst
    res_dct[f"h2_gxe_{h2_add}"] = h2gxe_lst
    res_dct[f"h2_nxe_{h2_add}"] = h2nxe_lst
    res_dct[f"h2_e_{h2_add}"] = h2e_lst

dfR = pd.DataFrame(res_dct)
dfR.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_GENIE_sum/Varying_h2_add2.csv", index=False, na_rep="NA")

