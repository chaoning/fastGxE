import pandas as pd
import numpy as np
import sys
import os

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_mom/")

# Varying h2_gxe
h2_add = 30
h2_nxe = 15
sample_size = 100000
res_dct = {}
for h2_gxe in [1, 5, 10]:
    res_h2_add_lst = []
    res_h2_add_noNxE_lst = []
    res_h2_gxe_lst = []
    res_h2_nxe_lst = []
    res_h2_gxe_noNxE_lst = []
    for rep in range(100):
        file = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_mom/h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_{sample_size}.rep_{rep+1}.var"
        try:
            arr = np.loadtxt(file)
            res_h2_add_lst.append(arr[0]/np.sum(arr))
            res_h2_gxe_lst.append(arr[1]/np.sum(arr))
            res_h2_nxe_lst.append(arr[2]/np.sum(arr))
        except Exception as e:
            print(f"Error processing {file}: {e}")
            res_h2_gxe_lst.append(np.nan)
            res_h2_nxe_lst.append(np.nan)
            res_h2_add_lst.append(np.nan)
        
        file = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_mom/h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_{sample_size}.rep_{rep+1}.noNxE.var"
        try:
            arr = np.loadtxt(file)
            res_h2_add_noNxE_lst.append(arr[0]/np.sum(arr))
            res_h2_gxe_noNxE_lst.append(arr[1]/np.sum(arr))
        except Exception as e:
            print(f"Error processing {file}: {e}")
            res_h2_gxe_noNxE_lst.append(np.nan)
            res_h2_add_noNxE_lst.append(np.nan)
    
    res_dct[f"h2_gxe_{h2_gxe}_h2_add_{h2_add}"] = res_h2_add_lst
    res_dct[f"h2_gxe_{h2_gxe}_h2_add_{h2_add}_noNxE"] = res_h2_add_noNxE_lst
    res_dct[f"h2_gxe_{h2_gxe}_h2_gxe"] = res_h2_gxe_lst
    res_dct[f"h2_gxe_{h2_gxe}_h2_gxe_noNxE"] = res_h2_gxe_noNxE_lst
    res_dct[f"h2_gxe_{h2_gxe}_h2_nxe"] = res_h2_nxe_lst

dfR = pd.DataFrame(res_dct)
dfR.to_csv(f"../res_mom_sum/Varying_h2_gxe.csv", index=False)

# Varying h2_nxe
h2_gxe = 5 
h2_add = 30
sample_size = 100000
res_dct = {}
for h2_nxe in [0, 5, 15, 25]:
    res_h2_add_lst = []
    res_h2_add_noNxE_lst = []
    res_h2_gxe_lst = []
    res_h2_nxe_lst = []
    res_h2_gxe_noNxE_lst = []
    for rep in range(100):
        file = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_mom/h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_{sample_size}.rep_{rep+1}.var"
        try:
            arr = np.loadtxt(file)
            res_h2_add_lst.append(arr[0]/np.sum(arr))
            res_h2_gxe_lst.append(arr[1]/np.sum(arr))
            res_h2_nxe_lst.append(arr[2]/np.sum(arr))
        except Exception as e:
            print(f"Error processing {file}: {e}")
            res_h2_add_lst.append(np.nan)
            res_h2_gxe_lst.append(np.nan)
            res_h2_nxe_lst.append(np.nan)
        
        file = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_mom/h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_{sample_size}.rep_{rep+1}.noNxE.var"
        try:
            arr = np.loadtxt(file)
            res_h2_add_noNxE_lst.append(arr[0]/np.sum(arr))
            res_h2_gxe_noNxE_lst.append(arr[1]/np.sum(arr))
        except Exception as e:
            print(f"Error processing {file}: {e}")
            res_h2_add_noNxE_lst.append(arr[0]/np.sum(arr))
            res_h2_gxe_noNxE_lst.append(np.nan)
    res_dct[f"h2_nxe_{h2_nxe}_h2_add_{h2_add}"] = res_h2_add_lst
    res_dct[f"h2_nxe_{h2_nxe}_h2_add_{h2_add}_noNxE"] = res_h2_add_noNxE_lst
    res_dct[f"h2_nxe_{h2_nxe}_h2_gxe"] = res_h2_gxe_lst
    res_dct[f"h2_nxe_{h2_nxe}_h2_nxe"] = res_h2_nxe_lst
    res_dct[f"h2_nxe_{h2_nxe}_h2_gxe_noNxE"] = res_h2_gxe_noNxE_lst

dfR = pd.DataFrame(res_dct)
dfR.to_csv(f"../res_mom_sum/Varying_h2_nxe.csv", index=False)
