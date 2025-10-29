import pandas as pd
import numpy as np
import sys
import os

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_mom/")


def calc_se(lst):
    """
    Calculate the Standard Error (SE) of a list or numpy array,
    automatically ignoring np.nan values.
    
    Parameters
    ----------
    lst : list or np.ndarray
        Input data containing numeric values (may include np.nan).
    
    Returns
    -------
    float
        Standard Error (SE). Returns np.nan if sample size <= 1.
    """
    arr = np.array(lst, dtype=float)
    
    # Remove NaN values
    valid = arr[~np.isnan(arr)]
    n = len(valid)
    
    # Cannot compute SE with <= 1 valid sample
    if n <= 1:
        return np.nan

    mean = np.mean(valid)
    # Sample standard deviation (ddof=1 gives unbiased estimate)
    std = np.std(valid, ddof=1)
    
    # Standard error formula: SE = SD / sqrt(n)
    se = std / np.sqrt(n)
    
    return mean, se

def get_h2gxe(h2_add, h2_gxe, h2_nxe, sample_size):
    res_h2_gxe_lst = []
    res_h2_gxe_noNxE_lst = []
    for rep in range(100):
        file = f"h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_{sample_size}.rep_{rep+1}.var"
        try:
            arr = np.loadtxt(file)
            res_h2_gxe_lst.append(arr[1]/np.sum(arr))
        except Exception as e:
            print(f"Error processing {file}: {e}")
            res_h2_gxe_lst.append(np.nan)
        file = f"h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_{sample_size}.rep_{rep+1}.noNxE.var"
        try:
            arr = np.loadtxt(file)
            res_h2_gxe_noNxE_lst.append(arr[1]/np.sum(arr))
        except Exception as e:
            print(f"Error processing {file}: {e}")
            res_h2_gxe_noNxE_lst.append(np.nan)
    return list(calc_se(res_h2_gxe_lst)) + list(calc_se(res_h2_gxe_noNxE_lst))


# Varying h2_nxe
h2_add = 30
h2_gxe = 5
h2_nxe = 15
sample_size = 100000
res_lst = []
for h2_nxe in [0, 5, 15, 25]:
    lst = [h2_add, h2_gxe, h2_nxe]
    lst.extend(get_h2gxe(h2_add, h2_gxe, h2_nxe, sample_size))
    res_lst.append(lst)


# Varying h2_g
h2_add = 30
h2_gxe = 5
h2_nxe = 15
sample_size = 100000
for h2_add in [5, 30, 50]:
    lst = [h2_add, h2_gxe, h2_nxe]
    lst.extend(get_h2gxe(h2_add, h2_gxe, h2_nxe, sample_size))
    res_lst.append(lst)

# Varying h2_gxe
h2_add = 30
h2_gxe = 5
h2_nxe = 15
sample_size = 100000
for h2_gxe in [1, 5, 10]:
    lst = [h2_add, h2_gxe, h2_nxe]
    lst.extend(get_h2gxe(h2_add, h2_gxe, h2_nxe, sample_size))
    res_lst.append(lst)


dfR = pd.DataFrame(res_lst, columns=["h2_add", "h2_gxe", "h2_nxe", "h2_gxe_mean", "h2_gxe_se", "h2_gxe_noNxE_mean", "h2_gxe_noNxE_se"])

dfR.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/h2GxE/5.2_mom_sum.csv", index=False)
