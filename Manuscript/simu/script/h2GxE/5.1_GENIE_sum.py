import os

import numpy as np
import pandas as pd
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_GENIE/")

import numpy as np

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

def get_h2_gxe(prefix):
    # prefix = f"h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_{sample_size}"
    h2gxe_lst = []
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
                h2gxe_lst.append(h2gxe / total)
        except Exception as e:
            print(f"Error processing {file}: {e}")
            h2gxe_lst.append(np.nan)
    return calc_se(h2gxe_lst)

# Varying h2_nxe
h2_add = 30
h2_gxe = 5
h2_nxe = 15
sample_size = 100000
res_lst = []
for h2_nxe in [0, 5, 15, 25]:
    prefix = f"h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_{sample_size}"
    mean, se = get_h2_gxe(prefix)
    res_lst.append([h2_add, h2_gxe, h2_nxe, mean, se])


# Varying h2_g
h2_add = 30
h2_gxe = 5
h2_nxe = 15
sample_size = 100000
for h2_add in [5, 30, 50]:
    prefix = f"h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_{sample_size}"
    mean, se = get_h2_gxe(prefix)
    res_lst.append([h2_add, h2_gxe, h2_nxe, mean, se])


# Varying h2_gxe
h2_add = 30
h2_gxe = 5
h2_nxe = 15
sample_size = 100000
for h2_gxe in [1, 5, 10]:
    prefix = f"h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_{sample_size}"
    mean, se = get_h2_gxe(prefix)
    res_lst.append([h2_add, h2_gxe, h2_nxe, mean, se])

dfR = pd.DataFrame(res_lst, columns=["h2_g", "h2_gxe", "h2_nxe", "mean_h2_gxe_est", "se_h2_gxe_est"])
dfR.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/h2GxE/5.1_GENIE_sum.csv", index=False)

