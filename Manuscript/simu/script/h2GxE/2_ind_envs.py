import numpy as np
import pandas as pd

dfE = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/envi.csv", dtype={"eid": str})

# Set parameters
n_samples = dfE.shape[0]  # Number of rows
n_variables = 40    # Number of independent variables

# Generate independent variables from standard normal distribution
np.random.seed(42)
data = np.random.randn(n_samples, n_variables)

dfE.loc[:, "age":"alcohol_frequency"] = data
dfE.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/envi_ind.csv", index=False)

corr_val = dfE.iloc[:, 1:].corr().values

# Extract non-diagonal elements
non_diag_vals = corr_val[~np.eye(corr_val.shape[0], dtype=bool)]

# Get max and min of non-diagonal values
max_corr = non_diag_vals.max()
min_corr = non_diag_vals.min()

print(max_corr, min_corr)