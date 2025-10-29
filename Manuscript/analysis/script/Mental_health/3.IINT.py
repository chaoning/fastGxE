import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm

sys.path.append("/net/zootopia/disk1/chaon/WORK/structLMM/")
from inverse_normal_transform import inverse_normal_transform

# Read input arguments
in_file = sys.argv[1]
out_file = sys.argv[2]

# Read the input file
df = pd.read_csv(in_file, sep="\s+")

# Extract and transform the response variable
y = df.iloc[:, 1].copy()

# Create the design matrices
X = df.loc[:, "Genotype_batch":"Sex_Age3"].copy()
X = sm.add_constant(X)

# Calculate beta coefficients and residuals
beta = np.linalg.lstsq(X, y, rcond=None)[0]
y_resi = y - np.dot(X, beta)

# Transform residuals
y_resi_trans = inverse_normal_transform(y_resi)

# Select columns and prepare the output DataFrame
envi_lst = list(df.loc[:, "Age":"Confide"].columns)
col_keep_lst = list(df.iloc[:, :2].columns)
col_keep_lst.extend(envi_lst)

dfR = df.loc[:, col_keep_lst].copy()
dfR.iloc[:, 1] = y_resi_trans
dfR.to_csv(out_file, sep=" ", index=False)
