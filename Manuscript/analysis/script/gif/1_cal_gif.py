import pandas as pd
import numpy as np
import sys

in_file = sys.argv[1]
out_file = sys.argv[2]

for i in range(10):
    file = in_file + f".10_{i+1}.res"
    df = pd.read_csv(file, sep=r"\s+", usecols=["p_gxe"])
    df = df.dropna()

quantile_arr = np.array([0.5, 0.05, 0.01, 0.001])
p_observed_arr = np.quantile(df["p_gxe"].values, quantile_arr)

dfR = pd.DataFrame({
    "quantile": quantile_arr,
    "gif": np.log10(p_observed_arr) / np.log10(quantile_arr)
})

dfR.to_csv(out_file, index=False)
