import os
import glob
import numpy as np
import pandas as pd

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/plink/")

file_lst = glob.glob("*.cov_singleE.csv")

res_lst = []
for file in file_lst:
    df = pd.read_csv(file)
    res_lst.append(df.iloc[:, 0].to_list())

res_df = pd.DataFrame(res_lst)
res_df["file"] = file_lst
mat = res_df.iloc[:, :-1].copy()
print(
    "min:", f"{np.min(mat):.3f}",
    "max:", f"{np.max(mat):.3f}",
    "mean:", f"{np.mean(mat):.3f}",
    "median:", f"{np.median(mat):.3f}",
    "25%:", f"{np.percentile(mat, 25):.3f}",
    "75%:", f"{np.percentile(mat, 75):.3f}"
)



res_df.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/covariantGxE/1.4_sum_singleE.csv", index=False)

import os
import glob
import numpy as np
import pandas as pd

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/plink/")

file_lst = glob.glob("*.cov_allE.csv")

res_lst = []
for file in file_lst:
    df = pd.read_csv(file)
    res_lst.append(df.iloc[0, :].to_list())

res_df = pd.DataFrame(res_lst)
res_df.columns = ["rGI", "add", "GxE", "num_snps"]
res_df["addToGxE"] = res_df.iloc[:, 1] / res_df.iloc[:, 2]
res_df["file"] = file_lst
mat = res_df.iloc[:, :-1].copy()
mat.iloc[:, [1,2]] = mat.iloc[:, [1,2]] * 100
print( "min:\n", mat.min(), 
      "\nmax:\n", mat.max(), 
      "\nmean:\n", mat.mean(), 
      "\nmedian:\n", mat.median(), 
      "\n25%:\n", mat.quantile(0.25), 
      "\n75%:\n", mat.quantile(0.75) )

res_df.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/covariantGxE/1.4_sum_cov_allE.csv", index=False)
