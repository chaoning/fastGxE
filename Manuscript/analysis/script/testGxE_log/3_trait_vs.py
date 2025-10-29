import pandas as pd
import numpy as np
import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


trait = sys.argv[1]
shortName = sys.argv[2]

# IINT
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/testGxE/")
df_lst = []
for i in range(10):
    file = f"{trait}.10_{i+1}.res"
    df = pd.read_csv(file, sep=r"\s+", usecols=["chrom", "SNP", "base", "allele1", "allele2", "af", "p_gxe"])
    df_lst.append(df)

dfIINT = pd.concat(df_lst, axis=0, ignore_index=True)

try:
    os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/log_trans/testGxE/")
    df_lst = []
    for i in range(10):
        file = f"{trait}.10_{i+1}.res"
        df = pd.read_csv(file, sep=r"\s+", usecols=["SNP", "p_gxe"])
        df_lst.append(df)
    df_log = pd.concat(df_lst, axis=0, ignore_index=True)
    df_log.rename(columns={"p_gxe": "p_gxe_log"}, inplace=True)
except:
    df_log = pd.DataFrame({"SNP": dfIINT["SNP"], "p_gxe_log": dfIINT["p_gxe"]})
    print(f"Log-transformed p-values not found, using IINT values: {trait}.")

df_res = pd.merge(dfIINT, df_log, on="SNP", how="inner")

df_res.to_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/log_trans/vs/{trait}.txt", sep="\t", na_rep="NA", index=False)

x = -np.log10(df_res["p_gxe"])
y = -np.log10(df_res["p_gxe_log"])

# Set global font to Arial and adjust font sizes
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,            # base font size
    "axes.labelsize": 10,      # axis titles
    "xtick.labelsize": 8,      # x-axis tick labels
    "ytick.labelsize": 8,      # y-axis tick labels
    "axes.titlesize": 10       # figure title
})

# Create figure with size in centimeters (converted to inches)
plt.figure(figsize=(6/2.54, 5/2.54))

# Scatter plot of -log10(p) from INT vs log transformation
plt.scatter(x, y, s=8, alpha=0.5, color="#376795")

# Determine maximum axis value (slightly expanded for margin)
maxval = max(x.max(), y.max()) * 1.05

# Add diagonal line (y = x) for reference
plt.plot([0, maxval], [0, maxval], color="black", linestyle="--", lw=1)

# Axis labels with LaTeX-style formatting
plt.xlabel(r'$-\mathrm{log}_{10}(p_{\mathrm{INT}})$')
plt.ylabel(r'$-\mathrm{log}_{10}(p_{\mathrm{log}})$')

# Set axis limits
plt.xlim(0, maxval)
plt.ylim(0, maxval)
plt.title(f"{shortName}")

# Remove top and right spines
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Calculate and annotate correlation
from scipy.stats import pearsonr
r = df_res.loc[:, ["p_gxe", "p_gxe_log"]].corr().iloc[0, 1]
plt.text(
    0.05 * maxval, 0.90 * maxval, 
    f"r = {r:.2f}", 
    fontsize=7, 
    ha="left", va="center"
)

# Optimize spacing so labels are not cut off
plt.tight_layout()

png_out = f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/log_trans/vs/{trait}.pval_compare.png"
plt.savefig(png_out, dpi=300)
plt.close()

df_res = df_res[df_res.loc[:, ["p_gxe", "p_gxe_log"]].min(axis=1) < 5e-8 / 32]
df_res.to_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/log_trans/vs/{trait}.filtered.csv", index=False)
