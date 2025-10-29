import matplotlib
matplotlib.use("Agg")  # Use non-interactive backend for headless environments

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import os

# Get trait name from command-line argument
trait = sys.argv[1]
dfT = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/trait32.txt", sep="\t")
trait_name = dfT[dfT["FieldID"] == trait]["ShortName"].values[0]

# === Step 1: Load p-values from model without NxE interaction ===
file_path = f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/Physical_measures/result/SGEM/SGEM-O/IINT/{trait}.gz"
df_no_nxe = pd.read_csv(
    file_path,
    sep=r"\s+",
    compression="gzip",
    usecols=["SNP", "p_SGEM_O"]
).rename(columns={"p_SGEM_O": "p_gxe_noNxE"})

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/testGxE/")
replicates = [
    pd.read_csv(f"{trait}.noNxE.10_{i+1}.res", sep=r"\s+", usecols=["SNP", "p_gxe"])
    for i in range(10)
]
df_no_nxe = pd.concat(replicates, ignore_index=True)
df_no_nxe.rename(columns={"p_gxe": "p_gxe_noNxE"}, inplace=True)


# === Step 2: Load p-values from model with NxE interaction (10 replicates) ===
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/testGxE/")
replicates = [
    pd.read_csv(f"{trait}.10_{i+1}.res", sep=r"\s+", usecols=["SNP", "p_gxe"])
    for i in range(10)
]
df_nxe = pd.concat(replicates, ignore_index=True)

# === Step 3: Merge and drop missing values ===
df_merged = pd.merge(df_no_nxe, df_nxe, on="SNP").dropna()

# === Step 4: Compute -log10(p) values ===
logp_no_nxe = -np.log10(df_merged["p_gxe_noNxE"])
logp_nxe = -np.log10(df_merged["p_gxe"])

# === Step 5: Plot scatter plot and save as PNG ===
plt.figure(figsize=(7/2.54, 6/2.54)) 
plt.scatter(logp_no_nxe, logp_nxe, s=5, alpha=0.5, color="#1e466e")

# Add y = x reference line
lims = [
    np.min([logp_no_nxe.min(), logp_nxe.min()]),
    np.max([logp_no_nxe.max(), logp_nxe.max()])
]
plt.plot(lims, lims, color='gray', linestyle='--', linewidth=1)

# Axis labels and title
plt.xlabel(r'$-\log_{10}(p)$ without NxE', fontsize=8)
plt.ylabel(r'$-\log_{10}(p)$ with NxE', fontsize=8)
plt.title(trait, fontsize=10)

# Grid and frame customization
plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.3)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

plt.tight_layout()

# Save figure
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/p_vs/")
plt.savefig(f"{trait_name}.comparison2.png", dpi=300)
plt.close()
