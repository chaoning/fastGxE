import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os 

# Set global font to Arial
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 7   # base font size

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/h2GxE/")

# ===== Load data =====
df1 = pd.read_csv("5.2_mom_sum.csv")
df2 = pd.read_csv("5.1_GENIE_sum.csv")

# Concatenate mom results (df1) with GENIE results (last two columns in df2)
df = pd.concat([df1, df2.iloc[:, -2:]], axis=1)

# Rename columns for clarity
df.columns = ["h2_g", "h2_gxe", "h2_nxe", 
              "fastGxE_mean", "fastGxE_se", 
              "fastGxE-noNxE_mean", "fastGxE-noNxE_se", 
              "GENIE_mean", "GENIE_se"]

df.iloc[:, :3] /= 100

# Define the three methods, columns, and colors
methods = {
    "fastGxE": ("fastGxE_mean", "fastGxE_se", "#376795"),
    "fastGxE-noNxE": ("fastGxE-noNxE_mean", "fastGxE-noNxE_se", "#ffd06f"),
    "GENIE": ("GENIE_mean", "GENIE_se", "#e76254")
}

# ===== Create figure with 1 row and 3 subplots =====
fig, axes = plt.subplots(1, 3, figsize=(18 / 2.54, 6 / 2.54), dpi=300)

# Define subsets of the data for each figure
subsets = [
    (df.iloc[[0, 2, 3], :], "h2_nxe", r"$\rho_{NxE}$"),
    (df.iloc[4:7, :], "h2_g",   r"$h^2_{G}$"),
    (df.iloc[7:10, :], "h2_gxe", r"$h^2_{GxE}$"),
]

# Loop over each subplot
for idx, (ax, (df_sub, factor, xlabel)) in enumerate(zip(axes, subsets)):
    x = np.arange(len(df_sub[factor]))
    width = 0.25
    
    # Plot mean ± 2*SE for each method
    for i, (method, (mean_col, se_col, color)) in enumerate(methods.items()):
        ax.bar(
            x + i*width,
            df_sub[mean_col],
            width,
            yerr=2*df_sub[se_col],
            capsize=3,
            label=method,
            color=color,
            error_kw={'elinewidth': 0.5}  # make error bars thinner
        )

    
    # Customize axis
    ax.set_xticks(x + width)
    ax.set_xticklabels(df_sub[factor])
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(r"Estimated $h^2_{GxE}$", fontsize=8)
    ax.tick_params(axis='both', labelsize=7)
    
    # Remove top and right spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # === Add horizontal reference lines ===
    if idx == 0:  # Figure 1: h2_nxe
        ax.axhline(y=0.05, linestyle="--", color="gray", linewidth=0.5)
    elif idx == 1:  # Figure 2: h2_g
        ax.axhline(y=0.05, linestyle="--", color="gray", linewidth=0.5)
    elif idx == 2:  # Figure 3: h2_gxe
        for y in [0.01, 0.05, 0.1]:
            ax.axhline(y=y, linestyle="--", color="gray", linewidth=0.5)

# Add subplot labels (a), (b), (c)
labels = ['a', 'b', 'c']
for idx, ax in enumerate(axes):
    ax.text(-0.2, 1.05, labels[idx], transform=ax.transAxes,
            fontsize=12, fontweight='normal', va='top', ha='right')
# Add legend below the second subplot
axes[1].legend(
    loc="upper center",
    bbox_to_anchor=(0.5, -0.4),
    ncol=3,
    frameon=False,
    fontsize=7
)

plt.tight_layout()
plt.subplots_adjust(wspace=0.4)
plt.savefig("Comparison_three_methods.png", dpi=300)
plt.close()
