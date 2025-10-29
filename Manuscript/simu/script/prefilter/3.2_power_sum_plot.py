import os
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for headless environments

import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

num_envi_power = sys.argv[1]

# ======== Set global font to Arial ========
plt.rcParams['font.family'] = 'Arial'

base_dir = "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power"
out_dir = f"{base_dir}/res_filter_sum"
os.chdir(out_dir)


# Read the CSV file
df = pd.read_csv(f"num_envi_{num_envi_power}.power_summary_all.csv")
df = df[df["add_gxe_ratio"] > 0.02]  # Filter out rows where add_gxe_ratio <= 0.02

# ======== Get unique add_gxe_ratio values ========
unique_ratios = sorted(df["add_gxe_ratio"].unique())

# ======== Create subplots ========
fig, axes = plt.subplots(
    1, len(unique_ratios),                 # 1 row, N columns
    figsize=(4 * len(unique_ratios) / 2.54, 6 / 2.54),    # Adjust figure width
    dpi=150,
    sharey=True                             # Share Y-axis for better comparison
)
# ======== Font size settings ========
title_fontsize = 8
label_fontsize = 8
tick_fontsize = 6
legend_fontsize = 6

# ======== Plot each ratio ========
for ax, ratio in zip(axes, unique_ratios):
    # Filter data
    sub_df = df[df["add_gxe_ratio"] == ratio].sort_values("rGI")
    
    # Plot with custom colors
    ax.plot(sub_df["rGI"], sub_df["power_gxe"], marker="o", color="#376795", label="fastGxE")
    ax.plot(sub_df["rGI"], sub_df["power_filter"], marker="s", color="#e76254", label="fastGxE-filter")
    
    # Format main effects value
    main_effect_val = 1.5 * ratio
    formatted_val = f"{main_effect_val:.4f}".rstrip('0').rstrip('.')
    
    # Title with larger font
    ax.set_title(f"SNP main effects {formatted_val} %", fontsize=title_fontsize)
    
    # X-axis label (GI subscript) with larger font
    ax.set_xlabel(r"$r_{GI}$", fontsize=label_fontsize)
    
    # Fixed X ticks and font size
    ax.set_xticks([-0.9, -0.3, 0, 0.3, 0.9])
    ax.tick_params(axis='x', labelsize=tick_fontsize)
    ax.tick_params(axis='y', labelsize=tick_fontsize)
    
    # Grid
    ax.grid(True, linestyle="--", alpha=0.5)
    
    # Y-axis label only for first subplot
    if ax == axes[0]:
        ax.set_ylabel("Power", fontsize=label_fontsize)
    
    # Legend only on last subplot
    if ax == axes[-1]:
        ax.legend(fontsize=legend_fontsize)

# ======== Adjust layout ========
plt.tight_layout()

# ======== Save figure ========
plt.savefig(f"power_compare_num_envi_{num_envi_power}.png")
