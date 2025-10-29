import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# Set font to Arial
plt.rcParams['font.family'] = 'Arial'

# Set working directory
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/GxE_signal_h2/")

# Read CSV file and convert h2gxe to percentage
df = pd.read_csv("2.1_real_data_5e-8.csv")
data_pct = df["h2gxe"] * 100  # convert to percentage

# 3-sigma rule: keep only values within mean ± 3*std
mean_val_all = np.mean(data_pct)
std_val_all = np.std(data_pct)
filtered_data_pct = data_pct[
    (data_pct >= mean_val_all - 3 * std_val_all) &
    (data_pct <= mean_val_all + 3 * std_val_all)
]

# Calculate summary statistics for filtered data
mean_val = np.mean(data_pct)
median_val = np.median(data_pct)
q25 = np.percentile(data_pct, 25)
q75 = np.percentile(data_pct, 75)

# Set small font sizes for publication style
plt.rcParams.update({
    'font.size': 6,         # base font size
    'axes.titlesize': 6,    # title font size
    'axes.labelsize': 6,    # x/y label font size
    'xtick.labelsize': 5,   # x-axis tick labels
    'ytick.labelsize': 5,   # y-axis tick labels
    'legend.fontsize': 5    # legend font size
})

# Create figure
fig, ax = plt.subplots(figsize=(6/2.54, 5/2.54))
ax.hist(filtered_data_pct, bins=30, edgecolor='black', alpha=0.7)

# Add vertical lines for summary statistics
ax.axvline(mean_val, color='red', linestyle='--', label=f"Mean = {mean_val:.3f}%")
ax.axvline(median_val, color='blue', linestyle='--', label=f"Median = {median_val:.3f}%")
ax.axvline(q25, color='green', linestyle=':', label=f"25% = {q25:.3f}%")
ax.axvline(q75, color='purple', linestyle=':', label=f"75% = {q75:.3f}%")

# Configure labels and limits
ax.set_xlabel(r"$h^2_{\mathrm{GxE}}\ (\%)$")
ax.set_ylabel("Frequency")

# Remove title, top & right border, and grid
ax.set_title("")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(False)

# Add legend inside plot
ax.legend(frameon=False, loc='upper right')

# Save figure
plt.tight_layout()
plt.savefig("2.1_real_data_5e-8_plot.png", dpi=300)
plt.close()

print("Plot completed and saved as 2.1_real_data_5e-8_plot.png")
