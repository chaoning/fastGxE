import matplotlib
matplotlib.use('Agg')  # Use the 'Agg' backend for headless mode
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/FTO_leadingSNP/")

# Load the CSV file
file_path = 'corr.csv'
corr_df = pd.read_csv(file_path, index_col=0)

# Create a custom colormap
colors = ["#e76254", "#f7aa58", "#ffe6b7"]  # Blue to white to red
n_bins = 100  # Increase this number for a smoother transition
cmap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=n_bins)


plt.figure(figsize=(10/2.54, 8/2.54))
sns.heatmap(corr_df, annot=True, cmap=cmap, fmt=".2f", cbar=True, annot_kws={"size": 8})
plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')
plt.savefig("corr.heatmap.pdf", format='pdf', bbox_inches='tight')
