import numpy as np
import matplotlib
matplotlib.use("Agg")  # Use Agg backend (no display required)
import matplotlib.pyplot as plt
import pandas as pd
import os
from tqdm import tqdm

def plot_three_hists(arrs, titles=None, bins=200, ylabel_name="Frequency", outfile="three_histograms.png"):
    """
    Plot three frequency histograms in one row and save as PNG.
    
    Parameters
    ----------
    arrs : list of array-like
        A list of three arrays/Series to plot.
    titles : list of str, optional
        Titles for each subplot (length must be 3). Default uses generic titles.
    bins : int or str, optional
        Number of bins (default=30). Can also use "auto".
    outfile : str, optional
        Filename of the saved PNG (default="three_histograms.png").
    """
    plt.rcParams.update({
        "font.size": 11,       # base font size
        "axes.titlesize": 10,  # title font size
        "axes.labelsize": 8,   # x/y label font size
        "xtick.labelsize": 6,  # x tick font size
        "ytick.labelsize": 6   # y tick font size
    })
    fig, axes = plt.subplots(1, 3, figsize=(18/2.54, 5/2.54), sharey=True)
    
    for ax, data, title in zip(axes, arrs, titles):
        data = np.asarray(data)
        data = data[~np.isnan(data)]  # remove NaN if any
        ax.hist(
            data,
            bins=bins,
            color="#376795",       # fill color
            edgecolor="#376795"    # edge color
        )
        if title is not None:
            ax.set_title(title)
        ax.set_xlabel("Value")
        ax.grid(False)
    axes[0].set_ylabel(ylabel_name)
    plt.tight_layout()
    
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"Saved figure to {outfile}")


# ===== Example usage =====
if __name__ == "__main__":
    os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/log_trans/pheno_dist/")
    dfT = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/trait32.txt", sep="\t")
    for trait, name in tqdm(zip(dfT.iloc[:, 0], dfT.iloc[:, 2])):
        yRaw = pd.read_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/pheno/Physical_measures/pheno.e42.{trait}.txt", sep=r"\s+", usecols=[trait]).iloc[:, 0]
        yINT = pd.read_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/pheno/Physical_measures/IINT/pheno.e42.{trait}.txt", sep=r"\s+", usecols=[trait]).iloc[:, 0]
        try:
            yLog = pd.read_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/log_trans/pheno/{trait}.txt", sep=r"\s+", usecols=[trait]).iloc[:, 0]
        except Exception as e:
            yLog = np.array([np.nan]*len(yRaw))
        plot_three_hists([yRaw, yINT, yLog], titles=[f"Raw", f"INT", f"Log"], outfile=f"{trait}.png", ylabel_name=name)
    dfT = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/trait67_shortName.txt", sep="\t", dtype=str)
    for trait, name in tqdm(zip(dfT.iloc[:, 0], dfT.iloc[:, -1])):
        yRaw = pd.read_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/pheno/Biological_samples/pheno.e42.{trait}.txt", sep=r"\s+", usecols=[trait]).iloc[:, 0]
        yINT = pd.read_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/pheno/Biological_samples/IINT/pheno.e42.{trait}.txt", sep=r"\s+", usecols=[trait]).iloc[:, 0]
        try:
            yLog = pd.read_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/log_trans/pheno/{trait}.txt", sep=r"\s+", usecols=[trait]).iloc[:, 0]
        except Exception as e:
            yLog = np.array([np.nan]*len(yRaw))
        plot_three_hists([yRaw, yINT, yLog], titles=[f"Raw", f"INT", f"Log"], outfile=f"{trait}.png", ylabel_name=name)
