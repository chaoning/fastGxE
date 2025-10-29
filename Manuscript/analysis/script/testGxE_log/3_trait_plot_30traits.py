# -*- coding: utf-8 -*-
"""
Make a 5x6 grid (30 traits) of scatter plots comparing -log10(p_INT) vs -log10(p_log).
Add vertical and horizontal threshold lines at -log10(5e-8/32).
No CSV outputs.

Author: Chao + ChatGPT
Date: 2025-09-11
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from tqdm import tqdm

# ----------------------------- Configs -----------------------------
TRAIT_FILE = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/trait32.txt"
DATA_DIR   = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/log_trans/vs/"
OUT_DIR    = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/log_trans/vs"

os.makedirs(OUT_DIR, exist_ok=True)

# 5 columns x 6 rows
N_COLS, N_ROWS = 5, 6
DPI = 300
FIG_W_CM, FIG_H_CM = 25.0, 24.0
FIGSIZE = (FIG_W_CM / 2.54, FIG_H_CM / 2.54)

# Threshold line
THRESH_P = 5e-8 / 32.0
THRESH_LOGP = -np.log10(THRESH_P)

# Matplotlib global style
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "savefig.dpi": DPI,
})

# ----------------------------- Helpers -----------------------------
def _neglog10_safe(pvals: np.ndarray) -> np.ndarray:
    """Safely compute -log10(p)."""
    arr = np.asarray(pvals, dtype=float)
    eps = np.finfo(float).tiny
    arr = np.where((arr > 0) & np.isfinite(arr), arr, eps)
    arr = np.clip(arr, eps, 1.0)
    return -np.log10(arr)

def _title_from_row(row: pd.Series) -> str:
    """Pick a nice title name from row info."""
    for k in ("ShortName", "shortName", "TraitName"):
        if k in row and pd.notna(row[k]) and str(row[k]).strip():
            return str(row[k])
    return str(row["FieldID"])

# ----------------------------- Load traits -----------------------------
dfT = pd.read_csv(TRAIT_FILE, sep="\t", header=0)
dfT = dfT[(dfT["FieldID"] != "78") & (dfT["FieldID"] != "WHRadjBMI")].copy()

if len(dfT) < 30:
    raise ValueError(f"Expected 30 traits after filtering, got {len(dfT)}.")
traits_df = dfT.iloc[:30].reset_index(drop=True)

id2title = {row["FieldID"]: _title_from_row(row) for _, row in traits_df.iterrows()}

# ----------------------------- Prepare figure -----------------------------
fig, axes = plt.subplots(nrows=N_ROWS, ncols=N_COLS, figsize=FIGSIZE)
axes = axes.flatten()

warnings.filterwarnings("ignore", category=RuntimeWarning)

r_val_lst = []
# ----------------------------- Plot each trait -----------------------------
for ax, field_id in tqdm(zip(axes, traits_df["FieldID"])):
    file_path = os.path.join(DATA_DIR, f"{field_id}.txt")
    df = pd.read_csv(file_path, sep=r"\s+", header=0)
    df.dropna(inplace=True)
    # Scatter data
    x = _neglog10_safe(df["p_gxe"].values)
    y = _neglog10_safe(df["p_gxe_log"].values)

    # Pearson correlation

    r_val, _ = pearsonr(df["p_gxe"].values, df["p_gxe_log"].values)
    r_val_lst.append(r_val)

    # Axis range
    maxval = float(np.nanmax([x.max(), y.max()])) if x.size and y.size else 1.0
    maxval = max(1.0, maxval) * 1.05

    # Scatter
    ax.scatter(x, y, s=6, alpha=0.45)

    # Reference diagonal
    ax.plot([0, maxval], [0, maxval], linestyle="--", linewidth=0.8, color="black")

    # Threshold lines
    ax.axvline(THRESH_LOGP, linestyle=":", linewidth=0.8, color="red")
    ax.axhline(THRESH_LOGP, linestyle=":", linewidth=0.8, color="red")

    # Style
    ax.set_xlim(0, maxval)
    ax.set_ylim(0, maxval)
    ax.set_title(id2title[field_id], pad=2)
    ax.set_xlabel(r'$-\mathrm{log}_{10}(p_{\mathrm{INT}})$')
    ax.set_ylabel(r'$-\mathrm{log}_{10}(p_{\mathrm{log}})$')

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    ax.text(0.05 * maxval, 0.90 * maxval,
            f"r = {r_val:.2f}" if np.isfinite(r_val) else "r = NA",
            ha="left", va="center", fontsize=7)

# Hide extra axes if any
for ax in axes[len(traits_df):]:
    ax.axis("off")

fig.tight_layout(pad=0.8)
out_png = os.path.join(OUT_DIR, "pval_compare_30traits_5x6.png")
fig.savefig(out_png, dpi=DPI, bbox_inches="tight")
plt.close(fig)

print(np.mean(r_val_lst), np.median(r_val_lst), np.min(r_val_lst), np.max(r_val_lst))
print(f"Saved figure to: {out_png}")
