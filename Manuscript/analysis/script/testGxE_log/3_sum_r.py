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
DATA_DIR   = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/testGxE_merge"


# ----------------------------- Load traits -----------------------------
dfT = pd.read_csv(TRAIT_FILE, sep="\t", header=0)
dfT = dfT[(dfT["FieldID"] != "78") & (dfT["FieldID"] != "WHRadjBMI")].copy()
print(dfT["FieldID"].tolist())
r_lst = []
for field_id in tqdm(dfT["FieldID"]):
    file = f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/testGxE_merge/{field_id}_correlation.csv"
    dfR = pd.read_csv(file, header=0, index_col=0)
    r_lst.append(dfR.loc["p_gxe_log", "p_gxe"])

print(np.min(r_lst), np.max(r_lst), np.mean(r_lst), np.median(r_lst))
