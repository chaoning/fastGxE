import sys
sys.path.append('/net/zootopia/disk1/chaon/WORK/structLMM/')
from qqplot import qqplot, qqplot_pure, qqplot_frame

import pandas as pd
import numpy as np
from tqdm import tqdm
import os


df = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_rGE_sum/power_snp_h2_add_0.2.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000.csv")

p_arr_lst = []
for R in [0, 0.001, 0.01, 0.1, 0.5]:
    colname = f"rGE_{R}"
    p_arr = df[colname].values
    p_arr_lst.append(p_arr)

color_lst = ["#1e466e", "#528fad", "#aadce0", "#d97c00", "#a94d00"]
label_lst = ["R=0", "R=0.001", "R=0.01", "R=0.1", "R=0.5"]
out_file = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_rGE_sum/power_snp_h2_add_0.2.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000.qqplot"
qqplot(p_arr_lst, color_lst, label_lst, out_file, output_format="png")
