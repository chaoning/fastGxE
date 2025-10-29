import sys
sys.path.append('/net/zootopia/disk1/chaon/WORK/structLMM/')
from qqplot import qqplot, qqplot_pure, qqplot_frame

import pandas as pd
import numpy as np
from tqdm import tqdm
import os


df = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_addToint_fastgxe_sum/sum.csv")

for unknown_h2SNP_h2E in [0.02, 0.08, 0.16]:
        p_arr_lst = []
        for R2 in [0, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9]:
            colname = f"h2_{unknown_h2SNP_h2E}.R2_{R2}.R2GE_0"
            p_arr = df[colname].values
            p_arr_lst.append(p_arr)     
        color_lst = ["#1e466e", "#528fad", "#aadce0", "#ffe6b7", "#f2b705", "#d97c00", "#a94d00"]
        label_lst = ["R2=0", "R2=0.01", "R2=0.1", "R2=0.3", "R2=0.5", "R2=0.7", "R2=0.9"]
        out_file = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_addToint_fastgxe_sum/qqplot_h2_{unknown_h2SNP_h2E}.R2.rGE_0.pdf"
        qqplot(p_arr_lst, color_lst, label_lst, out_file, output_format="png")
        for R2GE in [0, 0.01, 0.1, 0.3]:
            colname = f"h2_{unknown_h2SNP_h2E}.R2_{0.5}.R2GE_{R2GE}"
            p_arr = df[colname].values
            p_arr_lst.append(p_arr)
        color_lst = ["#1e466e", "#aadce0", "#ffe6b7", "#e76254"]
        label_lst = ["rGE=0", "rGE=0.01", "rGE=0.1", "rGE=0.3"]
        out_file = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_addToint_fastgxe_sum/qqplot_h2_{unknown_h2SNP_h2E}.R2_0.5.rGE.pdf"
        qqplot(p_arr_lst, color_lst, label_lst, out_file, output_format="png")
