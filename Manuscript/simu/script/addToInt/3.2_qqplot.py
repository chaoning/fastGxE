import sys
sys.path.append('/net/zootopia/disk1/chaon/WORK/structLMM/')
from qqplot import qqplot
import pandas as pd

for unknown_h2SNP_h2E in [0.02, 0.08, 0.16, 0.32, 0.64, 1]:
    p_arr_lst = []
    for R2 in [0, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9]:
        df = pd.read_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_addToint_lm/h2_{unknown_h2SNP_h2E}_R2_{R2}.csv")
        colname = "p_value_proxy"
        p_arr = df[colname].values
        p_arr_lst.append(p_arr)     
    color_lst = ["#1e466e", "#528fad", "#aadce0", "#ffe6b7", "#f2b705", "#d97c00", "#a94d00"]
    label_lst = ["R2=0", "R2=0.01", "R2=0.1", "R2=0.3", "R2=0.5", "R2=0.7", "R2=0.9"]
    out_file = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_addToint_lm_sum/qqplot_h2_{unknown_h2SNP_h2E}"
    qqplot(p_arr_lst, color_lst, label_lst, out_file, output_format="png")
