import sys
sys.path.append('/net/zootopia/disk1/chaon/WORK/structLMM/')
from qqplot import qqplot_frame, qqplot_pure, qqplot
import pandas as pd
import numpy as np
from tqdm import tqdm


# get trait list
dfT = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/gif/PT_trait4.txt",
                  sep=r"\s+")
trait_lst = list(dfT.iloc[:, 0])

def process_trait(trait):
    df = pd.read_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/merge/{trait}.gz", sep=r"\s+", usecols=["p_gxe"])
    df.dropna(inplace=True)
    p_arr = df["p_gxe"].values
    return p_arr

p_arr_lst = []
for trait in tqdm(trait_lst):
    p_arr = process_trait(trait)
    p_arr_lst.append(p_arr)


color_lst = [
    "#1e466e", "#72bcd5", "#ffd06f", "#e76254"
]


label_lst = list(dfT.iloc[:, 1])
out_file = f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/gif/QQ.traits4.PT"
qqplot_frame(p_arr_lst, color_lst, label_lst, out_file, output_format="pdf")
qqplot_pure(p_arr_lst, color_lst, label_lst, out_file, output_format="png")
# qqplot(p_arr_lst, color_lst, label_lst, out_file, output_format="png")
