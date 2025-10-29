import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed
from mmsusie import MMSuSiE
import logging
import sys
import random
sys.path.append('/net/zootopia/disk1/chaon/WORK/structLMM/')
from susie import susie


# Set up logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s') 

rep = int(sys.argv[1])
input_window_index = int(sys.argv[2])

start_pos = 16488635
end_pos = 51237712
window_size = 250_000

window_index = 1
for window_start in range(start_pos, end_pos, window_size):
    window_end = window_start + window_size
    if input_window_index == window_index:
        break
    window_index += 1

logging.info(f"Processing window {window_index} from {window_start} to {window_end}")

MS = MMSuSiE()

prefix = "h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000"

bedfile = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22"

data_file = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_fastgxe/{prefix}_part1.txt"
head_names = pd.read_csv(data_file, sep=r"\s+", nrows=0).columns.tolist()
env_int = head_names[1:41]

logging.info(f"Processing replicate {rep}")
MS = MMSuSiE()

MS.read_data(data_file, trait=f"trait{rep}", env_int=env_int)

grm_file = "/net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp"
MS.read_sp_grm(grm_file=grm_file)
E = MS.get_env_int()
y = MS.get_y(adjust=True)

file = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_fastgxe_testgxe/{prefix}.rep_{rep}.res"
df = pd.read_csv(file, sep=r"\s+")
df_window = df[(df["base"] >= window_start) & (df["base"] < window_end)]

np.random.seed(rep * 1000 + input_window_index)
random_idx = np.random.choice(len(df_window), size=1, replace=False)
lead_snp = df_window.iloc[random_idx[0], 2]

G = MS.get_genotype(bedfile, [lead_snp], scale=True)
y = y - G @ np.linalg.lstsq(G, y, rcond=None)[0]

varcom = np.loadtxt(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_fastgxe_testgxe/{prefix}.rep_{rep}.var", dtype=float)[:, 0]
MS.cal_spVi(varcom)

logging.info(f"Vi.shape = {MS.Vi.shape}")

logging.info(f"G.shape = {G.shape}, E.shape = {E.shape}, y.shape = {y.shape}")
EG = E * G
res = MS.mmsusie(EG, y, L=10)
MS.out_mmsusie(res, f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/AIE_split_mmsusie_random/{rep}.{input_window_index}.mmsusie_out")

with open(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/AIE_split_mmsusie_random/{rep}.{input_window_index}.mmsusie_out.lead_snp", "w") as f:
    f.write(lead_snp + "\n")
