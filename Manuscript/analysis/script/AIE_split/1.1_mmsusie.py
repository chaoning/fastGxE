import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed
from mmsusie import MMSuSiE
import logging
import sys
import random
sys.path.append('/net/zootopia/disk1/chaon/WORK/structLMM/')

trait_type = sys.argv[1]
trait = sys.argv[2]
snp = sys.argv[3]

# Set up logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s') 


MS = MMSuSiE()
bedfile = f"/net/zootopia/disk1/chaon/data/UKB/imp/ukb_imp_info.qc"

data_file = f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/AIE_split/pheno/{trait}.part1.txt"
head_names = pd.read_csv(data_file, sep=r"\s+", nrows=0).columns.tolist()
env_int = head_names[2:]


MS.read_data(data_file, trait=trait, env_int=env_int)

grm_file = "/net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp"
MS.read_sp_grm(grm_file=grm_file)
E = MS.get_env_int()
y = MS.get_y(adjust=True)

G = MS.get_genotype(bedfile, [snp], scale=True)
y = y - G @ np.linalg.lstsq(G, y, rcond=None)[0]

varcom = np.loadtxt(f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/{trait_type}/testGxE/{trait}.var", dtype=float)[:, 0]
MS.cal_spVi(varcom)

EG = E * G
res = MS.mmsusie(EG, y, L=10)
MS.out_mmsusie(res, f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/AIE_split/mmsusie/{trait}.{snp}")
