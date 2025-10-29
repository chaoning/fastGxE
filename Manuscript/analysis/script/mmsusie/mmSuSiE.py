import sys
import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed
from mmsusie import MMSuSiE
import logging
import sys
import os
sys.path.append('/net/zootopia/disk1/chaon/WORK/structLMM/')
from susie import susie

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/run/")

trait = sys.argv[1]
snp = sys.argv[2]
data_file = sys.argv[3]
varcom_file = sys.argv[4]

bedfile = "/net/zootopia/disk1/chaon/data/UKB/imp/ukb_imp_info.qc"
snp_on_disk = Bed(bedfile, count_A1=True)
nSNP = snp_on_disk.sid_count
dfFam = pd.read_csv(bedfile + ".fam", sep=r"\s+", header=None)


logging.info("Pheno file")

MS = MMSuSiE()

head_names = pd.read_csv(data_file, sep=r"\s+", nrows=0).columns.tolist()
env_int = head_names[2:]

MS.read_data(data_file, trait=f"{trait}", env_int=env_int)

grm_file = "/net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp"
MS.read_sp_grm(grm_file=grm_file)
E = MS.get_env_int()
y = MS.get_y(adjust=True)

G = MS.get_genotype(bedfile, [snp], scale=True)
effect = np.linalg.lstsq(G, y, rcond=None)[0]
y = y - G @ effect

np.savetxt(f"{trait}.{snp}.main_effect.txt", effect)


varcom = np.loadtxt(varcom_file, dtype=float)[:, 0]
MS.cal_spVi(varcom)

logging.info(f"G.shape = {G.shape}, E.shape = {E.shape}, y.shape = {y.shape}")
EG = E * G
res = MS.mmsusie(EG, y, L=10)
MS.out_mmsusie(res, f"{trait}.{snp}.mmsusie")

pip, alpha, mu, sets = susie(EG, y)

np.savetxt(f"{trait}.{snp}.susie.pip.txt", pip)
np.savetxt(f"{trait}.{snp}.susie.alpha.txt", alpha)
np.savetxt(f"{trait}.{snp}.susie.mu.txt", mu)
