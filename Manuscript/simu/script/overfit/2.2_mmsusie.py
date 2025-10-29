import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed
from mmsusie import MMSuSiE
import logging
import sys
import random
sys.path.append('/net/zootopia/disk1/chaon/WORK/structLMM/')
import os

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_overfit_leadsnp_mmsusie/")

prefix = sys.argv[1]
rep = int(sys.argv[2])

# Set up logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s') 


MS = MMSuSiE()


bedfile = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22"
data_file = f"../pheno_take_home_fastgxe/{prefix}.txt"
head_names = pd.read_csv(data_file, sep=r"\s+", nrows=0).columns.tolist()
env_int = head_names[1:41]

logging.info(f"Processing replicate {rep}")
MS = MMSuSiE()

MS.read_data(data_file, trait=f"trait{rep}", env_int=env_int)

grm_file = "/net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp"
MS.read_sp_grm(grm_file=grm_file)
E = MS.get_env_int()
y = MS.get_y(adjust=True)

lead_snp_file = f"../res_overfit_leadsnp/{prefix}.rep_{rep}.lead_snps.txt"
lead_snps = pd.read_csv(lead_snp_file, sep=r"\s+").iloc[:, 0].tolist()
if len(lead_snps) == 0:
    logging.warning(f"No lead SNPs found. Exiting.")
    sys.exit()
varcom = np.loadtxt(f"../res_overfit/{prefix}.rep_{rep}.var", dtype=float)[:, 0]
MS.cal_spVi(varcom)
logging.info(f"Vi.shape = {MS.Vi.shape}")

for lead_snp in lead_snps:
    G = MS.get_genotype(bedfile, [lead_snp], scale=True)
    y = y - G @ np.linalg.lstsq(G, y, rcond=None)[0]
    logging.info(f"G.shape = {G.shape}, E.shape = {E.shape}, y.shape = {y.shape}")
    EG = E * G
    logging.info(f"EG.shape = {EG.shape}, y.shape = {y.shape}")
    res = MS.mmsusie(EG, y, L=10)
    MS.out_mmsusie(res, f"{prefix}.rep_{rep}.{lead_snp}.mmsusie_out")
