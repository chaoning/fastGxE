import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed
from mmsusie import MMSuSiE
import logging
import sys
import os
sys.path.append('/net/zootopia/disk1/chaon/WORK/structLMM/')

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power")
# Set up logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s') 
power_snp_h2_gxe = sys.argv[1]
num_envi = sys.argv[2]
rep = int(sys.argv[3])

MS = MMSuSiE()

prefix = "h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000"

bedfile = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22"

data_file = f"./pheno_fastgxe/power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{prefix}.txt"
head_names = pd.read_csv(data_file, sep=r"\s+", nrows=0).columns.tolist()
env_int = head_names[1:41]

logging.info(f"Processing replicate {rep}")
MS = MMSuSiE()

MS.read_data(data_file, trait=f"trait{rep}", env_int=env_int)

grm_file = "/net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp"
MS.read_sp_grm(grm_file=grm_file)
E = MS.get_env_int()
y = MS.get_y(adjust=True)

snp_file = f"./pheno/power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{prefix}.rep_{rep}.power_snp.bim"
lead_snp = pd.read_csv(snp_file, sep=r"\s+", header=None).iloc[0, 1]

G = MS.get_genotype(bedfile, [lead_snp], scale=True)
y = y - G @ np.linalg.lstsq(G, y, rcond=None)[0]

varcom = np.loadtxt(f"./res_fastgxe_testgxe/power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{prefix}.rep_{rep}.var", dtype=float)[:, 0]
MS.cal_spVi(varcom)

logging.info(f"G.shape = {G.shape}, E.shape = {E.shape}, y.shape = {y.shape}")
EG = E * G
res = MS.mmsusie(EG, y, L=10)
MS.out_mmsusie(res, f"./res_mmSuSiE/power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi}.{rep}.mmsusie_out")
