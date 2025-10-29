import pandas as pd
import numpy as np
import sys
from pysnptools.snpreader import Bed
import os
import logging
import gc
from tqdm import tqdm

dfBMI = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/ldscore/data/BMI.baselineLD.results", sep=r"\s+")
dfBaseline = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/baseline_model_2.2/annot_baseline2.2.txt", sep=r"\s+")
def cal_baseline_sigma2(snp_lst, dfBMI, dfBaseline):
    coeff = dfBMI["Coefficient"].values
    dfSNP = pd.DataFrame({"SNP": snp_lst})
    dfM = pd.merge(dfSNP, dfBaseline, on="SNP", how="left")
    annot_mat = dfM.iloc[:, 1:].fillna(0.0).to_numpy(dtype=float)
    sigma2 = annot_mat.dot(coeff)
    matched = annot_mat.shape[0] - dfM.iloc[:, 1:].isna().all(axis=1).sum()
    logging.info("Matched SNPs with any annotation: {}/{}".format(matched, len(snp_lst)))
    sigma2[sigma2 < 0] = 1.0e-100
    logging.info("Minimum sigma2: {:.2e}, Maximum sigma2: {:.2e}, Mean sigma2: {:.2e}".format(np.min(sigma2), np.max(sigma2), np.mean(sigma2)))
    return sigma2


# Set working directory
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_baselineLD/")

# Read command line arguments
h2_envi = 0.05
h2_add = float(sys.argv[1]) / 100
h2_gxe = float(sys.argv[2]) / 100
h2_nxe = float(sys.argv[3]) / 100

sample_size = int(sys.argv[4])
rep = int(sys.argv[5])

num_add_SNP = 5000
num_gxe_SNP = 500

# Define file paths
bed_file = "/net/zootopia/disk1/chaon/data/UKB/imp/ukb_imp_info.qc"
id_file = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/random_{sample_size}.fam"
dfID = pd.read_csv(id_file, sep=r"\s+", header=None, names=["fid", "eid"], dtype=str)
dfE = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/envi.csv", dtype={"eid": str})
dfP = pd.merge(dfID, dfE, on="eid", how="left")
num_ID = dfP.shape[0]
dfBim = pd.read_csv(bed_file + ".bim", sep=r"\s+", header=None)
nSNP_exclude22 = dfBim[dfBim[0] != 22].shape[0]



# Set random seed
np.random.seed(rep)

logging.info("Generate environmental effects")
dataE = np.array(dfP.loc[:, "age":"alcohol_frequency"].copy())
num_envi = dataE.shape[1]
eff_E = np.random.normal(size=num_envi)
E_var = np.sum(np.var(dataE * eff_E.reshape(1, -1), axis=0))
eff_E *= np.sqrt(h2_envi / E_var)
E_sum = np.dot(dataE, eff_E)

# Read SNP data from BED file
dfFam = pd.read_csv(bed_file + ".fam", sep=r"\s+", header=None, dtype=str)
id_dct = dict(zip(dfFam.iloc[:, 1], range(dfFam.shape[0])))
order_usedID_in_bed = [id_dct[eid] for eid in dfP["eid"]]
snp_on_disk = Bed(bed_file, count_A1=False)


def process_snp_data(order_usedID_in_bed, snp_indices, snp_on_disk):
    snp_data = snp_on_disk[order_usedID_in_bed, snp_indices].read().val
    snp_data = pd.DataFrame(snp_data).fillna(pd.DataFrame(snp_data).mean())
    snp_data = (snp_data - snp_data.mean()) / snp_data.std()
    return snp_data


logging.info("Generate additive SNP effects")
add_snp_indices = np.random.choice(range(nSNP_exclude22), size=num_add_SNP, replace=False)
add_snp_indices = np.sort(add_snp_indices)
snp_data_add = process_snp_data(order_usedID_in_bed, add_snp_indices, snp_on_disk)
snp_lst = dfBim.iloc[add_snp_indices, 1].tolist()
sigma2_vec = cal_baseline_sigma2(snp_lst, dfBMI, dfBaseline)
eff_add = np.array([np.random.normal(0, np.sqrt(sigma2)) for sigma2 in tqdm(sigma2_vec)])

add_var = np.sum(np.var(snp_data_add * eff_add.reshape(1, -1), axis=0))
eff_add *= np.sqrt(h2_add / add_var)
add_sum = np.dot(snp_data_add, eff_add)
del snp_data_add  # Free memory
gc.collect()

logging.info("Generate GxE effects")
gxe_snp_indices = np.random.choice(range(nSNP_exclude22), size=num_gxe_SNP, replace=False)
gxe_snp_indices = np.sort(gxe_snp_indices)
snp_lst = dfBim.iloc[gxe_snp_indices, 1].tolist()
sigma2_vec = cal_baseline_sigma2(snp_lst, dfBMI, dfBaseline)
numbers = list(range(40))
weights = [10]*10 + [1]*30
lst1 = []
lst2 = []
eff_GxE_lst = []
for isnp in tqdm(range(num_gxe_SNP)):
    chosen_number = np.random.choice(numbers, p=np.array(weights) / sum(weights))
    randomE_index = np.random.choice(range(num_envi), chosen_number, replace=False)
    randomE_index = np.sort(randomE_index)
    lst1.extend([gxe_snp_indices[isnp]] * len(randomE_index))
    lst2.extend(randomE_index)
    eff_GxE_lst.extend(np.random.normal(0, np.sqrt(sigma2_vec[isnp]), chosen_number))

num_gxe_pair = len(lst1)
snp_data_gxe = process_snp_data(order_usedID_in_bed, lst1, snp_on_disk)
eff_GxE = np.array(eff_GxE_lst)
GxE_var = np.sum(np.var(snp_data_gxe * dataE[:, lst2] * eff_GxE.reshape(1, -1), axis=0))
eff_GxE *= np.sqrt(h2_gxe / GxE_var)
GxE_sum = np.dot(snp_data_gxe * dataE[:, lst2], eff_GxE)
del snp_data_gxe  # Free memory
gc.collect()

logging.info("Generate NxE effects")
nxe_sum = 0
for i in range(num_envi):
    dataEi = dataE[:, i].copy()
    nxei = dataEi * np.random.normal(size=num_ID) * np.sqrt(h2_nxe / num_envi)
    nxe_sum += nxei

error = np.random.normal(size=num_ID) * np.sqrt(1 - h2_envi - h2_add - h2_gxe - h2_nxe)
y = add_sum + GxE_sum + error + E_sum + nxe_sum

dfR = dfP.loc[:, ["fid", "eid"]].copy()
dfR["pheno"] = y

out_prefix = f"h2_add_{sys.argv[1]}.h2_gxe_{sys.argv[2]}.h2_nxe_{sys.argv[3]}.sample_{sample_size}.rep_{rep}"
dfR.to_csv(f"{out_prefix}.txt", sep=" ", index=False)
