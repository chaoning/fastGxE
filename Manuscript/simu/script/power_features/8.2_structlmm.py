import sys
sys.path.append('/net/zootopia/disk1/chaon/WORK/structLMM/')
import os
from structLMM import structLMM_null, integer_divide, structLMM_mme_fixed
import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

prefix = sys.argv[1]
sample_size = int(sys.argv[2])
rep = sys.argv[3]


os.chdir(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/pheno_baselineLD_fastgxe/")

logging.info("Plink file")
bed_file = f"../pheno_baselineLD/{prefix}.rep_{rep}.power_snps"
snp_on_disk = Bed(bed_file, count_A1=False)


dfFam = pd.read_csv(bed_file + ".fam", sep=r"\s+", header=None)
id_dct = dict(zip(dfFam[1], range(dfFam.shape[0])))
dfBim = pd.read_csv(bed_file + ".bim", sep=r"\s+", header=None)


logging.info("Pheno file")
pfile = f"{prefix}.txt"
dfp = pd.read_csv(pfile, sep=r"\s+")

logging.info("Remove related samples")
df_idKeep = pd.read_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/random_{sample_size}_unRel.agrm.id",
                        sep=r"\s+", header=None, names=["eid"], usecols=[1])

dfp = pd.merge(df_idKeep, dfp, on="eid")
logging.info(f"Keep {dfp.shape[0]} samples")


logging.info("Read geno")
id_index_keep_lst = []
for val in dfp["eid"]:
    id_index_keep_lst.append(id_dct[val])

logging.info("Data")
y = np.array(dfp.loc[:, [f"trait{rep}"]]).reshape(-1, 1)
X_sm = np.ones((y.size, 1))
beta = np.dot(np.linalg.inv(np.dot(X_sm.transpose(), X_sm)), np.dot(X_sm.transpose(), y))
yadj = y - np.dot(X_sm, beta)
E = np.array(dfp.loc[:, "age":"alcohol_frequency"])

varcom = structLMM_null(y, X_sm, E)
logging.info(f"Estimated variances: {varcom}")


print(len(id_index_keep_lst), "samples selected for analysis")

subset_on_disk = snp_on_disk[id_index_keep_lst, :]
print(snp_on_disk.sid_count, "SNPs in the file")
print(snp_on_disk.iid_count, "samples in the file")
print(subset_on_disk.sid_count, "SNPs selected for analysis")
print(subset_on_disk.iid_count, "samples selected for analysis")


snpData = pd.DataFrame(subset_on_disk.read().val)
mean_value = snpData.mean()
snpData.fillna(value=mean_value, inplace=True)
snpData = snpData - np.array(mean_value).reshape(1, -1)

num_snp = snpData.shape[1]
p_liu_arr = np.ones(num_snp)
p_sad_arr = np.ones(num_snp)
score_arr = np.ones(num_snp)

for i in range(num_snp):
    G = snpData.iloc[:, i].copy()
    _, score, p_liu, p_sad = structLMM_mme_fixed(yadj, None, G, E, varcom)
    score_arr[i] = score
    p_liu_arr[i] = p_liu
    p_sad_arr[i] = p_sad

dfR = pd.DataFrame({
    "score": score_arr,
    "p_liu": p_liu_arr,
    "p_sad": p_sad_arr
})


dfSNP = dfBim.iloc[:, [0, 1, 3]]
dfSNP.columns = ["chrom", "SNP", "base"]
dfR.index = range(dfR.shape[0])
dfSNP.index = range(dfSNP.shape[0])
dfM = pd.concat([dfSNP, dfR], axis=1)
dfM.to_csv(f"../res_baselineLD_structlmm/{prefix}.{rep}.csv", index=False)
