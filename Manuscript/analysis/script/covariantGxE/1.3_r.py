import pandas as pd
import numpy as np
import sys
import os
from pysnptools.snpreader import Bed

trait = sys.argv[1]
assoc_dir = sys.argv[2]
pheno_dir = sys.argv[3]

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/plink/")
file = f"{trait}.clumped.clumped"

df = pd.read_csv(file, sep=r"\s+")
dfSNP = df.loc[:, ["SNP"]].copy()
snp_lst = df["SNP"].tolist()

bed_file = "/net/zootopia/disk1/chaon/data/UKB/imp/ukb_imp_info.qc"
dfBim = pd.read_csv(f"{bed_file}.bim", sep=r"\s+", header=None)
dct = dict(zip(dfBim[1], range(len(dfBim))))

snp_index_lst = [dct[snp] for snp in snp_lst]


os.chdir(f"{assoc_dir}")
def read_assoc(input_file_prefix):
    # Read the input file
    df_lst = []
    for i in range(10):
        input_file = f"{input_file_prefix}.10_{i+1}.res"
        df = pd.read_csv(input_file, sep=r"\s+", usecols=["SNP", "beta_main"] + [f"beta{j+1}" for j in range(42)])
        df_lst.append(df)
    # Concatenate all DataFrames
    combined_df = pd.concat(df_lst, ignore_index=True)
    return combined_df

dfassoc = read_assoc(f"{trait}")

dfM = pd.merge(dfSNP, dfassoc, left_on="SNP", right_on="SNP", how="left")

beta_col_lst = ["beta_main"] + [f"beta{i+1}" for i in range(42)]

beta_mat = dfM[beta_col_lst].values

os.chdir(pheno_dir)

dfE = pd.read_csv(f"pheno.e42.{trait}.txt", sep=r"\s+")
envi_mat = dfE.iloc[:, 2:].values

iid_lst = dfE.iloc[:, 0].tolist()

dfFam = pd.read_csv(bed_file + ".fam", sep=r"\s+", header=None)
iid_dct = dict(zip(dfFam[1], range(len(dfFam))))
iid_index_lst = [iid_dct[iid] for iid in iid_lst]

snp_on_disk = Bed(bed_file, count_A1=True)

snp_mat = snp_on_disk[iid_index_lst, snp_index_lst].read().val
snp_mat = pd.DataFrame(snp_mat)
snp_mat.fillna(snp_mat.mean(), inplace=True)

G_sum = np.dot(snp_mat.values, beta_mat[:, 0])
var1_sum = np.sum(np.var(snp_mat.values * beta_mat[:, 0].reshape(1, -1), axis=0))

GxE_sum_allE = 0
num_snps = len(snp_index_lst)
var2_sum = 0
for i in range(envi_mat.shape[1]):
    GxE_sum =  np.dot(snp_mat.values, beta_mat[:, i + 1])
    GxE_sum_allE += GxE_sum
    r = np.corrcoef(G_sum, GxE_sum)[0, 1]
    var2_sum += np.sum(np.var(snp_mat.values * beta_mat[:, i + 1].reshape(1, -1), axis=0))


os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/plink/")

df_res_allE = pd.DataFrame([[np.corrcoef(G_sum, GxE_sum_allE)[0, 1], 
                             np.var(G_sum) / num_snps, 
                             np.var(GxE_sum_allE) / num_snps, 
                             num_snps,
                             var1_sum / num_snps,
                             var2_sum / num_snps]], 
                           columns=["r", "varG", "varGxE", "num_snps", "varG_sum", "varGxE_sum"])
df_res_allE.to_csv(f"{trait}.cov_allE.csv", index=False)

