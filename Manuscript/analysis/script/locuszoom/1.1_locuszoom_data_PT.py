import sys
import os
import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed
from tqdm import tqdm
import re

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/locuszoom/data/")
# os.chdir(sys.argv[1])

bed_file = "/net/zootopia/disk1/chaon/data/UKB/imp/ukb_imp_info.qc"
dfBim = pd.read_csv(bed_file + ".bim.hg38", sep="\s+", header=None)
dfBim = dfBim.iloc[:, [0, 1, 3]]
dfBim.columns = ["chrom", "SNP", "base"]
snp_on_disk = Bed(bed_file, count_A1=False)

dfD = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.addP_logP.csv")

k = int(sys.argv[1]) - 1
snp = dfD.loc[k, "trait_leading_snp"]
chr = dfD.loc[k, 'chrom_x']
loci = dfD.loc[k, 'GenomicLocus_hg38']
start = dfD.loc[k, 'start_hg38'] - 1000000
end = dfD.loc[k, 'end_hg38'] + 1000000
trait = dfD.loc[k, 'trait']
# FieldName = dfD.iloc[k, -1]
gwas_file = f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/merge/{trait}.gz"
out_file = f"PT.{loci}.{trait}.{snp}"
dfGWAS = pd.read_csv(gwas_file, sep="\s+", usecols=["SNP", "p_gxe"])
dfGWAS = pd.merge(dfBim, dfGWAS, on="SNP")
dfGWAS = dfGWAS[(dfGWAS["chrom"] == chr) & (dfGWAS["base"] >= start) & (dfGWAS["base"] <= end)]
dfGWAS.to_csv(f"{out_file}.GWAS.txt", sep=" ", index=False)


snp_lst = list(dfGWAS["SNP"])
snp_index_lst = []
snp_index_lst.extend(snp_on_disk.sid_to_index([snp]))
snp_index_lst.extend(snp_on_disk.sid_to_index(snp_lst))
snp_data = snp_on_disk[:, snp_index_lst].read().val
snp_data = pd.DataFrame(snp_data)
corr_lst = []
for j in tqdm(range(1, snp_data.shape[1])):
    r2 = (snp_data.iloc[:, [0, j]]).corr().iloc[0, 1]
    r2 = r2 * r2
    corr_lst.append(r2)
res_dct = {"SNP1": [snp] * len(corr_lst),
           "SNP2": snp_lst,
           "r2": corr_lst}
df = pd.DataFrame(res_dct)
df.to_csv(f"{out_file}.r2.txt", sep=" ", index=False)
