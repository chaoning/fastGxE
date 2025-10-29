import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed

# Load genotype data
bed_file = "/net/zootopia/disk1/chaon/data/UKB/imp/ukb_imp_info.qc.100000"
snp_on_disk = Bed(bed_file, count_A1=True)
bim_df = pd.read_csv(bed_file + ".bim", sep="\t", header=None)
bim_df.columns = ["chr", "snp", "cm", "pos", "a1", "a2"]

# Read SNP lists
GRL_new = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.addP.csv")
trait_lead_snp_new = set(GRL_new["trait_leading_snp"])

GRL_old = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/Biological_samples/result/SGEM/FUMA/IINT/GenomicRiskLoci.hg38.keygene.TraitSpecific.Methods.addP.csv")
trait_lead_snp_old = set(GRL_old["trait_leading_snp"])

# Map SNP to index in .bed file
snp_to_idx = {snp: i for i, snp in enumerate(bim_df["snp"])}

# Keep only SNPs present in the genotype file
new_indices = [snp_to_idx[s] for s in trait_lead_snp_new]
old_indices = [snp_to_idx[s] for s in trait_lead_snp_old]

# Read genotypes for these SNPs
geno_new = snp_on_disk[:, new_indices].read().val
geno_old = snp_on_disk[:, old_indices].read().val

# Center genotypes (mean=0) for correlation
geno_new -= np.nanmean(geno_new, axis=0)
geno_old -= np.nanmean(geno_old, axis=0)

# Replace NaN with 0 for correlation calculation
geno_new = np.nan_to_num(geno_new)
geno_old = np.nan_to_num(geno_old)

# Compute pairwise LD (r²)
ld_pairs = []
for i, idx_new in enumerate(new_indices):
    for j, idx_old in enumerate(old_indices):
        r = np.corrcoef(geno_new[:, i], geno_old[:, j])[0, 1]
        if r**2 > 0.1:
            ld_pairs.append((bim_df.loc[idx_new, "snp"], bim_df.loc[idx_old, "snp"], r**2))

# Convert to DataFrame
ld_df = pd.DataFrame(ld_pairs, columns=["SNP_new", "SNP_old", "r2"])
ld_df1 = (ld_df.sort_values("r2", ascending=False)
                      .drop_duplicates(subset=["SNP_new"], keep="first"))
ld_df2 = (ld_df.sort_values("r2", ascending=False)
                         .drop_duplicates(subset=["SNP_old"], keep="first"))
GRL_new = pd.merge(GRL_new, ld_df1, left_on="trait_leading_snp", right_on="SNP_new", how="left")
GRL_new.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.addP.LD.csv", index=False)
GRL_old = pd.merge(GRL_old, ld_df2, left_on="trait_leading_snp", right_on="SNP_old", how="left")
GRL_old.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.addP.LD.old.csv", index=False)
