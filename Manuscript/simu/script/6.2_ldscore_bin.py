import pandas as pd
import numpy as np

# Load LD file and name the columns for clarity
print("Loading LD file...")
ld_file = "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink_ld.ld"
df = pd.read_csv(ld_file, sep=r"\s+", header=0)

# Load SNP list from the BIM file (2nd column contains SNP IDs)
print("Loading SNP list from BIM file...")
bim_file = "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22.bim"
snp_lst = pd.read_csv(bim_file, sep=r"\s+", header=None)[1].tolist()

# Compute LD score by summing R2 values where SNP is used as SNP_A
# This avoids double-counting and gives one-directional LD score
print("Calculating LD scores...")
ldscore_single = df.groupby("SNP_A")["R2"].sum()

# Create result DataFrame and assign LD scores to all SNPs
# SNPs not found in the LD file get an LD score of 0
print("Creating result DataFrame...")
dfR = pd.DataFrame({"SNP": snp_lst})
dfR["ldscore"] = dfR["SNP"].map(ldscore_single).fillna(0)

# Divide SNPs into 5 bins based on LD score quantiles (equal-frequency groups)
print("Dividing SNPs into LD score bins...")
dfR["ldscore_bin"] = pd.qcut(dfR["ldscore"], 5, labels=False)

# Output directory to save results
print("Saving results...")
out_dir = "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/QQplot_ldscore"
dfR.to_csv(f"{out_dir}/ldscore.txt", sep="\t", index=False)

# Save each LD score bin to a separate file and print summary
for i in range(5):
    dfR_bin = dfR[dfR["ldscore_bin"] == i]
    min_val = dfR_bin["ldscore"].min()
    max_val = dfR_bin["ldscore"].max()
    dfR_bin.to_csv(f"{out_dir}/ldscore_bin_{i+1}.txt", sep="\t", index=False)
    print(f"Group {i+1}: {len(dfR_bin)} SNPs, LD score range: {min_val:.4f} - {max_val:.4f}")
