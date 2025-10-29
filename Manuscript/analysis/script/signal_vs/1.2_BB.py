import sys
import pandas as pd
import os
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

trait_index = int(sys.argv[1])  # Trait index from command line

# Directories
base_dir = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/"
merge_dir = os.path.join(base_dir, "testGxE_merge")

# Read traits
trait_file = os.path.join(base_dir, "trait67_shortName.txt")
traits = pd.read_csv(trait_file, sep="\t").iloc[:, 0].tolist()
print(f"The number of traits: {len(traits)}")

# Columns to read
cols_to_use = [
    "chrom", "SNP", "base", "allele1", "allele2", "af",
    "p_gxe", "p_gxe_liu",  "p_SGEM_O"
]
pval_cols = [c for c in cols_to_use if c.startswith("p_")]
print(pval_cols)

def process_trait(trait):
    file_path = os.path.join(merge_dir, f"{trait}.gz")
    df = pd.read_csv(file_path, sep=r"\s+", usecols=cols_to_use)
    df["trait"] = trait
    sig_df = df[df[pval_cols].min(axis=1) < 5e-8].copy()
    print("Processing", trait, len(sig_df))
    return sig_df

# Processing
trait = traits[trait_index]
res = process_trait(trait)

# Combine and save
output_file = f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/signal/{trait}.signal.txt"
final_df = pd.concat([res], ignore_index=True)
final_df.to_csv(output_file, sep="\t", index=False, na_rep="NA")
