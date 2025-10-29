import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed
import sys
from tqdm import tqdm


# === Load input files ===
# Summary data with trait-specific GxE beta estimates
dfS = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/AIE_split/3.2.3_sum.PT_trait_specific.csv", usecols=["GenomicLocus_hg38", "Genes", "trait", "trait_leading_snp"])
dfS.rename(columns={"trait_leading_snp": "SNP"}, inplace=True)

file = f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/signal/signal.txt"
summary_df = pd.read_csv(file, sep=r"\s+")

summary_df = pd.merge(dfS, summary_df, on=["trait", "SNP"])

print("# of GxE signals:", summary_df.shape[0])


# Load SNP information (bim file) and create a SNP-to-index dictionary
bed_prefix = "/net/zootopia/disk1/chaon/data/UKB/imp/ukb_imp_info.qc"
bim_df = pd.read_csv(bed_prefix + ".bim", sep=r"\s+", header=None)
snp_index_map = dict(zip(bim_df[1], range(bim_df.shape[0])))

# Load .bed genotype file using PySnpTools
snp_reader = Bed(bed_prefix, count_A1=True)

# sample size

file = f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/Physical_measures/result/sum/sumSampleAssoc/sampleSize.txt"

sample_size_df = pd.read_csv(file, sep=r"\s+")
dct_sample_size = dict(zip(sample_size_df["trait"], sample_size_df["sample_size"]))

# === Initialize result lists ===
h2gxe_list = []
sample_size_list = []

# Column names of beta estimates for 41 environments
beta_columns = [f"beta{i}" for i in range(1, 42)]

# === Loop through each trait-lead SNP pair ===
for idx, row in tqdm(summary_df.iterrows(), total=summary_df.shape[0]):
    trait = row["trait"]

    # Load individual-level phenotype data to get sample size
    sample_size = dct_sample_size[trait]
    sample_size_list.append(sample_size)

    # Get genotype index of the lead SNP
    snp_id = row["SNP"]
    snp_idx = snp_index_map[snp_id]

    # Read genotype for the SNP (n_samples x 1)
    genotype = snp_reader[:, [snp_idx]].read().val
    genotype_df = pd.DataFrame(genotype)

    # Calculate SNP variance using allele frequency: Var = 2p(1-p)
    allele_freq = genotype_df.mean().values[0] / 2
    snp_variance = 2 * allele_freq * (1 - allele_freq)

    # Extract beta estimates and compute GxE heritability
    beta_vector = row[beta_columns].values
    h2gxe = np.sum(snp_variance * beta_vector**2)
    h2gxe_list.append(h2gxe)

# === Add results back to dataframe and save ===
summary_df["h2gxe"] = h2gxe_list
summary_df["sample_size"] = sample_size_list

print(
    "h2gxe:",
    round(np.min(h2gxe_list) * 100, 3),
    round(np.max(h2gxe_list) * 100, 3),
    round(np.median(h2gxe_list) * 100, 3)
)

print(
    "sample_size:",
    np.min(sample_size_list),
    np.max(sample_size_list),
    np.median(sample_size_list)
)

output_path = f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/GxE_signal_h2/1.1_real_data_PT.csv"
summary_df.to_csv(output_path, index=False)
