import pandas as pd
import numpy as np
from typing import List
from pysnptools.snpreader import Bed
import sys
import os

def ld_clumping(snps: List[str], p_values: List[float], ld_matrix: np.ndarray, r2_threshold: float = 0.1) -> List[str]:
    """
    Perform LD clumping to select independent lead SNPs based on r^2 threshold.

    Parameters:
    - snps: List of SNP IDs.
    - p_values: List of corresponding p-values for SNPs.
    - ld_matrix: Symmetric LD matrix (r) of shape [n_snps, n_snps].
    - r2_threshold: Threshold for LD r^2 to consider SNPs as dependent (default: 0.1).

    Returns:
    - List of lead SNPs that are approximately independent.
    """
    if len(snps) != len(p_values) or ld_matrix.shape[0] != ld_matrix.shape[1] or ld_matrix.shape[0] != len(snps):
        raise ValueError("Input dimensions mismatch.")

    df = pd.DataFrame({'SNP': snps, 'P': p_values})
    df = df.sort_values('P').reset_index(drop=True)

    selected_snps = []
    excluded = set()

    for _, row in df.iterrows():
        snp_id = row['SNP']
        if snp_id in excluded:
            continue
        selected_snps.append(snp_id)
        idx = snps.index(snp_id)
        ld_r2 = ld_matrix[idx, :] ** 2
        for j, r2 in enumerate(ld_r2):
            if r2 >= r2_threshold:
                excluded.add(snps[j])
    return selected_snps

def cal_ld_matrix(snps: List[str], bed_file: str) -> np.ndarray:
    """
    Calculate the LD matrix from a BED file containing SNP data.
    Parameters:
    - snps: List of SNP IDs.
    - bed_file: Path to the BED file containing SNP data.
    Returns:
    - LD matrix (r) of shape [n_snps, n_snps].
    """
    snp_on_disk = Bed(bed_file, count_A1=True)
    dfBim = pd.read_csv(bed_file + ".bim", sep=r"\s+", header=None)
    dct = dict(zip(dfBim[1], range(len(dfBim))))
    snp_indices = [dct[snp] for snp in snps]
    snp_data = snp_on_disk[:, snp_indices].read().val

    ld_r = pd.DataFrame(snp_data, columns=snps).corr().values
    return ld_r


os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_overfit/")
prefix = sys.argv[1]

res_file = prefix + ".res"

if os.path.exists(res_file):
    df = pd.read_csv(res_file, sep=r"\s+")
    df.dropna(inplace=True)
    df = df[df['p_gxe'] < 5e-8]
    print(f"Found {len(df)} significant SNPs.")
    if df.empty:
        print("No significant SNPs found.")
        open(f"../res_overfit_leadsnp/{prefix}.lead_snps.txt", "w").close()
        sys.exit(0)
    snps = df['SNP'].tolist()
    p_values = df['p_gxe'].tolist()
    power_snp_name_file = f"../pheno_take_home/{prefix}.power_snp_name.txt"
    power_snp_name = pd.read_csv(power_snp_name_file, sep=r"\s+", header=None).iloc[0, 0]

    bed_file = "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22"
    ld_matrix = cal_ld_matrix(snps, bed_file)

    lead_snps = ld_clumping(snps, p_values, ld_matrix, r2_threshold=0.1)

    dfBim = pd.read_csv(bed_file + ".bim", sep=r"\s+", header=None)
    dct = dict(zip(dfBim[1], dfBim[3]))
    power_snp_base = dct[power_snp_name]
    base_diff = [dct[snp] - power_snp_base for snp in lead_snps if snp in dct]

    ld_matrix2 = cal_ld_matrix(lead_snps + [power_snp_name], bed_file)
    lead_snps_df = pd.DataFrame({'SNP': lead_snps, "ld_r2": ld_matrix2[-1, :-1]**2, "base_diff": base_diff})
    lead_snps_df.to_csv(f"../res_overfit_leadsnp/{prefix}.lead_snps.txt", sep="\t",
                        index=False, header=True)
else:
    print(f"Result file {res_file} does not exist.")
    sys.exit(1)
