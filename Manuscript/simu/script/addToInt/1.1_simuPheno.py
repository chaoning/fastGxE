import pandas as pd
import numpy as np
import sys
from pysnptools.snpreader import Bed
import os
import logging
import gc
from scipy.stats import pearsonr

def process_snp_data(order_usedID_in_bed, snp_indices, snp_on_disk):
    snp_data = snp_on_disk[order_usedID_in_bed, snp_indices].read().val
    snp_data = pd.DataFrame(snp_data).fillna(pd.DataFrame(snp_data).mean())
    snp_data = (snp_data - snp_data.mean()) / snp_data.std()
    snp_data = snp_data.values  # Convert DataFrame to numpy array
    return snp_data



def simulate_snp_data(snp_A, r_squared_target, n_iterations=1000000):
    """
    Simulates SNP A and SNP B with a target r² value.
    
    Parameters:
    - snp_A: Array-like, genotype data for SNP A.
    - r_squared_target: Target r² value for SNP A and SNP B.
    - n_iterations: Number of iterations to adjust SNP B to achieve target r².
    
    Returns:
    - Final r² value after simulation.
    - SNP B data.
    """
    
    # Convert SNP A to continuous data and standardize
    snp_A_continuous = snp_A.astype(float)
    snp_A_std = np.std(snp_A_continuous)
    snp_A_continuous = snp_A_continuous / snp_A_std  # Standardize SNP A
    
    r_squared_actual = r_squared_target  # Start with target r² value
    n_samples = len(snp_A)  # Get the number of samples

    for n in range(1, n_iterations + 1):
        # Calculate correlation from target r²
        correlation = np.sqrt(r_squared_actual)

        # Simulate SNP B with noise and correlation to SNP A
        snp_B_continuous = correlation * snp_A_continuous + np.random.normal(0, 1, n_samples) * np.sqrt(1 - r_squared_actual)

        # Convert SNP B to 0, 1, 2 genotypes
        snp_B = np.clip(np.round(snp_B_continuous), 0, 2).astype(int)

        # Calculate the Pearson correlation and compute r²
        corr, _ = pearsonr(snp_A, snp_B)
        r_squared_actual = corr ** 2

        print(f"Iteration {n}, r² = {r_squared_actual:.4f}")

        # Adjust r² if it diverges from the target
        if r_squared_actual < r_squared_target - 0.0001:
            r_squared_actual = r_squared_target + 0.0001 * n  # Increase r² if too low
        elif r_squared_actual > r_squared_target + 0.0001:
            r_squared_actual = r_squared_target - 0.0001 * n  # Decrease r² if too high
        else:
            break  # Stop if r² is within a small tolerance of the target
    
    maf = np.mean(snp_B) / 2  # Calculate minor allele frequency
    if maf > 0.5:
        snp_B = 2 - snp_B # Adjust if MAF is greater than 0.5 to ensure correct allele coding
    return r_squared_actual, snp_B

# Set working directory
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_addToint/")

# Read command line arguments
h2_envi = 0.05
unknown_h2SNP_h2E = float(sys.argv[1]) / 100
R2 = float(sys.argv[2])  # squared correlation coefficient between G and unknown G, E and unknown E
R2GE = float(sys.argv[3])  # squared correlation coefficient between G and E

h2_add = float(sys.argv[4]) / 100
h2_gxe = float(sys.argv[5]) / 100
h2_nxe = float(sys.argv[6]) / 100

sample_size = int(sys.argv[7])
rep = int(sys.argv[8])

if unknown_h2SNP_h2E > h2_add or unknown_h2SNP_h2E > h2_envi:
    logging.error("Unknown heritability cannot be greater than additive or environmental heritability.")
    sys.exit(1)

if h2_envi + h2_add + h2_gxe + h2_nxe >= 1:
    logging.error("Total heritability (additive + GxE + NxE) cannot be greater than 1.")
    sys.exit(1)



# Set random seed
np.random.seed(rep)

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
nSNP = dfBim.shape[0]
envi_names = dfP.loc[:, "age":"alcohol_frequency"].columns  # get environment names


# Read SNP data from BED file
dfFam = pd.read_csv(bed_file + ".fam", sep=r"\s+", header=None, dtype=str)
id_dct = dict(zip(dfFam.iloc[:, 1], range(dfFam.shape[0])))
order_usedID_in_bed = [id_dct[eid] for eid in dfP["eid"]]
snp_on_disk = Bed(bed_file, count_A1=False)


print("### additive effects for unknown G")

power_snp_index = np.random.randint(nSNP_exclude22, nSNP)

snp_data = process_snp_data(order_usedID_in_bed, [power_snp_index], snp_on_disk)
snp_data = snp_data[:, 0]

# Generate a random correlated G
snp_data_unknown = np.sqrt(R2) * snp_data + np.random.normal(0, 1, snp_data.size) * np.sqrt(1 - R2)
print(f"The correlation coefficient for G and the unknown G is {np.corrcoef(snp_data, snp_data_unknown)[0, 1]}")
snp_data_unknown = (snp_data_unknown - snp_data_unknown.mean()) / snp_data_unknown.std()
eff = np.sqrt(unknown_h2SNP_h2E)
add_power_sum = snp_data_unknown * eff


print("### Generate correlated environments")
randomEForSimu_index = np.random.randint(40)
known_envi_name = envi_names[randomEForSimu_index]
envi_values = dfP[known_envi_name].values
envi_values = (envi_values - envi_values.mean()) / envi_values.std()
# Generate a random correlated environment
correlated_envi_arr = np.sqrt(R2GE) * snp_data_unknown + np.sqrt(R2) * envi_values + np.sqrt(1 - R2 - R2GE) * np.random.normal(size=snp_data.size)
correlated_envi_arr = (correlated_envi_arr - correlated_envi_arr.mean()) / correlated_envi_arr.std()
dfcorr = pd.DataFrame({
    "known_G": snp_data,
    "unknown_G": snp_data_unknown,
    "known_E": envi_values,
    "unknown_E": correlated_envi_arr
})
print(f"The squared correlation coefficients: {dfcorr.corr() * dfcorr.corr()}")
eff = np.sqrt(unknown_h2SNP_h2E)
unknown_E_sum = correlated_envi_arr * eff

print("### Generate environmental effects")
dataE = np.array(dfP.loc[:, "age":"alcohol_frequency"].copy())
num_envi = dataE.shape[1]
eff_E = np.random.normal(size=num_envi)
E_var = np.sum(np.var(dataE * eff_E.reshape(1, -1), axis=0))
eff_E *= np.sqrt((h2_envi - unknown_h2SNP_h2E) / E_var)
E_sum = np.dot(dataE, eff_E)




print("Generate additive SNP effects")
add_snp_indices = np.random.choice(range(nSNP_exclude22), size=num_add_SNP, replace=False)
add_snp_indices = np.sort(add_snp_indices)
snp_data_add = process_snp_data(order_usedID_in_bed, add_snp_indices, snp_on_disk)
eff_add = np.random.normal(size=num_add_SNP)
add_var = np.sum(np.var(snp_data_add * eff_add.reshape(1, -1), axis=0))
eff_add *= np.sqrt((h2_add - unknown_h2SNP_h2E) / add_var)
add_sum = np.dot(snp_data_add, eff_add)
del snp_data_add  # Free memory
gc.collect()

print("Generate GxE effects")
gxe_snp_indices = np.random.choice(range(nSNP_exclude22), size=num_gxe_SNP, replace=False)
gxe_snp_indices = np.sort(gxe_snp_indices)
numbers = list(range(40))
weights = [10]*10 + [1]*30
lst1 = []
lst2 = []
for isnp in range(num_gxe_SNP):
    chosen_number = np.random.choice(numbers, p=np.array(weights) / sum(weights))
    randomE_index = np.random.choice(range(num_envi), chosen_number, replace=False)
    randomE_index = np.sort(randomE_index)
    lst1.extend([gxe_snp_indices[isnp]] * len(randomE_index))
    lst2.extend(randomE_index)

num_gxe_pair = len(lst1)
snp_data_gxe = process_snp_data(order_usedID_in_bed, lst1, snp_on_disk)
eff_GxE = np.random.normal(size=num_gxe_pair)
GxE_var = np.sum(np.var(snp_data_gxe * dataE[:, lst2] * eff_GxE.reshape(1, -1), axis=0))
eff_GxE *= np.sqrt(h2_gxe / GxE_var)
GxE_sum = np.dot(snp_data_gxe * dataE[:, lst2], eff_GxE)
del snp_data_gxe  # Free memory
gc.collect()

print("Generate NxE effects")
nxe_sum = 0
for i in range(num_envi):
    dataEi = dataE[:, i].copy()
    nxei = dataEi * np.random.normal(size=num_ID) * np.sqrt(h2_nxe / num_envi)
    nxe_sum += nxei

error = np.random.normal(size=num_ID) * np.sqrt(1 - h2_envi - h2_add - h2_gxe - h2_nxe)
y = add_power_sum + unknown_E_sum + add_sum + GxE_sum + error + E_sum + nxe_sum

dfR = dfP.loc[:, ["fid", "eid"]].copy()
dfR["pheno"] = y

print("Saving results")
out_prefix = f"h2_{sys.argv[1]}.R2_{sys.argv[2]}.R2GE_{sys.argv[3]}.h2_add_{sys.argv[4]}.h2_gxe_{sys.argv[5]}.h2_nxe_{sys.argv[6]}.sample_{sample_size}.rep_{rep}"
dfR.to_csv(f"{out_prefix}.txt", sep=" ", index=False)


np.savetxt(f"{out_prefix}.power_snp_index.txt", [power_snp_index - nSNP_exclude22], fmt="%d")
snpdata = Bed(bed_file)[order_usedID_in_bed, [power_snp_index]].read()
Bed.write(f"{out_prefix}.power_snp.bed",snpdata,count_A1=False)
dfBim.iloc[[power_snp_index], :].to_csv(f"{out_prefix}.power_snp.bim", sep=" ", header=False, index=False)
dfFam.iloc[order_usedID_in_bed, :].to_csv(f"{out_prefix}.power_snp.fam", sep=" ", header=False, index=False)

print("Results saved successfully.")
