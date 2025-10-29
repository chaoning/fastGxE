import pandas as pd
import numpy as np
import sys
from pysnptools.snpreader import Bed
import os
import logging
import gc


# Set working directory
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/pheno/")

# Read command line arguments
h2_envi = 0.05
power_snp_index = int(sys.argv[1])
power_snp_h2_gxe = float(sys.argv[2]) / 100
num_envi_power_snp_gxe = int(sys.argv[3])  # number of environments for GxE SNP

h2_add = float(sys.argv[4]) / 100
h2_gxe = float(sys.argv[5]) / 100
h2_nxe = float(sys.argv[6]) / 100

if power_snp_h2_gxe > h2_gxe:
    logging.error("SNP heritability for GxE cannot be greater than total GxE heritability.")
    sys.exit(1)

if h2_envi + h2_add + h2_gxe + h2_nxe >= 1:
    logging.error("Total heritability (additive + GxE + NxE) cannot be greater than 1.")
    sys.exit(1)

sample_size = int(sys.argv[7])
rep = int(sys.argv[8])

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


print("Generate environmental effects")
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
    snp_data = snp_data.values  # Convert DataFrame to numpy array
    return snp_data

print("SNP GxE effects for power analysis")
if power_snp_index < 0:
    power_snp_index = np.random.randint(nSNP_exclude22, nSNP)

snp_data = process_snp_data(order_usedID_in_bed, [power_snp_index], snp_on_disk)
power_envi_indices = np.random.choice(range(num_envi), num_envi_power_snp_gxe, replace=False)
power_envi_indices = np.sort(power_envi_indices)
dataE_power = dataE[:, power_envi_indices]
eff = np.random.normal(size=num_envi_power_snp_gxe)

eff_var = np.sum(np.var(snp_data * dataE_power * eff.reshape(1, -1), axis=0))
eff *= np.sqrt(power_snp_h2_gxe / eff_var)
GxE_power_sum = np.dot(snp_data * dataE_power, eff)

print("Generate additive SNP effects")
add_snp_indices = np.random.choice(range(nSNP_exclude22), size=num_add_SNP, replace=False)
add_snp_indices = np.sort(add_snp_indices)
snp_data_add = process_snp_data(order_usedID_in_bed, add_snp_indices, snp_on_disk)
eff_add = np.random.normal(size=num_add_SNP)
add_var = np.sum(np.var(snp_data_add * eff_add.reshape(1, -1), axis=0))
eff_add *= np.sqrt(h2_add / add_var)
add_sum = np.dot(snp_data_add, eff_add)
del snp_data_add  # Free memory
gc.collect()

print("Generate GxE effects")
h2_gxe = h2_gxe - power_snp_h2_gxe # Adjust GxE heritability after accounting for power SNP
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
y = GxE_power_sum + add_sum + GxE_sum + error + E_sum + nxe_sum

dfR = dfP.loc[:, ["fid", "eid"]].copy()
dfR["pheno"] = y

print("Saving results")
out_prefix = f"power_snp_h2_gxe_{sys.argv[2]}.num_envi_{sys.argv[3]}.h2_add_{sys.argv[4]}.h2_gxe_{sys.argv[5]}.h2_nxe_{sys.argv[6]}.sample_{sample_size}.rep_{rep}"
dfR.to_csv(f"{out_prefix}.txt", sep=" ", index=False)

dfPowerSNPEnvi = pd.DataFrame({
    "power_envi_indices": power_envi_indices,
    "power_envi_names": [envi_names[i] for i in power_envi_indices]
})
dfPowerSNPEnvi.to_csv(f"{out_prefix}.power_snp_envi.txt", sep=" ", index=False)

np.savetxt(f"{out_prefix}.power_snp_index.txt", [power_snp_index - nSNP_exclude22], fmt="%d")
snpdata = Bed(bed_file)[order_usedID_in_bed, [power_snp_index]].read()
Bed.write(f"{out_prefix}.power_snp.bed",snpdata,count_A1=False)
dfBim.iloc[[power_snp_index], :].to_csv(f"{out_prefix}.power_snp.bim", sep=" ", header=False, index=False)
dfFam.iloc[order_usedID_in_bed, :].to_csv(f"{out_prefix}.power_snp.fam", sep=" ", header=False, index=False)

print("Results saved successfully.")
