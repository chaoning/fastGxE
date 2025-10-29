import os
import sys
import pandas as pd
from tqdm import tqdm

# -------------------
# Input arguments
# -------------------
rGI = sys.argv[1]            # e.g., 0.2, 0.5, 0.8
add_gxe_ratio = sys.argv[2]  # e.g., 1.5
# prefix = "num_envi_30.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000"
num_envi_power = sys.argv[3]
prefix = f"num_envi_{num_envi_power}.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000"          # e.g., "num_envi_30.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000"
# -------------------
# Paths and constants
# -------------------
base_dir = "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power"
os.chdir(f"{base_dir}/res_filter_fastgxe/")


# Initialize counters
power_gxe = 0
power_filter = 0
total_snp = 0
res_lst = []

# -------------------
# Main loop over replicates
# -------------------
for rep in tqdm(range(1, 101), desc="Processing replicates"):
    try:
        # File paths
        file_res = f"rGI_{rGI}.add_gxe_ratio_{add_gxe_ratio}.{prefix}.rep_{rep}.res"
        file_snp = f"../pheno_filter/rGI_{rGI}.add_gxe_ratio_{add_gxe_ratio}.{prefix}.rep_{rep}.power_snp.bim"

        # Read result file
        df = pd.read_csv(file_res, sep=r"\s+")

        # Read SNP list
        dfSNP = pd.read_csv(
            file_snp, sep=r"\s+", usecols=[1], header=None, names=["SNP"]
        )

        # ---- Power for GxE without filtering ----
        dfM = pd.merge(dfSNP, df, on="SNP", how="left")
        power_gxe += (dfM["p_gxe"] < 5e-8).sum()

        # ---- Power for GxE with prefiltering ----
        df_filtered = df[df["p_main"] < 5e-8]
        if not df_filtered.empty:
            p_cut = 0.05 / df_filtered.shape[0]  # Bonferroni threshold
            dfM_filtered = pd.merge(dfSNP, df_filtered, on="SNP")
            power_filter += (dfM_filtered["p_gxe"] < p_cut).sum()

        total_snp += dfSNP.shape[0]  # Each replicate has 10 power SNPs

    except FileNotFoundError:
        print(f"[Warning] Missing file for replicate {rep}, skipping...")
        continue
    except Exception as e:
        print(f"[Error] Rep {rep} failed: {e}")
        continue

# -------------------
# Calculate proportions
# -------------------
if total_snp > 0:
    power_gxe /= total_snp
    power_filter /= total_snp

# -------------------
# Save results
# -------------------
print(
    f"rGI: {rGI}, add_gxe_ratio: {add_gxe_ratio}, "
    f"Power (GxE): {power_gxe:.4f}, Power (Filter): {power_filter:.4f}, Total SNP: {total_snp}"
)

res_lst.append((rGI, add_gxe_ratio, power_gxe, power_filter, total_snp))
df_res = pd.DataFrame(res_lst, columns=["rGI", "add_gxe_ratio", "power_gxe", "power_filter", "total_snp"])

out_dir = f"{base_dir}/res_filter_sum"
os.makedirs(out_dir, exist_ok=True)
out_file = f"{out_dir}/rGI_{rGI}.add_gxe_ratio_{add_gxe_ratio}.num_envi_{num_envi_power}.power_summary.csv"
df_res.to_csv(out_file, index=False)
