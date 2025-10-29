import pandas as pd
import os

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_overfit_leadsnp_mmsusie/")

h2_add=30
h2_gxe=5
h2_nxe=15
sample_size=100000

res_lst = []
for num_envi_power in [1, 2, 5, 10, 20, 30, 40]:
    for power_snp_h2_gxe in [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14]:
        total_snp_lst = []
        false_snp_lst = []
        for rep in range(1, 101):
            prefix=f"power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi_power}.h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_{sample_size}.rep_{rep}"
            lead_snp_file = f"../res_overfit_leadsnp/{prefix}.lead_snps.txt"
            try:
                df_lead_snps = pd.read_csv(lead_snp_file, sep=r"\s+")
                if df_lead_snps.empty:
                    print(f"File {lead_snp_file} is empty.")
                    continue
                total_snp_lst.extend(df_lead_snps['SNP'].tolist())
                false_snp_lst.extend(df_lead_snps[df_lead_snps['ld_r2'] < 0.1]["SNP"].tolist())
            except Exception as e:
                continue
        fdr = len(false_snp_lst) / len(total_snp_lst)
        res_lst.append((num_envi_power, power_snp_h2_gxe, len(total_snp_lst), len(false_snp_lst), fdr))

df_res = pd.DataFrame(res_lst, columns=["Num Envi Power", "Power SNP h2 GxE", "Total SNPs", "False SNPs", "FDR"])

df_res_aggregated = df_res.groupby('Num Envi Power').agg(
    Total_SNPs_sum=('Total SNPs', 'sum'),
    False_SNPs_sum=('False SNPs', 'sum')
).reset_index()

df_res_aggregated["FDR"] = df_res_aggregated["False_SNPs_sum"] / df_res_aggregated["Total_SNPs_sum"]
df_res_aggregated.to_csv("../res_overfit_sum/summary_results_AIE.csv", index=False)
