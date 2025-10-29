import numpy as np
import pandas as pd
import os
from tqdm import tqdm

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_overfit_leadsnp_mmsusie/")

h2_add=30
h2_gxe=5
h2_nxe=15
sample_size=100000
pip_cutoff=0.8

res_lst = []
for num_envi_power in [1, 2, 5, 10]:
    for power_snp_h2_gxe in [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14]:
        total_lst = []
        false_lst = []
        total_snps_num = 0
        false_snps_num = 0
        print(f"Processing num_envi_power={num_envi_power}, power_snp_h2_gxe={power_snp_h2_gxe}")
        for rep in tqdm(range(1, 101)):
            prefix=f"power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi_power}.h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_{sample_size}.rep_{rep}"
            lead_snp_file = f"../res_overfit_leadsnp/{prefix}.lead_snps.txt"
            power_snp_envi_file = f"../pheno_take_home/{prefix}.power_snp_envi.txt"
            power_snp_envi_lst = pd.read_csv(power_snp_envi_file, sep=r"\s+").iloc[:, 0].to_list()
            try:
                df_lead_snps = pd.read_csv(lead_snp_file, sep=r"\s+")
                if df_lead_snps.empty:
                    print(f"File {lead_snp_file} is empty.")
                    continue
                total_snps = df_lead_snps['SNP'].tolist()
                total_snps_num += len(total_snps)
                for snp in total_snps:
                    file = f"{prefix}.{snp}.mmsusie_out.pip.txt"
                    dfPIP = pd.read_csv(file, sep=r"\s+", header=None)
                    envi_lst = dfPIP[dfPIP.iloc[:, 0] >= pip_cutoff].index.to_list()
                    total_lst.extend(envi_lst)
                    false_lst.extend([snp for snp in envi_lst if snp not in power_snp_envi_lst])
                false_snps = df_lead_snps[df_lead_snps['ld_r2'] < 0.1]["SNP"].tolist()
                false_snps_num += len(false_snps)
                for snp in false_snps:
                    dfPIP = pd.read_csv(file, sep=r"\s+", header=None)
                    envi_lst = dfPIP[dfPIP.iloc[:, 0] >= pip_cutoff].index.to_list()
                    false_lst.extend(envi_lst)
            except Exception as e:
                continue
        res_lst.append([
            num_envi_power,
            power_snp_h2_gxe,
            total_snps_num,
            false_snps_num,
            len(set(total_lst)),
            len(set(false_lst)),
            len(set(false_lst)) / len(set(total_lst)) if len(set(total_lst)) > 0 else np.nan
        ])


df_res = pd.DataFrame(res_lst, columns=["Num Envi Power", "Power SNP h2 GxE", "total_snps_num", "false_snps_num", "Total", "False", "FDR"])
df_res.to_csv("../res_overfit_sum/summary_results_mmsusie.csv", index=False)

df_res_aggregated = df_res.groupby('Num Envi Power').agg(
    Total_sum=('Total', 'sum'),
    False_sum=('False', 'sum')
).reset_index()

df_res_aggregated["FDR"] = df_res_aggregated["False_sum"] / df_res_aggregated["Total_sum"]
df_res_aggregated.to_csv("../res_overfit_sum/summary_results_mmsusie_aggregated.csv", index=False)
