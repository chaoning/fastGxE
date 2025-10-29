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
dfPIP_lst = []
for num_envi_power in [1, 2, 5, 10]:
    for power_snp_h2_gxe in [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14]:
        total_lst = []
        false_lst = []
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
                for snp in total_snps:
                    file = f"{prefix}.{snp}.mmsusie_out.pip.txt"
                    dfPIP = pd.read_csv(file, sep=r"\s+", header=None)
                    dfPIP['in_lst'] = dfPIP.index.isin(power_snp_envi_lst).astype(int)
                false_snps = df_lead_snps[df_lead_snps['ld_r2'] < 0.1]["SNP"].tolist()
                dfPIP_lst.append(dfPIP)
                for snp in false_snps:
                    dfPIP = pd.read_csv(file, sep=r"\s+", header=None)
                    dfPIP['in_lst'] = 0
            except Exception as e:
                continue
dfPIP = pd.concat(dfPIP_lst, ignore_index=True)
dfPIP.columns = ["PIP", "in_lst"]
dfPIP["bin"] = pd.cut(dfPIP["PIP"], bins=10)
result = dfPIP.groupby("bin").agg(
    mean_PIP=("PIP", "mean"),
    std_PIP=("PIP", "std"),
    in_lst_sum=("in_lst", "sum"),
    bin_count=("in_lst", "count")
)
result["SE"] = result["std_PIP"] / np.sqrt(result["bin_count"])
result["in_lst_rate"] = result["in_lst_sum"] / result["bin_count"]
result = result.reset_index()
result.to_csv(f"../res_overfit_sum/mmsusie_PIP_summary.csv", index=False)
