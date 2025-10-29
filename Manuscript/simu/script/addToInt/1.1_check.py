import os
import pandas as pd
from tqdm import tqdm   
h2_add=30
h2_gxe=5
h2_nxe=15
sample_size=100000

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_addToint/")

k = 0
lst = []
for rep in tqdm(range(1, 101)):  # 1 to 100 inclusive
    for unknown_h2SNP_h2E in [0.02, 0.08, 0.16]:
        for R2 in [0, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9]:
            for R2GE in [0]:
                k += 1
                out_prefix = f"h2_{unknown_h2SNP_h2E}.R2_{R2}.R2GE_{R2GE}.h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_{sample_size}.rep_{rep}"
                # check if output file already exists
                try:
                    df = pd.read_csv(f"{out_prefix}.txt", sep=r"\s+")
                    if (df.shape[0] != sample_size) or (df.isnull().values.any()):
                        lst.append(k)
                        print(out_prefix)
                except Exception as e:
                    lst.append(k)
                    print(out_prefix)
                    print(f"Error reading {out_prefix}.txt: {e}")
        
        for R2 in [0.5]:
            for R2GE in [0, 0.01, 0.1, 0.3]:
                k += 1
                out_prefix = f"h2_{unknown_h2SNP_h2E}.R2_{R2}.R2GE_{R2GE}.h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_{sample_size}.rep_{rep}"
                # check if output file already exists
                try:
                    df = pd.read_csv(f"{out_prefix}.txt", sep=r"\s+")
                    if (df.shape[0] != sample_size) or (df.isnull().values.any()):
                        lst.append(k)
                        print(out_prefix)
                except Exception as e:
                    lst.append(k)
                    print(out_prefix)
                    print(f"Error reading {out_prefix}.txt: {e}")

print(k)
print(lst)

