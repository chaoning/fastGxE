import numpy as np
import pandas as pd
from scipy.stats import norm
import re
from scipy.stats import chi2



file = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/trait67_shortName.txt"
df = pd.read_csv(file, sep="\t", header=0)

k = 0
with open("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/variance.txt", "w") as fout:
    fout.write("Trait\th2_1\th2_2\th2_gxe1\th2_gxe2\th2_nxe\tLRT_stat\tp_nxe\n")
    for trait in df.iloc[:, 0]:
        file = f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/mom/{trait}.var"
        var1 = np.loadtxt(file)
        file = f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/mom/{trait}_noNxE.var"
        var2 = np.loadtxt(file)

        h2_1 = var1[0] / np.sum(var1)
        h2_2 = var2[0] / np.sum(var2)
        h2_gxe1 = var1[1] / np.sum(var1)
        h2_gxe2 = var2[1] / np.sum(var2)
        h2_nxe = var1[2] / np.sum(var1)

        k += 1
        file = f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Aout/fastgxe_biomarker.{k}.txt"
        LRT_stat = 0
        p_val = 0
        with open(file, "r") as f:
            content = f.read()
            matches = re.findall(r"-2logL in null model: ([\d\.]+)", content)
            values = [float(val) for val in matches]
            LRT_stat = values[1] - values[0]
            p_val = chi2.sf(LRT_stat, df=1)
        fout.write(f"{trait}\t{h2_1:.6f}\t{h2_2:.6f}\t{h2_gxe1:.6f}\t{h2_gxe2:.6f}\t{h2_nxe:.6f}\t{LRT_stat:.6f}\t{p_val:.6f}\n")
