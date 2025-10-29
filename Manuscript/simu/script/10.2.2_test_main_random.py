import pandas as pd
import numpy as np
import os
from scipy.stats import chi2
import glob
import re
import sys
from pysnptools.snpreader import Bed


rep = int(sys.argv[1])
input_window_index = int(sys.argv[2])

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/validate_mmsusie_random/")

lead_snp = pd.read_csv(f"{rep}.{input_window_index}.mmsusie_out.lead_snp", sep=r"\s+").iloc[0, 0]

bedfile = "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/random_100000"
dfBim = pd.read_csv(f"{bedfile}.bim", sep=r"\s+", header=None)

iid_lst = pd.read_csv(f"{bedfile}.fam", sep=r"\s+", header=None).iloc[:, 1].tolist()

dct = dict(zip(dfBim[1], range(len(dfBim[1]))))

lead_snp_index = dct[lead_snp]

snp_on_disk = Bed(bedfile, count_A1=True)

beta_lst = []
se_lst = []

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/validate_AIE_random/")

for group in range(1, 6):
    file = f"{rep}.{input_window_index}.score_group_{group}.txt"
    dfP = pd.read_csv(file, sep=r"\s+")

    snpdata = snp_on_disk[:, [lead_snp_index]].read().val
    df = pd.DataFrame(snpdata, columns=["SNP"])
    df.fillna(df.mean(), inplace=True)
    df["eid"] = iid_lst

    dfm = pd.merge(dfP, df, on="eid")
    y = dfm["trait"].values
    E = dfm.loc[:, "age":"alcohol_frequency"].values
    G = dfm["SNP"].values.reshape(-1, 1)

    # correct for covariates
    X = np.hstack((np.ones((y.shape[0], 1)), E))
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    y_resid = y - X @ beta

    # scale and center residuals and SNPs
    y_resid -= np.mean(y_resid)
    y_resid /= np.std(y_resid)
    G -= np.mean(G)
    G /= np.std(G)

    # linear regression to test for association
    beta = np.dot(G.T, y_resid) / np.sum(G**2)
    se = np.sqrt(1 / np.sum(G**2))

    beta_lst.append(beta)
    se_lst.append(se)


