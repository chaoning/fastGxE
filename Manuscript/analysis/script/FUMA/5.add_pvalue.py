import pandas as pd
import numpy as np
from tqdm import tqdm
import sys


dfSig = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/signal/signal.txt", sep=r"\s+")
dfSig = dfSig.iloc[:, 1:].copy()
dfSig.rename(columns={"SNP": "trait_leading_snp"}, inplace=True)

dfBim = pd.read_csv("/net/zootopia/disk1/chaon/data/UKB/imp/ukb_imp_info.qc.bim.hg38", sep=r"\s+", header=None)
dfBim = dfBim.iloc[:, [1, 3]]
dfBim.columns = ["trait_leading_snp", "base_hg38"]
dfSig = pd.merge(dfBim, dfSig, on="trait_leading_snp", how="right")

GL_file = sys.argv[1]
dfGLTrait = pd.read_csv(GL_file)

dfT = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/trait32.txt", sep="\t")
dfT = dfT.iloc[:, [0, 2]].copy()
dfT.columns = ["trait", "TraitName"]

dfR = pd.merge(dfGLTrait, dfT, on="trait", how="left")
dfR = pd.merge(dfR, dfSig, on=["trait", "trait_leading_snp"], how="left")

dfR.to_csv(sys.argv[2], index=False)

