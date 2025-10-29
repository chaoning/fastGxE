import pandas as pd
import numpy as np
from tqdm import tqdm
import sys


dfSig = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/signal/signal.txt", sep="\s+")
dfSig = dfSig.iloc[:, 1:].copy()
dfSig.rename(columns={"SNP": "trait_leading_snp"}, inplace=True)

GL_file = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.csv"
dfGLTrait = pd.read_csv(GL_file)
print(dfGLTrait.shape)
dfGLTrait2 = dfGLTrait.assign(trait_leading_snp=dfGLTrait['trait_leading_snp'].str.split(';')).explode('trait_leading_snp')
print(dfGLTrait2.shape)

dfT = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/trait67_shortName.txt", sep="\t")
dfT = dfT.iloc[:, [0, -1]].copy()
dfT.columns = ["trait", "TraitName"]

dfR = pd.merge(dfGLTrait2, dfT, on="trait", how="left")
dfR = pd.merge(dfR, dfSig, on=["trait", "trait_leading_snp"], how="left")

dfR.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.addP.csv",
           index=False)
