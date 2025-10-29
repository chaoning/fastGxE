import sys
import pandas as pd
import numpy as np
import os


file = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.addP_logP.csv"

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/run/")

df = pd.read_csv(file)

pip_lst = []
main_effect_lst = []
for i in range(df.shape[0]):
    trait = df.iloc[i, 13]
    snp = df.iloc[i, 14]
    file = "{}.{}.mmsusie.pip.txt".format(trait, snp)
    df2 = pd.read_csv(file, sep=r"\s+", header=None)
    pip = np.array(df2.iloc[:, 0])
    alpha = np.loadtxt("{}.{}.mmsusie.alpha.txt".format(trait, snp))
    mu = np.loadtxt("{}.{}.mmsusie.mu.txt".format(trait, snp))
    effect = np.sum(alpha * mu, axis=0)
    pip[effect < 0] = -pip[effect < 0]
    pip = np.array(pip).reshape(1, -1)
    pip_lst.append(pip)
    main_effect = np.loadtxt("{}.{}.main_effect.txt".format(trait, snp))
    main_effect_lst.append(main_effect)

head_df = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/envi42.txt",
                      header=0, sep=r"\s+")
head_lst = list(head_df.iloc[:, 0])
pip_df = pd.DataFrame(np.concatenate(pip_lst, axis=0), columns=head_lst)
pip_df.insert(0, "Main", main_effect_lst)
dfm = pd.concat([df, pip_df], axis=1)

dfm = dfm.sort_values(by=["chrom_x", "start_hg38", "end_hg38", "p_gxe"])

dfm.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/PT_PIP.csv", index=False)


file = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.addP_logP.csv"

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/run/")

df = pd.read_csv(file)

pip_lst = []
main_effect_lst = []
for i in range(df.shape[0]):
    trait = df.iloc[i, 13]
    snp = df.iloc[i, 14]
    file = "{}.{}.mmsusie.pip.txt".format(trait, snp)
    df2 = pd.read_csv(file, sep=r"\s+", header=None)
    pip = np.array(df2.iloc[:, 0])
    alpha = np.loadtxt("{}.{}.mmsusie.alpha.txt".format(trait, snp))
    mu = np.loadtxt("{}.{}.mmsusie.mu.txt".format(trait, snp))
    effect = np.sum(alpha * mu, axis=0)
    pip[effect < 0] = -pip[effect < 0]
    pip = np.array(pip).reshape(1, -1)
    pip_lst.append(pip)
    main_effect = np.loadtxt("{}.{}.main_effect.txt".format(trait, snp))
    main_effect_lst.append(main_effect)

head_df = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/envi42.txt",
                      header=0, sep=r"\s+")
head_lst = list(head_df.iloc[:, 0])
pip_df = pd.DataFrame(np.concatenate(pip_lst, axis=0), columns=head_lst)
pip_df.insert(0, "Main", main_effect_lst)
dfm = pd.concat([df, pip_df], axis=1)

dfm = dfm.sort_values(by=["chrom_x", "start_hg38", "end_hg38", "p_gxe"])

dfm.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/BB_PIP.csv", index=False)
