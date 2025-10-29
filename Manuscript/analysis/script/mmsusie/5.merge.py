import pandas as pd
df1 = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.addP_logP.csv")
df1 = df1.loc[:, ['chrom_x', 'start', 'end', 'GenomicLocus_hg38', 'Genes', 'trait', "TraitName", 'trait_leading_snp', 'p_main', 'p_gxe', 'p_gxe_log']]
df1.columns = ['chrom', 'start', 'end', 'GenomicLocus', 'Gene', 'FieldID', 'TraitName', 'LeadSNP', 'p_main', 'p_gxe', 'p_gxe_log']

dfT = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/trait32.txt", sep="\t", usecols=['FieldID', 'ShortName', ])

df1 = pd.merge(df1, dfT, on="FieldID", how="left")

df2 = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/PT.PIP.sum.loci.envi.csv")
df2.columns = ['FieldID', 'LeadSNP', 'Environment', 'PIP']

df = pd.merge(df1, df2, on=['FieldID', "LeadSNP"], how="right")
df["Direct"] = "+"

df.loc[df["PIP"] < 0 , "Direct"] = "-"
df["PIP"] = df["PIP"].abs()
df.sort_values(by=['chrom', 'start', 'end', 'p_gxe', 'PIP'], inplace=True)
df.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/PT_GL_envi.csv", index=False)


# Collapse Environment by trait+LeadSNP, keep merged string and count
df = (
    df.groupby([
        'chrom','start','end','GenomicLocus','Gene',
        'FieldID','TraitName', "ShortName", 'LeadSNP','p_main','p_gxe','p_gxe_log'
    ], dropna=False, as_index=False)
    .agg({
        'Environment': [
            lambda x: '; '.join(sorted(set(x.dropna()))),  # merged string
            lambda x: len(set(x.dropna()))                 # number of unique environments
        ]
    })
)

df.columns = [
    'chrom','start','end','GenomicLocus','Gene',
    'FieldID','TraitName', "ShortName", 'LeadSNP','p_main','p_gxe','p_gxe_log',
    'Environment','Environment_count'
]

# Sort
df.sort_values(by=['chrom','start','end','p_gxe'], inplace=True)

# Save
df.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/PT_GL_mergeE.csv", index=False)

df1 = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.addP_logP.csv")
df1 = df1.loc[:, ['chrom_x', 'start', 'end', 'GenomicLocus_hg38', 'Genes', 'trait', "TraitName", 'trait_leading_snp', 'p_main', 'p_gxe', 'p_gxe_log']]
df1.columns = ['chrom', 'start', 'end', 'GenomicLocus', 'Gene', 'FieldID', 'TraitName', 'LeadSNP', 'p_main', 'p_gxe', 'p_gxe_log']

dfT = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/trait67.txt", sep="\t", usecols=['FieldID', 'ShortName', ])

df1 = pd.merge(df1, dfT, on="FieldID", how="left")

df2 = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/BB.PIP.sum.loci.envi.csv")
df2.columns = ['FieldID', 'LeadSNP', 'Environment', 'PIP']

df = pd.merge(df1, df2, on=['FieldID', "LeadSNP"], how="outer")
df["Direct"] = "+"

df.loc[df["PIP"] < 0 , "Direct"] = "-"
df["PIP"] = df["PIP"].abs()
df.sort_values(by=['chrom', 'start', 'end', 'p_gxe', 'PIP'], inplace=True)
df.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/BB_GL_envi.csv", index=False, na_rep="NA")

# Collapse Environment by trait+LeadSNP, keep merged string and count
df = (
    df.groupby([
        'chrom','start','end','GenomicLocus','Gene',
        'FieldID','TraitName', "ShortName", 'LeadSNP','p_main','p_gxe','p_gxe_log'
    ], dropna=False, as_index=False)
    .agg({
        'Environment': [
            lambda x: '; '.join(sorted(set(x.dropna()))),  # merged string
            lambda x: len(set(x.dropna()))                 # number of unique environments
        ]
    })
)

df.columns = [
    'chrom','start','end','GenomicLocus','Gene',
    'FieldID','TraitName', "ShortName", 'LeadSNP','p_main','p_gxe','p_gxe_log',
    'Environment','Environment_count'
]

# Sort
df.sort_values(by=['chrom','start','end', 'TraitName', 'p_gxe'], inplace=True)

# Save
df.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/BB_GL_mergeE.csv", index=False, na_rep="NA")
