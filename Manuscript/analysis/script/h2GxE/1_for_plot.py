import pandas as pd

df1 = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/mom/variance.txt", sep="\t")

df1 = df1.loc[:, ["Trait", "h2_1", "h2_gxe1", "h2_nxe", "p_nxe"]]
df1.columns = ["FieldID", "h2G", "h2GxE", "h2NxE", "p_NxE"]
dfT = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/trait32.txt", sep="\t", usecols=["FieldID", "TraitName", "ShortName"])

df = pd.merge(dfT, df1, on="FieldID", how="left")
df.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mom/PT.variance_for_plot.csv", index=False)

df_sum = df.loc[:, ["h2G", "h2GxE", "h2NxE"]]

summary = df_sum.describe(percentiles=[0.25, 0.5, 0.75]).T[
    ["min", "25%", "50%", "75%", "max", "mean"]
]

summary = summary.rename(columns={
    "min": "Min",
    "25%": "Lower quartile",
    "50%": "Median",
    "75%": "Upper quartile",
    "max": "Max",
    "mean": "Mean"
})
summary_percent = summary.map(lambda x: f"{x*100:.2f}%")
summary_percent.to_csv(
    "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mom/PT.variance_summary.csv",
    index=True
)


df1 = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/mom/variance.txt", sep="\t")

df1 = df1.loc[:, ["Trait", "h2_1", "h2_gxe1", "h2_nxe", "p_nxe"]]
df1.columns = ["FieldID", "h2G", "h2GxE", "h2NxE", "p_NxE"]
cols = ["h2G", "h2GxE", "h2NxE"]
df1.loc[:, cols] = df1.loc[:, cols].clip(lower=0)

dfT = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/trait67.txt", sep="\t", usecols=["FieldID", "TraitName", "ShortName"])

df = pd.merge(dfT, df1, on="FieldID", how="left")
df.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mom/BB.variance_for_plot.csv", index=False)

df_sum = df.loc[:, ["h2G", "h2GxE", "h2NxE"]]

summary = df_sum.describe(percentiles=[0.25, 0.5, 0.75]).T[
    ["min", "25%", "50%", "75%", "max", "mean"]
]

summary = summary.rename(columns={
    "min": "Min",
    "25%": "Lower quartile",
    "50%": "Median",
    "75%": "Upper quartile",
    "max": "Max",
    "mean": "Mean"
})
summary_percent = summary.map(lambda x: f"{x*100:.2f}%")
summary_percent.to_csv(
    "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mom/BB.variance_summary.csv",
    index=True
)
