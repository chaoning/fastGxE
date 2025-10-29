import pandas as pd
import numpy as np
import os

dfT = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/trait32.txt",
                   sep="\t", header=0, usecols=["FieldID", "TraitName", "ShortName"])

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/gif/traits/")

df_lst = []
for trait in dfT.iloc[:, 0]:
    file = f"{trait}.csv"
    df = pd.read_csv(file, index_col=0)
    df_lst.append(df.T)

df = pd.concat(df_lst, axis=0, ignore_index=True)
df = pd.concat([dfT, df], axis=1)

df.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/gif/PT.csv", index=False)


df_sum = df.iloc[:, 1:]

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
summary_percent = summary.map(lambda x: f"{x:.3f}")
summary_percent.to_csv(
    "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/gif/PT_summary.csv",
    index=True
)

dfT = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/trait67.txt",
                   sep="\t", header=0, usecols=["FieldID", "TraitName", "ShortName"])

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/gif/traits/")

df_lst = []
for trait in dfT.iloc[:, 0]:
    file = f"{trait}.csv"
    df = pd.read_csv(file, index_col=0)
    df_lst.append(df.T)

df = pd.concat(df_lst, axis=0, ignore_index=True)
df = pd.concat([dfT, df], axis=1)


df.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/gif/BB.csv", index=False)


df_sum = df.iloc[:, 1:]

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
summary_percent = summary.map(lambda x: f"{x:.3f}")
summary_percent.to_csv(
    "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/gif/BB_summary.csv",
    index=True
)