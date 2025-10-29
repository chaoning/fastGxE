import pandas as pd

df = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/annovar/plink22.variant_function", sep=r"\s+", header=None, usecols=[0, 3], names=["variant_function", "pos"])

dfSNP= pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/annovar/plink22.bim", sep=r"\s+", header=None, usecols=[1, 3],  names=["SNP", "pos"])

dfm = pd.merge(df, dfSNP, on="pos", how="left")

print(dfm.shape, df.shape, dfSNP.shape)

variant_function_lst = dfm.variant_function.unique().tolist()

with open("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/QQplot_vf/variant_function.txt", "w") as f:
    for vf in variant_function_lst:
        df_vf = dfm[dfm.variant_function == vf].copy()
        df_vf = df_vf.drop_duplicates(subset="SNP")
        if df_vf.shape[0] > 100:
            df_vf.to_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/QQplot_vf/variant_function.{vf}.txt", sep="\t", index=False)
            f.write(f"{vf} {df_vf.shape[0]}\n")
