import pandas as pd
import numpy as np
df1 = pd.read_csv("/net/mulan/home/chaon/WORK/GxE/Analysis/Physical_measures/result/SGEM/FUMA/IINT/GenomicRiskLoci.hg38.keygene.TraitSpecific.addP.csv")

beta_columns = [f'beta{i}' for i in range(1, 43)]
p_columns = [f'p{i}' for i in range(1, 43)]

# Extracting the relevant columns
beta_data = df1[beta_columns]
beta_data.columns = range(42)
p_data = df1[p_columns]
p_data.columns = range(42)

# Melt the beta columns
beta_melted = beta_data.melt(var_name='envi', value_name='beta_GxE')
# Melt the p columns
p_melted = p_data.melt(var_name='envi', value_name='p_GxE')

# Concatenate the melted beta and p columns horizontally
melted_data = pd.concat([beta_melted, p_melted['p_GxE']], axis=1)

dfE = pd.read_csv("/net/mulan/home/chaon/WORK/GxE/Analysis/Physical_measures/result/SGEM/SGEM-C/e42.head", sep="\s+")
E_arr = dfE.iloc[:, 0].to_numpy()
melted_data["envi"] = E_arr[melted_data["envi"]]

# Combine beta and p-values with the repeated columns, excluding 'p_variable'
id_columns = list(df1.loc[:, :"p_main"].columns)
id_columns.extend(['p_min', 'p_SGEM_C', 'p_SGEM_V', 'p_SGEM_O'])

# Extract the columns that will be repeated
repeated_data = df1[id_columns]

# Repeat the other columns to match the number of rows in the melted data
repeated_data_expanded = pd.concat([repeated_data] * len(beta_columns), ignore_index=True)

# Combine all the data into one DataFrame
final_data = pd.concat([repeated_data_expanded, melted_data], axis=1)


df2 = pd.read_csv("/net/mulan/home/chaon/WORK/GxE/Analysis/Physical_measures/result/SGEM/susie/e42IINT_PIP/PIP.sum.loci.envi.addShortName.txt",
                  sep="\t")
df2 = df2.loc[:, :"pip"]
df2["pip"] = np.absolute(df2["pip"])
df2.rename(columns={"leadingSNP": "trait_leading_snp"}, inplace=True)
df = pd.merge(df2, final_data, on=["trait", "trait_leading_snp", "envi"], how="left")
df.to_csv("/net/mulan/home/chaon/WORK/GxE/Analysis/Physical_measures/result/SGEM/FUMA/IINT/GenomicRiskLoci.hg38.keygene.TraitSpecific.addP.addE.csv",
          index=False)

