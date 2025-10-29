import pandas as pd
import numpy as np
import math


df = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/description/trait.csv")

fam = pd.read_csv('/net/zootopia/disk1/chaon/data/UKB/imp/ukb_imp_info.qc.fam', sep=r'\s+',
                  names=['FID', 'eid', 'A', 'B', 'SEX', 'PHE'])
fam = fam.loc[:, "eid"]
df = pd.merge(df, fam, on="eid")

dfT = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/description/trait_Description.txt", sep="\t", dtype=str)
trait_lst = dfT.iloc[:, 0].to_list()
dct = dict(zip([f"{trait}-0.0" for trait in trait_lst], trait_lst))
df.rename(columns=dct, inplace=True)
categories_lst = list(dfT.iloc[:, -1])


import pandas as pd
import numpy as np

def cal_trait(df, trait, sigma=3):
    """
    Clean a single trait:
    1. Remove missing values
    2. Remove values < 0
    3. Remove outliers based on the sigma rule (default: 3σ)

    Args:
        df (pd.DataFrame): Input data containing 'eid' and the target trait
        trait (str): Column name of the trait to process
        sigma (float): Standard deviation threshold for outlier removal (default=3)

    Returns:
        pd.DataFrame: Cleaned data containing only 'eid' and the processed trait
    """
    df_new = df.loc[:, ["eid", trait]].copy()

    # Step 1: Remove missing values
    df_new = df_new.dropna(subset=[trait])
    print(f"Initial sample size: {df_new.shape[0]}")

    # Step 2: Remove values < 0
    neg_count = (df_new[trait] < 0).sum()
    df_new = df_new[df_new[trait] >= 0]
    print(f"Removed values < 0: {neg_count}, Remaining: {df_new.shape[0]}")

    # Step 3: Compute outlier bounds
    mean = df_new[trait].mean()
    std = df_new[trait].std()
    lower_bound = mean - sigma * std
    upper_bound = mean + sigma * std

    # Step 4: Remove outliers
    outlier_mask = (df_new[trait] < lower_bound) | (df_new[trait] > upper_bound)
    outlier_count = outlier_mask.sum()
    df_new = df_new[~outlier_mask]
    print(f"Removed outliers ({sigma}σ): {outlier_count}, Remaining: {df_new.shape[0]}")

    return df_new



import pandas as pd
import numpy as np

def envi_scale(df):
    """
    Standardize environmental variables and create derived interaction terms.

    Steps:
    1. Standardize all columns from 'Sex' onward (z-score normalization).
    2. Create Sleep_duration^2 and standardize it.
    3. Add polynomial and interaction terms for Age and Sex.

    Args:
        df (pd.DataFrame): Input DataFrame containing at least:
                           'Sex', 'Age', 'Sleep_duration', and other env variables.

    Returns:
        pd.DataFrame: Scaled and extended DataFrame.
    """
    df_scaled = df.copy()

    # Step 1: Standardize environmental variables from 'Sex' onward
    # Force float dtype once
    df_scaled = df.copy()
    df_scaled.loc[:, 'Sex':] = df_scaled.loc[:, 'Sex':].apply(pd.to_numeric, downcast=None).astype(float)
    # Now scaling works without warnings
    env_cols = df_scaled.loc[:, 'Sex':]
    
    df_scaled.loc[:, 'Sex':] = (env_cols - env_cols.mean()) / env_cols.std()

    # Step 2: Add Sleep_duration^2 (squared) and standardize it
    sleep_sq = df_scaled["Sleep_duration"] ** 2
    df_scaled["Sleep_duration2"] = (sleep_sq - sleep_sq.mean()) / sleep_sq.std()

    # Step 3: Create Age and Sex interaction/polynomial terms
    age = df_scaled["Age"].to_numpy()
    sex = df_scaled["Sex"].to_numpy()

    # Insert in reverse order so column positions match original intention
    df_scaled.insert(24, "Age2", age ** 2)
    df_scaled.insert(25, "Age3", age ** 3)
    df_scaled.insert(26, "Sex_Age", sex * age)
    df_scaled.insert(27, "Sex_Age2", sex * (age ** 2))
    df_scaled.insert(28, "Sex_Age3", sex * (age ** 3))

    return df_scaled



dfE = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/pheno/environments/envi42.txt", sep=r"\s+")

for i in range(len(trait_lst)):
    trait = trait_lst[i]
    print("########", trait)
    df_new = cal_trait(df, trait)
    df_res42 = pd.merge(df_new, dfE, on="eid")
    df_res42 = envi_scale(df_res42)
    print("The number of used iids: {}".format(df_res42.shape[0]))
    df_res42.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/pheno/{}.txt".format(trait), sep=" ",
                    index=False)
