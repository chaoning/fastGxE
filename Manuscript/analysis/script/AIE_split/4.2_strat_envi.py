import os
import pandas as pd
import numpy as np

# -----------------------------
# Paths & IO setup
# -----------------------------
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/AIE_split/split5/")

f_loci = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.addP.csv"
f_envs = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/PT_GL_envi.csv"
f_raw  = "/net/zootopia/disk1/chaon/WORK/GxE/pheno/environments/envi42.txt"

out_dir = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/AIE_split/split5_envi"
os.makedirs(out_dir, exist_ok=True)

N_GROUPS = 5  # score_group_1..5

# -----------------------------
# Load metadata once
# -----------------------------
dfL = pd.read_csv(f_loci)
dfE = pd.read_csv(f_envs)

# Load raw (unstandardized) environments once and index by eid
df_raw = pd.read_csv(f_raw, sep=r"\s+")
if "eid" not in df_raw.columns:
    raise ValueError("`eid` column is required in the raw environment file.")
df_raw = df_raw.set_index("eid")

# -----------------------------
# Helper: column-wise mean & SE
# -----------------------------
def col_mean_se(df: pd.DataFrame, cols):
    """
    Compute (mean, SE) for each column in `cols`.
    - Converts to numeric (coerce errors to NaN).
    - SE = sd(ddof=1) / sqrt(n_non_na); if n<=1 -> NaN.
    - Missing columns -> NaN.
    Returns: (means_list, ses_list)
    """
    means, ses = [], []
    for c in cols:
        if c in df.columns:
            s = pd.to_numeric(df[c], errors="coerce")
            n = s.notna().sum()
            mean = s.mean()
            sd = s.std(ddof=1)
            se = (sd / np.sqrt(n)) if n > 1 else np.nan
            means.append(mean)
            ses.append(se)
        else:
            means.append(np.nan)
            ses.append(np.nan)
    return means, ses

# -----------------------------
# Main loop over (trait, snp)
# -----------------------------
pairs = (
    dfL[["trait", "trait_leading_snp"]]
    .drop_duplicates()
    .itertuples(index=False, name=None)
)

for trait, snp in pairs:
    # Environments to summarise for this (trait, snp)
    envi_lst = dfE.loc[
        (dfE["FieldID"] == trait) & (dfE["LeadSNP"] == snp),
        "Environment"
    ].tolist()
    if not envi_lst:
        continue

    group_means, group_ses = [], []
    group_means_raw, group_ses_raw = [], []

    for g in range(1, N_GROUPS + 1):
        fn = f"{trait}.{snp}.score_group_{g}.txt"
        if not os.path.exists(fn):
            group_means.append([np.nan] * len(envi_lst))
            group_ses.append([np.nan] * len(envi_lst))
            group_means_raw.append([np.nan] * len(envi_lst))
            group_ses_raw.append([np.nan] * len(envi_lst))
            continue

        # Read current group file (whitespace separated)
        df = pd.read_csv(fn, sep=r"\s+")

        # Stats on group file (possibly standardized columns)
        m, se = col_mean_se(df, envi_lst)
        group_means.append(m)
        group_ses.append(se)

        # Stats on raw scale: align to eids present in this group
        if "eid" in df.columns:
            sub_raw = df_raw.reindex(df["eid"].values)
            m_raw, se_raw = col_mean_se(sub_raw, envi_lst)
        else:
            m_raw  = [np.nan] * len(envi_lst)
            se_raw = [np.nan] * len(envi_lst)

        group_means_raw.append(m_raw)
        group_ses_raw.append(se_raw)

    # Build wide tables (5 x E) and save
    idx = pd.Index(range(1, N_GROUPS + 1), name="Group")

    df_mean = pd.DataFrame(group_means,      index=idx, columns=envi_lst).reset_index()
    df_se   = pd.DataFrame(group_ses,        index=idx, columns=envi_lst).reset_index()
    df_mean_raw = pd.DataFrame(group_means_raw, index=idx, columns=envi_lst).reset_index()
    df_se_raw   = pd.DataFrame(group_ses_raw,   index=idx, columns=envi_lst).reset_index()

    base = f"{trait}.{snp}.score_group_envi"
    df_mean.to_csv(os.path.join(out_dir, f"{base}_mean.csv"), index=False)
    df_se.to_csv(os.path.join(out_dir,   f"{base}_se.csv"),   index=False)
    df_mean_raw.to_csv(os.path.join(out_dir, f"{base}_mean_raw.csv"), index=False)
    df_se_raw.to_csv(os.path.join(out_dir,   f"{base}_se_raw.csv"),   index=False)

print("Done!")
