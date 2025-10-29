import os
import re
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2
from numpy.linalg import LinAlgError


# ---------- Utilities ----------

def to_float64_matrix(a, name="array", make_2d=None):
    """
    Cast input to float64 ndarray and verify finiteness.
    make_2d: None (no reshape), "col" -> (n,1), "row" -> (1,n)
    """
    if hasattr(a, "to_numpy"):
        a = a.to_numpy()

    # robust numeric coercion (handles object arrays)
    if isinstance(a, np.ndarray) and a.dtype == object:
        a = pd.DataFrame(a).apply(pd.to_numeric, errors="coerce").to_numpy()

    a = np.asarray(a, dtype=np.float64)

    if make_2d == "col" and a.ndim == 1:
        a = a.reshape(-1, 1)
    elif make_2d == "row" and a.ndim == 1:
        a = a.reshape(1, -1)

    if not np.all(np.isfinite(a)):
        bad = np.argwhere(~np.isfinite(a))
        raise ValueError(f"{name} contains non-finite values at indices like {bad[:5].tolist()}")
    return a


def newest_match(directory: Path, pattern_regex: re.Pattern):
    """
    Return the newest file (by mtime) in `directory` whose name matches `pattern_regex`.
    Raise FileNotFoundError if none matches.
    """
    candidates = [p for p in directory.iterdir() if pattern_regex.search(p.name)]
    if not candidates:
        raise FileNotFoundError(f"No file matches regex in {directory}: {pattern_regex.pattern}")
    # pick the most recently modified file
    return max(candidates, key=lambda p: p.stat().st_mtime)


# ---------- Core ----------

def process_trait_specific(summary_csv: str, out_csv: str,
                           work_dir: str = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/AIE_split/test_main/",
                           verbose: bool = True,
                           use_cholesky: bool = False,
                           ridge_eps: float = 1e-8):
    """
    For each (trait, snp) pair from `summary_csv`, find 5 group result files,
    read beta/se/p, assemble Wald test for contrasts defined by `arr`,
    and write results to `out_csv`.

    - work_dir: directory where *.res files live
    - use_cholesky: if True and V is SPD, use cho_solve; else np.linalg.solve
    - ridge_eps: add eps*I to V when solve fails (ill-conditioned)
    """

    # 4x5 contrast matrix: tests equality of consecutive betas
    arr = np.array([[1, -1, 0,  0,  0],
                    [0,  1, -1, 0,  0],
                    [0,  0,  1, -1, 0],
                    [0,  0,  0,  1, -1]], dtype=np.float64)

    work_dir = Path(work_dir)

    # --- read summary and extract columns 14th & 15th (0-based: 13, 14) ---
    dfS = pd.read_csv(summary_csv)

    # Defensive: ensure the indices exist
    if dfS.shape[1] <= 14:
        raise ValueError(f"{summary_csv} must have at least 15 columns; got {dfS.shape[1]}")

    trait_lst = dfS.iloc[:, 13].astype(str).tolist()
    snp_lst   = dfS.iloc[:, 14].astype(str).tolist()

    # Collect rows: [trait, snp, beta1, se1, p1, ..., beta5, se5, p5]
    res_rows = []

    for trait, snp in zip(trait_lst, snp_lst):
        row = [trait, snp]
        for group in range(1, 6):
            # Example file name pattern: "{trait}.{snp}.{group}.<digits>_<digits>.res"
            regex = re.compile(rf"^{re.escape(trait)}\.{re.escape(snp)}\.{group}\.\d+_\d+\.res$")
            try:
                fpath = newest_match(work_dir, regex)
            except FileNotFoundError as e:
                # You can choose to raise or skip; here we raise to be explicit
                raise FileNotFoundError(str(e))

            df = pd.read_csv(fpath, sep=r"\s+", header=0)
            # Expect columns: beta, se, p
            if not {"beta", "se", "p"}.issubset(df.columns):
                raise KeyError(f"{fpath} must contain columns: beta, se, p")

            beta, se, pval = (pd.to_numeric(df.loc[0, "beta"], errors="coerce"),
                              pd.to_numeric(df.loc[0, "se"],   errors="coerce"),
                              pd.to_numeric(df.loc[0, "p"],    errors="coerce"))
            if not np.isfinite([beta, se, pval]).all():
                raise ValueError(f"Non-finite beta/se/p in {fpath}")

            row.extend([float(beta), float(se), float(pval)])
        res_rows.append(row)

    # Build result DataFrame
    colnames = ["trait", "snp",
                "beta1", "se1", "p1",
                "beta2", "se2", "p2",
                "beta3", "se3", "p3",
                "beta4", "se4", "p4",
                "beta5", "se5", "p5"]
    df_res = pd.DataFrame(res_rows, columns=colnames)

    # Enforce dtypes: strings for id cols, float64 for numeric cols
    df_res[["trait", "snp"]] = df_res[["trait", "snp"]].astype(str)
    num_cols = [c for c in df_res.columns if c not in ("trait", "snp")]
    df_res[num_cols] = df_res[num_cols].apply(pd.to_numeric, errors="coerce").astype("float64")

    # Attach columns from dfS (Genes etc.)
    # Adjust these source column names to what you need:
    attach_cols = []
    for name in ["GenomicLocus_hg38", "Genes", "TraitName"]:
        if name in dfS.columns:
            attach_cols.append(name)
    if attach_cols:
        df_res = pd.concat([dfS.loc[:, attach_cols].reset_index(drop=True), df_res], axis=1)

    # --- Wald test for each row ---
    pvals = np.empty(len(df_res), dtype=np.float64)

    for i in range(len(df_res)):
        betas = df_res.loc[i, ["beta1", "beta2", "beta3", "beta4", "beta5"]].to_numpy(dtype=np.float64)
        ses   = df_res.loc[i, ["se1",   "se2",   "se3",   "se4",   "se5"  ]].to_numpy(dtype=np.float64)

        # Convert to proper shapes
        betas = to_float64_matrix(betas, "betas", make_2d="col")  # (5,1)
        ses   = to_float64_matrix(ses,   "ses",   make_2d=None)   # (5,)

        # Covariance matrix under independence assumption
        cov = np.diag(ses ** 2).astype(np.float64)               # (5,5)
        A = arr                                                  # (4,5)
        x = A @ betas                                            # (4,1)
        V = A @ cov @ A.T                                        # (4,4)

        # Debug (optional)
        if verbose and i == 0:
            print(f"rank(V) = {np.linalg.matrix_rank(V)}; min diag = {np.min(np.diag(V)):.3e}")

        # Prefer solve over explicit inverse; ridge if ill-conditioned
        try:
            chi2_stat = float((x.T @ np.linalg.solve(V, x)).squeeze())
        except LinAlgError:
            V_ridge = V + ridge_eps * np.eye(V.shape[0])
            chi2_stat = float((x.T @ np.linalg.solve(V_ridge, x)).squeeze())

        pvals[i] = chi2.sf(chi2_stat, df=A.shape[0])

    df_res["p_wald"] = pvals

    print(f"[p < 5e-8:", (df_res["p_wald"] < 5e-8).sum())
    print(f"[p < 0.05 / N:", (df_res["p_wald"] < 0.05 / df_res.shape[0]).sum())
    print(f"[p < 0.05:", (df_res["p_wald"] < 0.05).sum())

    p_col = [f"p{i}" for i in range(1, 6)]
    df_res["p_main_cutoff1"] = (df_res[p_col] < (5e-8 / 5)).sum(axis=1)
    df_res["p_main_cutoff2"] = (df_res[p_col] < (0.05 / df_res.shape[0] / 5)).sum(axis=1)
    df_res["p_main_cutoff3"] = (df_res[p_col] < 0.05/5).sum(axis=1)

    # Write output
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    df_res.to_csv(out_csv, index=False)
    if verbose:
        print(f"Wrote: {out_csv} | rows={len(df_res)}")


# ---------- Run ----------

# PT
process_trait_specific(
    "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.addP.csv",
    "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/AIE_split/PT_trait_specific.csv"
)

# BB
process_trait_specific(
    "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.addP.csv",
    "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/AIE_split/BB_trait_specific.csv"
)
