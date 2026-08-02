"""Simulate a GxE phenotype for the fastGxE example.

The phenotype follows the model fastGxE / mmSuSiE assume,

    y = env_main + g_poly + gxe_background + gxe_focus + residual,

whose variance decomposition mirrors

    V = h2_g · K  +  h2_gxe · (K ∘ EE'/m)  +  h2_nxe · diag(‖e‖²/m)  +  h2_e · I,

built from *discrete causal variants* (not drawn from the covariance directly) so
that fastGxE has real SNPs to detect and mmSuSiE has a real locus to fine-map:

  - ``env_main``       — small environmental main effects (variance ``h2_env``).
  - ``g_poly``         — additive polygenic background from ``n_snp_g`` random SNPs
                         (variance ``h2_g``); approximates ``h2_g · K``.
  - ``gxe_background`` — polygenic GxE from ``n_snp_gxe`` random SNPs each interacting
                         with a few environments (variance ``h2_gxe_bg``); approximates
                         the ``h2_gxe · (K ∘ EE')`` background you want to *control for*.
  - ``gxe_focus``      — ONE focal SNP interacting with ``n_env_focus`` environments
                         (variance ``h2_gxe_focus``); this is the signal fastGxE should
                         flag and mmSuSiE should resolve to its driving environments.
  - ``nxe``            — OPTIONAL noise-by-environment (heteroscedastic) noise whose
                         per-individual variance scales with ``‖e_i‖²`` (variance
                         ``h2_nxe``, default 0); the ``h2_nxe · diag(‖e‖²/m)`` term.
  - ``residual``       — i.i.d. noise (the remaining variance, so the parts sum to 1).

The focal SNP defaults to ``rs550011`` to match the tutorial; every draw is seeded, and
the focal SNP / environments / effect sizes are written to a truth file for validation.

Note: ``h2_gxe_focus`` (default 0.20) is deliberately inflated so a *single* locus is
detectable in this small example (n≈427). Real single-SNP GxE effects are far smaller.

Usage:
    python simu.py                       # writes test_simu_pheno.txt + test_simu_truth.txt
    python simu.py --h2-gxe-focus 0.10 --focus-snp rs550011 --seed 2025
"""
import argparse
import os

import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _read_standardized(bed, indices):
    """Read the given SNP columns, mean-impute, and standardize to unit variance."""
    indices = [int(i) for i in np.atleast_1d(indices)]
    arr = np.asarray(bed[:, indices].read().val, dtype=float)
    col_mean = np.nanmean(arr, axis=0)
    nan_r, nan_c = np.where(np.isnan(arr))
    arr[nan_r, nan_c] = col_mean[nan_c]
    std = arr.std(axis=0, ddof=0)
    std[std < 1e-12] = 1.0
    return (arr - arr.mean(axis=0)) / std


def _scale_to_var(x, target_var):
    """Rescale ``x`` to have exactly ``target_var`` empirical variance (0 -> zeros)."""
    if target_var <= 0:
        return np.zeros_like(x, dtype=float)
    sd = np.std(x)
    if sd < 1e-12:
        return np.zeros_like(x, dtype=float)
    return x / sd * np.sqrt(target_var)


def _random_env_correlation(n_env, rng, max_offdiag=0.2):
    """A valid (PD, unit-diagonal) environmental correlation matrix."""
    corr = np.eye(n_env)
    iu = np.triu_indices(n_env, k=1)
    corr[iu] = rng.uniform(-max_offdiag, max_offdiag, size=iu[0].size)
    corr = corr + corr.T - np.diag(np.diag(corr))
    # Project to the nearest PD correlation matrix if needed.
    eigval, eigvec = np.linalg.eigh(corr)
    if eigval.min() < 1e-8:
        corr = (eigvec * np.clip(eigval, 1e-8, None)) @ eigvec.T
        d = np.sqrt(np.diag(corr))
        corr = corr / np.outer(d, d)
    corr = np.clip((corr + corr.T) / 2.0, -1.0, 1.0)
    np.fill_diagonal(corr, 1.0)
    return corr


# --------------------------------------------------------------------------- #
# simulation
# --------------------------------------------------------------------------- #
def simulate_gxe_phenotype(
    bed_file,
    n_env=40,
    h2_env=0.01,
    h2_g=0.30,
    h2_gxe_bg=0.10,
    h2_gxe_focus=0.20,
    n_snp_g=1000,
    n_snp_gxe=500,
    n_env_focus=2,
    focus_snp="rs550011",
    h2_nxe=0.0,
    env_max_corr=0.2,
    seed=2025,
):
    """Simulate a GxE phenotype; return ``(pheno_df, truth_df)``.

    ``pheno_df`` columns: ``iid``, ``pheno``, ``E1 … E{n_env}`` (standardized).
    ``focus_snp`` is a SNP id (matched in the .bim) or ``None`` to pick one at random.
    """
    h2_resid = 1.0 - (h2_env + h2_g + h2_gxe_bg + h2_gxe_focus + h2_nxe)
    if h2_resid < 0:
        raise ValueError(
            f"variance components sum to {1 - h2_resid:.3f} > 1; reduce the h2_* values."
        )

    rng = np.random.default_rng(seed)
    bim = pd.read_csv(f"{bed_file}.bim", sep=r"\s+", header=None, dtype=str)
    fam = pd.read_csv(f"{bed_file}.fam", sep=r"\s+", header=None, dtype=str)
    iids = fam.iloc[:, 1].tolist()
    sids = bim.iloc[:, 1].tolist()
    sid_to_idx = {s: i for i, s in enumerate(sids)}
    bed = Bed(bed_file, count_A1=True)
    n_iid, n_snp = bed.iid_count, bed.sid_count

    # Resolve the focal SNP first so it is excluded from the polygenic pools.
    if focus_snp is not None:
        if focus_snp not in sid_to_idx:
            raise ValueError(f"focus_snp {focus_snp!r} not found in {bed_file}.bim")
        focus_idx = sid_to_idx[focus_snp]
    else:
        focus_idx = int(rng.choice(n_snp))
        focus_snp = sids[focus_idx]

    # Disjoint SNP pools for polygenic G, polygenic GxE, and the focal SNP.
    pool = np.setdiff1d(np.arange(n_snp), [focus_idx], assume_unique=True)
    g_idx = rng.choice(pool, size=min(n_snp_g, pool.size), replace=False)
    pool = np.setdiff1d(pool, g_idx, assume_unique=True)
    gxe_idx = rng.choice(pool, size=min(n_snp_gxe, pool.size), replace=False)

    # ---- environments ----
    env_corr = _random_env_correlation(n_env, rng, env_max_corr)
    env = rng.multivariate_normal(np.zeros(n_env), env_corr, size=n_iid)
    env = (env - env.mean(0)) / env.std(0, ddof=0)          # standardized columns

    # ---- environmental main effect ----
    env_main = _scale_to_var(env @ rng.normal(0, 1, n_env), h2_env)

    # ---- additive polygenic background ~ h2_g · K ----
    Gg = _read_standardized(bed, g_idx)
    g_poly = _scale_to_var(Gg @ rng.normal(0, 1, g_idx.size), h2_g)

    # ---- polygenic GxE background ~ h2_gxe · (K ∘ EE') ----
    # Each background SNP interacts with a few environments (the first environments are
    # favoured, as real exposures are); accumulate SNP∘env products.
    Ggxe = _read_standardized(bed, gxe_idx)
    env_weights = np.array([10.0] * min(10, n_env) + [1.0] * max(0, n_env - 10))
    env_weights /= env_weights.sum()
    gxe_bg = np.zeros(n_iid)
    for j in range(gxe_idx.size):
        k = 1 + int(rng.integers(0, n_env)) if n_env > 1 else 1
        envs = rng.choice(n_env, size=min(k, n_env), replace=False, p=env_weights)
        beta = rng.normal(0, 1, envs.size)
        gxe_bg += (Ggxe[:, [j]] * env[:, envs]) @ beta
    gxe_bg = _scale_to_var(gxe_bg, h2_gxe_bg)

    # ---- focal GxE: one SNP × a few environments (the signal to fine-map) ----
    g_focus = _read_standardized(bed, [focus_idx])[:, 0]
    focus_envs = np.sort(rng.choice(n_env, size=min(n_env_focus, n_env), replace=False))
    focus_beta = rng.normal(0, 1, focus_envs.size)
    gxe_focus_raw = (g_focus[:, None] * env[:, focus_envs]) @ focus_beta
    gxe_focus = _scale_to_var(gxe_focus_raw, h2_gxe_focus)

    # ---- noise-by-environment (NxE): heteroscedastic residual whose per-individual
    #      variance scales with ‖e_i‖²/m — the fastGxE NxE component diag(‖e‖²/m).
    #      Off by default (drawn only when h2_nxe > 0, so the default output is
    #      unchanged). ----
    if h2_nxe > 0:
        nxe_scale = np.sum(env ** 2, axis=1) / n_env            # ‖e_i‖²/m per individual
        nxe = _scale_to_var(rng.normal(0, 1, n_iid) * np.sqrt(nxe_scale), h2_nxe)
    else:
        nxe = np.zeros(n_iid)

    # ---- residual ----
    residual = _scale_to_var(rng.normal(0, 1, n_iid), h2_resid)

    pheno = env_main + g_poly + gxe_bg + gxe_focus + nxe + residual

    pheno_df = pd.DataFrame({"iid": iids, "pheno": pheno})
    env_df = pd.DataFrame(env, columns=[f"E{i + 1}" for i in range(n_env)])
    pheno_df = pd.concat([pheno_df, env_df], axis=1)

    # ---- truth (targets + realized variances) ----
    truth_rows = [
        {"component": "env_main", "target_var": h2_env, "realized_var": float(np.var(env_main))},
        {"component": "g_poly", "target_var": h2_g, "realized_var": float(np.var(g_poly))},
        {"component": "gxe_background", "target_var": h2_gxe_bg, "realized_var": float(np.var(gxe_bg))},
        {"component": "gxe_focus", "target_var": h2_gxe_focus, "realized_var": float(np.var(gxe_focus))},
        {"component": "nxe", "target_var": h2_nxe, "realized_var": float(np.var(nxe))},
        {"component": "residual", "target_var": h2_resid, "realized_var": float(np.var(residual))},
    ]
    for e, b in zip(focus_envs, focus_beta):
        truth_rows.append({
            "component": "focus_gxe_env", "snp_id": focus_snp, "snp_index0": focus_idx,
            "env": f"E{e + 1}", "beta": float(b), "target_var": h2_gxe_focus,
        })
    truth_df = pd.DataFrame(truth_rows)
    return pheno_df, truth_df


def _build_parser():
    p = argparse.ArgumentParser(description="Simulate a GxE phenotype for the fastGxE example.")
    p.add_argument("--workdir", default=os.path.dirname(os.path.abspath(__file__)),
                   help="Directory with <bed>.bed/.bim/.fam; outputs are written here.")
    p.add_argument("--bed-file", default="test", help="PLINK bed prefix inside workdir.")
    p.add_argument("--out-file", default="test_simu_pheno.txt", help="Phenotype output file.")
    p.add_argument("--truth-file", default="test_simu_truth.txt", help="Truth output file.")
    p.add_argument("--n-env", type=int, default=40)
    p.add_argument("--h2-env", type=float, default=0.01)
    p.add_argument("--h2-g", type=float, default=0.30)
    p.add_argument("--h2-gxe-bg", type=float, default=0.10)
    p.add_argument("--h2-gxe-focus", type=float, default=0.20,
                   help="Focal single-SNP GxE variance (inflated for demo detectability).")
    p.add_argument("--n-snp-g", type=int, default=1000)
    p.add_argument("--n-snp-gxe", type=int, default=500)
    p.add_argument("--n-env-focus", type=int, default=2)
    p.add_argument("--h2-nxe", type=float, default=0.0,
                   help="Noise-by-environment (heteroscedastic) variance; 0 disables it.")
    p.add_argument("--focus-snp", default="rs550011",
                   help="Focal SNP id (matched in .bim), or 'random'.")
    p.add_argument("--seed", type=int, default=2025)
    return p


def main():
    args = _build_parser().parse_args()
    bed_path = os.path.join(args.workdir, args.bed_file)
    focus = None if args.focus_snp.lower() == "random" else args.focus_snp

    pheno_df, truth_df = simulate_gxe_phenotype(
        bed_path, n_env=args.n_env, h2_env=args.h2_env, h2_g=args.h2_g,
        h2_gxe_bg=args.h2_gxe_bg, h2_gxe_focus=args.h2_gxe_focus,
        n_snp_g=args.n_snp_g, n_snp_gxe=args.n_snp_gxe, n_env_focus=args.n_env_focus,
        h2_nxe=args.h2_nxe, focus_snp=focus, seed=args.seed,
    )

    out_pheno = os.path.join(args.workdir, args.out_file)
    out_truth = os.path.join(args.workdir, args.truth_file)
    pheno_df.to_csv(out_pheno, sep=" ", index=False)
    truth_df.to_csv(out_truth, sep="\t", index=False)

    focus_rows = truth_df[truth_df["component"] == "focus_gxe_env"]
    print(f"Saved phenotype -> {out_pheno}  ({pheno_df.shape[0]} individuals, {args.n_env} environments)")
    print(f"Saved truth     -> {out_truth}")
    print(f"Focal GxE SNP   : {focus_rows['snp_id'].iloc[0]} "
          f"(index {int(focus_rows['snp_index0'].iloc[0])}) "
          f"× {list(focus_rows['env'])}")
    print("Variance (target -> realized):")
    for _, r in truth_df[truth_df["target_var"].notna()].drop_duplicates("component").iterrows():
        if "realized_var" in r and pd.notna(r.get("realized_var")):
            print(f"  {r['component']:14s} {r['target_var']:.3f} -> {r['realized_var']:.3f}")


if __name__ == "__main__":
    main()
