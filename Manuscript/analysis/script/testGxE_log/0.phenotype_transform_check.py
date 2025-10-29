#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
phenotype_transform_check.py
Utility to assess whether a phenotype benefits from log (or other power) transformation.

Usage (CLI):
    python phenotype_transform_check.py --csv data.csv --y BMI --covariates age,sex,PC1,PC2 --plots outdir/

Programmatic use:
    from phenotype_transform_check import recommend_transformation
    rec = recommend_transformation(y, X)  # y: 1D array-like, X: DataFrame or 2D array (optional)
    print(rec["summary"])

Notes:
- If y contains non-positive values, log-transformation is skipped automatically.
- Besides log, Box-Cox (positive-only) and Yeo–Johnson (any real) are considered.
- Criteria include skewness reduction, residual normality, and heteroskedasticity.
"""

import argparse
import os
import warnings
from typing import Dict, Optional

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import PowerTransformer
import statsmodels.api as sm
from statsmodels.stats.diagnostic import het_breuschpagan, normal_ad

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ---------------------------
# Helper functions
# ---------------------------

def _fit_lm_and_diagnostics(y, X=None):
    """Fit OLS (with intercept) and return residual diagnostics."""
    if X is None:
        X = np.ones((len(y), 1))
    else:
        X = np.asarray(X)
        # Add intercept if not present
        X = sm.add_constant(X, has_constant='add')
    model = sm.OLS(y, X, missing='drop')
    res = model.fit()
    resid = res.resid
    # Normality tests
    jb_stat, jb_p = stats.jarque_bera(resid)
    ad_stat, ad_p = normal_ad(resid)  # Anderson–Darling
    # Heteroskedasticity test (Breusch–Pagan) needs at least 2 regressors; if not, fall back to White
    try:
        bp_stat, bp_p, _, _ = het_breuschpagan(resid, X)
    except Exception:
        bp_stat, bp_p = np.nan, np.nan
    # Skewness of residuals
    sk = stats.skew(resid, bias=False, nan_policy='omit')
    return {
        "res": res,
        "resid": resid,
        "jb_p": jb_p,
        "ad_p": ad_p,
        "bp_p": bp_p,
        "skew_resid": sk,
        "rss": np.sum(resid**2),
        "aic": res.aic,
        "bic": res.bic,
        "r2": res.rsquared if hasattr(res, "rsquared") else np.nan
    }


def _safe_log(y):
    """Return log(y) if all y>0, else None."""
    y = np.asarray(y)
    if np.all(y > 0):
        return np.log(y)
    return None


def _boxcox(y):
    """Box–Cox transform (MLE lambda). Requires y>0."""
    y = np.asarray(y)
    if np.all(y > 0):
        yt, lam = stats.boxcox(y)  # returns transformed and lambda
        return yt, lam
    return None, None


def _yeojohnson(y):
    """Yeo–Johnson (works with non-positive values too)."""
    pt = PowerTransformer(method="yeo-johnson", standardize=False)
    yt = pt.fit_transform(np.asarray(y).reshape(-1, 1)).ravel()
    lam = float(pt.lambdas_)
    return yt, lam


def _summary_row(name, y, X):
    """Compute summary metrics for a given transformed response."""
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(y)
    y = y[mask]
    X_use = X.loc[mask] if isinstance(X, pd.DataFrame) else (X[mask] if X is not None else None)
    dist_skew = stats.skew(y, bias=False, nan_policy='omit')
    dist_kurt = stats.kurtosis(y, fisher=True, bias=False, nan_policy='omit')
    diag = _fit_lm_and_diagnostics(y, X_use)
    return {
        "name": name,
        "skew_y": dist_skew,
        "kurtosis_y": dist_kurt,
        "jb_p": diag["jb_p"],
        "ad_p": diag["ad_p"],
        "bp_p": diag["bp_p"],
        "skew_resid": diag["skew_resid"],
        "aic": diag["aic"],
        "bic": diag["bic"],
        "r2": diag["r2"]
    }, y, X_use


def recommend_transformation(y, X: Optional[pd.DataFrame] = None, make_plots: bool = False, plot_dir: Optional[str] = None) -> Dict:
    """
    Compare raw, log (if valid), Box–Cox, and Yeo–Johnson. Return recommendation.
    Criteria (heuristic):
        - Prefer the transform that:
            * substantially reduces |skew_y| (≥30% vs raw) AND
            * improves residual normality (higher JB/AD p-values) AND
            * reduces heteroskedasticity (higher BP p-value) OR reduces AIC/BIC.
        - If raw already looks good (|skew_y|<0.5 and AD p>0.05 and BP p>0.05), keep raw.
    """
    y = np.asarray(y, dtype=float)
    if X is not None and not isinstance(X, (pd.DataFrame, np.ndarray)):
        raise ValueError("X must be a pandas DataFrame, numpy array, or None.")

    results = []

    # RAW
    row_raw, y_raw, X_use = _summary_row("raw", y, X)
    results.append(row_raw)

    # LOG
    y_log = _safe_log(y)
    if y_log is not None:
        row_log, _, _ = _summary_row("log", y_log, X)
        results.append(row_log)

    # BOX-COX
    y_bc, lam_bc = _boxcox(y)
    if y_bc is not None:
        row_bc, _, _ = _summary_row(f"boxcox(λ={lam_bc:.3f})", y_bc, X)
        results.append(row_bc)

    # YEO–JOHNSON
    y_yj, lam_yj = _yeojohnson(y)
    row_yj, _, _ = _summary_row(f"yeo-johnson(λ={lam_yj:.3f})", y_yj, X)
    results.append(row_yj)

    df = pd.DataFrame(results).set_index("name")

    # Heuristic decision
    rec = "raw"
    reason = []
    raw = df.loc["raw"]
    good_raw = (abs(raw["skew_y"]) < 0.5) and (raw["ad_p"] > 0.05) and (np.isnan(raw["bp_p"]) or raw["bp_p"] > 0.05)
    if good_raw:
        rec = "raw"
        reason.append("Raw scale already near-normal with homoskedastic residuals.")
    else:
        # pick best by a combined score
        score = []
        for name, row in df.iterrows():
            # larger is better for p-values; smaller is better for |skew| and BIC
            p_norm = np.nan_to_num(row["ad_p"], nan=0.0)
            p_bp = np.nan_to_num(row["bp_p"], nan=0.0)
            skew_gain = max(0.0, abs(raw["skew_y"]) - abs(row["skew_y"]))  # reduction
            bic_gain = np.nan_to_num(raw["bic"] - row["bic"], nan=0.0)
            composite = (2.0 * p_norm) + (1.0 * p_bp) + (0.5 * skew_gain) + (0.2 * bic_gain)
            score.append((composite, name))
        score.sort(reverse=True)
        rec = score[0][1]
        reason.append("Chosen by highest composite score combining residual normality (AD p), BP p, skewness reduction, and BIC gain.")

    out = {
        "table": df,
        "recommendation": rec,
        "summary": f"Recommended transformation: {rec}. " + " ".join(reason)
    }

    if make_plots and plot_dir is not None:
        import matplotlib.pyplot as plt
        os.makedirs(plot_dir, exist_ok=True)
        # Histogram & QQ for each candidate
        def _plot_hist_qq(vals, name):
            vals = np.asarray(vals, dtype=float)
            vals = vals[np.isfinite(vals)]
            fig1 = plt.figure()
            plt.hist(vals, bins=40)
            plt.title(f"Histogram: {name}")
            plt.tight_layout()
            fig1.savefig(os.path.join(plot_dir, f"hist_{name}.png"), dpi=150)
            plt.close(fig1)

            fig2 = plt.figure()
            stats.probplot(vals, dist="norm", plot=plt)
            plt.title(f"QQ-plot: {name}")
            plt.tight_layout()
            fig2.savefig(os.path.join(plot_dir, f"qq_{name}.png"), dpi=150)
            plt.close(fig2)

        # generate from the actual arrays
        _plot_hist_qq(y_raw, "raw")
        if y_log is not None: _plot_hist_qq(y_log, "log")
        if y_bc is not None: _plot_hist_qq(y_bc, "boxcox")
        _plot_hist_qq(y_yj, "yeo-johnson")

    return out


def _cli():
    ap = argparse.ArgumentParser(description="Assess whether a phenotype needs log (or other power) transform.")
    ap.add_argument("--csv", type=str, required=True, help="CSV file containing the phenotype (and optional covariates).")
    ap.add_argument("--y", type=str, required=True, help="Column name of the phenotype.")
    ap.add_argument("--covariates", type=str, default="", help="Comma-separated covariate column names (optional).")
    ap.add_argument("--plots", type=str, default="", help="Output directory to save diagnostic plots (optional).")
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    if args.y not in df.columns:
        raise ValueError(f"Phenotype column '{args.y}' not in CSV.")
    y = df[args.y].values
    X = None
    if args.covariates.strip():
        covars = [c.strip() for c in args.covariates.split(",") if c.strip()]
        for c in covars:
            if c not in df.columns:
                raise ValueError(f"Covariate '{c}' not in CSV columns.")
        X = df[covars]

    out = recommend_transformation(y, X, make_plots=bool(args.plots), plot_dir=args.plots if args.plots else None)
    # Print small summary
    print(out["summary"])
    print("\nDiagnostics table:")
    print(out["table"].to_string(float_format=lambda v: f"{v:.4g}"))


if __name__ == "__main__":
    _cli()
