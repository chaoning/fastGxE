import argparse
import numpy as np
import statsmodels.api as sm
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')  # headless environments (servers/CI)
import matplotlib.pyplot as plt
from tqdm import tqdm


# ------------------------------- Plotting ----------------------------------- #
def qqplot(p_arr_lst, color_lst, label_lst, out_file, dpi=300, figsize_cm=(6, 5),
           output_format='pdf', title=None):
    """
    Generate a Q-Q plot for multiple p-value datasets with a 95% CI band.

    Parameters
    ----------
    p_arr_lst : list[np.ndarray]
        List of 1D arrays; each array contains p-values for a dataset (e.g., an R² level).
        All arrays should be the same length (e.g., number of replicates).
    color_lst : list[str]
        Colors for each dataset.
    label_lst : list[str]
        Legend labels for each dataset (same length/order as p_arr_lst).
    out_file : str
        Output filename (without extension).
    dpi : int
        Figure DPI.
    figsize_cm : tuple(float, float)
        Figure size in centimeters (width, height).
    output_format : str
        File format for saving (e.g., 'pdf' or 'png').
    """
    assert len(p_arr_lst) == len(color_lst) == len(label_lst), "List lengths must match."

    # Convert cm → inches for matplotlib
    figsize_in = (figsize_cm[0] / 2.54, figsize_cm[1] / 2.54)

    plt.rcParams['font.family'] = ['Arial']
    fig, ax = plt.subplots(figsize=figsize_in, dpi=dpi)

    # Use the length of the first dataset for expected quantiles and CI band
    n = len(p_arr_lst[0])
    if any(len(arr) != n for arr in p_arr_lst):
        raise ValueError("All p-value arrays must have the same length for a shared QQ envelope.")

    # Expected -log10(p) under the null (uniform[0,1]) using (i - 0.5)/n
    ranks = np.arange(1, n + 1)
    expected = -np.log10((ranks - 0.5) / n)

    # 95% CI for order statistics of the uniform distribution via Beta quantiles
    p_up = beta.ppf(0.975, ranks, n - ranks + 1)
    p_low = beta.ppf(0.025, ranks, n - ranks + 1)
    ax.fill_between(expected,
                    -np.log10(p_up),
                    -np.log10(p_low),
                    color='gray', alpha=0.2, label='_nolegend_')

    # Scatter each dataset
    max_val = expected.max()
    for p_vals, color, label in zip(p_arr_lst, color_lst, label_lst):
        # Numerical safety: clip(1e-300, 1) to avoid log(0); sort for QQ
        pv = np.sort(np.clip(np.asarray(p_vals, dtype=float), 1e-300, 1.0))
        observed = -np.log10(pv)
        ax.scatter(expected, observed, s=6, color=color, marker='.', alpha=1, label=label, linewidths=0.5)

    # Reference y = x line
    ax.plot([0, max_val], [0, max_val], linestyle='--', color='black', linewidth=1)

    ax.set_xlabel('Expected -log10($p$)', fontsize=10)
    ax.set_ylabel('Observed -log10($p$)', fontsize=10)
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=6, handlelength=0.8, handletextpad=0.6, frameon=False)

    # Clean frame
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Add title if provided
    if title is not None:
        ax.set_title(title, fontsize=11, pad=8)

    fig.savefig(f"{out_file}.{output_format}", bbox_inches='tight', format=output_format)
    plt.close(fig)


# ------------------------------- Simulation --------------------------------- #
def simulate_pvalues_for_r2(r2, h2, sample_size, reps, rng):
    """
    For a given R² between true and proxy variables, run 'reps' simulations and
    return a vector of p-values for the GxE interaction term when fitting on proxies.

    Model:
        y = sqrt(h2)*G_true + sqrt(h2)*E_true + e,  Var(e) = 1 - 2*h2
    Design (fitted):
        y ~ 1 + G_proxy + E_proxy + (G_proxy * E_proxy)

    Parameters
    ----------
    r2 : float
        Target R² between (G_true, G_proxy) and similarly (E_true, E_proxy). We set corr = sqrt(R²).
    h2 : float
        Per-component heritability for both G and E contributions (0 < 2*h2 < 1).
    sample_size : int
        Number of individuals per replicate.
    reps : int
        Number of replicated datasets to simulate at this R² level.
    rng : np.random.Generator
        Numpy random generator for reproducibility.

    Returns
    -------
    np.ndarray
        1D array of length 'reps' with the p-values for the interaction term.
    """
    # Correlation between true and proxy is sqrt(R²)
    rho = np.sqrt(max(0.0, min(1.0, r2)))
    cov = np.array([[1.0, rho],
                    [rho, 1.0]])

    # Precompute noise std so Var(y) ≈ 1
    noise_sd = np.sqrt(max(1e-12, 1.0 - 2.0 * h2))

    pvals = np.empty(reps, dtype=float)

    for rep in range(reps):
        # Draw (G_true, G_proxy) and (E_true, E_proxy) jointly
        G_joint = rng.multivariate_normal(mean=[0.0, 0.0], cov=cov, size=sample_size)
        E_joint = rng.multivariate_normal(mean=[0.0, 0.0], cov=cov, size=sample_size)
        G_true, G_proxy = G_joint[:, 0], G_joint[:, 1]
        E_true, E_proxy = E_joint[:, 0], E_joint[:, 1]

        # Outcome: additive genetic + additive environment + noise
        y = np.sqrt(h2) * G_true + np.sqrt(h2) * E_true + rng.normal(loc=0.0, scale=noise_sd, size=sample_size)

        # Design matrix using proxies, include intercept
        GxE_proxy = G_proxy * E_proxy
        X = np.column_stack([G_proxy, E_proxy, GxE_proxy])
        X = sm.add_constant(X, has_constant='add')  # adds intercept as first column

        # OLS and p-value for the interaction coefficient (index 3: const, G, E, GxE)
        model = sm.OLS(y, X).fit()
        pvals[rep] = float(model.pvalues[3])

    # Numerical safety
    pvals = np.clip(pvals, 1e-300, 1.0)
    return pvals


# ---------------------------------- Main ------------------------------------ #
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="QQ plots of GxE p-values under proxy attenuation across R² levels.")
    parser.add_argument("--h2_pct", type=float, required=True,
                        help="Per-component h^2 as percentage (e.g., 5 means h2 = 0.05 for G and 0.05 for E). Require 2*h2 < 1.")
    parser.add_argument("--n", type=int, default=100_000, help="Sample size per replicate (default: 100000).")
    parser.add_argument("--reps", type=int, default=100, help="Number of replicates per R^2 level (default: 100).")
    parser.add_argument("--out", type=str, default=None, help="Output figure prefix (default: auto from h2).")
    parser.add_argument("--fmt", type=str, default="png", choices=["png", "pdf", "svg"], help="Output format (default: png).")
    parser.add_argument("--seed", type=int, default=123, help="Base RNG seed (default: 1).")
    args = parser.parse_args()

    h2 = args.h2_pct / 100.0
    if not (0.0 < 2.0 * h2 < 1.0):
        raise ValueError("Require 0 < 2*h2 < 1 to keep Var(y) positive. Try a smaller h2_pct.")

    sample_size = int(args.n)
    reps = int(args.reps)
    out_prefix = args.out or f"qqplot_h2_{str(args.h2_pct)}"
    output_format = args.fmt

    # Color & label palette (ordered by R² list)
    r2_levels = [0.0, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9]
    color_lst = ["#1e466e", "#528fad", "#aadce0", "#ffe6b7", "#f2b705", "#d97c00", "#a94d00"]
    label_lst = [f"R²={r2:g}" for r2 in r2_levels]

    # Single RNG for reproducibility; advance seed per R² to vary streams
    base_rng = np.random.default_rng(args.seed)

    # Collect p-values per R² (one array per dataset for the QQ function)
    p_arr_lst = []
    for i, r2 in enumerate(tqdm(r2_levels, desc="Simulating across R² levels")):
        rng = np.random.default_rng(base_rng.integers(1, 2**31 - 1))
        pvals = simulate_pvalues_for_r2(r2=r2, h2=h2, sample_size=sample_size, reps=reps, rng=rng)
        p_arr_lst.append(pvals)

    # Plot QQ
    qqplot(
        p_arr_lst=p_arr_lst,
        color_lst=color_lst,
        label_lst=label_lst,
        out_file=out_prefix,
        output_format=output_format,
        figsize_cm=(6, 5),
        dpi=300,
        title=f"h²={str(args.h2_pct)}%"
    )
