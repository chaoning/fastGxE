import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tqdm import tqdm


def _safe_max_minus_log10(series, clip_min=1e-323):
    """
    Return the finite maximum of -log10(p) from a pandas Series of p-values.
    Coerces to numeric, clips to avoid log10(0), and ignores non-finite values.
    Returns np.nan if no finite values remain.
    """
    arr = pd.to_numeric(series, errors="coerce").to_numpy()
    if arr.size == 0:
        return np.nan
    arr = np.clip(arr, clip_min, None)
    logp = -np.log10(arr)
    logp = logp[np.isfinite(logp)]
    if logp.size == 0:
        return np.nan
    return float(logp.max())


def manplot(
    file_lst,
    chro, pos, p,
    color_vec=("#ef8a47", "#376795"),
    scattersize=1.5,
    titleName=None,
    threshold=(1.0e-6, "black"),
    logp=True,
    figsize=(10, 5),
    dpi=600,
    font='Arial',
    xlabel=(None, 8),
    ylabel=('-log10($p$)', 8),
    ticksize=(8, 8),
    chro_show=None,
    out_file="test",
    ylim=(0, None),
    ax=None,
):
    """
    Draw a Manhattan plot from one or more summary-stat files.
    Supports plotting into an external matplotlib axis (ax).
    """

    # Set font family
    plt.rc('font', family=font)

    # Create new figure if no axis is given
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(figsize[0] / 2.54, figsize[1] / 2.54), dpi=dpi)
        created_fig = True
    else:
        fig = ax.figure

    # Remove top/right borders
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Whether x ticks/labels were set
    xticks_set = False

    # Keep global x-span for threshold lines
    xmin_global, xmax_global = None, None

    # Iterate over all files (supports overlays)
    for file in tqdm(file_lst, disable=(len(file_lst) == 1)):
        # Read necessary columns only
        df = pd.read_csv(file, sep=r"\s+", usecols=[chro, pos, p])

        # Ensure chromosome column is string type
        df[chro] = df[chro].astype(str)

        # Preserve observed chromosome order
        chro_set = list(pd.unique(df[chro]))
        df[chro] = df[chro].astype("category").cat.reorder_categories(chro_set, ordered=True)

        # Sort by chromosome and position
        df = df.sort_values([chro, pos])

        # Concatenate genome-wide x-positions
        pos_accum = [0.0]
        xtick_vec = []
        for i, chr_name in enumerate(chro_set):
            mask = (df[chro] == chr_name)
            posx = df.loc[mask, pos].to_numpy(dtype=float) + pos_accum[-1]
            pvaly = pd.to_numeric(df.loc[mask, p], errors="coerce").to_numpy(dtype=float)

            # Transform to -log10(p) if requested
            if logp:
                pvaly = -np.log10(np.clip(pvaly, 1e-323, None))

            # Keep only finite pairs to avoid NaNs/Infs sneaking in
            finite_mask = np.isfinite(posx) & np.isfinite(pvaly)
            posx = posx[finite_mask]
            pvaly = pvaly[finite_mask]

            if posx.size:
                ax.scatter(posx, pvaly, s=scattersize, c=color_vec[i % len(color_vec)], linewidths=0)

            # Update accumulators and xtick centers (use last valid posx if any)
            last_max = float(posx.max()) if posx.size else pos_accum[-1]
            pos_accum.append(last_max)
            xtick_vec.append((pos_accum[-1] + pos_accum[-2]) / 2.0)

        # Update global x-span
        xmin, xmax = pos_accum[0], pos_accum[-1]
        xmin_global = xmin if xmin_global is None else min(xmin_global, xmin)
        xmax_global = xmax if xmax_global is None else max(xmax_global, xmax)

        # Set xticks/labels only once
        if not xticks_set:
            xtick_vec_show, chro_vec_show = [], []
            if chro_show is None:
                # Auto subset chromosomes
                for i, chr_name in enumerate(chro_set):
                    if i <= 10 or (i % 2 == 1):
                        xtick_vec_show.append(xtick_vec[i])
                        chro_vec_show.append(str(chr_name))
            else:
                desired = [str(x) for x in chro_show]
                for chr_name in chro_set:
                    if str(chr_name) in desired:
                        idx = chro_set.index(chr_name)
                        xtick_vec_show.append(xtick_vec[idx])
                        chro_vec_show.append(str(chr_name))

            ax.set_xticks(xtick_vec_show)
            ax.set_xticklabels(chro_vec_show, fontsize=ticksize[0])
            xticks_set = True

    # Draw threshold lines
    if threshold[1] is not None and xmin_global is not None and xmax_global is not None:
        thr_val = -np.log10(threshold[0]) if logp else threshold[0]
        ax.hlines(thr_val, xmin_global, xmax_global, colors=threshold[1], linestyles='dashdot', linewidth=1)
        gwas_thr = -np.log10(5e-8) if logp else 5e-8
        ax.hlines(gwas_thr, xmin_global, xmax_global, colors=threshold[1], linestyles='dashdot', linewidth=1)

    # Set y-axis ticks and labels
    ax.tick_params(axis='y', labelsize=ticksize[1])
    if xlabel[0] is not None:
        ax.set_xlabel(xlabel[0], fontsize=xlabel[1])
    if ylabel is not None and ylabel[0] is not None:
        ax.set_ylabel(ylabel[0], fontsize=ylabel[1])

    # Apply y limits only if finite and valid
    ylow, yhigh = ylim
    if (ylow is not None or yhigh is not None):
        if (ylow is None or np.isfinite(ylow)) and (yhigh is None or np.isfinite(yhigh)):
            # If only one bound is given, matplotlib handles the other as auto
            ax.set_ylim(ylow, yhigh)

    # Title
    if titleName is not None:
        ax.set_title(str(titleName))

    # Save figure only if created inside
    if created_fig:
        fig.savefig(out_file + ".man.png", format='png', dpi=dpi)
        plt.close(fig)
    return ax


if __name__ == "__main__":
    trait = sys.argv[1]
    ShortName = sys.argv[2]
    chrom = "chrom"
    base = "base"

    file_lst = [f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/log_trans/vs/{trait}.txt"]
    pcut = 5e-8 / 32.0
    out_prefix = f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/log_trans/vs/{trait}"
    chro_show = list(range(1, 13)) + [14, 16, 18, 20]
    chro_show = [str(x) for x in chro_show]

    print(trait, ShortName, "p_gxe & p_gxe_log")

    # === Determine ylim_up dynamically (robust) ===
    df_check = pd.read_csv(file_lst[0], sep=r"\s+", usecols=[chrom, base, "p_gxe", "p_gxe_log"])
    max_val_1 = _safe_max_minus_log10(df_check["p_gxe"])
    max_val_2 = _safe_max_minus_log10(df_check["p_gxe_log"])

    if np.isnan(max_val_1) and np.isnan(max_val_2):
        print("[WARN] No finite p-values in both columns; falling back to ylim_up=10.0")
        ylim_up = 10.0
    else:
        ylim_up = np.nanmax([max_val_1, max_val_2]) * 1.05  # 5% headroom
        if not np.isfinite(ylim_up) or ylim_up <= 0:
            print("[WARN] Computed ylim_up is invalid; falling back to 10.0")
            ylim_up = 10.0

    # === Create side-by-side figure ===
    figsize_cm = (20, 5)
    dpi = 600
    font = 'Arial'

    fig, axes = plt.subplots(1, 2, figsize=(figsize_cm[0] / 2.54, figsize_cm[1] / 2.54),
                             dpi=dpi, sharey=True)

    # Left panel: p_gxe
    manplot(
        file_lst=file_lst,
        chro=chrom, pos=base, p="p_gxe",
        titleName=f"{ShortName} (INT)",
        threshold=(pcut, "black"),
        chro_show=chro_show,
        ylim=(0.0, ylim_up),
        ylabel=('-log10($p$)', 8),
        ax=axes[0],
        font=font
    )

    # Right panel: p_gxe_log
    manplot(
        file_lst=file_lst,
        chro=chrom, pos=base, p="p_gxe_log",
        titleName=f"{ShortName} (log)",
        threshold=(pcut, "black"),
        chro_show=chro_show,
        ylim=(0.0, ylim_up),
        ylabel=(None, 8),
        ax=axes[1],
        font=font
    )

    # Set common x label
    for ax in axes:
        ax.set_xlabel("Chromosome", fontsize=8)

    plt.tight_layout()
    fig.savefig(out_prefix + ".p_gxe_vs_p_gxe_log.dual.man.png", format='png', dpi=dpi)
    plt.close(fig)
