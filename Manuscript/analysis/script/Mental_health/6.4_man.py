import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tqdm import tqdm


def manplot(
    file_lst,
    chro, pos, p,
    color_vec=("#ef8a47", "#376795"),
    scattersize=1.5,
    titleName=None,
    threshold=(1.0e-6, "black"),
    logp=True,
    figsize=(13, 5),
    dpi=600,
    font='Arial',
    xlabel=(None, 8),
    ylabel=('-log10($p$)', 8),
    ticksize=(8, 8),
    chro_show=None,
    out_file="test",
    ylim=(0, None)
):
    """
    Draw a Manhattan plot from one or more summary-stat files.

    Parameters
    ----------
    file_lst : list[str]
        Paths to input files. Each file must contain columns `chro`, `pos`, and `p`.
    chro, pos, p : str
        Column names for chromosome, base position, and p-value.
    color_vec : tuple[str]
        Colors rotated by chromosome.
    scattersize : float
        Marker size for scatter points.
    titleName : str or None
        Figure title.
    threshold : (float, str or None)
        (p-value threshold, color). If color is None, skip drawing the line.
    logp : bool
        If True, plot -log10(p); otherwise plot raw p.
    figsize : (float, float)
        Figure size in centimeters.
    dpi : int
        Output DPI.
    font : str
        Font family used by Matplotlib.
    xlabel : (str or None, int)
        X label text and font size. If text is None, omit label.
    ylabel : (str, int)
        Y label text and font size.
    ticksize : (int, int)
        Font sizes for x and y ticks.
    chro_show : list or None
        Specific chromosome labels to show (as strings or ints). If None, auto-subsample.
    out_file : str
        Output file prefix (".man.png" will be appended).
    ylim : (float or None, float or None)
        Y-axis limits (min, max). Use None for auto.
    """
    # Set font globally for this figure
    plt.rc('font', family=font)

    # Create figure/axis
    fig, ax = plt.subplots(figsize=(figsize[0] / 2.54, figsize[1] / 2.54), dpi=dpi)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Whether we've already set x ticks/labels (do it once, from the first file)
    xticks_set = False

    # Loop over files to plot (supports overlay)
    for file in tqdm(file_lst):
        # Read minimal columns; use raw-string for regex separator
        df = pd.read_csv(file, sep=r"\s+", usecols=[chro, pos, p])

        # Force chromosome column to string to avoid category/index mismatches
        df[chro] = df[chro].astype(str)

        # Build ordered categories in the observed order
        chro_set = list(pd.unique(df[chro]))
        df[chro] = df[chro].astype("category").cat.reorder_categories(chro_set, ordered=True)

        # Sort by chromosome then base position
        df = df.sort_values([chro, pos])

        # Accumulate genome-wide x-positions by concatenating chromosomes
        pos_accum = [0.0]
        xtick_vec = []
        for i, chr_name in enumerate(chro_set):
            mask = (df[chro] == chr_name)
            posx = df.loc[mask, pos].to_numpy(dtype=float) + pos_accum[-1]
            pvaly = df.loc[mask, p].to_numpy(dtype=float)

            if logp:
                # Avoid log10(0); clip to tiny positive
                pvaly = -np.log10(np.clip(pvaly, 1e-323, None))

            ax.scatter(posx, pvaly, s=scattersize, c=color_vec[i % len(color_vec)], linewidths=0)

            # Update accumulators and x tick centers
            pos_accum.append(float(posx.max()) if posx.size else pos_accum[-1])
            xtick_vec.append((pos_accum[-1] + pos_accum[-2]) / 2.0)

        # Only set ticks/labels from the first file (shared layout for overlays)
        if not xticks_set:
            xtick_vec_show = []
            chro_vec_show = []

            if chro_show is None:
                # Auto-subset: show first ~11 chromosomes; then show every other to declutter
                for i, chr_name in enumerate(chro_set):
                    if i <= 10 or (i % 2 == 1):
                        xtick_vec_show.append(xtick_vec[i])
                        chro_vec_show.append(str(chr_name))
            else:
                # Respect the user-specified list (can be int or str)
                # Cast to str for consistent matching
                desired = [str(x) for x in chro_show]
                for chr_name in chro_set:
                    if str(chr_name) in desired:
                        idx = chro_set.index(chr_name)
                        xtick_vec_show.append(xtick_vec[idx])
                        chro_vec_show.append(str(chr_name))

            ax.set_xticks(xtick_vec_show)
            ax.set_xticklabels(chro_vec_show, fontsize=ticksize[0])
            xticks_set = True

        # Significance lines (draw for the full x-span of this file)
        if threshold[1] is not None:
            xmin, xmax = pos_accum[0], pos_accum[-1]
            thr_val = -np.log10(threshold[0]) if logp else threshold[0]
            ax.hlines(thr_val, xmin, xmax, colors=threshold[1], linestyles='dashdot', linewidth=1)
            # GWAS conventional threshold
            # gwas_thr = -np.log10(5e-8) if logp else 5e-8
            # ax.hlines(gwas_thr, xmin, xmax, colors=threshold[1], linestyles='dashdot', linewidth=1)

    # Y ticks/fontsize
    ax.tick_params(axis='y', labelsize=ticksize[1])

    # Axis labels
    if xlabel[0] is not None:
        ax.set_xlabel(xlabel[0], fontsize=xlabel[1])
    ax.set_ylabel(ylabel[0], fontsize=ylabel[1])

    # Y limits (set once, outside loop to avoid overrides in multi-file mode)
    if ylim[0] is not None or ylim[1] is not None:
        ax.set_ylim(ylim[0], ylim[1])

    # Title and save
    if titleName is not None:
        ax.set_title(str(titleName))
    fig.savefig(out_file + ".man.png", format='png', dpi=dpi)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/testGxE/")
    chrom = "chrom"
    base = "base"
    p = 'p_gxe'

    dfT = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/description/trait_Description.txt", sep="\t")
    trait_lst = dfT.iloc[:, 0].to_list()
    file_lst = [f"{trait}.res" for trait in trait_lst]
    pcut = 5e-8 / 36
    out_prefix = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Mental_health/6.4_man"
    chro_show = list(range(1, 13)) + [14, 16, 18, 20]
    chro_show = [str(x) for x in chro_show]
    manplot(
        file_lst=file_lst,
        chro=chrom, pos=base, p=p,
        threshold=(pcut, "black"),
        out_file=out_prefix,
        chro_show=chro_show
    )
