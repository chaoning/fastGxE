import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # non-interactive backend for headless save
import matplotlib.pyplot as plt
from tqdm import tqdm

# Lighter PDFs with selectable text (and better Illustrator compatibility)
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"]  = 42

# =========================
# Utilities
# =========================
CHR_MAP = {**{str(i): i for i in range(1, 23)}, "X": 23, "Y": 24, "XY": 25, "M": 26, "MT": 26}

def _normalize_chr(series: pd.Series):
    """
    Normalize chromosome IDs to numeric keys and pretty labels.
    Accepts 'chr1','1','X','Y','MT', etc.
    Returns:
      key  : int series (1..22, 23=X, 24=Y, 25=XY, 26=MT, others as increasing ints)
      label: str series ('1'..'22','X','Y','MT', or numeric fallback)
    """
    s = series.astype(str).str.replace("chr", "", regex=False).str.upper().str.strip()
    key = s.map(CHR_MAP)

    # Try numeric for unmapped (e.g., '23')
    mask = key.isna()
    if mask.any():
        num = pd.to_numeric(s[mask], errors="coerce")
        key[mask] = num

    # Fallback: for anything still NaN, assign stable new codes after max known
    if key.isna().any():
        base = int(np.nanmax(key.values.astype(float))) if np.isfinite(key.values.astype(float)).any() else 26
        rem_codes = s[key.isna()].astype("category").cat.codes + base + 1
        key.loc[key.isna()] = rem_codes

    key = key.astype(int)

    def _label_for(k):
        if 1 <= k <= 22: return str(k)
        return {23: "X", 24: "Y", 25: "XY", 26: "MT"}.get(k, str(k))

    label = key.map(_label_for)
    return key, label

def _read_concat_files(file_lst, chro, pos, p, p_floor=1e-60):
    """
    Read per-trait gz files: <trait>.gz with columns [chro, pos, p].
    Returns a normalized DataFrame with columns:
      chr_key, chr_label, bp, pval, trait
    """
    frames = []
    for file in tqdm(file_lst, desc="Reading files"):
        df = pd.read_csv(f"{file}.gz", sep=r"\s+", usecols=[chro, pos, p], low_memory=False)
        df.rename(columns={chro: "chr_raw", pos: "bp", p: "pval"}, inplace=True)

        # Normalize chromosome + sanitize numeric types
        chr_key, chr_label = _normalize_chr(df["chr_raw"])
        df["chr_key"]   = chr_key
        df["chr_label"] = chr_label
        df["bp"]   = pd.to_numeric(df["bp"], errors="coerce")
        df["pval"] = pd.to_numeric(df["pval"], errors="coerce")
        df.dropna(subset=["chr_key", "bp", "pval"], inplace=True)

        # Avoid -log10(0)
        df["pval"] = np.clip(df["pval"].to_numpy(), p_floor, 1.0)
        df["trait"] = os.path.basename(file)

        frames.append(df[["chr_key", "chr_label", "bp", "pval", "trait"]])

    if not frames:
        raise RuntimeError("No input frames were read. Check file list / paths.")

    return pd.concat(frames, ignore_index=True)

def _prepare_genome_axis(df, chro=None, pos=None, chro_show=None):
    """
    Compute genome-wide x positions, per-chromosome offsets (by max(bp)),
    tick centers and labels.

    Returns:
      axis_info dict with keys:
        chr_keys, offsets, tick_positions, tick_labels, x, x_max
    """
    chr_keys = sorted(df["chr_key"].unique())
    if not chr_keys:
        return {"chr_keys": [], "offsets": {}, "tick_positions": [], "tick_labels": [], "x": np.array([]), "x_max": 0.0}

    # Chromosome lengths by max coordinate
    lengths = df.groupby("chr_key")["bp"].max().reindex(chr_keys)
    ends = lengths.cumsum().to_numpy()
    starts = ends - lengths.to_numpy()
    offsets = dict(zip(chr_keys, starts))
    centers = (starts + lengths.to_numpy() / 2.0).tolist()

    # Which labels to show
    if chro_show is None:
        label_show_set = set(chr_keys)
    else:
        # Normalize input chro_show to keys
        tmp = pd.Series(chro_show).astype(str).str.replace("chr", "", regex=False).str.upper()
        tmp_key = tmp.map(CHR_MAP)
        tmp_key = pd.to_numeric(tmp_key.where(tmp_key.notna(), tmp), errors="coerce")
        label_show_set = set([int(k) for k in tmp_key.dropna().tolist() if int(k) in chr_keys])

    # Map key -> pretty label using first occurrence
    lab_map = df.drop_duplicates("chr_key").set_index("chr_key")["chr_label"].to_dict()
    tick_labels = [lab_map.get(k, str(k)) if k in label_show_set else "" for k in chr_keys]

    # Genome-wide x = bp + chromosome offset
    x = df["bp"].to_numpy() + df["chr_key"].map(offsets).to_numpy()

    axis_info = {
        "chr_keys": chr_keys,
        "offsets": offsets,
        "tick_positions": centers,
        "tick_labels": tick_labels,
        "x": x,
        "x_max": float(ends[-1]) if len(ends) else 0.0,
    }
    return axis_info

def _apply_axis(ax, axis_info, xlabel, ylabel, ticksize):
    ax.set_xticks(axis_info["tick_positions"])
    ax.set_xticklabels(axis_info["tick_labels"], fontsize=ticksize[0])
    ax.tick_params(axis='y', labelsize=ticksize[1])
    if xlabel[0]:
        ax.set_xlabel(xlabel[0], fontsize=xlabel[1])
    ax.set_ylabel(ylabel[0], fontsize=ylabel[1])
    ax.set_ylim(bottom=0)

def _sig_hline(ax, pcut, color, logp, x_end):
    sig_y = -np.log10(pcut) if logp else pcut
    ax.hlines(sig_y, 0, x_end, colors=color, linestyles='dashdot', linewidth=1)

def _chr_color_map(chr_keys, color_vec):
    return {c: color_vec[i % len(color_vec)] for i, c in enumerate(chr_keys)}

# =========================
# Plotters
# =========================
def manplot(file_lst, chro, pos, p, color_vec=("#e76254", "#1e466e"), scattersize=1.5, sig_color='red',
            threshold=(1.0e-6, "black"), logp=True,
            figsize=(13.5, 5), dpi=300, font='Arial', xlabel=(None, 8), ylabel=('-log10($p$)', 8),
            ticksize=(8, 8), chro_show=None, out_file="test", file_format="tiff", highlight_sig=True):

    plt.rc('font', family=font)
    fig, ax = plt.subplots(figsize=(figsize[0] / 2.54, figsize[1] / 2.54), dpi=dpi)
    for spine in ("right", "top"):
        ax.spines[spine].set_visible(False)

    df = _read_concat_files(file_lst, chro, pos, p)
    axis = _prepare_genome_axis(df, chro, pos, chro_show=chro_show)

    # Colors by chromosome
    chr_to_color = _chr_color_map(axis["chr_keys"], color_vec)

    # y values
    y = -np.log10(df["pval"].to_numpy()) if logp else df["pval"].to_numpy()

    # Background points
    ax.scatter(axis["x"], y, s=scattersize,
               c=[chr_to_color[c] for c in df["chr_key"]],
               rasterized=True, edgecolors="none")

    # Optional: highlight significant points on top
    if highlight_sig and threshold[1] is not None:
        sig_mask = (df["pval"].to_numpy() <= threshold[0])
        if np.any(sig_mask):
            ax.scatter(axis["x"][sig_mask], y[sig_mask], s=scattersize,
                       c=sig_color, rasterized=True, edgecolors="none", zorder=3)

    # Significance line
    if threshold[1] is not None:
        _sig_hline(ax, threshold[0], threshold[1], logp, axis["x_max"])

    _apply_axis(ax, axis, xlabel, ylabel, ticksize)

    out = f"{out_file}.{file_format}"
    plt.savefig(out, format=file_format, bbox_inches="tight")
    plt.close()

def manplot_frame(file_lst, chro, pos, p, color_vec=("#e76254", "#1e466e"), scattersize=1.5, sig_color='red',
                  threshold=(1.0e-6, "black"), logp=True,
                  figsize=(13.5, 5), dpi=300, font='Arial', xlabel=(None, 8), ylabel=('-log10($p$)', 8),
                  ticksize=(8, 8), chro_show=None, out_file="test", file_format="tiff"):

    plt.rc('font', family=font)
    fig, ax = plt.subplots(figsize=(figsize[0] / 2.54, figsize[1] / 2.54), dpi=dpi)
    for spine in ("right", "top"):
        ax.spines[spine].set_visible(False)

    df = _read_concat_files(file_lst, chro, pos, p)
    axis = _prepare_genome_axis(df, chro, pos, chro_show=chro_show)

    yvals = -np.log10(df["pval"].to_numpy()) if logp else df["pval"].to_numpy()
    df_plot = pd.DataFrame({"chr_key": df["chr_key"].values, "x": axis["x"], "y": yvals})

    for i, c in enumerate(axis["chr_keys"]):
        sub = df_plot[df_plot["chr_key"] == c].sort_values("x")
        if sub.empty:
            continue
        col = color_vec[i % len(color_vec)]
        # first & last
        ax.scatter(sub["x"].iloc[[0, -1]], sub["y"].iloc[[0, -1]],
                   s=scattersize, c=col, rasterized=True, edgecolors="none")
        # max y
        imax = sub["y"].idxmax()
        ax.scatter(sub.loc[imax, "x"], sub.loc[imax, "y"],
                   s=scattersize, c=col, rasterized=True, edgecolors="none")

    if threshold[1] is not None:
        _sig_hline(ax, threshold[0], threshold[1], logp, axis["x_max"])

    _apply_axis(ax, axis, xlabel, ylabel, ticksize)

    out = f"{out_file}.{file_format}"
    plt.savefig(out, format=file_format, bbox_inches="tight")
    plt.close()

def manplot_pure(file_lst, chro, pos, p, color_vec=("#e76254", "#1e466e"), scattersize=1.5, sig_color='red',
                 threshold=(1.0e-6, "black"), logp=True,
                 figsize=(13.5, 5), dpi=300, font='Arial', xlabel=(None, 8), ylabel=('-log10($p$)', 8),
                 ticksize=(8, 8), chro_show=None, out_file="test", file_format="tiff"):

    plt.rc('font', family=font)
    fig, ax = plt.subplots(figsize=(figsize[0] / 2.54, figsize[1] / 2.54), dpi=dpi)

    df = _read_concat_files(file_lst, chro, pos, p)
    axis = _prepare_genome_axis(df, chro, pos, chro_show=chro_show)

    y = -np.log10(df["pval"].to_numpy()) if logp else df["pval"].to_numpy()

    for i, c in enumerate(axis["chr_keys"]):
        mask = (df["chr_key"] == c).to_numpy()
        if not np.any(mask):
            continue
        col = color_vec[i % len(color_vec)]
        ax.scatter(axis["x"][mask], y[mask], s=scattersize, c=col,
                   rasterized=True, edgecolors="none")

    if threshold[1] is not None:
        _sig_hline(ax, threshold[0], threshold[1], logp, axis["x_max"])

    # Pure frame: hide all axes/spines
    ax.axis('off')

    out = f"{out_file}.{file_format}"
    plt.savefig(out, format=file_format, bbox_inches="tight")
    plt.close()

# =========================
# Main
# =========================
if __name__ == "__main__":
    # Paths / inputs
    analysis_dir = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/merge/"
    os.chdir(analysis_dir)

    # trait list file: one trait name per line/column; each trait resolves to <trait>.gz
    dfT = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/trait67.txt", sep="\t")
    trait_lst = dfT.iloc[:, 0].astype(str).tolist()
    print(f"Number of traits: {len(trait_lst)}")

    file_lst = [t for t in trait_lst]
    chrom, base, p = "chrom", "base", "p_gxe"

    # Bonferroni across traits
    pcut = 5e-8 / max(len(trait_lst), 1)

    out_prefix = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/Manhattan/BB"
    chro_show = [1,2,3,4,5,6,7,8,9,10,11,12,14,16,18,20]  # labels to show (others hidden)

    manplot(file_lst, chrom, base, p, threshold=(pcut, "black"),
            chro_show=chro_show, out_file=out_prefix, file_format="png")

    manplot_frame(file_lst, chrom, base, p, threshold=(pcut, "black"),
                  chro_show=chro_show, out_file=out_prefix + ".frame", file_format="pdf")

    manplot_pure(file_lst, chrom, base, p, threshold=(pcut, "black"),
                 chro_show=chro_show, out_file=out_prefix + ".pure", file_format="png")
