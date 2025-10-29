import matplotlib
matplotlib.use("Agg")
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# === Config ===
csv_path = "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_take_home_whichE_sum/summary_results_adjusted.csv"
out_png  = "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_take_home_whichE_sum/summary_results_adjusted.png"
dpi = 300

# === Load data ===
df = pd.read_csv(csv_path)

# Ensure numeric dtype for safety
for col in ["power_snp_h2_gxe", "sample_size", "num_envi_power", "topE", "power"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")
df = df.dropna(subset=["power_snp_h2_gxe", "sample_size", "num_envi_power", "topE", "power"])

# Restrict to selected h2 and numE (as in your code)
df = df[df["power_snp_h2_gxe"].isin([0.02, 0.08, 0.14])]
df = df[df["num_envi_power"].isin([2, 5, 10, 30])]

# Optional: stable legend order for sample sizes
desired_sample_order = [50_000, 100_000, 200_000, 300_000, 400_000]
present = [s for s in desired_sample_order if s in set(df["sample_size"])]
df["sample_size"] = pd.Categorical(df["sample_size"], categories=present, ordered=True)

# Sort for clean line connections
df = df.sort_values(["power_snp_h2_gxe", "num_envi_power", "sample_size", "topE"])

# === Equal-spaced categorical X (TopE as factor) ===
x_levels = [1, 2, 5, 10, 20, 30, 40]          # visible tick labels
x_map = {v: i for i, v in enumerate(x_levels)} # value -> position index

# === Facet grid ===
unique_h2 = sorted(df["power_snp_h2_gxe"].unique())
unique_num_envi = sorted(df["num_envi_power"].unique())
n_rows, n_cols = len(unique_h2), len(unique_num_envi)

# Size in inches (you used cm; convert cm -> inch by /2.54)
fig_w = 6 * n_cols / 2.54
fig_h = 5 * n_rows / 2.54
fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h), sharex=True, sharey=True)

def get_ax(i, j):
    if n_rows == 1 and n_cols == 1: return axes
    if n_rows == 1: return axes[j]
    if n_cols == 1: return axes[i]
    return axes[i, j]

legend_handles = None

for i, h2_val in enumerate(unique_h2):
    for j, num_envi_val in enumerate(unique_num_envi):
        ax = get_ax(i, j)

        # Subset for this facet
        facet = df[(df["power_snp_h2_gxe"] == h2_val) & (df["num_envi_power"] == num_envi_val)].copy()

        if facet.empty:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", fontsize=9, transform=ax.transAxes)
            ax.set_title(rf"$h^2_{{G\times E}}$={h2_val:g}, $N_{{E,\mathrm{{active}}}}$={int(num_envi_val)}", fontsize=9)
            continue

        handles, labels = [], []

        # Draw one line per sample size
        for sample_size, sub in facet.groupby("sample_size", sort=True):
            # Keep only desired x levels and map to equal-spaced positions
            sub = sub[sub["topE"].isin(x_levels)].copy()
            sub["x_pos"] = sub["topE"].map(x_map)
            sub = sub.sort_values("x_pos")

            h, = ax.plot(
                sub["x_pos"].values,      # equal-spaced x
                sub["power"].values,
                marker="o",
                linewidth=1.5,
                markersize=3,
                label=f"{int(sample_size):,}"
            )
            handles.append(h)
            labels.append(f"{int(sample_size):,}")

        ax.set_title(
            rf"$h^2_{{G\times E}}$={h2_val:g}, $N_{{E,\mathrm{{active}}}}$={int(num_envi_val)}",
            fontsize=9
        )

        legend_handles = (handles, labels)

# Set x ticks/labels and limits for ALL axes (shared x)
for ax in np.array(axes).ravel():
    ax.set_xticks(range(len(x_levels)))
    ax.set_xticklabels([str(v) for v in x_levels], fontsize=8)
    ax.set_xlim(-0.5, len(x_levels) - 0.5)

# === Global legend (top center, single row) ===
if legend_handles is not None:
    handles, labels = legend_handles
    fig.legend(handles, labels, title="Sample size",
               loc="upper center", bbox_to_anchor=(0.5, 0.98),
               ncol=len(labels), frameon=False, fontsize=8)

# === Global axis labels ===
fig.supxlabel("# of included envs", fontsize=10, y=0.01)  # lift/lower by tweaking y
fig.supylabel("Power", fontsize=10, x=0.01)

# Layout: leave top space for legend and left for y-label
plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.90])

plt.savefig(out_png, dpi=dpi)
print(f"Saved figure to: {os.path.abspath(out_png)}")
