# Load required packages
library(ggplot2)
library(ggpubr)   # for stat_cor and ggarrange

# Set working directory (adjust if needed)
setwd("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mom/")

# ---------- Helper to create a single panel ----------
make_panel <- function(df, title, ylim_max, label_y) {
  # Compute a safe left position for the correlation label
  label_x <- min(df$h2G, na.rm = TRUE) + 0.01
  
  ggplot(df, aes(x = h2G, y = h2GxE)) +
    geom_point(color = "#1e466e", size = 1) +
    geom_smooth(method = "lm", color = "#e76254", se = TRUE) +
    coord_cartesian(ylim = c(0, ylim_max)) +  # use coord_cartesian to avoid dropping points
    stat_cor(
      method = "pearson",
      label.x = label_x,
      label.y = label_y,
      size = 3, color = "black"
    ) +
    theme_bw(base_family = "Arial") +
    labs(
      x = "SNP heritability",
      y = "GxE heritability",
      title = title
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      text = element_text(size = 8),
      axis.title = element_text(size = 9),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      strip.text = element_text(size = 8)
    )
}

# ---------- Load data ----------
df_all <- read.csv("PT.variance_for_plot.csv", header = TRUE)

# Subset for "without height" (FieldID != 50)
df_rm50 <- subset(df_all, FieldID != 50)

# If you also need the blood biomarkers (BB), load here:
# (Replace filename if needed; your original code reused PT file for BB block)
df_bb <- read.csv("BB.variance_for_plot.csv", header = TRUE)

# ---------- Build panels ----------
p_pt     <- make_panel(df_all,  "Physical traits",                         ylim_max = 0.13, label_y = 0.12)
p_pt_rm  <- make_panel(df_rm50, "Physical traits (without height)",        ylim_max = 0.13, label_y = 0.12)
p_bb     <- make_panel(df_bb,   "Blood biomarkers",                        ylim_max = 0.10, label_y = 0.095)

# ---------- Combine into one figure (3 columns) ----------
p_combined <- ggarrange(p_pt, p_pt_rm, p_bb, ncol = 3, nrow = 1, align = "hv")

# ---------- Save combined figure (≈18 cm × 7 cm at 300 dpi) ----------
ggsave(
  filename = "variance.trend_correlation_plot.combined.png",
  plot = p_combined,
  width = 18 / 2.54,   # 18 cm
  height = 6 / 2.54,   # 7 cm
  dpi = 300
)
