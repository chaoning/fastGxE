#!/usr/bin/env Rscript

# ---- Load libraries ----
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

# ---- I/O paths ----
in_csv  <- "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/take_home.csv"
out_png1 <- "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/take_home.png"
out_png2 <- "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/take_home2.png"

# ---- Read data ----
data <- read.csv(in_csv, stringsAsFactors = FALSE)

# ---- Tidy/annotate data ----
# Map sample_size to a readable label with fixed ordering
data <- data %>%
  mutate(
    sample_size_label = factor(
      case_when(
        sample_size ==  50000 ~ "50K",
        sample_size == 100000 ~ "100K",
        sample_size == 200000 ~ "200K",
        sample_size == 300000 ~ "300K",
        sample_size == 400000 ~ "400K",
        TRUE ~ as.character(sample_size)
      ),
      levels = c("50K", "100K", "200K", "300K", "400K")
    )
  )

# Ensure GxE h^2 has both numeric (for x-axis) and factor (for legend) versions
# Use sorted unique values to keep legend in numeric order
gxe_levels <- sort(unique(as.numeric(data$power_snp_h2_gxe)))
data <- data %>%
  mutate(
    gxe_h2   = as.numeric(power_snp_h2_gxe),                    # numeric: for continuous x-axis
    gxe_h2_f = factor(gxe_h2, levels = gxe_levels)              # factor: for color legend
  )

# ---- Common theme and palettes ----
base_theme <- theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black"),
    text             = element_text(size = 10),
    axis.title       = element_text(size = 10),
    axis.text.x      = element_text(size = 8),
    axis.text.y      = element_text(size = 8),
    legend.title     = element_blank(),
    legend.text      = element_text(size = 9),
    strip.text       = element_text(size = 10)
  )

# 5-color palette for sample-size lines (first plot)
palette5 <- c("#1e466e", "#06948E", "#e76254", "#f2b134", "#7e6bc4")
# Keep only as many colors as you have levels
palette5 <- palette5[seq_len(nlevels(data$sample_size_label))]

# 7-color palette for GxE h^2 (second plot)
palette7 <- c("#1e466e", "#06948E", "#e76254", "#f2b134", "#7e6bc4", "#8b4513", "#ff69b4")
# Keep only as many colors as you have levels (safe if levels < 7)
palette7 <- palette7[seq_len(nlevels(data$gxe_h2_f))]

# Nicely spaced x breaks for gxe_h2 (use actual unique values if they are regular)
x_breaks_1 <- gxe_levels
# If you prefer fixed breaks like 0.02 to 0.14 by 0.02, uncomment:
# x_breaks_1 <- seq(0.02, 0.14, by = 0.02)

# ---- Figure 1: Power vs GxE h^2, colored by sample size ----
fig1 <- ggplot(
  data,
  aes(x = gxe_h2, y = power_val, color = sample_size_label)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  scale_x_continuous(breaks = x_breaks_1, labels = number_format(accuracy = 0.01)) +
  scale_color_manual(values = palette5) +
  labs(
    title = NULL,
    x = "SNP GxE heritability",
    y = "Power"
  ) +
  base_theme

ggsave(filename = out_png1, plot = fig1, width = 10/2.54, height = 7/2.54, dpi = 300)

# ---- Figure 2: Power vs sample size (labels as 50K, 100K, ...), colored by GxE h^2 ----
fig2 <- ggplot(
  data,
  aes(x = sample_size_label, y = power_val, color = gxe_h2_f)  # use label factor for x
) +
  geom_line(aes(group = gxe_h2_f), linewidth = 0.8) +  # group lines by GxE h^2
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  scale_color_manual(values = palette7) +
  labs(
    title = NULL,
    x = "Sample size",
    y = "Power",
    color = expression(paste("GxE ", h^2))
  ) +
  base_theme

ggsave(filename = out_png2, plot = fig2, width = 10/2.54, height = 7/2.54, dpi = 300)

# ---- Done ----
message("Saved: ", out_png1)
message("Saved: ", out_png2)
