#!/usr/bin/env Rscript

# ================================
# Plot power curves with auto palettes
# - Merge two CSVs (extend + base), enforce num_envi_power == 30
# - Fig1: Power vs GxE h^2 (x as factor), colored by sample size
# - Fig2: Power vs Sample size (labels 50K, 100K, ...), colored by GxE h^2
# - Color palettes chosen automatically with sensible fallbacks
# ================================

# ---- Load libraries ----
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

# ---- I/O paths ----
in_csv  <- "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/take_home_extend.csv"
in_csv2 <- "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/take_home.csv"
out_png1 <- "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/take_home_extend1.png"
out_png2 <- "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/take_home_extend2.png"

# ---- Helper: auto qualitative palette (colorblind/print friendly) ----
get_qual_palette <- function(n) {
  # Try ggsci → RColorBrewer → viridisLite → rainbow
  pal <- NULL
  if (requireNamespace("ggsci", quietly = TRUE)) {
    # Use journal palettes with good class separability; repeat if n is large
    gens <- list(ggsci::pal_npg(), ggsci::pal_aaas(), ggsci::pal_lancet())
    for (g in gens) {
      base <- g(min(max(8, n), 10))
      pal <- rep(base, length.out = n)
      if (length(pal) == n) break
    }
  }
  if (is.null(pal) && requireNamespace("RColorBrewer", quietly = TRUE)) {
    for (nm in c("Dark2", "Set2", "Paired")) {
      maxn <- RColorBrewer::brewer.pal.info[nm, "maxcolors"]
      if (n <= maxn) {
        pal <- RColorBrewer::brewer.pal(n, nm); break
      }
    }
    if (is.null(pal)) {
      base <- RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info["Paired", "maxcolors"], "Paired")
      pal <- rep(base, length.out = n)
    }
  }
  if (is.null(pal)) {
    if (requireNamespace("viridisLite", quietly = TRUE)) {
      pal <- viridisLite::viridis(n, option = "D")  # colorblind & printer friendly
    } else {
      pal <- grDevices::rainbow(n)                  # last resort
    }
  }
  pal
}

# ---- Read & harmonize data ----
# Extended grid: keep num_envi_power == 30 and large sample sizes only
data_ext <- read.csv(in_csv, stringsAsFactors = FALSE)
if (!all(c("num_envi_power","sample_size") %in% names(data_ext))) {
  stop("Columns 'num_envi_power' and/or 'sample_size' not found in extend CSV.")
}
data_ext <- data_ext %>%
  filter(num_envi_power == 30, sample_size > 400000)

# Base grid: insert num_envi_power=30 at column 3, then align columns by name
data_base <- read.csv(in_csv2, stringsAsFactors = FALSE)
data_base$num_envi_power <- 30
# Move to 3rd column (keeps names intact; later we bind by names safely)
data_base <- data_base[, append(1:2, c(ncol(data_base), 3:(ncol(data_base)-1)))]

# Bind by column names to avoid order issues
data <- dplyr::bind_rows(data_base, data_ext)

write.csv(data, file = "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/take_home_combined.csv", row.names = FALSE)

# ---- Basic checks ----
required_cols <- c("sample_size", "power_snp_h2_gxe", "power_val")
missing_cols <- setdiff(required_cols, names(data))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

# ---- Tidy/annotate data ----
# Human-readable sample-size labels with fixed order for x-axis in Fig2 and color legend in Fig1
data <- data %>%
  mutate(
    sample_size_label = factor(
      dplyr::case_when(
        sample_size ==   50000 ~ "50K",
        sample_size ==  100000 ~ "100K",
        sample_size ==  200000 ~ "200K",
        sample_size ==  300000 ~ "300K",
        sample_size ==  400000 ~ "400K",
        sample_size ==  800000 ~ "800K",
        sample_size == 1200000 ~ "1200K",
        sample_size == 1600000 ~ "1600K",
        sample_size == 2000000 ~ "2000K",
        TRUE ~ as.character(sample_size)
      ),
      levels = c("50K", "100K", "200K", "300K", "400K", "800K", "1200K", "1600K", "2000K")
    )
  )

# Make sure GxE h^2 is numeric, then create an ordered factor for discrete x in Fig1 and color legend in Fig2
gxe_numeric <- suppressWarnings(as.numeric(data$power_snp_h2_gxe))
if (anyNA(gxe_numeric)) {
  # If input had percent strings like "2%", strip and convert to proportion
  gxe_numeric <- suppressWarnings(as.numeric(gsub("%", "", data$power_snp_h2_gxe))) / 100
}
if (anyNA(gxe_numeric)) stop("Cannot parse 'power_snp_h2_gxe' as numeric.")

gxe_levels <- sort(unique(gxe_numeric))
data <- data %>%
  mutate(
    gxe_h2   = gxe_numeric,
    gxe_h2_f = factor(gxe_numeric, levels = gxe_levels)
  )

# ---- Theme ----
base_theme <- theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black"),
    text             = element_text(size = 8, family = "Arial"),
    axis.title       = element_text(size = 8),
    axis.text.x      = element_text(size = 6),
    axis.text.y      = element_text(size = 7),
    legend.title     = element_blank(),
    legend.text      = element_text(size = 8),
    strip.text       = element_text(size = 8)
  )

# ---- Palettes (auto) ----
pal_sample <- get_qual_palette(nlevels(data$sample_size_label))  # for sample-size legend
pal_gxe    <- get_qual_palette(nlevels(data$gxe_h2_f))           # for GxE h^2 legend

# ---- Figure 1: Power vs GxE h^2 (x as factor), colored by sample size ----
fig1 <- ggplot(
  data,
  aes(x = gxe_h2_f, y = power_val, color = sample_size_label, group = sample_size_label)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  # Show two decimal places for h^2 values on x-axis
  scale_x_discrete(labels = function(x) format(as.numeric(x), nsmall = 2)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  scale_color_manual(values = pal_sample) +
  labs(x = "SNP GxE heritability (%)", y = "Power") +
  base_theme

ggsave(filename = out_png1, plot = fig1, width = 10/2.54, height = 7/2.54, dpi = 300)

# ---- Figure 2: Power vs sample size (labels as 50K, 100K, ...), colored by GxE h^2 ----
fig2 <- ggplot(
  data,
  aes(x = sample_size_label, y = power_val, color = gxe_h2_f)
) +
  geom_line(aes(group = gxe_h2_f), linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  scale_color_manual(values = pal_gxe) +
  labs(x = "Sample size", y = "Power", color = expression(paste("GxE ", h^2))) +
  base_theme

ggsave(filename = out_png2, plot = fig2, width = 10/2.54, height = 7/2.54, dpi = 300)

message("Saved: ", out_png1)
message("Saved: ", out_png2)


# ---- Combine two figures side-by-side ----
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
library(patchwork)

# Combine with patchwork: fig1 | fig2 places them side by side
combined_plot <- fig1 | fig2

# Save combined figure
out_png_combined <- "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/take_home_extend_combined.png"
ggsave(
  filename = out_png_combined,
  plot = combined_plot,
  width = 20/2.54,   # total width ~ 10cm per plot
  height = 7/2.54,
  dpi = 300
)

message("Saved combined plot: ", out_png_combined)
