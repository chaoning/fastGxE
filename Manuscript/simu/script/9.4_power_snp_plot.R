library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)  # For combining plots

# ================================
# 1. Read and prepare first dataset (df_long1)
# ================================
setwd("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/")

# Read SNP-level GxE results
df1 <- read_csv("snp_h2_gxe.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000.csv")


# ================================
# 2. Create plot for df1
# ================================
p1 <- ggplot(df1, aes(
    x = factor(snp_h2_gxe),  # X-axis: GxE heritability values
    y = power,
    fill = factor(method, levels = c("fastGxE", "fastGxE-noNxE", "StructLMM", "fastGWA-GE"))
  )) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#376795", "#ffd06f", "#06948E", "#e76254")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.key.size = unit(0.3, "cm"),      # Smaller legend keys
    legend.text = element_text(size = 6),   # Smaller legend text
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 8)
  ) +
  labs(x = "GxE interaction effects", y = "Power")

# ================================
# 3. Read and prepare second dataset (df_long2)
# ================================
# Read number-of-environments results
df2 <- read_csv("nE.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000.csv")

# ================================
# 4. Create plot for df_long2
# ================================
p2 <- ggplot(df2, aes(
    x = factor(num_envi),  # X-axis: Number of environments
    y = power,
    fill = factor(method, levels = c("fastGxE", "fastGxE-noNxE", "StructLMM", "fastGWA-GE"))
  )) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#376795", "#ffd06f", "#06948E", "#e76254")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 6),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 8)
  ) +
  labs(x = "# of active environmental factors", y = "Power")

# ================================
# 5. Combine the two plots side by side
# ================================
p_combined <- p1 + p2 +
  plot_layout(ncol = 2, guides = "collect") &   # collect = shared legend
  theme(legend.position = "bottom")                # place shared legend at top


# ================================
# 6. Save the combined figure
# ================================
ggsave(
  "power_two_plots.png",
  p_combined,
  width = 14/2.54,  # Adjust width for two side-by-side plots
  height = 6/2.54,
  dpi = 300
)
