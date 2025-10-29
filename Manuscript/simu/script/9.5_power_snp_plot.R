library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(grid)   # for unit()

# ================================
# 1. Read and prepare datasets
# ================================
setwd("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/")

# --- Dataset 1: SNP-level GxE results ---
df1 <- read_csv("snp_h2_gxe.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000.csv") %>%
  filter(method != "fastGxE-noNxE")
df1$snp_h2_gxe <- factor(df1$snp_h2_gxe, levels = sort(unique(df1$snp_h2_gxe)))

# --- Dataset 2: Number of environments ---
df2 <- read_csv("nE.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000.csv") %>%
  filter(method != "fastGxE-noNxE")
df2$num_envi <- factor(df2$num_envi, levels = sort(unique(df2$num_envi)))

# --- Dataset 3: Masked environments ---
df3 <- read_csv("rmE.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000.csv") %>%
  filter(method != "fastGxE-noNxE")
df3$rmE <- factor(df3$rmE, levels = sort(unique(df3$rmE)))

# ================================
# 2. Common theme and settings
# ================================
common_theme <- theme_bw() +
  theme(
    # Remove background grid and border
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    
    # Add black axis lines
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    
    # Legend style
    legend.title = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 6),
    
    # Global text
    text = element_text(size = 10),
    
    # Axis titles (custom sizes)
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 10),
    
    # Axis tick labels
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 8)
  )


# Fill colors and method order
fill_vals <- c("fastGxE" = "#376795", "StructLMM" = "#06948E", "fastGWA-GE" = "#e76254")
method_levels <- c("fastGxE", "StructLMM", "fastGWA-GE")

# ================================
# 3. Individual plots
# ================================
p1 <- ggplot(df1, aes(x = factor(snp_h2_gxe),
                      y = power,
                      fill = factor(method, levels = method_levels))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = fill_vals) +
  labs(x = "GxE interaction effects", y = "Power") +
  common_theme

p2 <- ggplot(df2, aes(x = factor(num_envi),
                      y = power,
                      fill = factor(method, levels = method_levels))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = fill_vals) +
  labs(x = "# of active environmental factors", y = "Power") +
  common_theme

p3 <- ggplot(df3, aes(x = factor(rmE),
                      y = power,
                      fill = factor(method, levels = method_levels))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = fill_vals) +
  labs(x = "# of masked environmental factors", y = "Power") +
  common_theme

# ================================
# 4. Combine plots with one shared legend
# ================================
# (p1 | p2) on the first row, p3 on the second row
p1 <- p1 + geom_point(aes(color = factor(method, levels = method_levels)), y = 0, alpha = 0, show.legend = TRUE) + scale_color_manual(values = fill_vals)
p2 <- p2 + geom_point(aes(color = factor(method, levels = method_levels)), y = 0, alpha = 0, show.legend = TRUE) + scale_color_manual(values = fill_vals)
p3 <- p3 + geom_point(aes(color = factor(method, levels = method_levels)), y = 0, alpha = 0, show.legend = TRUE) + scale_color_manual(values = fill_vals)

p_combined <- (p1 | p2 | p3) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",   # 横向排列
    legend.direction = "horizontal"
  ) &
  guides(
    fill  = guide_legend(
      title = NULL, order = 1,
      override.aes = list(shape = 22, size = 4)   # 方块
    ),
    color = guide_legend(
      title = NULL, order = 2,
      override.aes = list(alpha = 1, shape = 19, size = 3) # 圆点
    )
  )



# ================================
# 5. Save figure
# ================================
ggsave(
  "power_three_plots.png",
  p_combined,
  width = 18/2.54, 
  height = 6/2.54, 
  dpi = 300
)

ggsave(
  "power_three_plots.pdf",
  p_combined,
  width = 18/2.54, 
  height = 6/2.54,
  useDingbats = FALSE
)
