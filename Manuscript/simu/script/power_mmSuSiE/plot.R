library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)  # For combining plots

# ================================
# 1. Read and prepare first dataset (df_long1)
# ================================
setwd("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_overfit_sum/")
df1 <- read_csv("mmsusie_PIP_summary.csv")


# ================================
# 2. Create plot for df1
# ================================
p1 <- ggplot(df1, aes(x = mean_PIP, y = in_lst_rate)) +
  geom_point(size = 3, color = "#e76254") +                             
  geom_errorbar(aes(ymin = in_lst_rate - 2 * SE, 
                    ymax = in_lst_rate + 2 * SE), 
                width = 0.01, color = "#e76254") +                    
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  
  labs(
    x = "Mean PIP",
    y = "Proportion of signals"
  ) +
  theme_minimal(base_size = 14) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
  theme_bw(base_family = "Arial") +  # Set font family here
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 10)
  ) 

# ================================
# 3. Read and prepare second dataset (df_long2)
# ================================
# Read number-of-environments results
df2 <- read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_mmSuSiE/power_summary.csv")
df2 <- df2[1:14,]
# ================================
# 4. Create plot for df_long2
# ================================
p2 <- ggplot(df2, aes(
    x = factor(power_snp_h2_gxe),  # X-axis: Number of environments
    y = Power,
    fill = factor(Method, levels = c("mmSuSiE", "fastGWA-GE"))
  )) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#376795", "#e76254")) +
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
  labs(x = "GxE interaction effects", y = "Power")

df3 <- read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_mmSuSiE/power_summary.csv")
df3 <- df3[15:28,]
# ================================
# 4. Create plot for df_long2
# ================================
p3 <- ggplot(df3, aes(
    x = factor(num_envi),  # X-axis: Number of environments
    y = Power,
    fill = factor(Method, levels = c("mmSuSiE", "fastGWA-GE"))
  )) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#376795", "#e76254")) +
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
    axis.title = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 8)
  ) +
  labs(x = "# of active environmental factors", y = "Power")
# ================================
# 5. Combine the two plots side by side
# ================================
p_combined <- p1 + p2 + p3
  plot_layout(ncol = 2, guides = "collect") &   # collect = shared legend
  theme(legend.position = "bottom")                # place shared legend at top


# ================================
# 6. Save the combined figure
# ================================
ggsave(
  "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_mmSuSiE/power_identifyE.png",
  p_combined,
  width = 19,
  height = 6,
  units = "cm",
  dpi = 300
)
