library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# Read data
setwd("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/")
df <- read_csv("baselineLD.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000.csv")
colnames(df) <- c("h2_gxe", "fastGxE", "structLMM", "fastGWA-GE")

# Convert to long format
df_long <- df %>%
  pivot_longer(cols = -h2_gxe, names_to = "Method", values_to = "Power")

# Plot grouped bar chart
p <- ggplot(df_long, aes(x = factor(h2_gxe), y = Power, fill = factor(Method, levels = c("fastGxE", "structLMM", "fastGWA-GE")))) +
  geom_bar(stat = "identity", position = "dodge") +
            scale_fill_manual(values = c("#376795", "#06948E", "#e76254")) +
            theme_bw() +
           theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.key.size = unit(0.3, "cm"),   # 缩小色块
    legend.text = element_text(size = 6), # 缩小文字
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 8)
) +
            labs(title = NULL, x = "GxE interaction effect size", y = "Power")


# Save figure
ggsave("baselineLD.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000.png", p,  width = 7/2.54, height = 5/2.54, dpi = 300)
