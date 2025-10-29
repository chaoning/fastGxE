# Load required libraries
library(ggplot2)
library(dplyr)

# Load power result CSV
data <- read.csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/take_home_rndE.csv")
data$num_envi_power <- factor(paste0("# of active envs: ", data$num_envi_power),
                              levels = paste0("# of active envs: ", c(2, 30)))
data <- data %>%
  mutate(sample_size_label = factor(case_when(
    sample_size == 50000 ~ "50K",
    sample_size == 100000 ~ "100K",
    sample_size == 200000 ~ "200K",
    sample_size == 300000 ~ "300K",
    sample_size == 400000 ~ "400K"
  ), levels = c("50K", "100K", "200K", "300K", "400K")))
  
# Plot power by SNP GxE heritability, colored by sample size
fig <- ggplot(data, aes(
  x = factor(rndE, levels = c(1, 2, 5, 10, 20, 30, 40)),
  y = power_val,
  color = sample_size_label,
  group = sample_size_label
)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  scale_color_manual(values = c("#1e466e", "#06948E", "#e76254", "#f2b134", "#7e6bc4")) +
  facet_wrap(~ num_envi_power, nrow = 1, scales = "fixed") +  # Split by condition: causal SNPs included/removed
  theme_bw() +
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
  ) +
  labs(
    title = NULL,
    x = "# of included envs",
    y = "Power"
  )

# Save the plot to PNG
out_file <- "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/take_home_rndE.png"
ggsave(filename = out_file, plot = fig, width = 15 / 2.54, height = 7 / 2.54, dpi = 300)
