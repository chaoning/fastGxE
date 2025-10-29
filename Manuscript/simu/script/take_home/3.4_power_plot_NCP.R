# Load required libraries
library(ggplot2)
library(dplyr)
library(patchwork)


# Load power result CSV
data <- read.csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/take_home_ncp.csv")
seq_vec <- c(0.005, 0.01, seq(0.02, 0.14, by = 0.02))

data <- data %>%
  mutate(
    sample_size_label = factor(case_when(
      sample_size == 50000 ~ "50K",
      sample_size == 100000 ~ "100K",
      sample_size == 200000 ~ "200K",
      sample_size == 300000 ~ "300K",
      sample_size == 400000 ~ "400K",
      sample_size == 800000 ~ "800K",
      sample_size == 1200000 ~ "1.2M",
      sample_size == 1600000 ~ "1.6M",
      sample_size == 2000000 ~ "2M"
    ), levels = c("50K", "100K", "200K", "300K", "400K", "800K", "1.2M", "1.6M", "2M")),
    
    # Convert x-axis to factor
    h2_factor = factor(power_snp_h2_gxe, levels = seq_vec)
  )

# Plot power by SNP GxE heritability, colored by sample size
fig1 <- ggplot(data, aes(
  x = h2_factor,
  y = power_val,
  color = sample_size_label,
  group = sample_size_label
)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme_bw(base_family = "ArialMT") +  # Set font family here
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    text = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    strip.text = element_text(size = 8)
  ) +
  
  labs(
    title = NULL,
    x = "SNP GxE heritability",
    y = "Power"
  )

fig2 <- ggplot(data, aes(
  x = sample_size_label,
  y = power_val,
  color = h2_factor,
  group = h2_factor
)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme_bw(base_family = "ArialMT") +  # Set font family here
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    text = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    strip.text = element_text(size = 8)
  ) +
  
  labs(
    title = NULL,
    x = "Sample size",
    y = "Power"
  )

fig_combined <- fig1 + fig2 + plot_layout(ncol = 2)
# Save the plot to PNG
out_file <- "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/take_home_ncp.png"
ggsave(filename = out_file, plot = fig_combined, width = 20 / 2.54, height = 7 / 2.54, dpi = 300)
