# Load required libraries
library(ggplot2)
library(dplyr)

# Load power result CSV
data <- read.csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/nonE.csv")


# Plot power by SNP GxE heritability, colored by nonE
fig <- ggplot(data, aes(
  x = nonE,
  y = power_val, color = "#1e466e"
)) +
  geom_line(linewidth = 0.8, color = "#1e466e") +
  geom_point(size = 1.5, color = "#1e466e") +
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
    x = "# of non-active environmental factors",
    y = "Power"
  )

# Save the plot to PNG
out_file <- "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/nonE.png"
ggsave(filename = out_file, plot = fig, width = 8 / 2.54, height = 7 / 2.54, dpi = 300)
