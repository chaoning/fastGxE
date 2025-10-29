# Load required packages
library(ggplot2)
library(ggpubr)  # For adding Pearson correlation coefficient

# Load input data
df <- read.table("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/variance.txt", header = TRUE)

# Assign variables
x <- df$h2_1           # SNP heritability
y <- df$h2_gxe1        # GxE heritability

# Create the plot
p <- ggplot(df, aes(x = x, y = y)) +
  geom_point(color = "#1e466e", size = 1) +
  geom_smooth(method = "lm", color = "#e76254", se = TRUE) +
  coord_cartesian(ylim = c(0, 0.13)) +     # Better than ylim() to avoid data clipping
  stat_cor(method = "pearson", 
           label.x = min(x, na.rm = TRUE) + 0.01, 
           label.y = 0.12,                # Manual y pos works better with coord_cartesian
           size = 4, color = "black") +
  theme_bw(base_family = "Arial") +
  labs(x = "SNP heritability", y = "GxE heritability", title = "Physical traits") +
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

# Save the plot to PNG (6cm × 5cm, dpi=300)
ggsave("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/variance.trend_correlation_plot.PT.png",
       plot = p, width = 8 / 2.54, height = 7 / 2.54, dpi = 300)

df <- subset(df, Trait != 50)
# Assign variables
x <- df$h2_1           # SNP heritability
y <- df$h2_gxe1        # GxE heritability

# Create the plot
p <- ggplot(df, aes(x = x, y = y)) +
  geom_point(color = "#1e466e", size = 1) +
  geom_smooth(method = "lm", color = "#e76254", se = TRUE) +
  coord_cartesian(ylim = c(0, 0.13)) +     # Better than ylim() to avoid data clipping
  stat_cor(method = "pearson", 
           label.x = min(x, na.rm = TRUE) + 0.01, 
           label.y = 0.12,                # Manual y pos works better with coord_cartesian
           size = 4, color = "black") +
  theme_bw(base_family = "Arial") +
  labs(x = "SNP heritability", y = "GxE heritability", title = "Physical traits (without height)") +
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

# Save the plot to PNG (6cm × 5cm, dpi=300)
ggsave("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/variance.trend_correlation_plot.PT.rm50.png",
       plot = p, width = 8 / 2.54, height = 7 / 2.54, dpi = 300)

