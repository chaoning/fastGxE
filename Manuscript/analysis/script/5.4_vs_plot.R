# Load required libraries
library(ggplot2)
library(dplyr)

# Load data
data <- read.table("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/variance.txt", header = TRUE)

# Prepare long-format data
data1 <- data.frame(h2_gxe = data$h2_gxe1, Model = "NxE included")
data2 <- data.frame(h2_gxe = data$h2_gxe2, Model = "NxE removed")
data_long <- rbind(data1, data2)
data_long$Trait <- rep(1:nrow(data), times = 2)  # Create pairing index

# Perform paired t-test
ttest_result <- t.test(data$h2_gxe1, data$h2_gxe2, paired = TRUE)
pval_text <- bquote("Paired t-test " ~ italic(p) == .(signif(ttest_result$p.value, 3)))

# Boxplot
fig <- ggplot(data_long, aes(x = factor(Model, levels = unique(Model)), y = h2_gxe, fill = Model)) +
  geom_boxplot(width = 0.5, outlier.size = 0.8, outlier.shape = 21, color = "black") +
  scale_fill_manual(values = c("#1e466e", "#e76254")) +
  theme_bw(base_family = "Arial") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",  # Remove legend
    strip.text = element_text(size = 10)
  ) +
  labs(title = "Blood biomarkers", x = NULL, y = expression(h[GxE]^2)) +
  annotate("text", x = 1.5, y = max(data_long$h2_gxe) * 1.05, label = pval_text, size = 3)

# Save plot
out_file <- "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/variance.GxEh2.BB.png"
ggsave(filename = out_file, plot = fig, width = 7 / 2.54, height = 5 / 2.54, dpi = 300)
