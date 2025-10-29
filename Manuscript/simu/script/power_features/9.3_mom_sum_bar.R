library(ggplot2)
library(dplyr)

# Set working directory
setwd("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_baselineLD_mom_sum/")

# Read CSV file
data <- read.csv("mom_sum.csv")

# Extract the second column (GxE heritability values)
gxe_vals <- data[[2]]

# Calculate sample size, mean, and standard error (SE)
n <- sum(!is.na(gxe_vals))
mean_val <- mean(gxe_vals, na.rm = TRUE)
se_val <- sd(gxe_vals, na.rm = TRUE) / sqrt(n)

# Create a summary data frame for plotting
summary_df <- tibble(
  h2   = "GxE heritability",         # x-axis label
  mean = mean_val,                   # bar height
  se   = se_val,                     # standard error
  lower = mean_val - 2 * se_val,     # lower bound (mean - 2*SE)
  upper = mean_val + 2 * se_val      # upper bound (mean + 2*SE)
)

# Build the bar plot
p <- ggplot(summary_df, aes(x = h2, y = mean)) +
  # Main bar
  geom_bar(stat = "identity", width = 0.6, color = "black", fill = "#376795") +
  # Error bar (± 2*SE)
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  # Reference lines at y = 0.05 and y = 0.3
  geom_hline(yintercept = c(0.05, 0.3), linetype = "dashed", 
             color = "black", linewidth = 0.4) +
  # Apply a clean theme
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "none"
  ) +
  labs(
    y = expression(h^2),   # y-axis label: h²
    x = ""                 # no x-axis label
  ) +
  # Set y-axis limits slightly above the maximum error bar
  scale_y_continuous(limits = c(0, 0.15))

# Save the plot to PNG file
png("h2gxe_only.png", width = 6, height = 5, units = "cm", res = 300)
print(p)
dev.off()

# Print the summary statistics
print(summary_df)
