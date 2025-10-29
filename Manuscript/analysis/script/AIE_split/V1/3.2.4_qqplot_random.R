# Load the required library
library(ggplot2)

setwd("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/AIE_split/")
# Read the data from CSV file
data <- read.table("3.2.3_test_main_py_random.txt", header = FALSE)

# Rename the column to 'pval' for convenience
colnames(data) <- c("pval")

# Number of p-values
n <- nrow(data)

# If the number of p-values exceeds 10,000, sample 10,000; otherwise use all
if (n > 10000) {
  set.seed(123)  # for reproducibility
  data <- data[sample(n, 10000), , drop = FALSE]
}
n <- nrow(data)
# Sort observed p-values and transform to -log10 scale
observed <- -log10(sort(data$pval))

# Expected uniform distribution quantiles, also in -log10 scale
expected <- -log10(ppoints(n))

# Set alpha for 95% confidence interval
alpha <- 0.05

# Calculate confidence interval bands using beta distribution
lower_ci <- -log10(qbeta(alpha / 2, 1:n, n:1))
upper_ci <- -log10(qbeta(1 - alpha / 2, 1:n, n:1))

# Create data frame for plotting
qq_df <- data.frame(
  expected = expected,
  observed = observed,
  lower_ci = lower_ci,
  upper_ci = upper_ci
)

# Save plot to a high-resolution PNG file
png("3.2.3_test_main_py_random.qq_plot_with_CI.png", width = 8, height = 7, res = 300, units = "cm")

# Generate the QQ plot
ggplot(qq_df, aes(x = expected, y = observed)) +
  # Add the confidence interval ribbon
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "gray80", alpha = 0.5) +
  # Add y = x reference line
  geom_abline(slope = 1, intercept = 0, color = "black") +
  # Add observed data points
  geom_point(size = 1.5, alpha = 0.6, color = "#e76254") +
  scale_x_continuous(breaks = seq(0, ceiling(max(qq_df$expected)), by = 1)) +
  scale_y_continuous(breaks = seq(0, ceiling(max(qq_df$observed)), by = 1)) +
  # Add plot title and axis labels
  labs(
  x = expression(Expected ~ -log[10](italic(p))),
  y = expression(Observed ~ -log[10](italic(p)))
) +
    theme_bw(base_family = "ArialMT") +  # Set font family here
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

# Close the PNG device to save the file
dev.off()
