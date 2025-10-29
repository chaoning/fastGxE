# Load the required library
library(ggplot2)

setwd("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/AIE_split/")
# Read the data from CSV file
data <- read.table("3.2.3_test_main_py_random.txt", header = FALSE)
# Number of p-values
n <- nrow(data)

# If the number of p-values exceeds 10,000, sample 10,000; otherwise use all
if (n > 10000) {
  set.seed(123)  # for reproducibility
  data <- data[sample(n, 10000), , drop = FALSE]
}

pval <- as.vector(as.matrix(data[, 1:(ncol(data)-1), drop = FALSE]))
main_pval <- data.frame(pval = pval)
pval <- as.vector(as.matrix(data[, ncol(data)-1]))
gxe_pval <- data.frame(pval = pval)



colnames(main_pval) <- "pval"
colnames(gxe_pval)  <- "pval"

# ==== Add type labels ====
main_pval$type <- "Genetic main test"
gxe_pval$type  <- "Stratified test"

# ==== Combine ====
data <- rbind(main_pval, gxe_pval)

# ==== Remove invalid p-values ====
data <- data[!is.na(data$pval) & data$pval > 0 & data$pval <= 1, , drop = FALSE]

# ==== Prepare QQ plot data ====
qq_df <- do.call(rbind, lapply(split(data, data$type), function(df) {
  n <- nrow(df)
  observed <- -log10(sort(df$pval))
  expected <- -log10(ppoints(n))
  alpha <- 0.05
  lower_ci <- -log10(qbeta(alpha / 2, 1:n, n:1))
  upper_ci <- -log10(qbeta(1 - alpha / 2, 1:n, n:1))
  data.frame(
    expected = expected,
    observed = observed,
    lower_ci = lower_ci,
    upper_ci = upper_ci,
    type = unique(df$type)
  )
}))

# ==== Plot and save ====
png("3.2.3_test_main_gxe.qq_plot_with_CI.png", width = 12, height = 6, res = 300, units = "cm")

ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "gray80", alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_point(size = 1.5, alpha = 0.6, color = "#e76254") +
  scale_x_continuous(breaks = seq(0, ceiling(max(qq_df$expected)), by = 1)) +
  scale_y_continuous(breaks = seq(0, ceiling(max(qq_df$observed)), by = 1)) +
  facet_wrap(~ type, nrow = 1) +
  labs(
    x = expression(Expected ~ -log[10](italic(p))),
    y = expression(Observed ~ -log[10](italic(p)))
  ) +
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
    strip.text = element_text(size = 10),
    legend.position = "none"
  )

dev.off()