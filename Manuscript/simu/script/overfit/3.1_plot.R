# Load required libraries
library(ggplot2)
library(dplyr)

data1 <- read.csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_overfit_sum/summary_results_AIE.csv")
data2 <- read.csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_overfit_sum/summary_results_mmsusie_aggregated.csv")

colnames(data1) <- c("nE", "Total", "False", "fdr")
colnames(data2) <- c("nE", "Total", "False", "fdr")
data1$Method <- "AIE-stratified analysis"
data2$Method <- "mmSuSiE"
data <- rbind(data2, data1)
data$Method <- factor(data$Method, levels = c("mmSuSiE", "AIE-stratified analysis"))

# Plot power by Variant function and method, split by causal SNPs inclusion
fig <- ggplot(data, aes(x = factor(nE, levels = unique(nE)), y = fdr)) +
  geom_bar(stat = "identity", position = "dodge", fill = "#1e466e") +  # Set color inside the bar
  facet_wrap(~ Method, nrow = 1) +  # Split by condition: causal SNPs included/removed
  theme_bw() + ylim(0, 0.05) + 
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
  labs(title = NULL, x = "# of active envs", y = "FDR")

# Save the plot to PNG
out_file="/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_overfit_sum/overfit.png"
ggsave(filename = out_file, plot = fig, width = 15 / 2.54, height = 6 / 2.54, dpi = 300)  # width doubled for two panels
