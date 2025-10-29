# Load required libraries
library(ggplot2)
library(dplyr)

# Load power result CSV
data <- read.csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/vf.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000.csv")

# Helper function to reshape power results for each method and scenario
reshape_power_data <- function(data, power_cols, methods, title_label) {
  bind_rows(lapply(seq_along(power_cols), function(i) {
    data %>%
      select(`Variant function` = 1, Power = all_of(power_cols[i])) %>%
      mutate(Method = methods[i], tit = title_label)
  }))
}

# Define column indices and corresponding method names
power_cols_with_causal <- c(2, 4, 6)
power_cols_no_causal   <- c(3, 5, 7)
methods <- c("fastGxE", "StructLMM", "fastGWA-GE")

# Reshape and tag datasets
dataM1 <- reshape_power_data(data, power_cols_with_causal, methods, "Causal SNPs Included")
dataM2 <- reshape_power_data(data, power_cols_no_causal,   methods, "Causal SNPs Removed")

# Combine both datasets
dataM <- bind_rows(dataM1, dataM2)

# Ensure consistent method label capitalization
dataM$Method <- factor(dataM$Method, levels = c("fastGxE", "StructLMM", "fastGWA-GE"))

# Plot power by Variant function and method, split by causal SNPs inclusion
fig <- ggplot(dataM, aes(x = factor(`Variant function`, levels = unique(`Variant function`)), y = Power, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#1e466e", "#06948E", "#e76254")) +
  facet_wrap(~ tit, nrow = 1) +  # Split by condition: causal SNPs included/removed
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8),
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 10)
  ) +
  labs(title = NULL, x = "Variant function", y = "Power")

out_file="/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/vf.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000.png"
# Save the plot to PDF
ggsave(filename = out_file, plot = fig, width = 15 / 2.54, height = 5 / 2.54, dpi = 300)  # width doubled for two panels
