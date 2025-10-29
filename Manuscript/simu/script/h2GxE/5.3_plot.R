library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

setwd("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/h2GxE/")
# ==== Step 1: Read and process fastGxE results ====
data_fastGxE <- read.csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_mom_sum/indE.Varying_h2_nxe.csv")

# Select relevant columns (with and without NxE effects)
data_fastGxE <- data_fastGxE[, c("h2_nxe_0_h2_gxe", "h2_nxe_0_h2_gxe_noNxE",
                                 "h2_nxe_5_h2_gxe", "h2_nxe_5_h2_gxe_noNxE",
                                 "h2_nxe_15_h2_gxe", "h2_nxe_15_h2_gxe_noNxE",
                                 "h2_nxe_25_h2_gxe", "h2_nxe_25_h2_gxe_noNxE")]

# Convert data to long format
data_long <- data_fastGxE %>%
  pivot_longer(cols = everything(), names_to = "label", values_to = "value") %>%
  mutate(
    NxE = str_extract(label, "(?<=h2_nxe_)\\d+"),
    Method = ifelse(str_detect(label, "noNxE"), "fastGxE-noNxE", "fastGxE")
  )

# Summarize mean and SE
data_fastGxE_summary <- data_long %>%
  group_by(Method, NxE) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    se = sd(value, na.rm = TRUE) / sqrt(n()),
    lower = mean - 2 * se,
    upper = mean + 2 * se,
    .groups = "drop"
  )

# ==== Step 2: Read and process GENIE results ====
data_GENIE <- read.csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_GENIE_sum/indE.Varying_h2_nxe2.csv")

# Rename columns for consistency
data_GENIE <- data_GENIE[, c("h2_gxe_0", "h2_gxe_5", "h2_gxe_15", "h2_gxe_25")]
colnames(data_GENIE) <- c("h2_gxe_0", "h2_gxe_5", "h2_gxe_15", "h2_gxe_25")

# Convert to long format and extract NxE values
data_GENIE_long <- data_GENIE %>%
  pivot_longer(cols = everything(), names_to = "label", values_to = "value") %>%
  mutate(
    NxE = str_extract(label, "(?<=h2_gxe_)\\d+"),
    Method = "GENIE"
  )

# Summarize mean and SE
data_GENIE_summary <- data_GENIE_long %>%
  group_by(Method, NxE) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    se = sd(value, na.rm = TRUE) / sqrt(n()),
    lower = mean - 2 * se,
    upper = mean + 2 * se,
    .groups = "drop"
  )
# ==== Step 3: Combine both summaries ====
data_combined <- bind_rows(data_fastGxE_summary, data_GENIE_summary)

# Set order of NxE levels
data_combined$NxE <- factor(paste0(data_combined$NxE, "%"), 
                            levels = c("0%", "5%", "15%", "25%"))

# ==== Step 4: Create plot ====
plot <- ggplot(data_combined, aes(x = NxE, y = mean, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, color = "black") +
           scale_fill_manual(values = c("#1e466e", "#06948E",  "#e76254")) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.7), width = 0.2) +
                  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black", linewidth = 0.4) +  # <--- added line here
  theme_bw(base_family = "Arial") +  # Set font family here
  theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line = element_line(color = "black"),
  text = element_text(size = 10, family = "Arial"),
  axis.title = element_text(size = 10),
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8),
  legend.title = element_blank(),
  legend.text = element_text(size = 7),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.key.size = unit(0.3, "cm"),
  legend.spacing.x = unit(0.2, "cm"),
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(0, 0, 0, 0),
  strip.text = element_text(size = 10)
) +
  labs(
    y = expression(h[GxE]^2),
    x = expression(h[NxE]^2)
  ) +
  scale_y_continuous(limits = c(0, max(data_combined$upper) * 1.1))


# ==== Step 5: Save plot to PNG ====
png("indE.h2gxe_vs_nxe2.png", width = 8, height = 6, units = "cm", res = 300)
print(plot)
dev.off()
