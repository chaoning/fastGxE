library(ggplot2)
# -----------------------------
# R ggplot2: grouped bars with ±2SE
# -----------------------------
library(tidyverse)
library(stringr)

# ---- 1) Load data ----
setwd("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_mom_sum/")
df <- read.csv("Varying_h2_nxe.csv")

# ---- 2) Select columns: *_h2_gxe and *_h2_gxe_noNxE ----
cols <- names(df) %>% str_subset("(?:_h2_gxe$|_h2_gxe_noNxE$)")
stopifnot(length(cols) > 0)

# ---- 3) Long format + parse NxE level and method ----
long <- df %>%
  select(all_of(cols)) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
  mutate(
    method    = if_else(str_detect(metric, "noNxE$"), "fastGxE-noNxE", "fastGxE"),
    NxE_level = str_match(metric, "h2_nxe_(\\d+)")[,2] %>% as.integer()
  )

# ---- 4) Summarize: mean, SE, ±2SE ----
summary_gxe <- long %>%
  group_by(method, NxE_level) %>%
  summarise(
    n    = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    se   = sd(value, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    lower = mean - 2*se,
    upper = mean + 2*se
  ) %>%
  arrange(method, NxE_level)

# ---- 5) Plot: grouped bars with ±2SE ----
p <- ggplot(summary_gxe, aes(x = factor(NxE_level), y = mean, fill = method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, color = "black") +
           scale_fill_manual(values = c("#376795", "#ffd06f")) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  scale_x_discrete(labels = function(x) paste0(x, "%")) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black", size = 0.6) +
  labs(
    x = expression(h[NxE]^2),
    y = expression("Estimated " * h[GE]^2),
    title = NULL,
    fill = "Method"
  ) +
  theme_bw(base_family = "ArialMT") +  # Set font family here
  theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line = element_line(color = "black"),
  text = element_text(size = 10, family = "ArialMT"),
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
) 

# If your values are proportions (0–1) and you want percentages, uncomment:
# library(scales)
# p <- p + scale_y_continuous(labels = percent_format(accuracy = 1))


# ---- 6) Save figure ----
outfile <- "Varying_h2_nxe_h2_gxe_and_noNxE_bar_2se_ggplot2.png"
ggsave(outfile, p, width = 6, height = 5, dpi = 300, units = "cm")
