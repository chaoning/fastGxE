library(ggplot2)
library(tidyr)
library(dplyr)
setwd("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mom/")

# Read the input data table
df <- read.csv(
  "PT.variance_for_plot.csv",
  header = TRUE)

# Reshape data from wide format to long format for ggplot
long <- df |>
  pivot_longer(cols = c(h2G, h2GxE, h2NxE),
               names_to = "component", values_to = "value") |>
  mutate(component = factor(component,
                            levels = c("h2G", "h2GxE", "h2NxE"),
                            labels = c("SNP Heritability", 
                                       "GxE Heritability", 
                                       "NxE effect level"))) |>
  group_by(ShortName) |>
  mutate(total = sum(value, na.rm = TRUE)) |>  # total per trait (for ordering)
  ungroup()

# Define colors for each component
cols <- c("SNP Heritability" = "#376795",
          "GxE Heritability" = "#06948E",
          "NxE effect level" = "#e76254")

# --- Horizontal stacked bar plot (recommended) ---
p_h <- ggplot(long, aes(x = reorder(ShortName, total), y = value, fill = component)) +
  geom_col(width = 0.8) +                                     # stacked bars
  coord_cartesian(ylim = c(0, 1), expand = FALSE) +           # keep y in [0,1]
  coord_flip() +                                              # flip axes (horizontal bars)
  scale_fill_manual(values = cols, breaks = names(cols)) +    # custom colors
  labs(x = "Trait", y = "Value", fill = NULL) +
  theme_bw(base_family = "Arial") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 9),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )

# Save the horizontal plot (PDF and PNG)
ggsave("PT_variance.pdf",
       p_h, width = 18, height = 14, units = "cm", dpi = 300, device = cairo_pdf)
ggsave("PT_variance.png",
       p_h, width = 18, height = 14, units = "cm", dpi = 300)



# Read the input data table
df <- read.csv(
  "BB.variance_for_plot.csv",
  header = TRUE)

# Reshape data from wide format to long format for ggplot
long <- df |>
  pivot_longer(cols = c(h2G, h2GxE, h2NxE),
               names_to = "component", values_to = "value") |>
  mutate(component = factor(component,
                            levels = c("h2G", "h2GxE", "h2NxE"),
                            labels = c("SNP Heritability", 
                                       "GxE Heritability", 
                                       "NxE effect level"))) |>
  group_by(ShortName) |>
  mutate(total = sum(value, na.rm = TRUE)) |>  # total per trait (for ordering)
  ungroup()

# Define colors for each component
cols <- c("SNP Heritability" = "#376795",
          "GxE Heritability" = "#06948E",
          "NxE effect level" = "#e76254")

# --- Horizontal stacked bar plot (recommended) ---
p_h <- ggplot(long, aes(x = reorder(ShortName, total), y = value, fill = component)) +
  geom_col(width = 0.8) +                                     # stacked bars
  coord_cartesian(ylim = c(0, 1), expand = FALSE) +           # keep y in [0,1]
  coord_flip() +                                              # flip axes (horizontal bars)
  scale_fill_manual(values = cols, breaks = names(cols)) +    # custom colors
  labs(x = "Trait", y = "Value", fill = NULL) +
  theme_bw(base_family = "Arial") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 9),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )

# Save the horizontal plot (PDF and PNG)
ggsave("BB_variance.pdf",
       p_h, width = 18, height = 20, units = "cm", dpi = 300, device = cairo_pdf)
ggsave("BB_variance.png",
       p_h, width = 18, height = 20, units = "cm", dpi = 300)
