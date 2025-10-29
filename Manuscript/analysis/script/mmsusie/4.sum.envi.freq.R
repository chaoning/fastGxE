# ===========================
# Packages
# ===========================
library(ggplot2)
library(dplyr)
library(readr)   # fast and safe CSV reader
library(scales)  # for percent_format()
library(grid)    # for unit()

# ===========================
# Global options
# ===========================
setwd("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/")

# Consistent colors for sign of PIP
COL_NEG <- "#e76254"   # negative category
COL_POS <- "#1e466e"   # positive category

# ===========================
# Helper: read, validate, aggregate
# ===========================
read_and_count <- function(csv_path) {
  # Read CSV; keep only needed columns and drop NA values early
  # Expect columns: 'envi' (environment label) and 'pip' (signed value)
  df <- read_csv(csv_path, show_col_types = FALSE) %>%
    select(envi, pip) %>%
    tidyr::drop_na(envi, pip)

  # Basic validation: check required columns exist
  required <- c("envi", "pip")
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing required columns in %s: %s",
      csv_path, paste(missing_cols, collapse = ", ")
    ))
  }

  # Tag category by sign of PIP
  df <- df %>%
    mutate(
      pip_category = if_else(pip > 0, "+", "-")
    )

  # Count by (environment, sign)
  count_data <- df %>%
    count(envi, pip_category, name = "count") %>%
    group_by(envi) %>%
    mutate(
      total = sum(count),
      percent = 100 * count / total
    ) %>%
    ungroup() %>%
    # Reorder x-axis by total count (descending)
    mutate(envi = reorder(envi, -total))

  return(count_data)
}

# ===========================
# Helper: plotting
# ===========================
plot_envi_pip_frequency <- function(count_data,
                                    outfile_base,          # base filename without extension
                                    width_cm = 8, height_cm = 6,
                                    use_percent = FALSE,
                                    show_values = TRUE,
                                    legend_pos = c(0.95, 0.95),
                                    png_dpi = 300, png_bg = "white") {
  # Select y-axis mapping: counts or percentages
  y_var <- if (use_percent) quo(percent) else quo(count)
  y_lab <- if (use_percent) "Percentage of GxE associations" else "# of GxE associations"

  # Define text labels shown inside stacked bars
  label_aes <- if (use_percent) aes(label = scales::percent(percent / 100, accuracy = 0.1))
               else aes(label = count)

  # Main ggplot object
  p <- ggplot(count_data, aes(x = envi, y = !!y_var, fill = pip_category)) +
    geom_bar(stat = "identity", width = 0.8) +
    # Optional numeric labels inside stacked bars
    { if (show_values) geom_text(label_aes, position = position_stack(vjust = 0.5), size = 2) } +
    scale_fill_manual(values = c("-" = COL_NEG, "+" = COL_POS), name = NULL) +
    # Axis labels (no title)
    labs(title = NULL, x = "Environmental factor", y = y_lab) +
    # Clean theme
    theme_bw(base_size = 8) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),                # remove full panel border
      axis.line.x = element_line(colour = "black"),  # keep bottom axis line
      axis.line.y = element_line(colour = "black"),  # keep left axis line
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1, colour = "black"),
      axis.text.y = element_text(size = 8, colour = "black"),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm"),
      legend.position = legend_pos,
      legend.justification = c("right", "top")
    )

  # --- Save plots in both PDF and PNG formats ---
  ggsave(filename = paste0(outfile_base, ".pdf"),
         plot = p, device = cairo_pdf,
         width = width_cm, height = height_cm, units = "cm")

  ggsave(filename = paste0(outfile_base, ".png"),
         plot = p,
         width = width_cm, height = height_cm, units = "cm",
         dpi = png_dpi, bg = png_bg)
}

# ===========================
# Run for both datasets
# ===========================
# 1) PT.PIP.sum.loci.envi.csv
count_PT <- read_and_count("PT.PIP.sum.loci.envi.csv")
plot_envi_pip_frequency(
  count_data   = count_PT,
  outfile_base = "PT_envi_pip_frequency",
  width_cm     = 8,
  height_cm    = 6,
  use_percent  = FALSE,   # set TRUE if you want percentages instead
  show_values  = TRUE
)

# 2) BB.PIP.sum.loci.envi.csv
count_BB <- read_and_count("BB.PIP.sum.loci.envi.csv")
plot_envi_pip_frequency(
  count_data   = count_BB,
  outfile_base = "BB_envi_pip_frequency",
  width_cm     = 10,
  height_cm    = 8,
  use_percent  = FALSE,   # set TRUE if you want percentages instead
  show_values  = TRUE
)
