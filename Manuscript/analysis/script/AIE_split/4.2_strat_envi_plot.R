# --- Packages ---
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)     # fast and robust read_*()
library(stringr)

# --- Config ---
setwd("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/AIE_split/split5_envi/")

# Output directory; will be created if it doesn't exist
out_dir <- "../split5_envi_fig"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Source table with traits/SNPs and locus info
dfS <- read.csv(
  "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.addP.csv",
  header = TRUE
)

# Color palette (expand or recycle automatically if needed)
color_lst <- c("#e76254", "#f7aa58", "#ffe6b7", "#aadce0", "#528fad", "#1e466e")

# --- Helper: build one plot safely for a (trait, snp) row ---
plot_one <- function(trait, snp, loci, color_vec = color_lst) {
  # Build input file names
  file_mean <- paste(trait, snp, "score_group_envi_mean.csv", sep = ".")
  file_se   <- paste(trait, snp, "score_group_envi_se.csv",   sep = ".")

  # Check existence
  if (!file.exists(file_mean) || !file.exists(file_se)) {
    message(sprintf("[Skip] Missing files for %s / %s", trait, snp))
    return(invisible(NULL))
  }

  # Read data
  # Expected wide format: columns = Group, Stair_climbing, Age, Alcohol, ...
  df_mean <- suppressMessages(readr::read_delim(file_mean, delim = ",", col_types = cols()))
  df_se   <- suppressMessages(readr::read_delim(file_se,   delim = ",", col_types = cols()))

  # Guard: ensure Group is present
  if (!("Group" %in% names(df_mean)) || !("Group" %in% names(df_se))) {
    message(sprintf("[Skip] 'Group' column missing for %s / %s", trait, snp))
    return(invisible(NULL))
  }

  # Long format
  df_mean_long <- df_mean |>
    pivot_longer(cols = -Group, names_to = "Variable", values_to = "Mean")
  df_se_long <- df_se |>
    pivot_longer(cols = -Group, names_to = "Variable", values_to = "SE")

  # Merge mean + SE
  df_plot <- df_mean_long |>
    left_join(df_se_long, by = c("Group", "Variable")) |>
    mutate(
      # Order quintiles Q1..Q5
      Group = factor(Group, levels = 1:5, labels = paste0("Q", 1:5)),
      Variable = factor(Variable)  # keep stable ordering
    ) |>
    arrange(Variable, Group)

  # Expand/recycle colors to match the number of variables
  n_vars <- nlevels(df_plot$Variable)
  if (length(color_vec) < n_vars) {
    color_vec <- rep(color_vec, length.out = n_vars)
  }

  # Plot
  p <- ggplot(df_plot, aes(x = Group, y = Mean, fill = Variable)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                  position = position_dodge(width = 0.8), width = 0.25) +
    scale_fill_manual(values = color_vec, name = NULL) +
    labs(
      x = "AIE Quintile",
      y = "Standardized value",
      title = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.25, "cm"),
      axis.title = element_text(size = 9),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8)
    )

  # Save (size in cm; 6x5 cm works well for multi-panel layouts)
  out_file <- file.path(out_dir, paste(loci, trait, snp, "ScoreStratEnvi.bar.pdf", sep = "."))
  ggsave(filename = out_file, plot = p, width = 6, height = 5, units = "cm", device = cairo_pdf)
  message(sprintf("[Saved] %s", out_file))
  invisible(p)
}

# --- Main loop ---
# If your CSV has columns named exactly like this:
#   trait, trait_leading_snp, GenomicLocus_hg38
# adapt if names differ.
needed_cols <- c("trait", "trait_leading_snp", "GenomicLocus_hg38")
if (!all(needed_cols %in% names(dfS))) {
  stop(sprintf("Required columns %s not found in dfS.", toString(needed_cols)))
}

for (i in seq_len(nrow(dfS))) {
  trait <- dfS$trait[i]
  snp   <- dfS$trait_leading_snp[i]
  loci  <- dfS$GenomicLocus_hg38[i]
  try(plot_one(trait, snp, loci))
}

# --- Optional: Facet version (per-variable panels) ---
# To switch to facets, replace geom_col/geom_errorbar block with:
#   p <- ggplot(df_plot, aes(x = Group, y = Mean)) +
#     geom_col(fill = "#528fad", width = 0.7) +
#     geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.25) +
#     facet_wrap(~ Variable, nrow = 1, scales = "free_y") + ...
