# ===========================
# Packages
# ===========================
library(ggplot2)
library(dplyr)
library(readr)      # robust CSV reader
library(stringr)    # for safe file names

# ===========================
# Paths
# ===========================
out_dir <- "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/AIE_split/plot_main/"
in_file <- file.path(out_dir, "..", "PT_trait_specific.csv")  # ../PT_trait_specific.csv

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ===========================
# Load data
# ===========================
# Use readr::read_csv for speed and type safety
dfS <- read_csv(in_file, show_col_types = FALSE)

# Ensure numeric columns are numeric (avoid "character" issues)
num_cols <- c("beta1","beta2","beta3","beta4","beta5",
              "se1","se2","se3","se4","se5",
              "p1","p2","p3","p4","p5")
present <- intersect(num_cols, names(dfS))
dfS <- dfS %>%
  mutate(across(all_of(present), ~ suppressWarnings(as.numeric(.))))

# ===========================
# Helper: safe filename
# ===========================
safe_name <- function(x) str_replace_all(x, "[^A-Za-z0-9._-]", "_")

# ===========================
# Plot per row
# ===========================
parts <- 5
x_labels <- paste0("Q", seq_len(parts))

for (i in seq_len(nrow(dfS))) {

  loci  <- dfS$GenomicLocus_hg38[i]
  trait <- dfS$trait[i]
  snp   <- dfS$snp[i]

  beta_vec <- c(dfS$beta1[i], dfS$beta2[i], dfS$beta3[i], dfS$beta4[i], dfS$beta5[i])
  se_vec   <- c(dfS$se1[i],   dfS$se2[i],   dfS$se3[i],   dfS$se4[i],   dfS$se5[i])
  p_vec    <- c(dfS$p1[i],    dfS$p2[i],    dfS$p3[i],    dfS$p4[i],    dfS$p5[i])

  # Build row-wise plotting data
  df <- tibble::tibble(
    x  = factor(seq_len(parts), labels = x_labels),
    beta = beta_vec,
    se   = se_vec,
    p    = p_vec,
    p_lab = sprintf("%.2e", p_vec)  # scientific notation for labels
  )

  # y-limits with padding; always include 0 for geom_col
  y_low  <- df$beta - df$se
  y_high <- df$beta + df$se
  rng <- range(c(0, y_low, y_high), na.rm = TRUE)  # <-- include 0
  pad <- 0.10 * diff(rng)
  if (!is.finite(pad) || pad == 0) pad <- 1
  y_lim <- c(rng[1] - pad, rng[2] + pad)


  # Main figure
  fig <- ggplot(df, aes(x = x, y = beta)) +
    geom_col(fill = "#1e466e", width = 0.6, na.rm = TRUE) +                                 # bars
    geom_errorbar(aes(ymin = beta - se, ymax = beta + se), width = 0.12, size = 0.35,
                  color = "black", na.rm = TRUE) +                                           # error bars
    geom_text(aes(label = p_lab), vjust = -0.6, size = 2, na.rm = TRUE) +                    # p-value labels
	geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_y_continuous(labels = function(y) sprintf("%.3f", y), limits = y_lim) +            # y ticks & limits
    labs(x = "AIE quintile", y = "SNP main effects") +
    theme_classic(base_family = "Arial") +
    theme(
      text = element_text(size = 7),
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      legend.position = "none"
    )

  # Output base name (make it filesystem-safe)
  base <- file.path(out_dir, paste0(safe_name(loci), ".", safe_name(trait), ".", safe_name(snp), ".beta"))

  # Save PDF and PNG (6x5 cm; PNG at 300 dpi)
  ggsave(paste0(base, ".pdf"), plot = fig, width = 6, height = 5, units = "cm", device = cairo_pdf)
}
