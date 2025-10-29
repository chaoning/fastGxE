library(ggplot2)
library(readr)

setwd("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_overfit_sum/")
df <- read_csv("mmsusie_PIP_summary.csv")


png("mmsusie_PIP_calibration_plot.png", width = 8, height = 7, res = 300, units = "cm")


ggplot(df, aes(x = mean_PIP, y = in_lst_rate)) +
  geom_point(size = 3, color = "#e76254") +                             
  geom_errorbar(aes(ymin = in_lst_rate - 2 * SE, 
                    ymax = in_lst_rate + 2 * SE), 
                width = 0.01, color = "#e76254") +                    
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  
  labs(
    x = "Mean PIP",
    y = "Proportion of signals"
  ) +
  theme_minimal(base_size = 14) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
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

dev.off()
