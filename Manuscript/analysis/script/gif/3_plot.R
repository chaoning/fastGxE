# 加载必要的库
library(ggplot2)
library(reshape2) # 用于 melt
setwd("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/gif/")

# 读取数据
data <- read.csv("BB.csv", header = TRUE)

# 选择相关列并转化为长格式
relevant_columns <- c('Trait', 'X0.5', 'X0.05', 'X0.01', 'X0.001')
tidy_data <- melt(data[, relevant_columns], 
                  id.vars = 'Trait', 
                  variable.name = 'Quantile', 
                  value.name = 'Value')

# 创建箱线图
p <- ggplot(tidy_data, aes(x = Quantile, y = Value, fill = Quantile)) +
  geom_boxplot(fill = "#e76254", color = "#1e466e", size = 0.2, outlier.size = 0.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.5) + # 添加 y = 1 的虚线
  scale_x_discrete(labels = c("X0.5" = "0.5", "X0.05" = "0.05", "X0.01" = "0.01", "X0.001" = "0.001")) +
  labs(x = "p-values quantile",
       y = "Genomic inflation factor") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(), # 移除边框
    axis.line.x = element_line(colour = "black"), # 保留下边框
    axis.line.y = element_line(colour = "black"), # 保留左边框
    legend.position = "none", # 移除图例
    text = element_text(size = 6),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8, colour = "black"),
    axis.text.y = element_text(size = 8, colour = "black")
  ) +
  coord_cartesian(ylim = c(0.0, 1.9)) # 避免裁剪超出范围的点

# 保存到 PDF 文件
pdf("BB_Quantile_Boxplots.pdf", width = 6 / 2.54, height = 5 / 2.54) # 设置文件名和尺寸（厘米为单位）
print(p) # 将图形输出到 PDF
dev.off() # 关闭 PDF 文件
