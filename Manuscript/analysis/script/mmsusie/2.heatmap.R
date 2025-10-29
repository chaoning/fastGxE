library(ComplexHeatmap)

# args <- commandArgs(trailingOnly = TRUE)


setwd("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/")
file = "PT_PIP_plot.csv"
data <- read.csv(file,header=T, row.names=1)

keep <- c(TRUE, colSums(data[, -1] != 0, na.rm = TRUE) > 0)
data <- data[, keep, drop = FALSE]

print(dim(data))
data[,1] = 0
vec = rep("Y", dim(data)[1])
vec[9]= "N"
vec[10]= "N"

#data_row2 <- read.table("/net/mulan/home/chaon/WORK/GxE/Analysis/Physical_measures/sum/3.man_pheatmap.20M.row_order", sep="\t")
#data <- data[match(data_row2$V1, rownames(data)), ]
#data <- data[1:33,]

#file = paste0("heatmap_", eee, ".png")
#png(file=file,width=13,height=15,units="cm",res=600)
file = paste0("PT_PIP_heatmap.pdf")
cairo_pdf(file=file,width=8/2.54,height=10/2.54)
ht<-Heatmap(data,
        name = " ",
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = T,
        show_row_names = T,
        row_names_side =  'left',
        column_title = NULL,
        column_names_gp = gpar(fontsize = 6.5, fontfamily = "Arial"),
        row_names_gp = gpar(fontsize = 6, fontfamily = "Arial"),
        col = c('#e76254', "#ffd06f", '#FEFCFC', "#72bcd5", '#1e466e'),
        border = 'black',
        rect_gp = gpar(col = "grey", lwd = 0.5),
        heatmap_legend_param = list(title = "PIP",  
        title_gp = gpar(fontsize = 8, fontface = "plain"),
        legend_width = unit(1, "cm"),
        legend_height = unit(2, "cm"),
        grid_width = unit(0.1, "cm"),
        grid_height = unit(0.1, "cm"),
        legend_padding = unit(c(0.5, 0.5, 0.5, 0.5), "mm"),
        width=unit(0.1, "cm"), height=unit(0.1, "cm"),
        color_bar = "discrete", border = "black", 
        labels_gp = gpar(fontfamily = "Arial", fontsize = 8), 
        labels=c(">0.9 +", ">0.5 +", "0~0.5", ">0.5 -", ">0.9 -") , legend_position = "top" # 这里添加了legend_position参数
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (j == 1) {
                grid.text(sprintf("%s", vec[i]), x, y, gp = gpar(fontsize = 8, fontface = "plain"))
            }
        }
        )
draw(ht)
dev.off()



setwd("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/mmsusie/")
file = "BB_PIP_plot.csv"
data <- read.csv(file,header=T, row.names=1)

keep <- c(TRUE, colSums(data[, -1] != 0, na.rm = TRUE) > 0)
data <- data[, keep, drop = FALSE]

print(dim(data))
data[,1] = 0
vec = rep("Y", dim(data)[1])
vec[28]= "N"
vec[36]= "N"

#data_row2 <- read.table("/net/mulan/home/chaon/WORK/GxE/Analysis/Physical_measures/sum/3.man_pheatmap.20M.row_order", sep="\t")
#data <- data[match(data_row2$V1, rownames(data)), ]
#data <- data[1:33,]

#file = paste0("heatmap_", eee, ".png")
#png(file=file,width=13,height=15,units="cm",res=600)
file = paste0("BB_PIP_heatmap.pdf")
cairo_pdf(file=file,width=15/2.54,height=20/2.54)
ht<-Heatmap(data,
        name = " ",
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = T,
        show_row_names = T,
        row_names_side =  'left',
        column_title = NULL,
        column_names_gp = gpar(fontsize = 6.5, fontfamily = "Arial"),
        row_names_gp = gpar(fontsize = 6, fontfamily = "Arial"),
        col = c('#e76254', "#ffd06f", '#FEFCFC', "#72bcd5", '#1e466e'),
        border = 'black',
        rect_gp = gpar(col = "grey", lwd = 0.5),
        heatmap_legend_param = list(title = "PIP",  
        title_gp = gpar(fontsize = 8, fontface = "plain"),
        legend_width = unit(1, "cm"),
        legend_height = unit(2, "cm"),
        grid_width = unit(0.1, "cm"),
        grid_height = unit(0.1, "cm"),
        legend_padding = unit(c(0.5, 0.5, 0.5, 0.5), "mm"),
        width=unit(0.1, "cm"), height=unit(0.1, "cm"),
        color_bar = "discrete", border = "black", 
        labels_gp = gpar(fontfamily = "Arial", fontsize = 8), 
        labels=c(">0.9 +", ">0.5 +", "0~0.5", ">0.5 -", ">0.9 -") , legend_position = "top" # 这里添加了legend_position参数
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (j == 1) {
                grid.text(sprintf("%s", vec[i]), x, y, gp = gpar(fontsize = 8, fontface = "plain"))
            }
        }
        )
draw(ht)
dev.off()
