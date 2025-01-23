rm(list = ls())
gc()

setwd("/mnt/data00/minghui/project_NMF/04_Test_NMF_Program_0120/")
library(Seurat)
library(dplyr)
library(cowplot)
library(harmony)
library(readr) 
library(tidyr)
library(gplots)
library(reshape2)
library(paletteer)
library(AUCell)
library(pheatmap)
library(gridExtra)
suppressPackageStartupMessages(library(NMF))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(scibet))
suppressMessages(library(viridis))
suppressMessages(library(ggsci))

d_palettes<- palettes_d_names
mycol<-paletteer_d( "ggsci::default_igv",n=51)

source('/mnt/data00/minghui/Fuctions/Markers_kno.R')
source('/mnt/data00/minghui/Fuctions/seurat_functions.R')

library(ComplexHeatmap)
library(circlize)

library(scater)


library(htmlwidgets)

library(GSVA)
library(ggplot2)
library(ggsignif)
library(ggbeeswarm)

cluster_cols <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")


my36colors <- c('#53A85F', '#F1BB72', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175', '#E5D2DD', '#E0D4CA', '#F3B1A0' )

col.list <- list(Tissuetype = c('BE-IM' = my36colors[1], 
                                'BSCJ' = my36colors[2], 
                                'CAG' = my36colors[3],
                                'CIM' = my36colors[4],
                                'Colon' = my36colors[5],
                                'E-GM' = my36colors[6],
                                'GIM' = my36colors[7],
                                'Ileum' = my36colors[8],
                                'NAG' = my36colors[9],
                                'ND' = my36colors[10],
                                'NGB' = my36colors[11],
                                'NGC' = my36colors[12],
                                'NSCJ' = my36colors[13],
                                'Rectum' = my36colors[14]
))

AllTissue_merge <- readRDS('./01_AllTissue_merge.rds')

filtered_modules_list <- readRDS('01_filtered_modules_list.rds')

############################################################

#-----------------#AUCell -----------------

############################################################


DefaultAssay(AllTissue_merge) <- 'RNA'

cells_rankings <- AUCell_buildRankings(AllTissue_merge@assays$RNA@data,splitByBlocks=TRUE) 

cells_AUC <- AUCell_calcAUC(filtered_modules_list, cells_rankings, 
                            aucMaxRank=nrow(cells_rankings)*0.1)


moduleScore <- getAUC(cells_AUC)

names(filtered_modules_list)

for (geneSet in names(filtered_modules_list)) {
  
  AUCell_auc <- as.numeric(getAUC(cells_AUC)[geneSet, ])
  
  
  AllTissue_merge@meta.data[, geneSet] <- AUCell_auc
}


############################################################

#-----------------#Violin plot -----------------

############################################################

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}




vplot1 <- ggplot(AllTissue_merge@meta.data, aes_string(x = 'Tissue_in_paper', y = names(filtered_modules_list)[1], fill = 'Detailed_Cell_Type')) +
  geom_split_violin(trim = TRUE, colour = NA) +
  geom_point(stat = 'summary', fun = mean, position = position_dodge(width = 0.9)) +
  stat_summary(fun.min = function(x) { quantile(x)[2] },
               fun.max = function(x) { quantile(x)[4] },
               geom = 'errorbar', color = 'black', width = 0.01, size = 0.5,
               position = position_dodge(width = 0.9)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 17),  
        axis.text.y = element_text(size = 16), 
        axis.title = element_text(size = 16), 
        legend.text = element_text(size = 10),  # Reduce legend text size
        legend.title = element_text(size = 12),  # Reduce legend title size
        axis.line = element_line(size = 0.7), 
        panel.border = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_manual(values = cluster_cols) +
  stat_compare_means(aes(group = Detailed_Cell_Type), data = AllTissue_merge@meta.data, label = "p.signif", label.x = 1.5, size = 6)

print(vplot1)


# 
# plot_list <- list()
# 
# for (i in seq_along(filtered_modules_list)) {
#   module_name <- names(filtered_modules_list)[i]
#   plot <- ggplot(AllTissue_merge@meta.data, aes_string(x = 'Tissue_in_paper', y = module_name, fill = 'Detailed_Cell_Type')) +
#     ggtitle(module_name) + 
#     geom_split_violin(trim = TRUE, colour = NA) +
#     geom_point(stat = 'summary', fun = mean, position = position_dodge(width = 0.9)) +
#     stat_summary(fun.min = function(x) { quantile(x)[2] },
#                  fun.max = function(x) { quantile(x)[4] },
#                  geom = 'errorbar', color = 'black', width = 0.01, size = 0.5,
#                  position = position_dodge(width = 0.9)) +
#     theme_bw() +
#     theme(axis.text.x = element_text(size = 17),  
#           axis.text.y = element_text(size = 16), 
#           axis.title = element_text(size = 16), 
#           legend.text = element_text(size = 15),
#           legend.title = element_text(size = 15),
#           axis.line = element_line(size = 0.7), 
#           panel.border = element_blank(),
#           panel.grid = element_blank()) +
#     scale_fill_manual(values = cluster_cols) +
#     stat_compare_means(aes(group = Detailed_Cell_Type), data = AllTissue_merge@meta.data, label = "p.signif", label.x = 1.5, size = 6)
#   
#   plot_list[[i]] <- plot
# }


for (i in seq_along(filtered_modules_list)) {
  module_name <- names(filtered_modules_list)[i]
  plot <- ggplot(AllTissue_merge@meta.data, aes_string(x = 'Tissue_in_paper', y = module_name, fill = 'Detailed_Cell_Type')) +
    ggtitle(module_name) + 
    geom_split_violin(trim = TRUE, colour = NA) +
    geom_point(stat = 'summary', fun = mean, position = position_dodge(width = 0.9)) +
    stat_summary(fun.min = function(x) { quantile(x)[2] },
                 fun.max = function(x) { quantile(x)[4] },
                 geom = 'errorbar', color = 'black', width = 0.01, size = 0.5,
                 position = position_dodge(width = 0.9)) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 17),  
          axis.text.y = element_text(size = 16), 
          axis.title = element_text(size = 16), 
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          axis.line = element_line(size = 0.7), 
          panel.border = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none") +  # This line removes the legend
    scale_fill_manual(values = cluster_cols) +
    stat_compare_means(aes(group = Detailed_Cell_Type), data = AllTissue_merge@meta.data, label = "p.signif", label.x = 1.5, size = 6)
  
  plot_list[[i]] <- plot
}

print(plot_list)

pdf('./05_violin_combined_plots.pdf',width = 28,height = 4*length(plot_list)/2)
do.call(grid.arrange, c(plot_list, ncol=2))
dev.off()


pdf('./05_legend1.pdf',width = 14,height =4)
vplot1
dev.off()


saveRDS(AllTissue_merge,'./05_AllTissue_merge.rds')


