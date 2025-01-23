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
library(reshape2)
library(paletteer)
library(corrplot)
library(cowplot)
library(psych)
library(corrr)

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

AllTissue_merge <- readRDS('./05_AllTissue_merge.rds')

filtered_modules_list <- readRDS('01_filtered_modules_list.rds')

############################################################

#-----------------#Correlation between modules-----------------

############################################################

dim(AllTissue_merge@meta.data[names(filtered_modules_list)])

module_data <- data.frame(AllTissue_merge@meta.data[names(filtered_modules_list)])

cor_matrix <- cor(module_data , method = "spearman")
p_matrix <- cor.mtest(module_data , method = "spearman")$p

threshod <- 0.6
cor_matrix[abs(cor_matrix) <= threshold | p_matrix >= 0.05] <- 0

my_color = rev(paletteer_d("RColorBrewer::RdYlBu"))
my_color = colorRampPalette(my_color)(10) 

pdf('./07_Heatmap.pdf',width = 15,height = 15)
corrplot(cor_matrix, type = "upper", 
         method = "pie",
         order = "hclust", 
         col = my_color,
         tl.col = "black", 
         tl.srt = 45)

dev.off()

rownames(cor_matrix)
ph1 <- pheatmap(cor_matrix,color = my_color)
ph1


rownames <- rownames(cor_matrix)

pdf('./06_Corr_modules_Heatmap.pdf',height =9, width = 13,onefile = F)

# Create a vector of labels
labels <- ifelse(grepl("BE_", rownames), "BE", 
                 ifelse(grepl("intestine", rownames), "intestine", 
                        ifelse(grepl("stomach", rownames), "stomach", "others")))

annotation_row <- data.frame(Group = labels)
annotation_colors <- list(Group = c(BE = "firebrick3", intestine = "limegreen", stomach = "cornflowerblue"))
rownames(annotation_row) <- rownames

pheatmap(cor_matrix, annotation_row = annotation_row,show_colnames = F,annotation_colors = annotation_colors)
dev.off()
# 
# 
# # Create a vector of labels
# labels <- ifelse(grepl("BE_", rownames), "BE", 
#                  ifelse(grepl("intestine", rownames), "intestine", 
#                         ifelse(grepl("stomach", rownames), "stomach", "others")))
# 
# annotation_row <- data.frame(Group = labels)
# annotation_colors <- list(Group = c(BE = "firebrick3", intestine = "limegreen", stomach = "cornflowerblue"))
# rownames(annotation_row) <- rownames
png('./06_Corr_modules_Heatmap.png',height =9, width = 13,units = 'in',res = 600)

pheatmap(cor_matrix, annotation_row = annotation_row,show_colnames = F,annotation_colors = annotation_colors)
dev.off()


# > filtered_modules_list$BE_modules5
# [1] "ADH1C"    "C15orf48" "CAPS"     "CLDN15"   "GSTA1"    "MT1E"     "MT1F"     "MT1G"     "MT1H"     "MT1M"     "MT1X"     "MT2A"    
# [13] "NKX6.2"   "NUPR1"    "SCGB2A1"  "TAGLN2"   "UGT2B15"
# 
# 
# cor_matrix <- module_data %>%
#   correlate(method = 'spearman') %>%
#   rearrange()
# 
# # Create network plot with smaller label font and custom colors
# network_plot(cor_matrix, min_cor = 0 ,colours= c(  "#2D3FA6","#F2E1AE","#F21313"), ) +
#   theme(axis.text = element_text(size = 10),  # Smaller label font size
#         legend.text = element_text(size = 10),  # Smaller legend text size
#         legend.title = element_text(size = 12))
# # +  # Smaller legend title size
#   # scale_color_manual(values = c("red", "blue", "green"))  # Custom colors
# 
# 
# 




