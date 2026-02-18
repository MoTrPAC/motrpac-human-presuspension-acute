#NFIC phospho heatmap

#load packages
library(MotrpacHumanPreSuspension)
library(MotrpacHumanPreSuspensionData)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(cowplot)
library(tidyverse)

#load DA and filter for NFIC phosphosites
da <- MUSCLE_PROT_PH_DA
da <- da %>% dplyr::filter(contrast_category=="EE-CON" | contrast_category=="RE-CON")
da <- left_join(da,HUMAN_FEATURE_TO_GENE,by="feature_id")

NFIC <- da %>% dplyr::filter(gene_symbol=="NFIC")
NFIC_filter <- NFIC %>% dplyr::select("feature_id","contrast_short","logFC","p_value","adj_p_value")
NFIC_wide <- pivot_wider(NFIC_filter, names_from=contrast_short, values_from=c("logFC","p_value","adj_p_value")) 
NFIC_wide <- NFIC_wide[,grepl("feature_id|delta-delta",colnames(NFIC_wide))]
rownames(NFIC_wide) <- NFIC_wide$feature_id

#separate out FC and adj p-value for the heatmap
fc <- NFIC_wide %>% dplyr::select(contains("logFC")) %>% as.matrix()
pval <- NFIC_wide %>% dplyr::select(contains("adj_p_value")) %>% as.matrix()
rownames(fc) = rownames(pval) = NFIC_wide$feature_id

#set heatmap annotations
cdesc <- data.frame(group=c("ADUEndur","ADUEndur","ADUEndur","ADUResist","ADUResist","ADUResist"),
                    time=c("post_15_30_45_min","post_3.5_4_hr","post_24_hr","post_15_30_45_min","post_3.5_4_hr","post_24_hr"))
ha <- HeatmapAnnotation(
  Modality = cdesc$group,
  Timepoint = cdesc$time,
  col = list(
    Modality = HUMAN_EXERCISE_GROUP_COLORS,
    Timepoint = HUMAN_ACUTE_TIMEPOINT_COLORS
  ))

#Heatmap dimensions
a <- unit(3, "mm")
width <- ncol(fc) * a
height <- nrow(fc) * a

#color scale
 col_fun <- circlize::colorRamp2(
   breaks = c(-0.5,0,0.5),
   colors = c("blue", "white", "red"))

#plot heatmap
ht <- Heatmap(fc,
              heatmap_legend_param = list(
                title = "logFC",
                legend_direction = "vertical",
                legend_width = unit(50, "mm")
              ),
              width = width,
              height = height,
              top_annotation = ha,
              #left_annotation = row_ha,
              cluster_columns = F,
              cluster_rows=T,
              show_column_names = F,
              show_row_names = T,
              col = col_fun,
              column_title_rot = 90,
              row_title_rot = 0,
              column_names_gp = grid::gpar(fontsize = 8),
              row_names_gp = grid::gpar(fontsize = 8),
              
              cell_fun = function(j, i, x, y, width, height, fill) {
              #Point size
              if(pval[i, j] < 0.05) {
                gb = textGrob("*")
                gb_w = convertWidth(grobWidth(gb), "mm")
                gb_h = convertHeight(grobHeight(gb), "mm")
                grid.text("*", x, y - gb_h*0.5 + gb_w*0.4)
              }
              }
)

#save to PDF
pdf("NFIC_phospho_heatmap.pdf",height=8,width=8)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     legend_gap = unit(10, "mm"))
dev.off()
