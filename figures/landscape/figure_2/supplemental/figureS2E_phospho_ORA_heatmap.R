rm(list=ls())
session_packages = c(
  "colorRamp2",
  "ComplexHeatmap",
  "ComplexUpset",
  "conflicted",
  "dplyr",
  "ggplot2",
  "ggpubr",
  "hrbrthemes",
  "magrittr",
  "MotrpacHumanPreSuspensionAnalysis",
  "purrr",
  "stringr",
  "tidyverse",
  "TMSig",
  "here"
)
for (lib_name in session_packages){
  tryCatch({library(lib_name,character.only = TRUE)}, error = function(e) {
    print(paste("Cannot load",lib_name,", please install"))
  })
}

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("union","base")
conflict_prefer("intersect","base")
conflict_prefer("setdiff", "base")

######
# phospho heatmap
#
# phospho_curated_terms <- read.csv("precovid_landscape_figureS2D_phospho_curated_terms.csv",
#                                   header = T, stringsAsFactors = F, check.names = F)
# ptmsigdb_ora_full <- read.csv("figure2_PTMSigDB_ORA_results_oct2025.csv",
#                               header = T, stringsAsFactors = F, check.names = F) %>%
#   mutate(statistic_column=-log10(p_value))
# ptmsigdb_ora_filt = ptmsigdb_ora_full %>%
#   filter(set_id %in% phospho_curated_terms$set_id)
# write.csv(ptmsigdb_ora_filt, file = "curated_S2E.csv", row.names = FALSE)

phospho_plot_df = readxl::read_excel(path = file.path(here(), "figures", "landscape", "curated_pathways.xlsx"),
                              sheet = "S2E") %>%
  mutate(set_short=gsub("PTMSIGDB_","",set_short)) %>%
  mutate(contrast = factor(contrast,
                           levels = c(
                             "Muscle_only",
                             "Adipose_only",
                             "Adipose_Muscle"))) %>%
  arrange(contrast)
# Nicer labels for plot
levels(phospho_plot_df$contrast) <- gsub("_", " ",
                                         gsub("_only","",
                                              levels(phospho_plot_df$contrast)))

pdf(file = "landscape_figureS2E_heatmap_newDA_newORA.pdf",
    width = 8.5*1.25,height=7.5*1.25)
TMSig::enrichmap(phospho_plot_df,
                 set_column = "set_short",
                 n_top=30,
                 statistic_column = "statistic_column",
                 contrast_column = "contrast",
                 padj_column = "adj_p_value",
                 padj_legend_title = "adj p \n(background)",
                 padj_fill = "grey80",
                 colors = c("white", "#543483"),
                 heatmap_args = list(name = "-log(p)",
                                     na_col="grey90",
                                     cluster_columns=FALSE,
                                     cluster_rows=TRUE,
                                     border = T,
                                     rect_gp = gpar(col = "grey70",
                                                    lwd = 2)))
dev.off()
