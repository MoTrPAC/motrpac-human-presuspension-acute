rm(list=ls())
session_packages = c(
  "colorRamp2",
  "ComplexHeatmap",
  "ComplexUpset",
  "conflicted",
  "dplyr",
  "ggplot2",
  "ggpubr",
  "magrittr",
  "MotrpacHumanPreSuspensionAnalysis",
  "purrr",
  "stringr",
  "tidyverse",
  "TMSig"
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

################################################################################
# ora1 <- read.csv(file.path(here(), "figures", "landscape", "precovid_ORA_cross_tissue_newDA_newORA.csv"),
#                  header = T,
#                  stringsAsFactors = F,
#                  check.names = F) %>%
#   filter(contrast!="Adipose_Blood_Muscle")
# shared_pathways_list <- read.csv(file.path(here(), "figures", "landscape", "precovid_landscape_figure2_shared_curated_terms.csv"),
#                                  header = T, stringsAsFactors = F,
#                                  check.names = F) %>%
#   select(set_short) %>%
#   pull()
#
# plot_df1 <- ora1 %>%
#   filter(set_short %in% shared_pathways_list)

plot_df1 = readxl::read_excel(path = file.path(here(), "figures", "landscape", "curated_pathways.xlsx"),
                              sheet = "2F") %>%
  select(contrast, set_short, statistic_column, adj_p_value) %>%
  pivot_wider(id_cols = set_short,
              names_from = contrast,
              values_from = c(adj_p_value, statistic_column),
              names_sep = "@") %>%
  tidyr::pivot_longer(cols = -set_short,
                      names_to = c(".value", "contrast"),
                      names_sep = "@") %>%
  mutate(contrast = factor(contrast,
                           levels = c("Adipose_only",
                                      "Muscle_only",
                                      "Blood_only",
                                      "Adipose_and_Muscle_only",
                                      "Blood_and_Muscle_only"))) %>%
  arrange(contrast)
levels(plot_df1$contrast) <- gsub(" only","",
                                  gsub("_", " ", levels(plot_df1$contrast)))

# jpeg(file = "landscape_figure2C_ORA_newDA_newORA.jpg",
#      width = 8.5*1.25,height=7.5*1.25,units = "in",res=1200)
pdf(file = "landscape_figure2F_heatmap_newDA_newORA.pdf",
     width = 8.5*1.25,height=7.5*1.25)
TMSig::enrichmap(plot_df1,
                 set_column = "set_short",
                 n_top=35,
                 statistic_column = "statistic_column",
                 contrast_column = "contrast",
                 padj_column = "adj_p_value",
                 padj_legend_title = "adj p \n(background)",
                 padj_fill = "grey70",
                 colors = c("white","#800000"),
                 #colors = c("white", "#543483"),
                 heatmap_args = list(name = "-log(p)",
                                     na_col="grey90",
                                     cluster_columns=FALSE,
                                     cluster_rows=TRUE,
                                     border = T,
                                     rect_gp = gpar(col = "grey70",
                                                    lwd = 2)))
dev.off()
