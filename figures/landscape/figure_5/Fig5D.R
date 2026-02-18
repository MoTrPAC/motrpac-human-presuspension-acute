library(MotrpacHumanPreSuspensionAnalysis)
library(ggpubr)
library(grid)

hsf1 = SET_TO_ID %>% filter(set == "REACTOME_HSF1_ACTIVATION") %>% pull(set_id)

MotrpacHumanPreSuspensionAnalysis::plot_feature_heatmap(set_id = hsf1,
                                                        selected_ome = "transcript-rna-seq")
