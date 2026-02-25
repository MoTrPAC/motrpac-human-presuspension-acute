# Figure Naming and Organization

This document aims to give relevant background to support replicatability. The code for each of the figures is found in the locations described in each of the corresponding files, to allow others to replicate and reinterpret the code and values as easily as possible. In some cases, the figures were then post-processed to keep naming, font, color, etc. consistent, purely for aesthetic purposes. No values have been changed.

Note: several different plots made with `plot_single_feature` are 
combined into `single_feature_plots_for_manuscript.Rmd`. Each of these are also highlighted within each of the panels specifically.

Similarly, several different plots made for timewise comparisons or fuzzy c-means are found in the `figures/landscape/CAMERA-PR/` or `figures/landscape/fuzzy_c-means/` respectively. Whenever code found in one of those files is relevant, that is specified below as well.

For each of the plots, if we wanted to highlight selected pathways (all of which are FDR < 0.05, if not otherwise specified), they are highlighted in: "figures/landscape/curated_pathways.xlsx" in the sheet according to the pathway. The rationale for highlighting specific pathways is: (1) to avoid situations where all highlighted pathways have redudant sounding names (2) to not purely prioritize based on effect size/variance. In every case, the Benjamini Hotchberg adjustment was done across all pathways tested, so we are NOT selecting pathways prior to multiple testing. 

For CAMERA-PR Plots: `figures/landscape/CAMERA-PR/data/Curated_multi-tissue...`

For all other panels: `figures/landscape/curated_pathways.xlsx`

## Figure 1:


A: Illustrator, no accompanying code.

B: Illustrator, no accompanying code.

C: `figures/landscape/figure_1/updated_available_data_heatmap.Rmd/` -> Requires `MotrpacHumanPreSuspensionData` access to fully recreate.

### S1:

S1A: `figures/landscape/figure_1/supplemental/thousand_genomes_pca/` -> Requires WGS (dbGaP only) access to fully recreate.

S1B: `figures/landscape/figure_1/supplemental/omic_overlap_visualizations_S1B.Rmd` -> Requires `MotrpacHumanPreSuspensionData` access to fully recreate.

S1C: `figures/landscape/figure_1/supplemental/pca_var_exp_plots` -> Requires `MotrpacHumanPreSuspensionData` access to fully recreate.

S1D: `figures/landscape/figure_1/supplemental/pca_var_exp_plots` -> Requires `MotrpacHumanPreSuspensionData` access to fully recreate.


### Tables S1: 

Participants per group: Directly taken from the `pheno` object. -> Requires `MotrpacHumanPreSuspensionData` access to fully recreate.

Participants per group/sex: Directly taken from the `pheno` object. -> Requires `MotrpacHumanPreSuspensionData` access to fully recreate.

Participants per tissue/omics assay: `figures/landscape/figures_1/supplemental/omic_overlap_S1B.Rmd` -> Requires `MotrpacHumanPreSuspensionData` access to fully recreate.

Outliers: Each of the outliers were first identified from `QC/ relevant ome QC vignette` as described in the methods.

The full list of outliers is saved in `MotrpacHumanPreSuspensionAnalysis::OUTLIERS`

## Figure 2:

A: `figures/landscape/figure_2/DA_feature_matrix.Rmd`

B: `figures/landscape/figure_2/figure2B_upsetPlot_barPlot.R`

C: `figures/landscape/figure_2/figure2C_heatmap.R`

D: `figures/landscape/single_feature_plots_for_manuscript.Rmd`

E:`figures/landscape/figure_2/figure2E_heatmap.R`

F: `figures/landscape/figure_2/figure2F_ORA_heatmap.R`

### Tables S2: 

COVARIATES: Each of the covariates were first identified from `QC/ relevant ome QC vignette` as described in the methods.
The full list of covariates is saved in `MotrpacHumanPreSuspensionAnalysis::COVARIATES_FILE`

Total number of features: `figures/landscape/supplemental/features_database.R`

Number of DA features: `figures/landscape/figure_2/DA_feature_matrix.Rmd` 

DA_PercentOverlap: `figures/landscape/figure_2/supplemental/figureS2C_upsetPlots.R`

CrossTissue_ORA: `figures/landscape/figure_2/figure2B_upsetPlot_barPlot.R`




### S2: 

S2A: `figures/landscape/figure_2/supplemental/DA_hit_type_matrix.Rmd`

S2B: `figures/landscape/single_feature_plots_for_manuscript.Rmd`

S2C: `figures/landscape/figure_2/supplemental/figureS2C_upsetPlots.R`

S2D: `figures/landscape/single_feature_plots_for_manuscript.Rmd`

S2E: `figures/landscape/figure_2/supplemental/figureS2E_phospho_ORA_heatmap.R`

S2F: `figures/landscape/figure_2/supplemental/figureS2F_upsetPlot_by_ome.R`

## Figure 3: 

3A-B: `figures/landscape/CAMERA-PR/curated_CAMERA-PR_heatmaps.R`

3C-J: `figures/landscape/figure_3/c-means_main_figure.Rmd`

### Tables S3:

FCM_ORA_filtered: Taken directly from `MotrpacHumanPreSuspensionAnalysis::FCM_ORA`. See `MotrpacHumanPreSuspensionAnalysis/data-raw/` for the implementation of the package functions themselves. 

### S3:

S3A: `figures/landscape/CAMERA-PR/curated_CAMERA-PR_heatmaps.R`

S3B-D: `figures/landscape/figure_3/c-means_main_figure.Rmd`

## Figure 4:

A: Built directly in illustrator

B-E: `figures/landscape/figure_4/cross_tissue_rna_chris.R` -> Requires `MotrpacHumanPreSuspensionData` access to fully recreate.

F-H: `figures/landscape/single_feature_plots_for_manuscript.Rmd`

J-L: `figures/landscape/figure_4/cross_tissue_metab_chris.R` -> Requires `MotrpacHumanPreSuspensionData` access to fully recreate.

### Tables S4:

Cross Tissue RNAseq PLIER: `figures/landscape/figure_4/cross_tissue_rna_chris.R` -> Requires `MotrpacHumanPreSuspensionData` access to fully recreate.

Transcript Affinity: `figures/landscape/figure_4/cross_tissue_rna_chris.R` -> Requires `MotrpacHumanPreSuspensionData` access to fully recreate.

Cross Tissue Metabolomics PLIER: `figures/landscape/figure_4/cross_tissue_metab_chris.R` -> Requires `MotrpacHumanPreSuspensionData` access to fully recreate.

Metabolite Affinity: `figures/landscape/figure_4/cross_tissue_metab_chris.R` -> Requires `MotrpacHumanPreSuspensionData` access to fully recreate.

### S4

A-H: `figures/landscape/figure_4/cross_tissue_rna_chris.R` -> Requires `MotrpacHumanPreSuspensionData` access to fully recreate.

I-O: `figures/landscape/figure_4/cross_tissue_metab_chris.R` -> Requires `MotrpacHumanPreSuspensionData` access to fully recreate.

## Figure 5: 

A-C: `figures/landscape/figure_5/re_vs_ee_visualizations.Rmd`

D: `figures/landscape/figure_5/Fig5D.R`

E: `figures/landscape/figure_5/re_vs_ee_visualizations.Rmd` + Illustrator for the icons

F: `figures/landscape/single_feature_plots_for_manuscript.Rmd`

### Tables S5:

EEmuscleonly_ORA_features: `figures/landscape/figure_5/re_vs_ee_visualizations.Rmd`

REmuscleonly_ORA_features: `figures/landscape/figure_5/re_vs_ee_visualizations.Rmd`

all_muscle_re_ee_only_ORA: `figures/landscape/figure_5/re_vs_ee_visualizations.Rmd`

triangle_selected_pathway: `figures/landscape/figure_5/re_vs_ee_visualizations.Rmd`

da_features_with_EE_RE_opposite: `figures/landscape/figure_5/re_vs_ee_visualizations.Rmd`

### S5: 

A: `figures/landscape/supplemental/upset_plot_S5A.R`

B: `figures/landscape/single_feature_plots_for_manuscript.Rmd`

C: `figures/landscape/figure_5/re_vs_ee_visualizations.Rmd` + Illustrator

D: `figures/landscape/single_feature_plots_for_manuscript.Rmd`

## Figure 6: 
For just about all panels/tables for fig 6 -> Requires `MotrpacHumanPreSuspensionData` access to fully recreate.

A - H: `figures/landscape/figure_6/homer_analysis_chris.R`

### Tables S6: 

Transcription factor enrichment and expression: `figures/landscape/figure_6/homer_analysis_chris.R` 

### S6: 

A-B: `figures/landscape/figure_6/homer_analysis_chris.R`

## Figure 7: 

For just about all panels/tables for fig 7 -> Requires `MotrpacHumanPreSuspensionData` access to fully recreate.

A: Illustrator

B: Cytoscape. See also: `figures/landscape/figure_7/Figure7_SCION.ai`

C-E: `figures/landscape/figure_7/scion_figures.Rmd`.  -> Requires `MotrpacHumanPreSuspensionData` access to fully recreate.
See also: `figures/landscape/figure_7/Figure_7_script.R` -> Requires `MotrpacHumanPreSuspensionData` access to fully recreate.

F: `figures/landscape/single_feature_plots_for_manuscript.Rmd`

G: Illustrator

### Tables S7:

Combined EE+RE Edge Table: Cytoscape + `figures/landscape/figure_7/scion_figures.Rmd` 

Combined EE+RE Node Table: Cytoscape + `figures/landscape/figure_7/scion_figures.Rmd`

TFEB ChIP Targets: `figures/landscape/figure_7/scion_figures.Rmd`

EE_3.5_4_hr_up: `figures/landscape/figure_7/scion_figures.Rmd`

EE_3.5_4_hr_down: `figures/landscape/figure_7/scion_figures.Rmd`

RE_3.5_4_hr_up: `figures/landscape/figure_7/scion_figures.Rmd`

RE_3.5_4_hr_down: `figures/landscape/figure_7/scion_figures.Rmd`


### S7:

A: `figures/landscape/figure_7/scion_figures.Rmd` 

B: `figures/landscape/single_feature_plots_for_manuscript.Rmd`

C: `figures/landscape/figure_7/scion_figures.Rmd`

D: `figures/landscape/single_feature_plots_for_manuscript.Rmd`

E: `figures/landscape/figure_7/scion_figures.Rmd`

## Figure 8 

For just about all panels/tables for fig 8 -> Requires `MotrpacHumanPreSuspensionData` access to fully recreate.


A: `figures/landscape/figure_8/secretome_final.R`

B: `figures/landscape/figure_8/single_feature_plots_for_manuscript.Rmd`

C: `figures/landscape/figure_8/single_feature_plots_for_manuscript.Rmd`

D: `figures/landscape/figure_8/landscape_ccn1.R`

E: `figures/landscape/figure_8/secretome_final.R`


### Tables S8:

Exercine Candidate List: `figures/landscape/figure_8/secretome_final.R`

### S8:

A-B: `figures/landscape/figure_8/secretome_final.R`

C: `figures/landscape/figure_8/single_feature_plots_for_manuscript.Rmd`

D: `figures/landscape/figure_8/secretome_final.R`

