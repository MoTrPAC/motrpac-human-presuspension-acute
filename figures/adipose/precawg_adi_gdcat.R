### Gene-Derived Correlation Across Tissue (GD-CAT) applying to Motrpac's rat training study.

#see: https://elifesciences.org/reviewed-preprints/88863/reviews for gd cat info.
library(fgsea)
library(msigdbr)
library(MotrpacRatTraining6mo)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(TMSig)

output_folder = ""

gene_sets_human <- msigdbr(species = "Homo sapiens", category = "C5")

# custom function assigning gene names to pass1b bicor data
assign_row_and_col_names_to_bicor <- function(dat, row_names, col_names) {
  # Loop through each condition (e.g., control, training_1w, etc.)
  for (condition in names(dat)) {
    # Check if bicor matrix exists in the current condition
    if (!is.null(dat[[condition]]$bicor)) {
      # Assign row names to the bicor matrix
      rownames(dat[[condition]]$bicor) <- row_names
      # Assign column names to the bicor matrix
      colnames(dat[[condition]]$bicor) <- col_names
    }
  }
  return(dat)
}

generate_fgsea_heatmap <- function(Direction, target_gene, manual = NULL) {

  # Embedded gene sets and conditions
  gene_sets_human <- gene_sets_human
  condition_colors <- c(
    "control" = "#000000B3",
    "training_1w" = "#EDF8B1",
    "training_2w" = "#7FCDBB",
    "training_4w" = "#1D91C0",
    "training_8w" = "#0C2C84"
  )
  contrast_order <- c("training_1w", "training_2w", "training_4w", "training_8w", "control")

  # Function to extract bicor values for the target gene across conditions
  extract_bicor_values <- function(Direction, target_gene) {
    bicor_values <- list()

    for (condition in names(Direction)) {
      bicor_matrix <- Direction[[condition]]$bicor
      if (target_gene %in% colnames(bicor_matrix)) {
        bicor_values[[condition]] <- bicor_matrix[, target_gene]
      } else {
        warning(paste("Gene", target_gene, "not found in condition", condition))
        bicor_values[[condition]] <- NULL
      }
    }
    return(bicor_values)
  }

  # Extract bicor values for the target gene across conditions
  bicor_values <- extract_bicor_values(Direction, target_gene)

  # Convert Ensembl Rat IDs to Human gene symbols
  for (timepoint in names(bicor_values)) {
    rat_ensembl_ids <- names(bicor_values[[timepoint]])
    human_gene_symbols <- RAT_TO_HUMAN_GENE$HUMAN_ORTHOLOG_SYMBOL[match(rat_ensembl_ids, RAT_TO_HUMAN_GENE$RAT_ENSEMBL_ID)]
    names(bicor_values[[timepoint]]) <- human_gene_symbols
  }

  bicor_values <- lapply(bicor_values, function(bicor) {
    bicor[!is.na(bicor)]
  })

  # Initialize a list to store FGSEA results
  fgsea_results_list <- list()

  contrasts <- names(Direction)

  for (contrast in contrasts) {
    control_bicor <- bicor_values[[contrast]]
    ranked_genes_df <- data.frame(gene = names(control_bicor), score = control_bicor)

    ranked_genes_df <- ranked_genes_df[!is.na(ranked_genes_df$gene) & ranked_genes_df$gene != "unknown", ]
    ranked_genes_df <- ranked_genes_df[order(ranked_genes_df$score, decreasing = TRUE),]

    names(ranked_genes_df$score) <- ranked_genes_df$gene

    ranked_genes_df <- ranked_genes_df %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(score = ifelse(all(score > 0), max(score),
                                      ifelse(all(score < 0), min(score),
                                             ifelse(abs(max(score)) > abs(min(score)), max(score), min(score))))) %>%
      ungroup()

    stats <- ranked_genes_df$score
    names(stats) <- ranked_genes_df$gene

    gene_sets <- gene_sets_human %>%
      dplyr::group_by(gs_name) %>%
      dplyr::summarise(genes = list(gene_symbol)) %>%
      ungroup()

    pathways <- setNames(gene_sets$genes, gene_sets$gs_name)

    fgsea_res <- fgseaMultilevel(
      pathways = pathways,
      stats = stats,
      minSize = 15,
      maxSize = 500
    )

    fgsea_res$contrast <- contrast
    fgsea_results_list[[contrast]] <- fgsea_res
  }

  all_fgsea_results <- do.call(rbind, lapply(names(fgsea_results_list), function(contrast) {
    res <- fgsea_results_list[[contrast]]
    res$contrast <- contrast
    return(res)
  }))

  all_fgsea_results <- all_fgsea_results[!grepl("^HP", all_fgsea_results$pathway), ]

  # Select the top 5 pathways for each contrast
  top_pathways <- all_fgsea_results %>%
    group_by(contrast) %>%
    arrange(contrast, pval) %>%
    slice_head(n = 2) %>%
    pull(pathway)

  # Identify significantly enriched pathways in training_8w but NOT in control
  training_8w_sig <- all_fgsea_results %>%
    filter(contrast == "training_8w" & padj < 0.05) %>%
    pull(pathway)

  control_sig <- all_fgsea_results %>%
    filter(contrast == "control" & padj < 0.05) %>%
    pull(pathway)

  extra_training_8w <- setdiff(training_8w_sig, control_sig) %>% head(3)  # Pick top 2 unique to training_8w
  extra_control <- setdiff(control_sig, training_8w_sig) %>% head(3)  # Pick top 2 unique to control

  # Combine all selected pathways
  final_pathways <- unique(c(top_pathways, extra_training_8w, extra_control, manual)) # Include manual pathways

  contrasts_in_top_pathways <- unique(all_fgsea_results$contrast[all_fgsea_results$pathway %in% final_pathways])

  # Ensure the contrasts in `columnAnnotation` match the ones available in the data
  col_anno <- columnAnnotation(
    Condition = factor(contrasts_in_top_pathways, levels = contrasts_in_top_pathways),
    col = list(
      Condition = condition_colors
    ),
    annotation_legend_param = list(
      Condition = list(title = "Condition")
    )
  )

  # Define the output file path with direction and target gene in the file name
  output_file <- paste0(output_folder, deparse(substitute(Direction)), "_", target_gene, ".pdf")

  # Create the heatmap
  all_fgsea_results %>%
    filter(pathway %in% final_pathways) %>%
    enrichmap(n_top = Inf,
              set_column = "pathway",
              statistic_column = "NES",
              contrast_column = "contrast",
              padj_column = "padj",
              padj_legend_title = "BH Adjusted\nP-Value",
              heatmap_args = list(
                top_annotation = col_anno,
                column_order = contrast_order,
                show_column_names = FALSE
              ),
              filename = output_file,
              height = 10,
              width = 15)

  return(output_file)
}
# Figure 7G: Igfbp7 from watsc to skm-vl
wat_svl <- readRDS("/Volumes/Rapid/PASS1B_QENIE/bicor_results_sex/bicor_results_sexadj/WAT-SC/bicor_SKM-VL.rds")
svl_t <- TRNSCRPT_SKMVL_NORM_DATA # from MotrpacRatTraining6mo
wat_t <- TRNSCRPT_WATSC_NORM_DATA # from MotrpacRatTraining6mo

wat_svl <- assign_row_and_col_names_to_bicor(wat_svl, svl_t$feature_ID, wat_t$feature_ID)
# Example Usage: wat to svl:
generate_fgsea_heatmap(Direction = wat_svl, target_gene = "ENSRNOG00000002050") #igfbp7
rm(wat_svl) # clear R space

# Figure 7H: Igfbp7 from watsc to liver
wat_liv <- readRDS("/Volumes/Rapid/PASS1B_QENIE/bicor_results_sex/bicor_results_sexadj/WAT-SC/bicor_LIVER.rds")
liv_t <- TRNSCRPT_LIVER_NORM_DATA # from MotrpacRatTraining6mo

wat_liv <- assign_row_and_col_names_to_bicor(wat_liv, liv_t$feature_ID, wat_t$feature_ID)

generate_fgsea_heatmap(Direction = wat_liv, target_gene = "ENSRNOG00000002050") #igfbp7
rm(wat_liv) # clear R space
