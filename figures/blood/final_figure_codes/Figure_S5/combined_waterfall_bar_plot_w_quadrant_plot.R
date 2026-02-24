generate_combined_quadrant_waterfall_plot <- function(
    repo_local_dir,
    selected_tissue = "blood",
    selected_ome = "prot-ol",
    timepoint = "post_10_min",
    threshold = 0.05,
    x_axis_label = "Protein Assay",
    output_prefix = "Blood_proteomics",
    figure_dir,
    output_format = "pdf"
) {
  # Load packages
  required_libs <- c(
    "dplyr", "ggplot2", "tidyr", "glue", "cowplot", "ggrepel", "stringr",
    "MotrpacHumanPreSuspension", "MotrpacHumanPreSuspensionData"
  )
  missing_libs <- required_libs[!sapply(required_libs, requireNamespace, quietly = TRUE)]
  if (length(missing_libs) > 0) {
    stop(paste0("Missing packages: ", paste(missing_libs, collapse = ", "),
                ". Please install them first."))
  }
  suppressPackageStartupMessages(lapply(required_libs, library, character.only = TRUE))
  
  if (!dir.exists(figure_dir)) {
    dir.create(figure_dir, recursive = TRUE)
    message(glue::glue("📂 Created Plot output directory: {figure_dir}"))
  }
  
  color_exercise <- MotrpacHumanPreSuspension::HUMAN_EXERCISE_GROUP_COLORS
  color_resist <- color_exercise[["ADUResist"]]
  color_endur  <- color_exercise[["ADUEndur"]]
  color_control <- color_exercise[["ADUControl"]]
  
  create_quadrant_plot_clean <- function(end_estimates, res_estimates, 
                                         end_p_adj, res_p_adj, compound_ids) {
    
    plot_data <- data.frame(
      End_estimates = end_estimates, 
      Res_estimates = res_estimates, 
      End_p_adj = end_p_adj, 
      Res_p_adj = res_p_adj, 
      Compound_ID = compound_ids
    )
    
    # Assign significance category
    plot_data$Significance <- ifelse(
      plot_data$End_p_adj < 0.05 & plot_data$Res_p_adj < 0.05, "Both Significant",
      ifelse(plot_data$End_p_adj > 0.05 & plot_data$Res_p_adj > 0.05, "Not Significant",
             ifelse(plot_data$End_p_adj < 0.05 & plot_data$Res_p_adj > 0.05, "EE Significant",
                    ifelse(plot_data$End_p_adj > 0.05 & plot_data$Res_p_adj < 0.05, "RE Significant", NA)
             )
      )
    )
    
    top_5_ee_pos <- plot_data$Compound_ID[order(plot_data$End_estimates, decreasing = TRUE)][1:2]
    top_5_ee_neg <- plot_data$Compound_ID[order(plot_data$End_estimates, decreasing = FALSE)][1:2]
    top_5_re_pos <- plot_data$Compound_ID[order(plot_data$Res_estimates, decreasing = TRUE)][1:2]
    top_5_re_neg <- plot_data$Compound_ID[order(plot_data$Res_estimates, decreasing = FALSE)][1:2]
    
    top_ids <- unique(c(top_5_ee_pos, top_5_ee_neg, top_5_re_pos, top_5_re_neg))
    plot_data_filtered <- plot_data[!(plot_data$Compound_ID %in% top_ids), ]
    
    x_limits <- c(-max(abs(plot_data_filtered$End_estimates)), max(abs(plot_data_filtered$End_estimates)))
    y_limits <- c(-max(abs(plot_data_filtered$Res_estimates)), max(abs(plot_data_filtered$Res_estimates)))
    
    x_margin <- 0.1 * diff(x_limits)
    y_margin <- 0.1 * diff(y_limits)
    
    adjusted_x_limits <- c(x_limits[1] - x_margin, x_limits[2] + x_margin)
    adjusted_y_limits <- c(y_limits[1] - y_margin, y_limits[2] + y_margin)
    
    plot_data = plot_data_filtered
    
    # Quadrant scatter plot
    p <- ggplot(plot_data, aes(x = End_estimates, y = Res_estimates, color = Significance)) +
      geom_point() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      scale_color_manual(
        name = NULL,
        values = c(
          "Both Significant" = "black",
          "EE Significant" = color_endur,
          "RE Significant" = color_resist,
          "Not Significant" = "grey"
        ),
        breaks = c("Both Significant", "EE Significant", "RE Significant", "Not Significant")
      ) +
      guides(color = guide_legend(nrow = 1)) + 
      xlim(adjusted_x_limits) +
      ylim(adjusted_y_limits) +
      theme_minimal() +
      theme(legend.position = "bottom",
            plot.title = element_blank(),
            legend.text = element_text(size = 6.5),
            legend.title = element_text(size = 8)) +
      # guides(color = guide_legend(title = "Significance"))+
      labs(
        x = "EE logFC",
        y = "RE logFC"
      )
    
    return(p)
  }
  
  
  create_waterfall_plot_clean <- function(df_resistance, df_endurance, timepoint, threshold, x_axis_label, timepoint_label = NULL) {
    
    df_resistance <- df_resistance[df_resistance$Timepoint == timepoint, ]
    df_endurance  <- df_endurance[df_endurance$Timepoint == timepoint, ]
    
    significant_resistance <- df_resistance %>% filter(adj_p_value < threshold)
    significant_endurance  <- df_endurance %>% filter(adj_p_value < threshold)
    
    list_a <- setdiff(significant_resistance$feature_id, significant_endurance$feature_id)
    list_b <- setdiff(significant_endurance$feature_id, significant_resistance$feature_id)
    combined_list <- union(list_a, list_b)
    
    plot_data_resistance <- df_resistance %>%
      filter(feature_id %in% combined_list) %>%
      mutate(Group = "RE") %>%
      dplyr::select(feature_id, gene_symbol, logFC, Group, sig) %>%
      arrange(desc(logFC))
    
    plot_data_endurance <- df_endurance %>%
      filter(feature_id %in% combined_list) %>%
      mutate(Group = "EE") %>%
      dplyr::select(feature_id, gene_symbol, logFC, Group, sig) %>%
      arrange(logFC)
    
    plot_data <- bind_rows(plot_data_resistance, plot_data_endurance)
    
    assays_in_list_a <- plot_data_resistance$gene_symbol[plot_data_resistance$feature_id %in% list_a]
    assays_in_list_b <- plot_data_endurance$gene_symbol[plot_data_endurance$feature_id %in% list_b]
    
    plot_data <- plot_data %>%
      mutate(gene_symbol = factor(gene_symbol, levels = unique(c(assays_in_list_a, assays_in_list_b))))
    
    assay_colors <- ifelse(levels(plot_data$gene_symbol) %in% assays_in_list_a, color_resist, color_endur)
    
    max_y <- max(plot_data$logFC, na.rm = TRUE) * 1.2
    mid_x <- length(list_a) + 0.5
    
    arrow_y <- min(plot_data$logFC, na.rm = TRUE) - (0.1 * max_y)
    label_offset <- 0.08 * max_y   # adjust spacing above arrow
    significant_assays <- plot_data %>%
      group_by(gene_symbol) %>%
      summarise(is_significant = any(sig == TRUE), .groups = "drop") %>%
      filter(is_significant) %>%
      pull(gene_symbol)
    
    # Barplot
    p <- ggplot(plot_data, aes(x = gene_symbol, y = logFC, fill = Group)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
      scale_fill_manual(values = c("EE" = color_endur, "RE" = color_resist)) +
      labs(y = "logFC", x = x_axis_label, fill = "Group") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 16,
                                       face = "bold", color = assay_colors),
            legend.position = "bottom",
            legend.title = element_text(face = "bold"),
            legend.text = element_text(size = 10)) +
      geom_text(
        data = subset(plot_data, gene_symbol %in% significant_assays),
        aes(label = "*", y = max_y * 0.9),
        position = position_dodge(width = 0.8),
        size = 6, fontface = "bold", color = "black"
      ) +
      
      # Add RE and EE text + arrows BELOW the plot
      annotate("text", x = length(list_a) / 2, y = arrow_y+ label_offset, label = "RE", size = 8, fontface = "bold", color = color_resist) +
      annotate("text", x = length(list_a) + length(list_b) / 2 + 0.5, y = arrow_y+ label_offset, label = "EE", size = 8, fontface = "bold", color = color_endur) +
      
      # annotate("segment", x = 0.5, xend = mid_x - 0.05, y = arrow_y , yend = arrow_y,
      #          color = color_resist, arrow = arrow(length = unit(0.2, "cm"), ends = "first")) +
      # annotate("segment", x = length(list_a) + 0.55, xend = length(combined_list) + 0.5, y = arrow_y, yend = arrow_y,
      #          color = color_endur, arrow = arrow(length = unit(0.2, "cm"), ends = "last")) 
      
      annotate("segment",
               x = 0.5, xend = mid_x - 0.05,
               y = arrow_y, yend = arrow_y,
               color = color_resist, size = 0.8) +
      annotate("segment",
               x = length(list_a) + 0.55, xend = length(combined_list) + 0.5,
               y = arrow_y, yend = arrow_y,
               color = color_endur, size = 0.8)
    
    
    
    # 🔹 Add timepoint label at top-left inside bar plot
    if (!is.null(timepoint_label)) {
      p <- p + annotate("text", x = 0.5, y = max_y, 
                        label = paste("Time point:", timepoint_label),
                        hjust = 0, vjust = 1, fontface = "bold", size = 5)
    }
    
    return(p)
  }
  
  # ---- Master Function: Combine Waterfall + Quadrant ----
  create_combined_plot <- function(df_resistance, df_endurance, data,
                                   timepoint,
                                   threshold = 0.05,
                                   x_axis_label = "Proteins",
                                   timepoint_label = NULL) {
    library(ggplot2)
    library(cowplot)
    library(dplyr)
    
    # Choose label column
    label_col <- if ("assay" %in% names(data)) {
      "assay"
    } else if ("gene_symbol" %in% names(data)) {
      "gene_symbol"
    } else {
      "feature_id"
    }
    
    # Column names based on timepoint
    end_estimate_col <- paste0("Endur_logFC_", timepoint)
    res_estimate_col <- paste0("Resist_logFC_", timepoint)
    end_p_adj_col    <- paste0("Endur_adj_p_value_", timepoint)
    res_p_adj_col    <- paste0("Resist_adj_p_value_", timepoint)
    
    # Pretty timepoint label if not supplied
    if (is.null(timepoint_label)) {
      timepoint_label <- case_when(
        timepoint == "post_10_min"       ~ "Post 10 Min",
        timepoint == "post_15_30_45_min" ~ "Post 30 Min",
        timepoint == "post_3.5_4_hr"     ~ "Post 3.5 Hour",
        timepoint == "post_24_hr"        ~ "Post 24 Hour",
        TRUE ~ timepoint
      )
    }
    
    # Extract vectors from master data
    end_estimates <- data[[end_estimate_col]]
    res_estimates <- data[[res_estimate_col]]
    end_p_adj     <- data[[end_p_adj_col]]
    res_p_adj     <- data[[res_p_adj_col]]
    compound_ids  <- data[[label_col]]
    
    # Quadrant plot (clean, no title)
    quad_clean <- create_quadrant_plot_clean(
      end_estimates, res_estimates,
      end_p_adj, res_p_adj,
      compound_ids
    ) + theme(plot.title = element_blank())
    
    # Waterfall plot (clean, no EE/RE arrows, no title)
    waterfall_clean <- create_waterfall_plot_clean(
      df_resistance, df_endurance,
      timepoint, threshold,
      x_axis_label = x_axis_label,
      timepoint_label = timepoint_label
    ) + ggtitle("")
    
    if (timepoint == "post_3.5_4_hr") {
      inset_x <- 0.70
      inset_y <- 0.78
      inset_w <- 0.25
      inset_h <- 0.22
    } else {
      inset_x <- 0.65
      inset_y <- 0.55
      inset_w <- 0.35
      inset_h <- 0.45
    }
    
    # Combine plots using the dynamic layout
    final_plot <- ggdraw() +
      draw_plot(waterfall_clean) +
      draw_plot(quad_clean, x = inset_x, y = inset_y, width = inset_w, height = inset_h)
    
    return(final_plot)
  }
  
  
  # ---- Load differential analysis ----
  df_da_full <- load_differential_analysis(
    repo_local_dir = repo_local_dir,
    selected_omes = "all",
    selected_tissues = selected_tissue,
    single_matrix = TRUE,
    combine_with_featgene = TRUE,
    epigen = FALSE
  ) %>%
    dplyr::mutate(
      assay = dplyr::case_when(
        stringr::str_detect(assay, "metab") ~ "metab",
        TRUE ~ assay
      ),
      gene_symbol = dplyr::if_else(
        assay == "metab",
        feature_id,
        gene_symbol
      )
    )
  
  df_da_full_resist <- df_da_full %>%
    filter(contrast_type == "exercise_with_controls", randomGroupCode == "ADUResist")
  
  df_da_full_endur <- df_da_full %>%
    filter(contrast_type == "exercise_with_controls", randomGroupCode == "ADUEndur")
  
  # ---- Prepare wide format ----
  df_da_full_resist_wide <- df_da_full_resist %>%
    dplyr::select(assay, feature_id, logFC, p_value, adj_p_value, Timepoint, tissue) %>%
    pivot_wider(names_from = Timepoint, values_from = c(logFC, p_value, adj_p_value))
  
  df_da_full_endur_wide <- df_da_full_endur %>%
    dplyr::select(assay, feature_id, logFC, p_value, adj_p_value, Timepoint, tissue) %>%
    pivot_wider(names_from = Timepoint, values_from = c(logFC, p_value, adj_p_value))
  
  colnames(df_da_full_resist_wide)[4:ncol(df_da_full_resist_wide)] <-
    paste("Resist", colnames(df_da_full_resist_wide)[4:ncol(df_da_full_resist_wide)], sep = "_")
  colnames(df_da_full_endur_wide)[4:ncol(df_da_full_endur_wide)] <-
    paste("Endur", colnames(df_da_full_endur_wide)[4:ncol(df_da_full_endur_wide)], sep = "_")
  
  merged_endur_resist_da <- merge(
    df_da_full_resist_wide, df_da_full_endur_wide,
    by = c("assay", "feature_id", "tissue")
  )
  
  df_unique_blood <- merged_endur_resist_da[merged_endur_resist_da$tissue == selected_tissue, ]
  df_unique_blood_prot <- df_unique_blood[df_unique_blood$assay == selected_ome, ] %>%
    dplyr::select(-assay, -tissue)
  
  # ---- Mapping ----
  mapping <- MotrpacHumanPreSuspensionData::HUMAN_FEATURE_TO_GENE %>%
    filter(assay == selected_ome) %>%
    dplyr::select(feature_id, uniprot, gene_symbol)
  
  df_unique_blood_prot <- merge(mapping, df_unique_blood_prot, by = "feature_id")
  
  # ---- Endur vs Resist p-values ----
  df_da_blood_end_vs_res_prot <- df_da_full %>%
    filter(assay == selected_ome, contrast_type == "Endur_vs_Resist", Timepoint != "pre_exercise")
  
  df_da_blood_end_vs_res_prot_p_adj <- df_da_blood_end_vs_res_prot %>%
    dplyr::select(feature_id, uniprot, gene_symbol, Timepoint, contains("adj_p_value")) %>%
    mutate(sig = adj_p_value <= 0.05) %>%
    dplyr::select(-adj_p_value)
  
  # ---- Merge for plotting ----
  df_da_blood_prot <- df_da_full %>%
    filter(assay == selected_ome, contrast_type == "exercise_with_controls") %>%
    dplyr::select(feature_id, uniprot, gene_symbol, logFC, p_value, adj_p_value, randomGroupCode, Timepoint)
  
  df_da_blood_prot_plot <- merge(
    df_da_blood_prot, df_da_blood_end_vs_res_prot_p_adj,
    by = c("feature_id", "uniprot", "gene_symbol", "Timepoint")
  )
  
  df_da_blood_prot_plot_endur <- df_da_blood_prot_plot %>%
    filter(randomGroupCode == "ADUEndur")
  
  df_da_blood_prot_plot_resist <- df_da_blood_prot_plot %>%
    filter(randomGroupCode == "ADUResist")
  
  # ---- Create combined plot ----
  final_plot <- create_combined_plot(
    df_da_blood_prot_plot_resist,
    df_da_blood_prot_plot_endur,
    df_unique_blood_prot,
    timepoint = timepoint,
    threshold = threshold,
    x_axis_label = x_axis_label,
    timepoint_label = NULL
  )
  
  # ---- Save plot ----
  output_file <- file.path(
    figure_dir,
    paste0(output_prefix, "_", selected_tissue, "_", selected_ome, "_", timepoint, ".", output_format)
  )
  
  ggsave(output_file, plot = final_plot, width = 16, height = 12, units = "in", dpi = 2500)
  message(glue::glue("✅ Saved combined plot to: {output_file}"))
  
  return(invisible(final_plot))
}