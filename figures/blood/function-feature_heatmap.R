#' @title Create a heatmap of select features
#'
#' @description Create a heatmap for user-specified features for a given tissue
#'   and ome combination.
#'
#' @param DA_list nested list of DA results. The names at the top level should
#'   be tissues, while the names at the second level should be omes. Each
#'   element of the nested list is a \code{data.table}.
#' @param feature_ids character or \code{NULL}; vector of feature IDs that will
#'   appear in the heatmap. Must be a subset of
#'   \code{HUMAN_FEATURE_TO_ID[["feature_id"]]}.
#' @param set_id character of \code{NULL}; if not \code{NULL}, the feature IDs
#'   in the set will be used to filter rows for the heatmap, but those IDs will
#'   not be used for the heatmap row labels.
#' @param platforms character or \code{NULL}; the platform(s) used to filter
#'   metabolites. If \code{NULL} (default), metabolites in \code{feature_ids}
#'   will be selected from all available platforms. If metabolites appear in
#'   more than one platform, the platform will appear before the metabolite name
#'   in the row names of the heatmap.
#' @param selected_tissue character; the tissue that will be used to create the
#'   heatmap.
#' @param selected_ome character; the ome that will be used to create the
#'   heatmap.
#' @param filename character; optional file name used to save the heatmap. If
#'   provided, the heatmap will not be drawn.
#' @param post_min numeric; for experiments, timepoints could be either 15, 30, or
#'   45 minutes depending on the analysis, this argument allows users to specify
#'   which of those 3 values to use, by default the value is NULL and will provide
#'   a generic label of post 15/30/45 min
#' @param post_hr numeric; for experiments, timepoints could be either 3.5 or 4 hours
#'   depending on the analysis, this argument allows users to specify
#'   which of those 2 values to use, by default the value is NULL and will provide
#'   a generic label of post 3.5/4 hr
#' @param full_modality_names logical; if TRUE the modality values are set as Endurance
#'   Exercise and Resistance Exercise but if FALSE the modality values are set as EE
#'   and RE (by default the value is FALSE)
#'
#' @returns Nothing. A heatmap is drawn or saved to a file if \code{filename} is
#'   provided.
#'
#' @export feature_heatmap
#'
#' @author Tyler Sagendorf, Damon Leach
#'
#' @import ComplexHeatmap
#' @importFrom dplyr %>% filter mutate bind_rows
#' @importFrom grDevices dev.off
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr pivot_wider


## Example
# feature_heatmap(set_id = "11725",
#                 DA_list = DA_list,
#                 contrast_type = "exercise_with_controls",
#                 selected_tissue = "muscle",
#                 selected_ome = "prot-ph",
#                 filename = "sandbox/test_feature_heatmap.pdf")

library(MotrpacHumanPreSuspension)
library(ComplexHeatmap)
library(gridtext)
library(dplyr)
library(grDevices)
library(tibble)
library(tidyr)
library(TMSig) # enrichmap

feature_heatmap <- function(feature_ids = NULL,
                            set_id = NULL, # if provided, feature_ids is ignored
                            platforms = NULL, # metabolomics only
                            selected_tissue = c("adipose", "blood", "muscle"),
                            selected_ome = c("transcript-rna-seq", "prot-pr",
                                             "prot-ph", "metab", "prot-ol"),
                            contrast_type = "exercise_with_controls",
                            column_title = "",
                            filename,
                            post_min = NULL,
                            post_hr = NULL,
                            full_modality_names = FALSE)
{
  DA_list <- load_differential_analysis(combine_with_featgene = TRUE)

  # check that post_min is formatted correctly
  if(!is.null(post_min)){
    if(length(post_min) > 1){
      stop("If not NULL, post_min can only have length 1")
    }
    if(!is.numeric(post_min)){
      stop("If not NULL, post_min must be numeric")
    }
    if(!post_min %in% c(15,30,45)){
      stop("If not NULL, post_min can only take on the numeric values of 15, 30, or 45")
    }
  }

  # check that post_hr is formatted correctly
  if(!is.null(post_hr)){
    if(length(post_hr) > 1){
      stop("If not NULL, post_hr can only have length 1")
    }
    if(!is.numeric(post_hr)){
      stop("If not NULL, post_hr must be numeric")
    }
    if(!post_hr %in% c(3.5,4)){
      stop("If not NULL, post_hr can only take on the numeric values of 3.5 or 4")
    }
  }

  # check that full_modality_names is formatted correctly
  if(length(full_modality_names) > 1){
    stop("full_modality_names can only have length 1")
  }
  if(!is.logical(full_modality_names)){
    stop("full_modality_names must be a logical argument")
  }

  # If set_id is provided, select the genes/features from the set
  if (!is.null(set_id)) {
    if (!is.character(set_id) || length(set_id) != 1L) {
      stop("`set_id` must be a character string specifying a set ID ",
           "from the `SET_TO_ID` object.")
    }
    pathway <- SET_TO_ID %>%
      filter(set_id == !!set_id) %>%
      pull(set)

    pathways <- unlist(structure(MOLECULAR_SIGNATURES, names = NULL),
                       recursive = FALSE)
    feature_ids <- pathways[[pathway]]
  }

  if (!is.vector(feature_ids, mode = "character")) {
    stop("If `set_id` is not provided, `feature_ids` must be a character ",
         "vector of feature IDs that will appear in the heatmap.")
  }

  selected_tissue <- match.arg(selected_tissue,
                               choices = c("adipose", "blood", "muscle"),
                               several.ok = TRUE)

  selected_ome <- match.arg(selected_ome,
                            choices = c("transcript-rna-seq", "prot-pr",
                                        "prot-ph", "metab", "prot-ol"))

  contrast_type <- match.arg(contrast_type,
                             choices = c("exercise_with_controls",
                                         "exercise_no_controls",
                                         "Endur_vs_Resist",
                                         "baseline",
                                         "control_only"))

  ## Prepare DA results ----
  # DA_res <- DA_list[[selected_tissue]][[selected_ome]] %>%
  #   filter(contrast_type == !!contrast_type) %>%
  #   # Keeping these columns as factors will cause problems with filtering and
  #   # modifying the feature_id column.
  #   mutate(across(.cols = any_of(c("feature_id", "gene_symbol",
  #                                  "platform", "flanking_sequence")),
  #                 .fns = as.character)) %>%
  #   droplevels.data.frame()

  DA_res_list <- purrr::map(selected_tissue, function(st){
    filt_dat = DA_list[[st]][[selected_ome]] %>%
      filter(contrast_type == !!contrast_type) %>%
      # Keeping these columns as factors will cause problems with filtering and
      # modifying the feature_id column.
      mutate(across(.cols = any_of(c("feature_id", "gene_symbol",
                                     "platform", "flanking_sequence")),
                    .fns = as.character)) %>%
      droplevels.data.frame()
    filt_dat
  })
  names(DA_res_list) <- selected_tissue

  DA_res = dplyr::bind_rows(DA_res_list)

  if (is.null(set_id)) {
    # Filter before modifying IDs
    DA_res <- filter(DA_res, feature_id %in% feature_ids)
  }

  if (!nrow(DA_res)) {
    stop("`feature_ids` do not match feature IDs in the ",
         sprintf("%s %s differential analysis results.",
                 selected_tissue, selected_ome))
  }

  if (selected_ome %in% c("prot-pr", "transcript-rna-seq", "prot-ol")) {
    DA_res <- DA_res %>%
      mutate(feature_id = ifelse(!is.na(gene_symbol),
                                 gene_symbol,
                                 feature_id))
  } else if (selected_ome == "prot-ph") {
    DA_res <- DA_res %>%
      mutate(feature_id = gsub("[sty]", ";", feature_id),
             feature_id = sub("(.*);$", "\\1", feature_id),
             feature_id = sub(".*_", "", feature_id),
             feature_id = paste(gene_symbol, feature_id))

  } else  if (selected_ome == "metab" && !is.null(platforms)) {
    platforms <- match.arg(arg = platforms,
                           choices = paste0("metab-",
                                            c("t-amines", "t-conv", "t-imm-crt",
                                              "t-imm-glc", "t-imm-ins",
                                              "t-oxylipneg", "t-tca",
                                              "u-hilicpos", "u-ionpneg",
                                              "u-lrpneg", "u-lrppos",
                                              "u-rpneg", "u-rppos")),
                           several.ok = TRUE)

    DA_res <- DA_res %>%
      filter(platform %in% platforms)

    if (!nrow(DA_res)) {
      stop("Metabolites in `feature_ids` not found in the provided platforms.")
    }
  }

  if (!is.null(set_id)) {
    if (selected_ome %in% c("prot-pr", "transcript-rna-seq")) {
      DA_res <- filter(DA_res, feature_id %in% feature_ids)
    } else if (selected_ome == "prot-ph") {
      DA_res <- filter(DA_res, flanking_sequence %in% feature_ids)
    } else if (selected_ome == "metab") {
      DA_res <- filter(DA_res, feature_id %in% feature_ids)
    }
  }

  if (nrow(DA_res) == 0L) {
    stop("The combination of `set_id`, `selected_tissue`, and `selected_ome` ",
         "is not valid.")
  }

  # if (selected_ome == "metab") {
  #   DA_res <- DA_res %>%
  #     # Include platform in the row names of the heatmap
  #     mutate(feature_id = paste(sub("^metab-", "", platform), feature_id))
  # }

  x <- DA_res %>%
    dplyr::mutate(contrast2 = paste(tissue, contrast),
                  contrast2 = factor(contrast2, levels = unique(contrast2))) %>%
    arrange(contrast2, abs(z.std), feature_id) %>%
    filter(.by = contrast2,
           !duplicated(feature_id)) %>%
    select(feature_id, contrast, contrast2, z.std, adj_p_value,tissue) %>%
    droplevels.data.frame()

  # Better contrast labels
  contrast_df <- .add_contrast_labels() %>%
    filter(contrast %in% levels(x$contrast)) %>%
    droplevels.data.frame()

  contrast_colors_vector <- .contrast_colors()

  if(!is.null(post_min)){
    if(!post_min %in% c(15,30,45)){
      stop("post_min can only take on the values 15, 30, or 45")
    } else {
      post_min_str = paste0("post ", post_min, " min")
      # update level information
      new_levels = levels(contrast_df$contrast_labels)
      new_levels[new_levels == "post 15/30/45 min"] = post_min_str
      # update values in df
      contrast_df$contrast_labels <- as.character(contrast_df$contrast_labels)
      contrast_df$contrast_labels[contrast_df$contrast_labels == "post 15/30/45 min"] = post_min_str
      # reset to factor with new levels
      contrast_df$contrast_labels <- factor(contrast_df$contrast_labels, levels = new_levels)
      # update color vectors
      names(contrast_colors_vector)[names(contrast_colors_vector) == "post 15/30/45 min"] = post_min_str
    }
  }

  if(!is.null(post_hr)){
    if(!post_hr %in% c(3.5,4)){
      stop("post_hr can only take on the values 3.5 or 4")
    } else {
      post_hr_str = paste0("post ", post_hr, " hr")
      # update level information
      new_levels = levels(contrast_df$contrast_labels)
      new_levels[new_levels == "post 3.5/4 hr"] = post_hr_str
      # update values in df
      contrast_df$contrast_labels <- as.character(contrast_df$contrast_labels)
      contrast_df$contrast_labels[contrast_df$contrast_labels == "post 3.5/4 hr"] = post_hr_str
      # reset to factor with new levels
      contrast_df$contrast_labels <- factor(contrast_df$contrast_labels, levels = new_levels)
      # update color vectors
      names(contrast_colors_vector)[names(contrast_colors_vector) == "post 3.5/4 hr"] = post_hr_str

    }
  }

  ## Create heatmap ----
  column_df <- distinct(x,tissue,contrast) %>%
    mutate(tissue = factor(tissue, levels = selected_tissue)) %>%
    arrange(contrast,tissue) %>%
    left_join(contrast_df,
              by = "contrast") %>%
    rename(modality = anno_group) %>%
    mutate(modality = factor(modality,
                             levels = c("EE", "RE")))

  if(full_modality_names == TRUE){
    column_df$modality <- as.character(column_df$modality)
    column_df$modality[column_df$modality == "EE"] = "Endurance Exercise"
    column_df$modality[column_df$modality == "RE"] = "Resistance Exercise"
    column_df$modality <- factor(column_df$modality, levels = c("Endurance Exercise","Resistance Exercise"))
  }

  # order the annotations by tissue -> modality -> contrast_labels
  anno_df <- dplyr::select(column_df,
                           Tissue = tissue,
                           Modality = modality,
                           Timepoint = contrast_labels) %>%
    dplyr::arrange(Tissue,Modality,Timepoint)

  if(full_modality_names == TRUE){
    anno_col <- list(
      "Tissue" = HUMAN_TISSUE_COLORS[selected_tissue],
      "Modality" = setNames(c("#d95f02", "#1b9e77"),
                            c("Endurance Exercise", "Resistance Exercise")),
      "Timepoint" = contrast_colors_vector[levels(anno_df$Timepoint)]
    )
  } else {
    anno_col <- list(
      "Tissue" = HUMAN_TISSUE_COLORS[selected_tissue],
      "Modality" = setNames(c("#d95f02", "#1b9e77"),
                            c("EE", "RE")),
      "Timepoint" = contrast_colors_vector[levels(anno_df$Timepoint)]
    )
  }


  # contrasts from different tissues must be treated as distinct
  column_df <- column_df %>%
    dplyr::mutate(contrast2 = paste(tissue, contrast),
                  contrast2 = factor(contrast2, levels = unique(contrast2)))
  # If there is a single tissue, remove the tissue annotation
  if (length(selected_tissue) == 1L) {
    anno_df$Tissue <- NULL
    anno_col["Tissue"] <- NULL
  }

  show_column_names <- FALSE

  if (contrast_type == "baseline") {
    anno_df$Timepoint <- NULL
    anno_col["Timepoint"] <- NULL

    show_column_names <- TRUE
  }

  if (!contrast_type %in% c("exercise_with_controls", "exercise_no_controls")) {
    anno_df$Modality <- NULL
    anno_col["Modality"] <- NULL
  }

  if (length(anno_col)) {
    if (!is.null(anno_df[["Tissue"]]) && !is.null(anno_df[["Modality"]])) {
      column_split <- anno_df %>%
        mutate(column_split = paste(Tissue, Modality),
               row_order = 1:n()) %>%
        arrange(Tissue, Modality) %>%
        mutate(column_split = factor(column_split,
                                     levels = unique(column_split))) %>%
        arrange(row_order) %>%
        pull(column_split)
    } else if (!is.null(anno_df[["Tissue"]])) {
      column_split <- anno_df[["Tissue"]]
    } else if (!is.null(anno_df[["Modality"]])) {
      column_split <- anno_df[["Modality"]]
    }

    top_annotation <- HeatmapAnnotation(
      df = anno_df,
      col = anno_col,
      which = "column",
      border = TRUE,
      gap = unit(2, "pt"),
      annotation_name_gp = gpar(fontsize = 0.9 * 14),
      annotation_legend_param = list(
        border = TRUE,
        title_gp = gpar(fontsize = 0.9 * unit(14, "pt"),
                        fontface = "bold"),
        labels_gp = gpar(fontsize = 0.9 * unit(14, "pt"))
      )
    )
  } else {
    top_annotation <- column_split <- NULL
  }

  n_features <- length(unique(x[["feature_id"]]))

  # Dynamic heatmap height. Each cell of the heatmap is 14 pt, by default
  height <- convertUnit(n_features * unit(14, "pt"), "in")
  height <- max(as.numeric(height), 4.5) + 1.5 + show_column_names

  # Dynamic heatmap width
  row_label_width <- ComplexHeatmap::max_text_width(
    text = x[["feature_id"]],
    gp = gpar(fontsize = 0.9 * 14)
  )
  row_label_width <- convertUnit(row_label_width, "in")

  width <- convertUnit(nlevels(x[["contrast"]]) * unit(14, "pt"), "in")
  width_extra = ifelse(length(selected_tissue) > 1, 5, 3)
  width <- as.numeric(width + row_label_width) + width_extra

  draw_args <- list(
    annotation_legend_list = NULL,
    merge_legends = TRUE
  )

  if (n_features <= 10L) {
    draw_args[["padding"]] <- unit(c(80, 0, 0, 0), "pt")
  }

  extended_range <- TMSig::extendRangeNum(x[["z.std"]], nearest = 0.1)

  if (all(extended_range <= 0)) {
    breaks <- c(extended_range[1], 0)
  } else if (all(extended_range >= 0)) {
    breaks <- c(0, extended_range[2])
  } else {
    max_r <- max(abs(extended_range))
    breaks <- c(extended_range[1], 0, extended_range[2])
  }

  # if we are comparing multiple tissues, likely to have
  # missing values which likely will cause issues with clustering
  clust_row_info <- ifelse(length(selected_tissue) > 1, FALSE, TRUE)

  TMSig::enrichmap(
    x = x,
    n_top = Inf,
    set_column = "feature_id",
    statistic_column = "z.std",
    contrast_column = "contrast2",
    padj_column = "adj_p_value",
    plot_sig_only = FALSE,
    #filename = filename,
    height = height,
    width = width,
    heatmap_color_fun = .feature_color_function,
    heatmap_args = list(
      layer_fun = .feature_layer_fun,
      cluster_rows = clust_row_info,
      column_split = column_split,
      column_labels = contrast_df$specific,
      show_column_names = show_column_names,
      column_names_side = "top",
      column_title = gt_render(column_title,
                               padding = unit(c(0, 0, 0, 0.8), "in")),
      column_title_gp = gpar(fontsize = 12),
      top_annotation = top_annotation,
      na_col = "grey80",
      heatmap_legend_param = list(
        title = "Z-Score",
        at = breaks,
        labels = breaks
      )
    ),
    draw_args = draw_args
  )
}



## Helper functions ------------------------------------------------------------

.add_contrast_labels <- function() {
  contrast_labels <- c(
    # exercise_with_controls and exercise_no_controls
    rep(c("during 20 min", "during 40 min", "post 10 min", "post 15/30/45 min",
          "post 3.5/4 hr", "post 24 hr", "post 10 min", "post 15/30/45 min",
          "post 3.5/4 hr", "post 24 hr"), times = 2L),
    # Endur vs. Resist
    c("post 10 min", "post 15/30/45 min", "post 3.5/4 hr", "post 24 hr"),
    # baseline
    c("Endur - Resist", "Endur - Control", "Resist - Control"),
    # control_only
    c("during 20 min", "during 40 min", "post 10 min", "post 15/30/45 min",
      "post 3.5/4 hr", "post 24 hr")
  )

  contrast_labels <- factor(
    x = contrast_labels,
    levels = c("Endur - Resist", "Endur - Control", "Resist - Control",
               "during 20 min", "during 40 min", "post 10 min",
               "post 15/30/45 min", "post 3.5/4 hr", "post 24 hr")
  )

  anno_group <- c(
    # exercise_with_controls and exercise_no_controls
    rep(rep(c("EE", "RE"), c(6L, 4L)), times = 2L),
    # Endur_vs_Resist and baseline, and control_only
    rep(NA, 13L)
  )

  out <- cbind(CONTRAST_CONVERTER, contrast_labels, anno_group)

  return(out)
}

.contrast_colors <- function() {
  structure(
    c("#fde725",
      "#bad071",
      "#d1bbd7",
      "#ae76a3",
      "#882e72",
      "#61194f"),
    names = c("during 20 min",
              "during 40 min",
              "post 10 min",
              "post 15/30/45 min",
              "post 3.5/4 hr",
              "post 24 hr")
  )
}

.feature_layer_fun <- function(j, i, x, y, w, h, f) {
  grid.rect(x = x, y = y, width = w, height = h,
            gp = gpar(col = heatmap_args[["rect_gp"]][["col"]],
                      fill = f)
  )

  grid.text("*",
            x = x, y = y,
            #r = pindex(dmat, i, j) / 2 * cell_size,
            # Significant bubbles get a black outline to separate from padj_fill
            gp = gpar(col = ifelse(pindex(padj_mat, i, j) < padj_cutoff,
                                   "black", NA))
  )
}


.feature_color_function <- function(statistics,
                                    colors = c("#3366ff", "darkred"))
{
  r <- range(statistics, na.rm = TRUE)

  # Extend range of values out to the nearest tenth
  extended_range <- TMSig::extendRangeNum(r, nearest = 0.1)

  max_r <- max(abs(extended_range))

  if (all(r >= 0)) {
    breaks <- c(0, +1) * max_r
    colors <- c("white", colors[2])
  } else if (all(r <= 0)) {
    breaks <- c(-1, 0) * max_r
    colors <- c(colors[1], "white")
  } else {
    breaks <- c(-1, 0, 1) * max_r
    colors <- c(colors[1], "white", colors[2])
  }

  return(list(breaks = breaks, colors = colors))
}


