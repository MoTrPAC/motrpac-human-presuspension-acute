#!/usr/bin/env Rscript
# Figure 2 — Differential abundance landscape
#
# Panels:  FIG2A   percent of features differentially abundant, ome x tissue x
#                  exercise group                            (also writes ST2c)
#          FIG2B   differentially abundant features shared across tissues
#                                                            (also writes ST2e)
#          FIG2Bii enriched pathways among the features differential in all
#                  three tissues
#          FIG2C   features differentially abundant in all three tissues
#          FIG2D   MYC and AREG transcript trajectories
#          FIG2E   phosphosites differentially abundant in adipose and muscle
#          FIG2F   pathway enrichment across the tissue-membership sets
# Tables:  ST2c    number of DA features                     (written by FIG2A)
#          ST2e    cross-tissue ORA                          (written by FIG2B)
#
# Needs consortium data access. Six of the seven panels read one
# load_differential_analysis(epigen = TRUE) pass; the ATAC and methylCap DA
# tables are not in the data package and are downloaded into EPIGEN_QC_DIR on
# first use, about 7.7 GB, then reused.
#
#   Rscript figures/landscape/FIG2.R          every panel
#   Rscript figures/landscape/FIG2.R FIG2D    one panel

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(ComplexUpset)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(MotrpacHumanPreSuspensionAnalysis)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(here, "..", "lib", "table_export.R"))
source(file.path(here, "..", "lib", "supp_table_helpers.R"))
source(file.path(here, "..", "lib", "highlights.R"))
source(file.path(here, "..", "lib", "da_overlap_helpers.R"))
source(file.path(here, "..", "lib", "single_feature_helpers.R"))
source(file.path(here, "..", "single_feature_plots.R"))

# ---- shared ----------------------------------------------------------------

# Set order down the UpSet matrix. Stated rather than left to sort_sets, which
# orders by set size and therefore by the data: the legacy passed a fixed vector
# of stripe colours alongside sort_sets = "ascending", so the tissue a stripe
# coloured depended on which tissue happened to have the most hits that release.
# Four of the ported panels carried four different hardcoded orders.
#
# It is the presence matrix's column order too, which FIG2B reads and the other
# panels do not: they index the matrix by name.
TISSUE_SETS <- c("Muscle", "Blood", "Adipose")

# One load_differential_analysis(epigen = TRUE) pass behind six of the seven
# panels. Lazy and memoised rather than loaded at the top of the script: 7.7 GB
# should not be read to build FIG2D, which is the one panel that does not need
# it.
#
# The epigenomics DA tables are not shipped in the package; they are downloaded
# into EPIGEN_QC_DIR, which the package requires as its cache. The bucket is the
# package default, the staging prefix the current precovid-repro cycle writes:
# ATAC DA is v2.0 there, against v1.2 on the published collection.
differential_analysis <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- cross_tissue_da(epigen = TRUE)
    }
    cache
  }
})

# The FDR < 0.05 subset with the display columns attached, and the feature x
# tissue presence matrix built off it. Five panels start from one or the other
# and none of them may re-derive it: reading the DA tables twice in one run is
# what this script exists to stop.
significant_features <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- significant_da(differential_analysis())
    }
    cache
  }
})

tissue_presence_matrix <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- feature_presence_matrix(significant_features(), by = "Tissue",
                                        sets = TISSUE_SETS)
    }
    cache
  }
})

# The over-representation analysis of the six membership sets. FIG2B writes it
# as ST2e, FIG2Bii draws the bar plot of the three-tissue set and FIG2F draws
# the heatmap of all six: three readings of one computation. One script per
# panel would have to recompute it three times, because a panel cannot depend
# on another panel having run; here the three are one run and the
# memoised accessor is what stands in for that.
membership_set_ora <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- cross_tissue_ora(differential_analysis(),
                                 significant_features(),
                                 tissue_presence_matrix())
    }
    cache
  }
})

# ---- FIG2A — percent of features differentially abundant -------------------

# Each assayed feature is reduced to its smallest FDR-adjusted p-value across
# timepoints within an exercise group and tissue, using only the
# exercise-with-controls contrasts. Omes are collapsed to display names and the
# share of features reaching FDR < 0.05 becomes a heatmap cell, labelled with
# the underlying count of significant features.
#
# The panel also writes ST2c, the same DA counts kept per timepoint instead of
# collapsed over them. Both come off the one all_da object. The panel reduces
# each feature to its most significant timepoint, so a feature counts once; the
# table counts it in every timepoint it is significant at, which is why the two
# do not sum to each other.

fig2a <- function() {
  panel_init("FIG2A")
  all_da <- differential_analysis()

  write_st2c(all_da)

  # One row per feature within each group x tissue x ome: the timepoint at which
  # that feature reached its smallest adjusted p-value.
  all_da_sum <- all_da %>%
    filter(contrast_type == "exercise_with_controls") %>%
    mutate(Ome = assay) %>%
    select(Ome, randomGroupCode, tissue, adj_p_value, feature_id) %>%
    group_by(randomGroupCode, tissue, feature_id, Ome) %>%
    slice_min(adj_p_value) %>%
    ungroup()

  # Not lib/da_overlap_helpers.R's OME_LEVELS: this panel is a different figure
  # making a different point, and its row order is stated here.
  ome_levels <- c(
    "Chromatin Accessibility (ATAC)",
    "Transcriptomics",
    "Proteomics",
    "Phosphoproteomics",
    "Metabolomics",
    "Methylation"
  )

  # Percent and count of DA features per assay x tissue x group. The Olink and
  # mass-spec proteomics assays are merged into one Proteomics row after the
  # counts are taken, so each keeps its own denominator.
  all_da_slice <- all_da_sum %>%
    group_by(Ome, tissue, randomGroupCode) %>%
    summarize(
      num_feat = n(),
      percent_sig_fdr = 100 * sum(adj_p_value < 0.05) / n(),
      num_sig_fdr005 = sum(adj_p_value < 0.05),
      .groups = "drop"
    ) %>%
    mutate(
      Group = if_else(randomGroupCode == "ADUEndur", "EE", "RE"),
      Group = factor(Group, levels = c("EE", "RE"))
    ) %>%
    mutate(Ome = if_else((Ome == "prot-ol" | Ome == "prot-pr"), "prot-pr", Ome)) %>%
    mutate(Ome = case_when(
      str_detect(Ome, "metab") ~ "Metabolomics",
      str_detect(Ome, "prot-pr") ~ "Proteomics",
      str_detect(Ome, "prot-ph") ~ "Phosphoproteomics",
      str_detect(Ome, "rna") ~ "Transcriptomics",
      str_detect(Ome, "atac") ~ "Chromatin Accessibility (ATAC)",
      str_detect(Ome, "methyl") ~ "Methylation",
      .default = Ome
    )) %>%
    mutate(Ome = factor(Ome, levels = ome_levels))

  # Tissue order across the columns, within each exercise group. Stated rather
  # than left to whatever order pivot_wider() happens to emit, which follows the
  # sort of the summarise() upstream and would move if a grouping column changed.
  tissue_order <- c("adipose", "blood", "muscle")

  # Omes down the rows, the three EE tissue columns then the three RE ones across.
  # The fill matrix and the label matrix must be laid out the same way or the
  # counts land on the wrong cells, so both come from here.
  ome_by_group_tissue <- function(value_col) {
    wide <- all_da_slice %>%
      select(Ome, tissue, all_of(value_col), Group) %>%
      pivot_wider(values_from = all_of(value_col),
                  names_from = c("Group", "tissue")) %>%
      arrange(Ome) %>%
      column_to_rownames(var = "Ome") %>%
      as.matrix()

    # "<Group>_<tissue>", every group's tissues in tissue_order. A name the pivot
    # did not produce is an error: a silently missing column would shift every
    # count after it onto the wrong cell.
    col_order <- paste(rep(levels(all_da_slice$Group), each = length(tissue_order)),
                       tissue_order, sep = "_")
    missing <- setdiff(col_order, colnames(wide))
    if (length(missing) > 0) {
      stop("no column for ", paste(missing, collapse = ", "),
           " — the group x tissue grid is not complete", call. = FALSE)
    }
    wide[, col_order, drop = FALSE]
  }

  da_feat_perc <- ome_by_group_tissue("percent_sig_fdr")
  da_feat_num <- ome_by_group_tissue("num_sig_fdr005")

  # EE/RE/CON are the endurance, resistance and control arms under the short
  # labels this panel annotates with.
  group_colors <- MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS[
    c("ADUEndur", "ADUResist", "ADUControl")
  ]
  names(group_colors) <- c("EE", "RE", "CON")

  ann_colors <- list(
    Tissue = MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS,
    Group = group_colors,
    Ome = MotrpacHumanPreSuspensionAnalysis::HUMAN_OME_COLORS
  )

  annotation_data <- data.frame(
    Ome = factor(rownames(da_feat_perc), levels = ome_levels)
  )

  # Column names are "<Group>_<tissue>"; split them back into the two annotation
  # tracks, Group above Tissue.
  annotation_col <- data.frame(
    Group = sub("_.*$", "", colnames(da_feat_perc)),
    Tissue = sub("^[^_]*_", "", colnames(da_feat_perc))
  )
  rownames(annotation_col) <- colnames(da_feat_perc)

  ht <- ComplexHeatmap::pheatmap(
    da_feat_perc,
    heatmap_legend_param = list(title = "% features DA"),
    border_color = "gray3",
    scale = "none",
    color = c("white", "#9e9ac8"),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    annotation_row = annotation_data,
    annotation_colors = ann_colors,
    annotation_col = annotation_col,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_names_col = FALSE,
    display_numbers = as.matrix(da_feat_num),
    number_color = "black",
    fontsize_number = 15,
    angle_col = "90",
    gaps_col = 3,
    na_col = "white"
  )

  export_panel(
    function() {
      ComplexHeatmap::draw(
        ht,
        merge_legends = TRUE,
        heatmap_legend_side = "right",
        annotation_legend_side = "right"
      )
    },
    "FIG2A"
  )
}

# ---- FIG2B — differentially abundant features shared across tissues --------

# Every feature reaching FDR < 0.05 in at least one tissue, counted once, in the
# UpSet intersection of the three tissues it was found in. Bars are stacked by
# ome, so the panel reads both as how much overlap there is and as which
# platforms supply it.
#
# The panel also writes ST2e, the over-representation analysis of the six
# membership sets. FIG2Bii draws the bar plot of the three-tissue set and FIG2F
# draws the heatmap of all six.

fig2b <- function() {
  panel_init("FIG2B")
  presence <- tissue_presence_matrix()

  # Proteomics is the one ome whose cross-tissue identity is not its feature_id.
  # Blood ships Olink and muscle and adipose ship mass spec, so every protein
  # overlap in this panel runs through the OID-to-UniProt mapping that
  # overlap_identity() applies; significant_da() has already done it, and this
  # panel's sets and the ORA below are both keyed on overlap_id.
  #
  # The mapped identity is kept in its own column rather than written over
  # feature_id, so a prot-ol row stays joinable to the rest of the differential
  # analysis on the id it was measured under.
  write_st2e(membership_set_ora())

  # The stripe behind each set row, keyed by name so the colour follows the set
  # wherever the row order puts it.
  stripe_colors <- unname(
    MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS[tolower(TISSUE_SETS)]
  )
  if (anyNA(stripe_colors)) {
    stop("tissue not in HUMAN_TISSUE_COLORS: ",
         paste(TISSUE_SETS[is.na(stripe_colors)], collapse = ", "), call. = FALSE)
  }

  upset_plot <- ComplexUpset::upset(
    presence[, c(TISSUE_SETS, "Ome")],
    TISSUE_SETS,
    encode_sets = FALSE,
    width_ratio = 0.2,
    height_ratio = 1,
    set_sizes = FALSE,
    base_annotations = list(
      "Number of \ndifferential features" = intersection_size(
        counts = TRUE,
        mapping = aes(fill = Ome),
        text = list(size = 3)
      ) +
        scale_fill_manual(
          values = MotrpacHumanPreSuspensionAnalysis::HUMAN_OME_COLORS[OME_LEVELS] |>
            stats::setNames(OME_LEVELS)
        )
    ),
    themes = upset_modify_themes(
      list(intersections_matrix = theme(axis.title.x = element_blank()))
    ),
    sort_sets = FALSE,
    wrap = TRUE,
    stripes = stripe_colors
  )

  export_panel(function() print(upset_plot), "FIG2B")
}

# ---- FIG2Bii — enriched pathways in the three-tissue set -------------------

# The curated pathways from the over-representation analysis of the
# Adipose_Blood_Muscle membership set — the features called differential in
# adipose, blood and muscle alike — as a horizontal bar of -log10(p), labelled
# inside the bar.
#
# Numbered FIG2Bii rather than given a letter of its own: the manuscript labels
# both this and the UpSet as panel 2B, and Figure 2's letters A-F are all
# spoken for. The two are separate PDFs because they were separate artifacts in
# the legacy and are placed independently in the figure.

fig2bii <- function() {
  panel_init("FIG2Bii")
  ora <- membership_set_ora()

  # The curated selection, from config/highlights.json. The enrichment statistics
  # plotted are this run's.
  curated_sets <- highlight_pathways("FIG2Bii")$set

  plot_df <- ora |>
    dplyr::filter(.data$contrast == "Adipose_Blood_Muscle",
                  .data$set %in% curated_sets)

  missing <- setdiff(curated_sets, plot_df$set)
  if (length(missing) > 0) {
    stop(length(missing), " curated pathway(s) absent from this run's ORA: ",
         paste(utils::head(sort(missing), 5), collapse = ", "),
         "\n  The FIG2Bii selection in config/highlights.json was made against a different ",
         "ORA than the one this build produced.", call. = FALSE)
  }

  # Bars ascend, so the strongest enrichment sits at the top after the flip. The
  # factor is built from the sorted frame, which is what fixes the order.
  plot_df <- plot_df[order(plot_df$statistic_column), ]
  plot_df$label <- toupper(gsub("_", " ", sub("^[^_]*_", "", plot_df$set_short)))
  plot_df$label <- factor(plot_df$label, levels = plot_df$label)

  p <- ggplot(plot_df, aes(x = .data$label, y = .data$statistic_column)) +
    geom_col(width = 0.9, fill = "bisque3") +
    geom_text(aes(label = .data$label), hjust = 1.1, color = "grey20", size = 3) +
    coord_flip() +
    labs(x = NULL, y = "-log10(p)") +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  export_panel(p, "FIG2Bii")
}

# ---- FIG2C — features differentially abundant in all three tissues ---------

# The features that reach FDR < 0.05 in adipose, blood and muscle alike, one row
# each, with every exercise-with-controls contrast across the columns: tissue,
# then exercise arm, then timepoint. The cell carries that contrast's z
# statistic and is dotted where the contrast is itself significant, so a row
# shows both where a feature responds and how consistently the three tissues
# agree about it.
#
# Rows are clustered; columns are not. The column order is the tissue x arm x
# timepoint order stated below, not whatever the pivot happened to emit.

fig2c <- function() {
  panel_init("FIG2C")
  all_da <- differential_analysis()
  presence <- tissue_presence_matrix()

  # The tissue order across the heatmap columns. Separate from TISSUE_SETS,
  # which is the UpSet's row order: the presence matrix is indexed by name here,
  # so the two need not agree and this one is the column grid.
  TISSUE_COLUMN_ORDER <- c("Adipose", "Blood", "Muscle")

  shared_ids <- presence$overlap_id[
    presence$Adipose == 1 & presence$Blood == 1 & presence$Muscle == 1
  ]
  if (length(shared_ids) == 0) {
    stop("no feature is differentially abundant in all three tissues", call. = FALSE)
  }

  # The full DA record for those features, significant contrasts and not. The
  # heatmap shows every timepoint, so it is drawn from all_da rather than from the
  # significant subset the overlap was selected on.
  #
  # overlap_id is recomputed here rather than carried over from sig: all_da holds
  # the rows sig filtered away, and those rows need the same identity. It is a
  # pure function of assay, feature_id and uniprot, so the two agree by
  # construction.
  all_da$overlap_id <- overlap_identity(all_da)

  shared_da <- all_da |>
    filter(.data$overlap_id %in% shared_ids)

  # Tissue, then arm, then timepoint. Stated rather than sorted: "post_3.5_4_hr"
  # sorts before "post_10_min" alphabetically, and the columns have to run
  # forwards in time.
  TIMEPOINT_ORDER <- c(
    "during_20_min", "during_40_min", "post_10_min",
    "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr"
  )
  TIMEPOINT_LABELS <- c(
    "during_20_min"     = "during 20 min",
    "during_40_min"     = "during 40 min",
    "post_10_min"       = "post 10 min",
    "post_15_30_45_min" = "post 15/30/45 min",
    "post_3.5_4_hr"     = "post 3.5/4 hr",
    "post_24_hr"        = "post 24 hr"
  )

  # contrast_short is "<arm>.<timepoint>"; split it on the first dot only. The
  # timepoint itself contains dots ("post_3.5_4_hr"), which is why this is a
  # sub() on an anchored prefix and not a fixed-width split.
  shared_da <- shared_da |>
    mutate(
      Modality = unname(c("Endur" = "EE", "Resist" = "RE")[
        sub("\\..*$", "", .data$contrast_short)]),
      tp = sub("^[^.]*\\.", "", .data$contrast_short)
    )

  unknown_tp <- setdiff(unique(shared_da$tp), TIMEPOINT_ORDER)
  if (length(unknown_tp) > 0) {
    stop("timepoint not in the stated column order: ",
         paste(sort(unknown_tp), collapse = ", "),
         "\n  Add it to TIMEPOINT_ORDER and TIMEPOINT_LABELS.", call. = FALSE)
  }

  # Tissue comes off the full DA table's lowercase `tissue`: shared_da is a subset
  # of all_da, and the display columns are only attached to the significant frame.
  shared_da <- shared_da |>
    mutate(Tissue = tissue_display(.data$tissue),
           column = paste(.data$Tissue, .data$Modality, .data$tp, sep = "|"))

  # Every column that exists in the data, in the stated order. Built from the grid
  # rather than from unique(shared_da$column) so the order cannot follow the data.
  column_grid <- expand.grid(
    tp = TIMEPOINT_ORDER,
    Modality = c("EE", "RE"),
    Tissue = TISSUE_COLUMN_ORDER,
    stringsAsFactors = FALSE
  )
  column_grid$column <- with(column_grid, paste(Tissue, Modality, tp, sep = "|"))
  column_grid <- column_grid[column_grid$column %in% shared_da$column, ]

  # One row per feature, one column per contrast. The two matrices are built by
  # the same pivot over the same frame and then indexed to the same column order,
  # so the dot layer and the fill layer cannot come apart — cell_fun() addresses
  # adj_p_mat by the fill matrix's own (i, j).
  to_matrix <- function(value_col, fill) {
    wide <- shared_da |>
      select("overlap_id", "column", all_of(value_col)) |>
      pivot_wider(names_from = "column", values_from = all_of(value_col),
                  values_fn = function(x) x[1], values_fill = fill)

    dup <- shared_da |>
      count(.data$overlap_id, .data$column) |>
      filter(.data$n > 1)
    if (nrow(dup) > 0) {
      stop(nrow(dup), " feature x contrast cell(s) carry more than one value, ",
           "e.g. ", dup$overlap_id[1], " at ", dup$column[1],
           "\n  Two assays have been merged onto one overlap id within one tissue.",
           call. = FALSE)
    }

    m <- as.matrix(wide[, -1, drop = FALSE])
    rownames(m) <- wide$overlap_id
    missing <- setdiff(column_grid$column, colnames(m))
    for (col in missing) m <- cbind(m, stats::setNames(rep(fill, nrow(m)), NULL))
    colnames(m)[seq.int(ncol(m) - length(missing) + 1L, length.out = length(missing))] <-
      missing
    m[, column_grid$column, drop = FALSE]
  }

  # A contrast a feature was never tested in reads as no effect and no
  # significance, which is what an absent cell means here: the assay ran in that
  # tissue, so a gap is a feature the platform did not measure at that timepoint.
  z_mat <- to_matrix("z.std", 0)
  adj_p_mat <- to_matrix("adj_p_value", 1)

  # Row labels are gene symbols where the feature has one. Taken from the DA table
  # rather than re-joined to HUMAN_FEATURE_TO_GENE: combine_with_featgene = TRUE
  # already put it there, and a second join on a key the two tables spell
  # differently is how the legacy lost its Olink rows.
  row_labels <- shared_da |>
    distinct(.data$overlap_id, .keep_all = TRUE) |>
    mutate(
      # as.character() first: gene_symbol arrives as a factor, which nzchar() and
      # if_else()'s type check both reject.
      gene_symbol = as.character(.data$gene_symbol),
      label = if_else(is.na(.data$gene_symbol) | !nzchar(.data$gene_symbol),
                      as.character(.data$overlap_id), .data$gene_symbol))
  rownames(z_mat) <- row_labels$label[match(rownames(z_mat), row_labels$overlap_id)]

  column_ha <- ComplexHeatmap::columnAnnotation(
    Tissue = column_grid$Tissue,
    Modality = column_grid$Modality,
    Timepoint = factor(unname(TIMEPOINT_LABELS[column_grid$tp]),
                       levels = unname(TIMEPOINT_LABELS)),
    col = list(
      Tissue = MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS[
        c("adipose", "blood", "muscle")] |>
        stats::setNames(c("Adipose", "Blood", "Muscle")),
      Modality = MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS[
        c("ADUEndur", "ADUResist")] |>
        stats::setNames(c("EE", "RE")),
      Timepoint = MotrpacHumanPreSuspensionAnalysis::HUMAN_ACUTE_TIMEPOINT_COLORS[
        TIMEPOINT_ORDER] |>
        stats::setNames(unname(TIMEPOINT_LABELS))
    )
  )

  ht <- ComplexHeatmap::Heatmap(
    z_mat,
    name = "z.std",
    top_annotation = column_ha,
    border = TRUE,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    row_names_gp = grid::gpar(fontsize = 12),
    column_names_gp = grid::gpar(fontsize = 10),
    rect_gp = grid::gpar(col = "#2F4F4F", lwd = 1),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (adj_p_mat[i, j] < 0.05) {
        grid::grid.points(x, y, pch = 16, size = grid::unit(5, "pt"),
                          gp = grid::gpar(col = "grey20"))
      }
    }
  )

  export_panel(function() ComplexHeatmap::draw(ht), "FIG2C")
}

# ---- FIG2D — MYC and AREG transcript trajectories --------------------------

# Two features, each drawn as the exercise-group mean with a 95% confidence
# interval across the acute timepoints, one facet per tissue. A point is filled
# black where that timepoint's adjusted p value is below 0.05. The two rows sit
# under one collected legend.
#
# The plot is one entry of the single-feature catalog in single_feature_plots.R.

fig2d <- function() {
  panel_init("FIG2D")
  export_panel(single_feature_plot("FIG2D"), "FIG2D")
}

# ---- FIG2E — phosphosites differential in adipose and muscle ---------------

# The phosphosites reaching FDR < 0.05 in adipose and in muscle alike, one row
# each, labelled by gene symbol and site. Columns are the exercise-with-controls
# contrasts, split by tissue; the cell carries the z statistic and is dotted
# where that contrast is itself significant.
#
# Blood is excluded rather than absent: there is no blood phosphoproteomics, so
# the two-tissue set is the whole of the phospho overlap.

fig2e <- function() {
  panel_init("FIG2E")

  # The shared epigen = TRUE load, filtered to prot-ph. As its own script this
  # panel called cross_tissue_da(epigen = FALSE), because the ATAC and methylCap
  # tables are ~7.7 GB that nothing here reads; in one script per figure the
  # cheaper load would be a second read of the DA tables, and the prot-ph rows
  # are the same rows either way.
  all_da <- differential_analysis()

  phospho_sig <- significant_features() |>
    filter(.data$assay == "prot-ph")

  presence <- feature_presence_matrix(phospho_sig, by = "Tissue",
                                      sets = c("Adipose", "Muscle"))

  shared_ids <- presence$overlap_id[presence$Adipose == 1 & presence$Muscle == 1]
  if (length(shared_ids) == 0) {
    stop("no phosphosite is differentially abundant in both adipose and muscle",
         call. = FALSE)
  }

  # overlap_id and feature_id are the same thing here: overlap_identity() rewrites
  # prot-ol rows only, and this panel is prot-ph. Filtering on feature_id is
  # therefore exact rather than approximately right, but the equality is asserted
  # so that a change to overlap_identity() cannot silently make it approximate.
  stopifnot(identical(shared_ids, as.character(shared_ids)))
  shared_da <- all_da |>
    filter(.data$assay == "prot-ph", .data$feature_id %in% shared_ids)

  if (length(unique(shared_da$feature_id)) != length(shared_ids)) {
    stop("phosphosite overlap ids do not round-trip to feature ids: ",
         length(shared_ids), " selected, ",
         length(unique(shared_da$feature_id)), " found in the DA table",
         call. = FALSE)
  }

  # A prot-ph feature_id ends in the site list: "<protein>_<site><residue>", and a
  # multi-site peptide carries several. The legacy stripped exactly one trailing
  # character with sub(".*_(.+).{1}$", "\\1", x), which is right for a single site
  # ("..._S575s" -> "S575") and wrong for every multi-site one
  # ("..._S663sS665s" -> "S663sS665"). 5,002 of the 31,705 prot-ph features in
  # HUMAN_FEATURE_TO_GENE carry more than one site.
  #
  # Taken instead as everything after the last underscore, with the lowercase
  # residue markers dropped, so a multi-site label reads "S663 S665".
  phosphosite_label <- function(feature_id) {
    sites <- sub("^.*_", "", feature_id)
    sites <- gsub("([A-Z][0-9]+)[a-z]", "\\1 ", sites)
    trimws(sites)
  }

  row_meta <- shared_da |>
    distinct(.data$feature_id, .keep_all = TRUE) |>
    mutate(
      # as.character() first: gene_symbol and feature_id arrive as factors, and
      # both nzchar() and if_else()'s type check reject those.
      gene_symbol = as.character(.data$gene_symbol),
      feature_id = as.character(.data$feature_id),
      symbol = if_else(is.na(.data$gene_symbol) | !nzchar(.data$gene_symbol),
                       .data$feature_id, .data$gene_symbol),
      label = paste(.data$symbol, phosphosite_label(.data$feature_id))
    )

  if (anyDuplicated(row_meta$label) > 0) {
    dupes <- row_meta$label[duplicated(row_meta$label)]
    stop(length(unique(dupes)), " phosphosite label(s) are not unique, e.g. ",
         paste(utils::head(sort(unique(dupes)), 3), collapse = ", "),
         "\n  Two features would land on one heatmap row.", call. = FALSE)
  }

  TIMEPOINT_LABELS <- c(
    "post_10_min"       = "post 10 min",
    "post_15_30_45_min" = "post 15/30/45 min",
    "post_3.5_4_hr"     = "post 3.5/4 hr",
    "post_24_hr"        = "post 24 hr"
  )

  shared_da <- shared_da |>
    mutate(
      Modality = unname(c("Endur" = "EE", "Resist" = "RE")[
        sub("\\..*$", "", .data$contrast_short)]),
      tp = sub("^[^.]*\\.", "", .data$contrast_short)
    )

  unknown_tp <- setdiff(unique(shared_da$tp), names(TIMEPOINT_LABELS))
  if (length(unknown_tp) > 0) {
    stop("timepoint not in the stated column order: ",
         paste(sort(unknown_tp), collapse = ", "),
         "\n  Add it to TIMEPOINT_LABELS.", call. = FALSE)
  }

  # Tissue comes off the full DA table's lowercase `tissue`: shared_da is a subset
  # of all_da, and the display columns are only attached to the significant frame.
  shared_da$Tissue <- tissue_display(shared_da$tissue)
  shared_da$column <- paste(shared_da$Tissue, shared_da$Modality, shared_da$tp,
                            sep = "|")

  # Tissue, then arm, then timepoint, from the grid rather than from the data.
  # The legacy hardcoded the split as c(rep("Adipose", 2), rep("Muscle", 6)); this
  # derives both the order and the split from the columns that actually exist, so
  # a timepoint appearing or disappearing cannot silently shift the tissue
  # boundary.
  column_grid <- expand.grid(
    tp = names(TIMEPOINT_LABELS),
    Modality = c("EE", "RE"),
    Tissue = c("Adipose", "Muscle"),
    stringsAsFactors = FALSE
  )
  column_grid$column <- with(column_grid, paste(Tissue, Modality, tp, sep = "|"))
  column_grid <- column_grid[column_grid$column %in% shared_da$column, ]

  to_matrix <- function(value_col, fill) {
    dup <- shared_da |>
      count(.data$feature_id, .data$column) |>
      filter(.data$n > 1)
    if (nrow(dup) > 0) {
      stop(nrow(dup), " site x contrast cell(s) carry more than one value, e.g. ",
           dup$feature_id[1], " at ", dup$column[1], call. = FALSE)
    }

    wide <- shared_da |>
      select("feature_id", "column", all_of(value_col)) |>
      pivot_wider(names_from = "column", values_from = all_of(value_col))

    m <- as.matrix(wide[, -1, drop = FALSE])
    rownames(m) <- wide$feature_id
    m[is.na(m)] <- fill
    absent <- setdiff(column_grid$column, colnames(m))
    if (length(absent) > 0) {
      pad <- matrix(fill, nrow = nrow(m), ncol = length(absent),
                    dimnames = list(rownames(m), absent))
      m <- cbind(m, pad)
    }
    m[, column_grid$column, drop = FALSE]
  }

  z_mat <- to_matrix("z.std", 0)
  adj_p_mat <- to_matrix("adj_p_value", 1)
  rownames(z_mat) <- row_meta$label[match(rownames(z_mat), row_meta$feature_id)]

  column_ha <- ComplexHeatmap::columnAnnotation(
    Tissue = column_grid$Tissue,
    Modality = column_grid$Modality,
    Timepoint = factor(unname(TIMEPOINT_LABELS[column_grid$tp]),
                       levels = unname(TIMEPOINT_LABELS)),
    col = list(
      Tissue = MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS[
        c("adipose", "muscle")] |>
        stats::setNames(c("Adipose", "Muscle")),
      Modality = MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS[
        c("ADUEndur", "ADUResist")] |>
        stats::setNames(c("EE", "RE")),
      Timepoint = MotrpacHumanPreSuspensionAnalysis::HUMAN_ACUTE_TIMEPOINT_COLORS[
        names(TIMEPOINT_LABELS)] |>
        stats::setNames(unname(TIMEPOINT_LABELS))
    )
  )

  ht <- ComplexHeatmap::Heatmap(
    z_mat,
    name = "z.std",
    top_annotation = column_ha,
    border = TRUE,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    column_split = factor(column_grid$Tissue, levels = c("Adipose", "Muscle")),
    column_gap = grid::unit(2, "mm"),
    show_column_names = FALSE,
    row_names_gp = grid::gpar(fontsize = 12),
    column_names_gp = grid::gpar(fontsize = 10),
    rect_gp = grid::gpar(col = "#2F4F4F", lwd = 1),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (adj_p_mat[i, j] < 0.05) {
        grid::grid.points(x, y, pch = 16, size = grid::unit(5, "pt"),
                          gp = grid::gpar(col = "grey20"))
      }
    }
  )

  export_panel(function() ComplexHeatmap::draw(ht), "FIG2E")
}

# ---- FIG2F — pathway enrichment across the tissue-membership sets ----------

# The curated pathways from the same over-representation analysis FIG2B writes
# as ST2e, one column per membership set and one row per pathway, coloured by
# -log10(p). A cell whose enrichment does not clear the adjusted-p cut is filled
# grey rather than left to read as a weak result.
#
# Columns run Adipose, Muscle, Blood, then Adipose and Muscle, Blood and Muscle,
# then the three-tissue set. Rows cluster; columns do not.

fig2f <- function() {
  panel_init("FIG2F")
  ora <- membership_set_ora()

  # The curated selection, from config/highlights.json. The enrichment statistics
  # drawn are this run's.
  curated_sets <- highlight_pathways("FIG2F")$set

  missing <- setdiff(curated_sets, ora$set)
  if (length(missing) > 0) {
    stop(length(missing), " curated pathway(s) absent from this run's ORA: ",
         paste(utils::head(sort(missing), 5), collapse = ", "),
         "\n  The FIG2F selection in config/highlights.json was made against a different ORA ",
         "than the one this build produced.", call. = FALSE)
  }

  plot_df <- ora |>
    dplyr::filter(.data$set %in% curated_sets) |>
    dplyr::select("set_short", "statistic_column", "adj_p_value", "contrast")

  # Every set x pathway cell, including the ones the ORA returned no row for.
  # enrichmap() needs the full grid: without it a pathway tested in four of the
  # six sets is drawn against a legend that says the other two were not
  # significant, when in fact they were never plotted. Widening and re-lengthening
  # is what materialises the absent cells as NA.
  plot_df <- plot_df |>
    tidyr::pivot_wider(
      id_cols = "set_short",
      names_from = "contrast",
      values_from = c("adj_p_value", "statistic_column"),
      names_sep = "@"
    ) |>
    tidyr::pivot_longer(
      cols = -"set_short",
      names_to = c(".value", "contrast"),
      names_sep = "@"
    )

  # The column order, and the labels under it. Applied to levels() rather than to
  # the values so the two cannot come apart.
  column_order <- c("Adipose_only", "Muscle_only", "Blood_only",
                    "Adipose_and_Muscle_only", "Blood_and_Muscle_only",
                    "Adipose_Blood_Muscle")
  stopifnot(setequal(column_order, names(CROSS_TISSUE_ORA_SETS)))
  plot_df$contrast <- factor(plot_df$contrast, levels = column_order)
  levels(plot_df$contrast) <- gsub(" only", "", gsub("_", " ",
                                                     levels(plot_df$contrast)))
  plot_df <- plot_df[order(plot_df$contrast), ]

  # enrichmap() finishes with ComplexHeatmap::draw(), whose newpage default starts
  # a fresh page on the device export_panel() already opened — leaving page 1
  # blank and the heatmap on page 2, the page Illustrator discards. draw_args
  # carries newpage = FALSE through, the same way FIG3B does.
  heatmap <- function() {
    TMSig::enrichmap(
      as.data.frame(plot_df),
      set_column = "set_short",
      n_top = 35,
      statistic_column = "statistic_column",
      contrast_column = "contrast",
      padj_column = "adj_p_value",
      padj_legend_title = "adj p \n(background)",
      padj_fill = "grey70",
      colors = c("white", "#543483"),
      heatmap_args = list(
        name = "-log10(p)",
        na_col = "grey90",
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        border = TRUE,
        rect_gp = grid::gpar(col = "grey70", lwd = 2)
      ),
      draw_args = list(newpage = FALSE)
    )
  }

  export_panel(heatmap, "FIG2F")
}

# ---- ST2c — DA features per timepoint --------------------------------------

# Counted off all_da directly rather than off FIG2A's all_da_sum: a feature
# significant at three timepoints is three rows here, where the panel keeps only
# its most significant one. The table's legend says so explicitly.
write_st2c <- function(all_da) {
  st2c_long <- all_da %>%
    filter(.data$contrast_type == "exercise_with_controls",
           .data$adj_p_value < 0.05) %>%
    count(.data$assay, .data$randomGroupCode, .data$tissue, .data$Timepoint,
          name = "n") %>%
    mutate(column = paste0(str_to_title(.data$tissue), "_", .data$Timepoint))

  st2c_declared <- table_spec("ST2c")$columns
  st2c_value_cols <- setdiff(st2c_declared, c("assay", "randomGroupCode"))
  st2c_unexpected <- setdiff(unique(st2c_long$column), st2c_value_cols)
  if (length(st2c_unexpected) > 0) {
    stop("ST2c: DA results for tissue x timepoint combinations the manifest does ",
         "not declare: ", paste(sort(st2c_unexpected), collapse = ", "),
         "\n  Add them to `columns` for ST2c in config/table_map.json.",
         call. = FALSE)
  }

  # No roster here, unlike ST1c and ST2b. The rows are the assay x group
  # combinations the differential analysis actually produced, and a combination
  # with no significant feature anywhere is absent rather than blank — there is no
  # fixed set of DA results the way there is a fixed set of assays.
  st2c <- st2c_long %>%
    select("assay", "randomGroupCode", "column", "n") %>%
    tidyr::pivot_wider(names_from = "column", values_from = "n") %>%
    as.data.frame()
  for (missing_col in setdiff(st2c_value_cols, names(st2c))) {
    st2c[[missing_col]] <- NA_integer_
  }
  st2c <- order_by_assay(st2c, "assay", st2c[["randomGroupCode"]])[, st2c_declared]

  export_table(st2c, "ST2c")
}

# ---- ST2e — over-representation of each membership set ---------------------

write_st2e <- function(ora) {
  export_table(
    as.data.frame(ora)[, table_spec("ST2e")$columns, drop = FALSE],
    "ST2e"
  )
}

run_panels(list(FIG2A = fig2a, FIG2B = fig2b, FIG2Bii = fig2bii, FIG2C = fig2c,
                FIG2D = fig2d, FIG2E = fig2e, FIG2F = fig2f))
