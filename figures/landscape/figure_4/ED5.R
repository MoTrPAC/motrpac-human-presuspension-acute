#!/usr/bin/env Rscript
# Extended Data 5 — fuzzy c-means supplement
#
# Panels:  ED5A  RefMet classes enriched in the acute metabolomics response
#          ED5B  cluster ome composition and enrichment across the curated
#                pathways
#          ED5C  overlap of VEGF signalling features across enriched clusters
#          ED5D  VEGF features shared by the early muscle and adipose clusters
#          ED5E  VEGF features shared by the middle muscle and blood clusters
#
# Every panel here reads the CAMERA-PR results or the fuzzy c-means objects
# shipped in MotrpacHumanPreSuspensionAnalysis; nothing is refitted. The hand-set
# cluster numbers are in figure_4/cmeans.env and the curated pathways in
# config/highlights.json.
#
#   Rscript figures/landscape/ED5.R          every panel
#   Rscript figures/landscape/ED5.R ED5C     one panel

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(dplyr)
  library(ggnewscale)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(MotrpacHumanPreSuspensionAnalysis)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(here, "..", "lib", "highlights.R"))
source(file.path(here, "FIG4_ED5_helpers.R"))

# ---- shared ----------------------------------------------------------------

# Tissue order across the c-means panels, as in FIG4.R. The same three names as
# the helper's CMEANS_TISSUES, which is the heatmap strip order; this is the
# panel layout order.
TISSUE_ORDER <- c("adipose", "blood", "muscle")

# There is no shared loader here. ED5B, ED5C and the two cross-tissue heatmaps
# all start from the c-means objects, but ED5D and ED5E call
# cmeans_heatmap_features() with different tissues and different clusters, so
# their selections cannot be shared. The one call every panel makes with the same
# arguments is cmeans_plot_frame(), and that is memoised in
# helpers/FIG4_ED5.R so FIG4 and ED5 alike read FCM_CLUSTERS once.

# ---- ED5A — RefMet classes enriched in the metabolomics response -----------

# FIG4A and FIG4B over the metabolomics differential analysis. The rows are
# RefMet chemical classes rather than pathways, because REFMET is the only
# collection CAMERA-PR is run against for metab.
#
# THE SELECTION IS COMPUTED, NOT CURATED. Where FIG4A and FIG4B draw a hand
# list from config/highlights.json, this panel takes the six most significant
# classes from each tissue and contrast and draws their union — n_top = 6 in
# camera_enrich_heatmap(), which is the upstream default and what the legacy
# call used. Nothing about it is hand-set, so it has no highlights.json entry.
#
# Ordering is on the raw p-value rather than the adjusted one, to avoid ties
# among classes that share an adjusted p.

ed5a <- function() {
  panel_init("ED5A")

  # ---- plot ----

  heatmap <- camera_enrich_heatmap(
    set_ids = NULL,
    n_top = 6L,
    selected_ome = "metab",
    selected_tissues = c("adipose", "blood", "muscle"),
    contrast_type = "exercise_with_controls"
  )

  message(sprintf("        ED5A: drawing %d RefMet class(es) on a %.1f x %.1f in page",
                  heatmap$n_sets, heatmap$width, heatmap$height))

  export_panel(heatmap$draw, "ED5A",
               width = heatmap$width, height = heatmap$height)
}

# ---- ED5B — cluster composition and enrichment, all curated pathways -------

# The same three registers as FIG4D — cluster number, ome composition, pathway
# enrichment — over every curated pathway rather than VEGF alone. FIG4D is the
# one-row reading of this panel; this is the whole grid behind it.
#
# Row labels are coloured from CMEANS_SET_PALETTE, one hue per pathway, and
# each row carries a guide line in that colour: at fifteen rows the eye needs
# help tracking a row across 37 columns.
#
# Clusters are ordered tissue by tissue and, within a tissue, by when they peak
# and in which direction — the same grouping that colours FIG4C. Reading left to
# right is reading the response in time.
#
# The pathway rows are ordered by clustering the enrichment matrix, so pathways
# that light up in the same clusters sit together. Adipose is excluded from that
# clustering: it enriches for far fewer sets than the other two and would
# otherwise dominate the distance.

ed5b <- function() {
  panel_init("ED5B")

  BASE_SIZE <- 6
  TEXT_SIZE <- BASE_SIZE / 3.5
  POINT_MAX <- BASE_SIZE / 1.5      # the largest enrichment dot
  CIRCLE_SIZE <- 3                  # the cluster circles
  KEY_PT <- 8

  # ---- the x axis: tissue x cluster, ordered by when the cluster peaks ----

  centroids <- cmeans_plot_frame()

  cluster_meta <- unique(centroids[, c("tissue", "cluster", "peak_timepoint",
                                       "peak_sign", "peak_group")])
  cluster_meta$tissue <- factor(cluster_meta$tissue, levels = TISSUE_ORDER)

  # Peak group order: every timepoint, positive before negative, so the columns
  # run forwards in time within each tissue.
  PEAK_ORDER <- c("P15-45M", "P3.5/4H", "P24H")
  group_levels <- as.vector(t(outer(PEAK_ORDER, c("Positive", "Negative"),
                                    paste, sep = "_")))
  cluster_meta$peak_group <- factor(cluster_meta$peak_group, levels = group_levels)

  cluster_meta <- cluster_meta[order(cluster_meta$tissue, cluster_meta$peak_group,
                                     as.integer(cluster_meta$cluster)), ,
                               drop = FALSE]
  cluster_meta <- cluster_meta[!is.na(cluster_meta$peak_group), , drop = FALSE]
  cluster_meta$label <- paste(cluster_meta$cluster, cluster_meta$tissue, sep = " - ")
  cluster_meta$label <- factor(cluster_meta$label, levels = cluster_meta$label)

  # One colour per tissue x cluster, keyed by the same label the x axis uses, so
  # the circle above a column is the colour that cluster carries in FIG4C.
  colour_lookup <- unique(centroids[, c("tissue", "cluster", "cluster_color")])
  colour_lookup$label <- paste(colour_lookup$cluster, colour_lookup$tissue,
                               sep = " - ")
  cluster_colors <- stats::setNames(
    colour_lookup$cluster_color[match(levels(cluster_meta$label),
                                      colour_lookup$label)],
    levels(cluster_meta$label)
  )
  if (anyNA(cluster_colors)) {
    stop("cluster with no colour: ",
         paste(names(cluster_colors)[is.na(cluster_colors)], collapse = ", "),
         call. = FALSE)
  }

  # ---- register 2: how many features of each ome per cluster ----

  assignments <- cmeans_hard_assignments()
  assignments <- assignments[!is.na(assignments$cluster), , drop = FALSE]

  # The ome's display name. metab is set explicitly: assay_codes carries a row per
  # metabolomics PLATFORM, so joining on it alone would split one cluster's
  # metabolites across a dozen legend entries.
  ome_text <- MotrpacBicQC::assay_codes[, c("assay_code", "omics_text")]
  ome_text <- unique(ome_text[!is.na(ome_text$assay_code), ])

  composition <- assignments |>
    count(.data$tissue, .data$cluster, .data$assay, name = "Count") |>
    left_join(ome_text, by = c("assay" = "assay_code"),
              relationship = "many-to-many") |>
    mutate(omics_text = if_else(.data$assay == "metab", "Metabolomics",
                                .data$omics_text)) |>
    filter(!is.na(.data$omics_text)) |>
    distinct(.data$tissue, .data$cluster, .data$assay, .data$omics_text,
             .data$Count) |>
    mutate(label = paste(.data$cluster, .data$tissue, sep = " - ")) |>
    filter(.data$label %in% levels(cluster_meta$label)) |>
    mutate(label = factor(.data$label, levels = levels(cluster_meta$label)))

  # ---- register 3: which curated pathways each cluster is enriched for ----

  # The curated pathways and their figure labels, from config/highlights.json.
  # FIG4D narrows the same machinery to the single cmeans selection.
  curated <- highlight_pathways("ED5B")
  if (anyNA(curated$label)) {
    stop("ED5B selection in config/highlights.json has a set with no label",
         call. = FALSE)
  }
  set_labels <- stats::setNames(curated$label, curated$set)

  hits <- cmeans_ora_hits()
  hits$label <- paste(as.character(as.numeric(hits$cluster)), hits$tissue,
                      sep = " - ")

  enrichment <- hits |>
    filter(.data$set %in% names(set_labels),
           .data$label %in% levels(cluster_meta$label)) |>
    left_join(ome_text, by = c("assay" = "assay_code"),
              relationship = "many-to-many") |>
    mutate(omics_text = if_else(.data$assay == "metab", "Metabolomics",
                                .data$omics_text)) |>
    # The share of the cluster that belongs to the pathway.
    mutate(included_ratio = .data$set_size_in_cluster / .data$cluster_size) |>
    distinct(.data$set, .data$label, .data$omics_text, .data$included_ratio) |>
    mutate(label = factor(.data$label, levels = levels(cluster_meta$label)))

  if (nrow(enrichment) == 0) {
    stop("none of the ", length(set_labels), " curated ED5B pathways is ",
         "enriched in any cluster", call. = FALSE)
  }
  absent <- setdiff(names(set_labels), unique(enrichment$set))
  if (length(absent) > 0) {
    message(sprintf("        ED5B: %d curated pathway(s) enriched nowhere and not drawn: %s",
                    length(absent), paste(unname(set_labels[absent]), collapse = ", ")))
  }

  # ---- pathway row order, by clustering the enrichment matrix ----

  wide <- enrichment |>
    select("set", "label", "included_ratio") |>
    group_by(.data$set, .data$label) |>
    summarise(included_ratio = max(.data$included_ratio), .groups = "drop") |>
    pivot_wider(names_from = "label", values_from = "included_ratio",
                values_fill = 0)
  matrix_rows <- as.matrix(wide[, -1, drop = FALSE])
  rownames(matrix_rows) <- wide$set
  # Rows enter the clustering in the curated order: hclust breaks ties and lays
  # out leaves by input order, and the legacy fed it the sheet's order.
  matrix_rows <- matrix_rows[order(match(rownames(matrix_rows), names(set_labels))), ,
                             drop = FALSE]

  non_adipose <- !grepl("adipose", colnames(matrix_rows), ignore.case = TRUE)
  if (sum(non_adipose) >= 2 && nrow(matrix_rows) >= 2) {
    ordering <- stats::hclust(stats::dist(matrix_rows[, non_adipose, drop = FALSE],
                                          method = "euclidean"),
                              method = "complete")$order
    set_order <- rev(rownames(matrix_rows)[ordering])
  } else {
    set_order <- rownames(matrix_rows)
  }
  enrichment$set <- factor(enrichment$set, levels = set_order)

  set_colors <- stats::setNames(rep_len(CMEANS_SET_PALETTE, length(set_order)),
                                set_order)

  # ---- the tissue bands behind the dot plot ----

  bands <- do.call(rbind, lapply(split(cluster_meta, cluster_meta$tissue), function(g) {
    if (nrow(g) == 0) return(NULL)
    idx <- match(as.character(g$label), levels(cluster_meta$label))
    data.frame(tissue = as.character(g$tissue[1]),
               xmin = min(idx) - 0.5, xmax = max(idx) + 0.5,
               stringsAsFactors = FALSE)
  }))

  # ---- plot ----

  # Legend order follows the stack. With the axis reversed the stack is drawn from
  # the top down, so reverse = TRUE keeps the key in the order the bars read.
  circles <- ggplot(cluster_meta, aes(x = .data$label, y = 1)) +
    geom_point(aes(fill = .data$label), shape = 21, size = CIRCLE_SIZE,
               color = "black", show.legend = FALSE) +
    geom_text(aes(label = .data$cluster), size = TEXT_SIZE) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = cluster_colors) +
    theme_void(base_size = BASE_SIZE) +
    theme(plot.margin = margin(t = 2, unit = "pt"))

  # fill is mapped INSIDE geom_col, not at the top level. ggnewscale gives the new
  # scale to layers added after new_scale_fill(), and a top-level aes is bound
  # before it — which silently left every bar on the default grey.
  bars <- ggplot(composition, aes(x = .data$label, y = .data$Count)) +
    # The tissue bands run behind this register too, so a column belongs to a
    # tissue for its whole height rather than only where the dots are.
    geom_rect(data = bands, aes(xmin = .data$xmin, xmax = .data$xmax,
                                fill = .data$tissue),
              ymin = -Inf, ymax = Inf, alpha = 0.2, inherit.aes = FALSE) +
    scale_fill_manual(values = cmeans_tissue_colors(TISSUE_ORDER)) +
    guides(fill = "none") +
    new_scale_fill() +
    geom_col(aes(fill = .data$omics_text), width = 0.8) +
    scale_x_discrete(drop = FALSE) +
    # Bars hang DOWN from the cluster circles rather than rising from a baseline,
    # so the three registers read top to bottom in one direction: which cluster,
    # then what it is made of, then what it is enriched for.
    scale_y_reverse(expand = expansion(mult = c(0.05, 0))) +
    scale_fill_manual(name = "Omics Layer",
                      values = MotrpacHumanPreSuspensionAnalysis::HUMAN_OME_COLORS) +
    labs(x = NULL, y = "Number of Features") +
    theme_minimal(base_size = BASE_SIZE) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_line(),
      panel.grid.major.y = element_line(color = NA),
      panel.grid.minor.y = element_line(color = NA),
      plot.margin = margin(t = 0.1, b = 0.1, unit = "pt"),
      # Legend to the left, outside the registers, so it does not eat plot width
      # on the right where the muscle columns already run to the edge.
      legend.position = "left",
      legend.title = element_text(face = "plain")
    ) +
    guides(fill = guide_legend(reverse = TRUE,
                               keywidth = unit(KEY_PT, "pt"),
                               keyheight = unit(KEY_PT, "pt")))

  dots <- ggplot(enrichment, aes(x = .data$label, y = .data$set)) +
    geom_rect(data = bands, aes(xmin = .data$xmin, xmax = .data$xmax,
                                fill = .data$tissue),
              ymin = -Inf, ymax = Inf, alpha = 0.2, inherit.aes = FALSE) +
    scale_fill_manual(values = cmeans_tissue_colors(TISSUE_ORDER)) +
    guides(fill = "none") +
    # A second fill scale on one plot: the bands are filled by tissue and the dots
    # by ome, and ggplot allows one scale per aesthetic without this.
    new_scale_fill() +
    geom_hline(aes(yintercept = .data$set, color = .data$set), alpha = 0.6) +
    geom_point(aes(size = .data$included_ratio, fill = .data$omics_text),
               shape = 21, color = "black") +
    scale_color_manual(values = set_colors) +
    scale_fill_manual(values = MotrpacHumanPreSuspensionAnalysis::HUMAN_OME_COLORS,
                      na.translate = FALSE) +
    # drop = FALSE is load-bearing. A discrete scale shows only the levels present
    # in ITS OWN data, so a cluster with no enrichment was dropped from this
    # register alone — compressing 37 columns to 21 and sliding every dot out from
    # under the bar it belongs to. The registers must share one x.
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(labels = set_labels,
                     expand = expansion(mult = c(0.04, 0.07))) +
    scale_size_continuous(range = c(1, POINT_MAX)) +
    labs(x = NULL, y = "Pathway Set", size = "Included Ratio") +
    theme_minimal(base_size = BASE_SIZE) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = unname(set_colors[set_order])),
      axis.ticks.y = element_line(),
      panel.grid.major.x = element_line(color = "grey90"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.margin = margin(t = -5, unit = "pt"),
      legend.box.margin = margin(t = -5, unit = "pt"),
      plot.margin = margin(t = 0.1, unit = "pt")
    ) +
    guides(color = "none", fill = "none")

  combined <- circles / bars / dots +
    patchwork::plot_layout(heights = c(0.2, 1, 1.8))

  export_panel(function() print(combined), "ED5B")
}

# ---- ED5C — overlap of VEGF features across the enriched clusters ----------

# Each cluster enriched for the pathway contributes the set of its own features
# that belong to that pathway, and the UpSet shows which features those sets
# share. A feature in one bar only is specific to that cluster; a feature in an
# intersection responds the same way in two tissues at once.
#
# Two pairs of clusters are combined into one set before the intersections are
# taken — blood 7 with blood 9, muscle 6 with muscle 10. They are the same
# trajectory split in two by the clustering rather than two findings, and
# leaving them apart produces intersections that say only that. Which pairs
# merge is CMEANS_UPSET_<TISSUE>_COMBINED in figure_4/cmeans.env.
#
# The right-hand bars keep the ORIGINAL clusters: a merged set is drawn as a
# stack of its constituents in their own colours, so the merge is visible rather
# than hidden.

ed5c <- function() {
  spec <- panel_init("ED5C")

  # ---- data ----

  pathway <- cmeans_pathway()
  ome <- cmeans_setting("OME")

  enriched <- cmeans_enriched_clusters(pathway, ome)

  # The pathway's features, restricted to the ome the enrichment was computed on.
  feature_map <- MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE
  pathway_genes <- cmeans_pathway_features(pathway)
  pathway_features <- unique(as.character(
    feature_map$feature_id[as.character(feature_map$gene_symbol) %in% pathway_genes &
                             as.character(feature_map$assay) == ome]
  ))
  if (length(pathway_features) == 0) {
    stop("'", pathway, "' maps to no ", ome, " feature", call. = FALSE)
  }

  assignments <- cmeans_hard_assignments()
  members <- assignments |>
    filter(!is.na(.data$cluster),
           .data$assay == ome,
           .data$feature_id %in% pathway_features) |>
    mutate(group_id = paste0(.data$tissue, "-", .data$cluster)) |>
    filter(.data$group_id %in% enriched$group_id)

  if (nrow(members) == 0) {
    stop("no ", ome, " feature of '", pathway, "' falls in an enriched cluster",
         call. = FALSE)
  }

  # One set per enriched cluster, before merging.
  granular <- split(members$feature_id, members$group_id)

  # ---- merge the split trajectories ----

  merges <- list()
  for (tissue in c("adipose", "blood", "muscle")) {
    var <- paste0("UPSET_", toupper(tissue), "_COMBINED")
    if (!nzchar(Sys.getenv(paste0("CMEANS_", var), unset = ""))) next
    clusters <- cmeans_clusters(var)
    ids <- paste0(tissue, "-", clusters)
    present <- intersect(ids, names(granular))
    if (length(present) < 2) next
    name <- paste0(tissue, "-", paste(clusters, collapse = "/"), " (Combined)")
    merges[[name]] <- present
  }

  merged <- granular
  for (name in names(merges)) {
    merged[[name]] <- unique(unlist(merged[merges[[name]]], use.names = FALSE))
    merged[setdiff(merges[[name]], name)] <- NULL
  }
  message(sprintf("        ED5C: %d enriched cluster(s) -> %d set(s) after merging",
                  length(granular), length(merged)))

  # ---- the right-hand bars, stacked by original cluster ----

  # A merged set's bar is its constituents stacked; an unmerged set's bar is
  # itself. Built as a matrix of set x original cluster so ComplexHeatmap can
  # stack it.
  originals <- names(granular)
  counts <- matrix(0, nrow = length(merged), ncol = length(originals),
                   dimnames = list(names(merged), originals))
  for (set_name in names(merged)) {
    constituents <- if (set_name %in% names(merges)) {
      merges[[set_name]]
    } else {
      intersect(set_name, originals)
    }
    for (original in constituents) {
      counts[set_name, original] <- length(granular[[original]])
    }
  }

  # Each original cluster keeps the colour it carries in FIG4C and FIG4D.
  colour_lookup <- unique(cmeans_plot_frame()[, c("tissue", "cluster",
                                                  "cluster_color")])
  colour_lookup$group_id <- paste0(colour_lookup$tissue, "-", colour_lookup$cluster)
  bar_colors <- colour_lookup$cluster_color[match(originals, colour_lookup$group_id)]
  bar_colors[is.na(bar_colors)] <- "grey50"

  # ---- plot ----

  combinations <- ComplexHeatmap::make_comb_mat(merged)

  # Everything is sized in inches against the manifest page, not in points: at
  # 2.3 x 1.68 in a default 12pt label is most of the panel.
  TOP_H <- spec$height_in * 0.40
  RIGHT_W <- spec$width_in * 0.20
  FS <- 3.5

  upset <- ComplexHeatmap::UpSet(
    combinations,
    pt_size = grid::unit(1.5, "mm"),
    lwd = 1,
    comb_order = order(ComplexHeatmap::comb_size(combinations), decreasing = TRUE),
    set_order = order(ComplexHeatmap::set_size(combinations), decreasing = TRUE),
    top_annotation = ComplexHeatmap::upset_top_annotation(
      combinations,
      add_numbers = TRUE,
      numbers_offset = grid::unit(1, "mm"),
      height = grid::unit(TOP_H, "inch"),
      numbers_rot = 0,
      numbers_gp = grid::gpar(fontsize = FS * 0.8),
      annotation_name_gp = grid::gpar(fontsize = FS),
      axis_param = list(gp = grid::gpar(fontsize = FS))
    ),
    right_annotation = ComplexHeatmap::rowAnnotation(
      "Set Size" = ComplexHeatmap::anno_barplot(
        counts[ComplexHeatmap::set_name(combinations), , drop = FALSE],
        gp = grid::gpar(fill = bar_colors, col = NA),
        border = FALSE,
        axis_param = list(side = "bottom", gp = grid::gpar(fontsize = FS))
      ),
      width = grid::unit(RIGHT_W, "inch"),
      annotation_name_gp = grid::gpar(fontsize = FS)
    ),
    row_names_gp = grid::gpar(fontsize = FS),
    column_title = paste0("Feature Overlap: ", pathway),
    column_title_gp = grid::gpar(fontsize = FS, fontface = "bold")
  )

  # draw() with newpage = FALSE: export_panel() has already opened the device, and
  # the default would lay a blank page down in front of the plot.
  export_panel(function() ComplexHeatmap::draw(upset, newpage = FALSE), "ED5C")
}

# ---- ED5D — VEGF features shared by early muscle and adipose clusters ------

# One row per feature, one column per exercise-with-controls contrast in every
# tissue, adipose, blood and muscle. Rows are the features measured in all
# three, clustered. Three strips on the right name each feature's cluster per
# tissue in the per-cluster colours FIG4C and FIG4D use.
#
# Narrowed to features that ALSO land in the named adipose cluster(s). That
# intersection is the point of the panel: these features move together in muscle
# and in adipose, not merely in muscle.
#
# The pathway comes from config/highlights.json; the ome and clusters from
# figure_4/cmeans.env.
#
# ComplexHeatmap is attached, not qualified, at the top of this script: the
# Analysis package's cell-drawing function calls textGrob() unqualified, and
# only finds it because attaching ComplexHeatmap attaches grid. Without that the
# panel dies with "could not find function textGrob".

ed5d <- function() {
  panel_init("ED5D")

  # ---- data ----

  pathway <- cmeans_pathway()
  ome <- cmeans_setting("OME")

  selection <- cmeans_heatmap_features(
    pathway, ome,
    tissue = "muscle",
    clusters = cmeans_clusters("HEATMAP_MUSCLE_ADIPOSE_EARLY"),
    with_tissue = "adipose",
    with_clusters = cmeans_clusters("VEGF_ADIPOSE")
  )
  message(sprintf("        ED5D: %d %s feature(s) of %s",
                  length(selection$feature_ids), ome, cmeans_pathway_label()))

  # ---- plot ----

  # The heatmap draws itself onto the device export_panel() opens rather than
  # returning an object; see cmeans_feature_heatmap() for why.
  draw_heatmap <- function() {
    cmeans_feature_heatmap(
      annotation = selection$annotation,
      feature_ids = selection$feature_ids,
      selected_tissue = CMEANS_TISSUES,
      selected_ome = ome,
      contrast_type = "exercise_with_controls",
      verbose = FALSE
    )
  }

  export_panel(draw_heatmap, "ED5D")
}

# ---- ED5E — VEGF features shared by middle muscle and blood clusters -------

# ED5D against blood instead of adipose: the same heatmap, the same three
# cluster strips, narrowed to the features that also land in the named blood
# cluster(s). These features move together in muscle and in blood, not merely in
# muscle.
#
# The pathway comes from config/highlights.json; the ome and clusters from
# figure_4/cmeans.env.

ed5e <- function() {
  panel_init("ED5E")

  # ---- data ----

  pathway <- cmeans_pathway()
  ome <- cmeans_setting("OME")

  selection <- cmeans_heatmap_features(
    pathway, ome,
    tissue = "muscle",
    clusters = cmeans_clusters("HEATMAP_MUSCLE_BLOOD_MIDDLE"),
    with_tissue = "blood",
    with_clusters = cmeans_clusters("VEGF_BLOOD")
  )
  message(sprintf("        ED5E: %d %s feature(s) of %s",
                  length(selection$feature_ids), ome, cmeans_pathway_label()))

  # ---- plot ----

  # The heatmap draws itself onto the device export_panel() opens rather than
  # returning an object; see cmeans_feature_heatmap() for why.
  draw_heatmap <- function() {
    cmeans_feature_heatmap(
      annotation = selection$annotation,
      feature_ids = selection$feature_ids,
      selected_tissue = CMEANS_TISSUES,
      selected_ome = ome,
      contrast_type = "exercise_with_controls",
      verbose = FALSE
    )
  }

  export_panel(draw_heatmap, "ED5E")
}

run_panels(list(
  ED5A = ed5a,
  ED5B = ed5b,
  ED5C = ed5c,
  ED5D = ed5d,
  ED5E = ed5e
))
