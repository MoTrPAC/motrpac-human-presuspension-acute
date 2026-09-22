#!/usr/bin/env Rscript
# Figure 4 — Clustering and PLIER latent variables
#
# Panels:  FIG4A  curated transcriptomics pathways enriched in the acute
#                 response (CAMERA-PR)
#          FIG4B  the same over global proteomics
#          FIG4C  fuzzy c-means cluster centroid trajectories, labelled
#          FIG4D  cluster ome composition and top pathway enrichment
#          FIG4E  cluster centroid trajectories, peak-spaced labels
#          FIG4F  VEGF signalling features in the late muscle cluster
#          FIG4G  cross-tissue RNA latent variables against the acute response
#          FIG4H  top pathways of those latent variables
#          FIG4I  TAMALIN transcript trajectory in adipose, blood and muscle
#
# FIG4G and FIG4H read the cross-tissue PLIER fit, which analysis/04_plier.R
# produces and no figure script refits. Until it has been run those two report
# SKIPPED with the command to run; the other seven panels build without it.
#
#   Rscript figures/landscape/analysis/04_plier.R   once, hours on a cold start
#   Rscript figures/landscape/FIG4.R                every panel
#   Rscript figures/landscape/FIG4.R FIG4D          one panel

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(dplyr)
  library(ggnewscale)
  library(ggplot2)
  library(gplots)
  library(patchwork)
  library(pheatmap)
  library(tidyr)
  library(MotrpacHumanPreSuspensionAnalysis)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(here, "..", "lib", "highlights.R"))
source(file.path(here, "plier_helpers.R"))
source(file.path(here, "..", "lib", "single_feature_helpers.R"))
source(file.path(here, "..", "single_feature_plots.R"))
source(file.path(here, "FIG4_ED5_helpers.R"))

# ---- shared ----------------------------------------------------------------

# Tissue order across the c-means panels. The same three names as the helper's
# CMEANS_TISSUES, which is the heatmap strip order; this is the panel layout
# order, and FIG4C, FIG4D and FIG4E all read left to right in it.
TISSUE_ORDER <- c("adipose", "blood", "muscle")

# FIG4A and FIG4B do NOT share a CAMERA load: they call camera_enrich_heatmap()
# with different omes and different curated selections, so there is one result
# per panel and nothing to memoise.
#
# The c-means frame IS shared, and is memoised in helpers/FIG4_ED5.R rather than
# here: cmeans_plot_frame() takes no arguments, and ED5's panels and the
# cross-tissue heatmaps' row annotation want the same frame.

# One plier_selection("rna") pass behind FIG4G and FIG4H. Lazy and memoised: the
# fit is hundreds of megabytes and is not read at all when neither panel is
# selected, nor when analysis/04_plier.R has not been run.
rna_plier <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- plier_selection("rna")
    }
    cache
  }
})

# ---- FIG4A — curated transcriptomics pathways ------------------------------

# CAMERA-PR over the transcriptomics differential analysis, drawn as a bubble
# heatmap: one row per pathway, one column per tissue x modality x timepoint
# contrast, colour the z-score and the bubble the significance. Only the
# exercise-with-controls contrasts are drawn, which is the delta-delta against
# the control arm rather than the raw post-versus-pre change.
#
# The twelve pathways are a hand selection, in config/highlights.json. The test
# itself is not recomputed here: CAMERA_RESULTS in the Analysis package is the
# published run over all three tissues, and this panel selects rows from it.
#
# Adipose, blood and muscle at once. camera_enrich_heatmap() keeps only sets
# tested in at least two tissues, so a pathway measured in one is dropped rather
# than drawn as a row of blanks.

fig4a <- function() {
  panel_init("FIG4A")

  # ---- the curated selection ----

  curated <- highlight_pathways("FIG4A")

  message(sprintf("        FIG4A: %d curated pathway(s) requested",
                  nrow(curated)))

  # ---- plot ----

  heatmap <- camera_enrich_heatmap(
    set_ids = curated$set_id,
    selected_ome = "transcript-rna-seq",
    selected_tissues = c("adipose", "blood", "muscle"),
    contrast_type = "exercise_with_controls"
  )

  message(sprintf("        FIG4A: drawing %d pathway(s) on a %.1f x %.1f in page",
                  heatmap$n_sets, heatmap$width, heatmap$height))

  export_panel(heatmap$draw, "FIG4A",
               width = heatmap$width, height = heatmap$height)
}

# ---- FIG4B — curated proteomics pathways -----------------------------------

# FIG4A over the global proteomics differential analysis instead of the
# transcriptomics one: the same CAMERA-PR results, the same bubble heatmap, the
# same exercise-with-controls contrasts, eight curated pathways rather than
# twelve.
#
# prot-pr is assayed in fewer tissues than transcriptomics, so this panel is
# narrower than FIG4A. The tissue set is not restricted here — the filter to
# sets tested in at least two tissues decides what is left, and which tissues
# those are is a property of the assay rather than of this figure.

fig4b <- function() {
  panel_init("FIG4B")

  # ---- the curated selection ----

  curated <- highlight_pathways("FIG4B")

  message(sprintf("        FIG4B: %d curated pathway(s) requested",
                  nrow(curated)))

  # ---- plot ----

  heatmap <- camera_enrich_heatmap(
    set_ids = curated$set_id,
    selected_ome = "prot-pr",
    selected_tissues = c("adipose", "blood", "muscle"),
    contrast_type = "exercise_with_controls"
  )

  message(sprintf("        FIG4B: drawing %d pathway(s) on a %.1f x %.1f in page",
                  heatmap$n_sets, heatmap$width, heatmap$height))

  export_panel(heatmap$draw, "FIG4B",
               width = heatmap$width, height = heatmap$height)
}

# ---- FIG4C — cluster centroid trajectories, labelled -----------------------

# Every cluster in every tissue, one sparkline each, grouped into columns by the
# timepoint that cluster peaks at. Reading across the panel is reading when a
# response happens: the early column, then the middle, then late, within each
# tissue in turn.
#
# Two header strips run above — tissue, then peak timepoint — coloured from the
# package palettes so the column a sparkline sits in is legible without axis
# text. The sparklines themselves carry none: at half an inch per column the
# trajectories are shape, not measurement, and every panel shares one y scale so
# heights are comparable across the figure.
#
# Unlike FIG4E this labels EVERY cluster. Where the spread runs out of
# timepoints the leftovers fall back to their own peak and share one, which is
# the price of naming them all.

fig4c <- function() {
  panel_init("FIG4C")

  # ---- data ----

  centroids <- cmeans_plot_frame()

  # The columns: one per tissue x peak timepoint that any cluster actually falls
  # into. A peak group with no cluster in it draws no column rather than an empty
  # one, which is why this comes from the data and the ORDER comes from here.
  PEAK_ORDER <- c("P15-45M", "P3.5/4H", "P24H")

  groups <- unique(centroids[, c("tissue", "peak_timepoint")])
  groups <- groups[groups$tissue %in% TISSUE_ORDER &
                     groups$peak_timepoint %in% PEAK_ORDER, , drop = FALSE]
  groups <- groups[order(match(groups$tissue, TISSUE_ORDER),
                         match(groups$peak_timepoint, PEAK_ORDER)), , drop = FALSE]
  groups$column <- seq_len(nrow(groups))

  if (nrow(groups) == 0) {
    stop("no cluster peaks at any of the declared timepoints: ",
         paste(PEAK_ORDER, collapse = ", "), call. = FALSE)
  }

  # One y scale for the whole figure, headroom included, so a tall trajectory in
  # muscle and a flat one in adipose are drawn on the same ruler.
  LABEL_OFFSET <- 0.3
  y_limits <- c(min(centroids$value, na.rm = TRUE),
                max(centroids$value + LABEL_OFFSET, na.rm = TRUE))

  # ---- one column ----

  sparkline <- function(tissue, peak, first) {
    rows <- centroids[centroids$tissue == tissue &
                        centroids$peak_timepoint == peak, , drop = FALSE]
    rows$cluster <- as.character(rows$cluster)
    # See FIG4E: the axis factor spans every tissue, but only blood has P10M.
    rows$timepoint_short <- droplevels(rows$timepoint_short)

    # Ranked on distance from zero in either direction: a column of strongly
    # negative clusters would otherwise place every label at whichever timepoint
    # happened to be least negative.
    placed <- cmeans_label_positions(rows, fallback = TRUE, score = "abs")
    placed <- cmeans_apply_label_overrides(placed, tissue)

    labels <- merge(rows, placed, by = c("cluster", "timepoint_short"))
    labels$y <- labels$value + LABEL_OFFSET

    colors <- unique(rows[, c("cluster", "cluster_color")])
    colors <- stats::setNames(colors$cluster_color, colors$cluster)

    p <- ggplot(rows, aes(x = .data$timepoint_short, y = .data$value,
                          group = .data$cluster)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
                 linewidth = 0.25) +
      geom_line(color = "black", linewidth = 0.85, alpha = 0.7) +
      geom_line(aes(color = .data$cluster), linewidth = 0.75, alpha = 0.9) +
      geom_point(data = labels, aes(y = .data$y, fill = .data$cluster),
                 shape = 21, color = "black", size = 1.2, stroke = 0.2) +
      geom_text(data = labels, aes(y = .data$y, label = .data$cluster),
                color = "black", size = 0.8, fontface = "bold") +
      scale_color_manual(values = colors, guide = "none") +
      scale_fill_manual(values = colors, guide = "none") +
      scale_y_continuous(limits = y_limits,
                         expand = expansion(mult = c(0.06, 0.14))) +
      coord_cartesian(clip = "off") +
      labs(x = NULL, y = NULL) +
      theme_bw(base_size = 8) +
      theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0, "lines"),
        plot.margin = margin(1, 1, 1, 1, "pt")
      )

    # Only the leftmost column names the arms.
    if (first) {
      p + facet_grid(arm ~ ., switch = "y") +
        theme(strip.text.y = element_text(size = 5, angle = 0),
              strip.background = element_blank())
    } else {
      p + facet_grid(arm ~ .) + theme(strip.text.y = element_blank())
    }
  }

  columns <- lapply(seq_len(nrow(groups)), function(i) {
    sparkline(groups$tissue[i], groups$peak_timepoint[i], first = i == 1L)
  })

  # ---- header strips ----

  # Blood's column reads "P10-45M" because its P10M peaks were folded in; the
  # other tissues have no P10M to fold. Label only — the grouping already merged.
  groups$display <- ifelse(groups$tissue == "blood" &
                             groups$peak_timepoint == "P15-45M",
                           "P10-45M", groups$peak_timepoint)
  groups$timepoint_key <- names(CMEANS_TIMEPOINTS)[
    match(groups$peak_timepoint, unname(CMEANS_TIMEPOINTS))]

  timepoint_colors <- MotrpacHumanPreSuspensionAnalysis::HUMAN_ACUTE_TIMEPOINT_COLORS[
    unique(groups$timepoint_key)]
  if (anyNA(timepoint_colors)) {
    stop("timepoint not in HUMAN_ACUTE_TIMEPOINT_COLORS: ",
         paste(unique(groups$timepoint_key)[is.na(timepoint_colors)],
               collapse = ", "), call. = FALSE)
  }

  strip_x <- scale_x_continuous(limits = c(0.5, nrow(groups) + 0.5),
                                expand = c(0, 0))

  peak_strip <- ggplot(groups) +
    geom_rect(aes(xmin = .data$column - 0.5, xmax = .data$column + 0.5,
                  ymin = 0, ymax = 1, fill = .data$timepoint_key), color = NA) +
    geom_text(aes(x = .data$column, y = 0.5, label = .data$display),
              color = "black", fontface = "bold", size = 1.3) +
    scale_fill_manual(values = timepoint_colors, guide = "none") +
    strip_x + scale_y_continuous(expand = c(0, 0)) + theme_void()

  tissue_spans <- do.call(rbind, lapply(split(groups, groups$tissue), function(g) {
    data.frame(tissue = g$tissue[1], centre = mean(g$column),
               xmin = min(g$column) - 0.5, xmax = max(g$column) + 0.5,
               stringsAsFactors = FALSE)
  }))

  tissue_strip <- ggplot(tissue_spans) +
    geom_rect(aes(xmin = .data$xmin, xmax = .data$xmax, ymin = 0, ymax = 1,
                  fill = .data$tissue), color = NA) +
    geom_text(aes(x = .data$centre, y = 0.5, label = .data$tissue),
              color = "black", fontface = "bold", size = 1.5) +
    scale_fill_manual(values = cmeans_tissue_colors(TISSUE_ORDER), guide = "none") +
    strip_x + scale_y_continuous(expand = c(0, 0)) + theme_void()

  # ---- plot ----

  combined <- tissue_strip / peak_strip / patchwork::wrap_plots(columns, nrow = 1) +
    patchwork::plot_layout(heights = c(0.1, 0.08, 1))

  export_panel(function() print(combined), "FIG4C")
}

# ---- FIG4D — cluster ome composition and top pathway enrichment ------------

# Three registers over one shared x axis of tissue x cluster: the cluster number,
# then how many features of each ome fall in it, then which curated pathways it
# is enriched for and how much of each pathway it accounts for.
#
# Clusters are ordered tissue by tissue and, within a tissue, by when they peak
# and in which direction — the same grouping that colours FIG4C. Reading left to
# right is reading the response in time.
#
# The pathway rows are ordered by clustering the enrichment matrix, so pathways
# that light up in the same clusters sit together. Adipose is excluded from that
# clustering: it enriches for far fewer sets than the other two and would
# otherwise dominate the distance.

fig4d <- function() {
  panel_init("FIG4D")

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

  # One pathway, the cmeans selection in config/highlights.json. The panel is
  # about where the VEGF signature lands across the clusters; ED5B draws the
  # fifteen curated sets.
  pathway <- cmeans_pathway()
  pathway_label <- cmeans_pathway_label()

  hits <- cmeans_ora_hits()
  hits$label <- paste(as.character(as.numeric(hits$cluster)), hits$tissue,
                      sep = " - ")

  enrichment <- hits |>
    filter(.data$set == pathway | as.character(.data$set_id) == pathway,
           .data$label %in% levels(cluster_meta$label)) |>
    left_join(ome_text, by = c("assay" = "assay_code"),
              relationship = "many-to-many") |>
    mutate(omics_text = if_else(.data$assay == "metab", "Metabolomics",
                                .data$omics_text)) |>
    # The share of the cluster that belongs to the pathway.
    mutate(included_ratio = .data$set_size_in_cluster / .data$cluster_size) |>
    distinct(.data$label, .data$omics_text, .data$included_ratio) |>
    mutate(label = factor(.data$label, levels = levels(cluster_meta$label)),
           set = pathway_label)

  # Dots are sized on ED5B's scale, the range over every curated pathway, so the
  # VEGF dots here are the size they are in the full grid.
  curated_sets <- highlight_pathways("ED5B")$set
  size_limits <- range(hits$set_size_in_cluster[hits$set %in% curated_sets] /
                         hits$cluster_size[hits$set %in% curated_sets])

  if (nrow(enrichment) == 0) {
    stop("'", pathway, "' is not enriched in any cluster.",
         "\n  The cmeans selection in config/highlights.json names it.", call. = FALSE)
  }
  message(sprintf("        FIG4D: %s enriched in %d of %d clusters",
                  pathway_label, length(unique(enrichment$label)),
                  nlevels(cluster_meta$label)))

  # One row, so there is nothing to order and nothing to cluster. The legacy ran
  # a hierarchical clustering here to group fifteen pathways that behaved alike.
  #
  # ED5B colours those fifteen by position in that clustering, recycling
  # CMEANS_SET_PALETTE down the rows. VEGF lands second and is #006400. One row
  # cannot compute a position in an ordering it does not build, so the hue is a
  # literal here; re-derive it after a re-clustering by building ED5B and reading
  # the VEGF row.
  VEGF_ROW_COLOR <- "#006400"

  set_colors <- stats::setNames(VEGF_ROW_COLOR, pathway_label)
  set_order <- pathway_label

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
    scale_y_discrete(expand = expansion(mult = c(0.6, 0.6))) +
    scale_size_continuous(range = c(1, POINT_MAX), limits = size_limits) +
    # No y title. With one pathway the row label already names it, and patchwork
    # stacks this title into the same gutter as "Number of Features" above, where
    # the two overprint.
    labs(x = NULL, y = NULL, size = "Included Ratio") +
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
    # 1.8 for the dot register was sized for fifteen pathway rows. One row needs a
    # sliver, and giving it more just floats the dots in white space.
    patchwork::plot_layout(heights = c(0.2, 1, 0.35))

  export_panel(function() print(combined), "FIG4D")
}

# ---- FIG4E — cluster centroid trajectories, peak-spaced labels -------------

# One sparkline per selected cluster, faceted exercise arm by tissue. Each
# cluster carries its number once, placed at the timepoint where that cluster
# peaks — and no two clusters in a tissue take the same timepoint, so the labels
# spread along the x axis instead of stacking on whichever timepoint happens to
# be tallest.
#
# The panel is deliberately tiny: half an inch per tissue, drawn to sit inside
# Figure 4 rather than to be read on its own. Axis text is suppressed for the
# same reason; the trajectories are shape, not measurement.
#
# Which clusters appear comes from figure_4/cmeans.env (CMEANS_SPARKLINE_<TISSUE>),
# not from the enrichment. See helpers/FIG4_ED5.R.

fig4e <- function() {
  panel_init("FIG4E")

  # ---- data ----

  centroids <- cmeans_plot_frame()

  # The clusters drawn in each tissue.
  selections <- lapply(TISSUE_ORDER, function(tissue) {
    cmeans_clusters(paste0("SPARKLINE_", toupper(tissue)))
  })
  names(selections) <- TISSUE_ORDER

  for (tissue in TISSUE_ORDER) {
    have <- unique(centroids$cluster[centroids$tissue == tissue])
    missing <- setdiff(selections[[tissue]], have)
    if (length(missing) > 0) {
      stop(tissue, ": CMEANS_SPARKLINE_", toupper(tissue), " names cluster(s) the ",
           "clustering does not have: ", paste(missing, collapse = ", "),
           "\n  ", tissue, " is clustered at k = ", length(have),
           ". Re-clustering renumbers these; see figure_4/cmeans.env.", call. = FALSE)
    }
  }

  # ---- plot ----

  sparkline <- function(tissue, first) {
    rows <- centroids[centroids$tissue == tissue &
                        centroids$cluster %in% selections[[tissue]], , drop = FALSE]
    rows$cluster <- factor(rows$cluster, levels = selections[[tissue]])
    # droplevels: the axis factor is built across ALL tissues so the levels are in
    # one order, but only blood has P10M. Left in, adipose and muscle draw an
    # empty leading slot and their trajectories start a quarter of the way across
    # the panel — a timepoint that looks missing rather than absent.
    rows$timepoint_short <- droplevels(rows$timepoint_short)

    # The label sits slightly above the trajectory it names.
    # No fallback: FIG4E keeps one label per timepoint and reports what that
    # leaves out. FIG4C takes the same placement WITH a fallback.
    placed <- cmeans_label_positions(rows, fallback = FALSE, score = "signed")
    # The algorithm retires a timepoint once used, so it can place at most one
    # label per timepoint: a tissue with more selected clusters than timepoints
    # leaves the weakest unlabelled. That is the legacy's behaviour and the price
    # of the spread; it is reported rather than left to be noticed.
    unlabelled <- setdiff(selections[[tissue]], placed$cluster)
    if (length(unlabelled) > 0) {
      message(sprintf("        FIG4E: %s cluster(s) %s unlabelled - only %d timepoint(s) to place %d cluster(s)",
                      tissue, paste(unlabelled, collapse = ", "),
                      length(unique(rows$timepoint_short)),
                      length(selections[[tissue]])))
    }
    labels <- merge(rows, placed, by = c("cluster", "timepoint_short"))
    labels$y <- labels$value + 0.3

    colors <- unique(rows[, c("cluster", "cluster_color")])
    colors <- stats::setNames(colors$cluster_color, as.character(colors$cluster))

    strip_fill <- unname(
      MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS[tissue]
    )

    p <- ggplot(rows, aes(x = .data$timepoint_short, y = .data$value,
                          group = .data$cluster)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
                 linewidth = 0.25) +
      # Black underlay then the cluster colour, so a pale trajectory still reads
      # against the panel.
      geom_line(color = "black", linewidth = 0.85, alpha = 0.7) +
      geom_line(aes(color = .data$cluster), linewidth = 0.75, alpha = 0.9) +
      geom_point(data = labels, aes(y = .data$y, fill = .data$cluster),
                 shape = 21, color = "black", size = 1.2, stroke = 0.2) +
      geom_text(data = labels, aes(y = .data$y, label = .data$cluster),
                color = "black", size = 0.8, fontface = "bold") +
      scale_color_manual(values = colors, guide = "none") +
      scale_fill_manual(values = colors, guide = "none") +
      coord_cartesian(ylim = c(-0.2, 2.55), clip = "off") +
      labs(x = NULL, y = NULL) +
      theme_bw(base_size = 8) +
      theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0.01, "cm"),
        plot.margin = margin(0, 5, 0, 0),
        strip.background.x = element_rect(fill = strip_fill, color = NA),
        strip.text.x = element_text(face = "bold", size = 2,
                                    margin = margin(0.01, 0, 0.01, 0, "cm"))
      )

    # Only the leftmost tissue carries the arm labels down the side; repeating
    # them on all three would triple the width of a panel this narrow.
    if (first) {
      p + facet_grid(arm ~ tissue, switch = "y") +
        theme(strip.text.y = element_text(size = 5, angle = 180),
              strip.background.y = element_blank())
    } else {
      p + facet_grid(arm ~ tissue) +
        theme(strip.text.y = element_blank(),
              strip.background.y = element_blank())
    }
  }

  panels <- lapply(seq_along(TISSUE_ORDER), function(i) {
    sparkline(TISSUE_ORDER[i], first = i == 1L)
  })

  export_panel(function() print(patchwork::wrap_plots(panels, nrow = 1)), "FIG4E")
}

# ---- FIG4F — VEGF signalling features in the late muscle cluster -----------

# One row per feature, one column per exercise-with-controls contrast in every
# tissue, adipose, blood and muscle. Rows are the features measured in all
# three, clustered. Three strips on the right name each feature's cluster per
# tissue in the per-cluster colours FIG4C and FIG4D use.
#
# No second tissue: a feature is here because it lands in the named muscle
# cluster and for no other reason.
#
# The pathway comes from config/highlights.json; the ome and clusters from
# figure_4/cmeans.env.
#
# ComplexHeatmap is attached, not qualified, at the top of this script: the
# Analysis package's cell-drawing function calls textGrob() unqualified, and
# only finds it because attaching ComplexHeatmap attaches grid. Without that the
# panel dies with "could not find function textGrob".

fig4f <- function() {
  panel_init("FIG4F")

  # ---- data ----

  pathway <- cmeans_pathway()
  ome <- cmeans_setting("OME")

  selection <- cmeans_heatmap_features(
    pathway, ome,
    tissue = "muscle",
    clusters = cmeans_clusters("HEATMAP_MUSCLE_LATE")
  )
  message(sprintf("        FIG4F: %d %s feature(s) of %s",
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

  export_panel(draw_heatmap, "FIG4F")
}

# ---- FIG4G — cross-tissue RNA latent variables -----------------------------

# One row per selected latent variable, one column per tissue x arm x timepoint.
# A bubble's colour is the median difference between that arm and the controls in
# the LV's within-participant change from baseline; its size is -log10 of the
# adjusted p of that comparison, capped at 3; it is outlined where the adjusted p
# clears 0.05.
#
# Rows run in the cross-tissue agreement order every LV panel of these two
# figures uses — ascending blood-vs-muscle correlation — so this panel, ED6B and
# ED6C read left to right the same way.
#
# The response grid itself is analysis/04_plier.R's, not this panel's: ED6E draws
# the same grid over the metabolomics fit, and one computation means the two
# cannot disagree about what a comparison is.

fig4g <- function() {
  skip <- skip_if_no_fit("FIG4G")
  if (!is.null(skip)) return(skip)

  panel_init("FIG4G")

  selection <- rna_plier()
  plier_report_selection(selection)

  response <- selection$response
  rows <- selection$lv_names

  bubble <- plier_bubble_heatmap(
    effect = response$effect[rows, , drop = FALSE],
    size = response$neglog10_padj[rows, , drop = FALSE],
    comparisons = response$comparisons,
    cell_border = "black",
    radius_divisor = 1
  )

  draw <- function() {
    ComplexHeatmap::draw(bubble$heatmap,
                         heatmap_legend_list = bubble$legends,
                         merge_legend = FALSE)
  }

  export_panel(draw, "FIG4G")
}

# ---- FIG4H — what those latent variables are made of -----------------------

# Rows are the union of the strongest pathways of each selected latent variable,
# columns those latent variables in the same order FIG4G draws them, and a cell
# is PLIER's AUC for that pathway in that LV.
#
# Ranked on U, coloured by Uauc, which is what the legacy did. The two are not
# the same order: U is a prior set's loading on the LV, Uauc the AUC of that
# association, and where an LV is associated with more than PLIER_PATHWAYS_PER_LV
# sets the choice decides which make the cut. On this fit that is LV13 and LV33,
# and it moves five rows.
#
# ED6Aii ranks the same pathways by Uauc. That split is the legacy's — its bar
# plots read Uauc and this heatmap read U — and it is reproduced rather than
# unified.
#
# Two things the legacy stated by hand are derived here:
#
#   - How many pathways each LV contributes. It took 5 from most LVs, 4 from one
#     and 1 from two others, as literal `[1:5]` / `[1:4]` / `[1:1]` subscripts.
#     Those are min(5, sets the LV is associated with at all): LV15 has 4, LV54
#     and LV84 have 1. Dropping the unassociated entries is that rule stated once
#     rather than three numbers a refit invalidates.
#   - The row labels. They were a vector of 23 display strings positionally
#     matched to the union, so a refit that reordered one pathway relabelled
#     every row below it. They come off the pathway ids now.

fig4h <- function() {
  skip <- skip_if_no_fit("FIG4H")
  if (!is.null(skip)) return(skip)

  panel_init("FIG4H")

  PATHWAYS_PER_LV <- plier_env_int("PLIER_PATHWAYS_PER_LV")

  selection <- rna_plier()
  plier_report_selection(selection)

  loading <- selection$fit$U
  auc <- selection$fit$Uauc
  columns <- plier_lv_index(selection$lv_names)

  top_pathways <- unique(unlist(lapply(columns, function(i) {
    ranked <- rownames(loading)[order(-loading[, i])][seq_len(PATHWAYS_PER_LV)]
    ranked[loading[ranked, i] > 0]
  })))
  message(sprintf("        %d pathways, from the top %d by loading of each of %d latent variables",
                  length(top_pathways), PATHWAYS_PER_LV, length(columns)))

  auc_mat <- auc[top_pathways, columns, drop = FALSE]
  colnames(auc_mat) <- selection$lv_names

  export_panel(
    plier_auc_heatmap(auc_mat, plier_pathway_label(rownames(auc_mat))),
    "FIG4H"
  )
}

# ---- FIG4I — TAMALIN transcript trajectories -------------------------------

# The same transcript in three tissues, one sub-plot each, under one collected
# legend. Exercise-group means with 95% confidence intervals across the acute
# timepoints.
#
# The plot is one entry of the single-feature catalog in single_feature_plots.R.

fig4i <- function() {
  panel_init("FIG4I")
  export_panel(single_feature_plot("FIG4I"), "FIG4I")
}

run_panels(list(
  FIG4A = fig4a,
  FIG4B = fig4b,
  FIG4C = fig4c,
  FIG4D = fig4d,
  FIG4E = fig4e,
  FIG4F = fig4f,
  FIG4G = fig4g,
  FIG4H = fig4h,
  FIG4I = fig4i
))
