#!/usr/bin/env Rscript
# Extended Data 6 — cross-tissue PLIER latent variables
#
# Panels:  ED6A   each selected RNA latent variable across every sample
#          ED6Aii top pathways of the three latent variables ED6A highlights
#          ED6B   cross-tissue agreement of each selected RNA latent variable
#          ED6C   transcripts belonging to each selected RNA latent variable
#          ED6D   strongest transcripts of one RNA latent variable
#          ED6E   metabolomics latent variables against the acute response
#          ED6F   RefMet classes of the selected metabolomics latent variables
#          ED6G   cross-tissue agreement of each selected metabolomics LV
#          ED6H   metabolites belonging to each selected metabolomics LV
#          ED6I   strongest metabolites of one metabolomics latent variable
#
# ED6Aii is a manufactured id, not a manuscript panel letter: the legacy drew
# those bar plots as three separate graphics to be placed beside the ED6A
# schematic, and ED6's letters are all spoken for. FIG2Bii is the precedent.
#
# Every panel here reads a cross-tissue PLIER fit, which analysis/04_plier.R
# produces and no figure script refits. Until it has been run every panel
# reports SKIPPED with the command to run:
#
#   Rscript figures/landscape/analysis/04_plier.R   once, hours on a cold start
#
# Which latent variables each panel draws, the affinity cutoffs and the two
# marker-heatmap LVs are figure_4/plier_lvs.env's. A refit renumbers latent
# variables; that file's header says how to re-decide them.
#
#   Rscript figures/landscape/ED6.R          every panel
#   Rscript figures/landscape/ED6.R ED6E     one panel

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(ggplot2)
  library(gplots)
  library(grid)
  library(patchwork)
  library(pheatmap)
  library(MotrpacHumanPreSuspensionAnalysis)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(here, "plier_helpers.R"))

# ---- shared ----------------------------------------------------------------

# Two accessors, not one. The RNA and the metabolomics halves of this figure are
# two different PLIER fits over two different matrices with their own latent
# variables, so there is nothing to share between them; within each half the five
# panels call plier_selection() with the same argument and share everything.
#
# Lazy and memoised: a fit is hundreds of megabytes and is not read at all when
# only the other half's panels are selected, nor when analysis/04_plier.R has not
# been run.
rna_plier <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- plier_selection("rna")
    }
    cache
  }
})

metab_plier <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- plier_selection("metab")
    }
    cache
  }
})

# ---- the three graphics that exist once per ome ----------------------------

# ED6B/ED6G, ED6C/ED6H and ED6D/ED6I are the same three graphics over the two
# PLIER arms: the cross-tissue agreement of each selected latent variable, how
# many features belong to it, and one LV's top features across every sample.
#
# Six panels, three functions. What differs between the arms — the affinity
# cutoff, which LV gets a heatmap, how a feature is labelled — is an argument or
# comes from figure_4/plier_lvs.env, so nothing here knows which ome it is drawing.
#
# They live in this script rather than in lib/ because only ED6 draws them;
# lib/plier_helpers.R holds what ED6 and Figure 4 share.

# The three tissue-pair colours. Okabe-Ito orange, sky blue and reddish purple,
# which is what the legacy hardcoded in both scripts. They are pair colours, not
# tissue colours, so HUMAN_TISSUE_COLORS is not what they should be: no one
# tissue owns "Adipose vs Blood".
ED6_PAIR_COLORS <- c(
  "Adipose vs Blood"  = "#E69F00",
  "Adipose vs Muscle" = "#56B4E9",
  "Blood vs Muscle"   = "#CC79A7"
)

#' ED6B / ED6G — each selected LV's correlation between each pair of tissues.
ed6_correlation_plot <- function(selection) {
  correlation <- selection$correlation
  keep <- correlation$LV %in% selection$lv_names
  long <- data.frame(
    LV = rep(correlation$LV[keep], 3),
    Comparison = rep(names(ED6_PAIR_COLORS), each = sum(keep)),
    Correlation = c(correlation$Adipose_v_Blood[keep],
                    correlation$Adipose_v_Muscle[keep],
                    correlation$Blood_v_Muscle[keep])
  )
  long$LV <- factor(long$LV, levels = selection$lv_names)

  ggplot(long) +
    geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
    geom_point(aes(x = .data$LV, y = .data$Correlation, color = .data$Comparison),
               size = 3) +
    scale_color_manual(values = ED6_PAIR_COLORS) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      legend.position = "top",
      text = element_text(size = 10),
      legend.key.spacing = unit(0, "cm")
    ) +
    guides(col = guide_legend(ncol = 1))
}

#' ED6C / ED6H — how many features are markers of each selected LV.
#'
#' The legacy built a frame with one row per marker and let geom_bar count them.
#' Same bars off the counts themselves, which is also what the panel reports.
ed6_membership_plot <- function(selection) {
  counts <- vapply(selection$markers[selection$lv_names], length, integer(1))
  df <- data.frame(LV = factor(names(counts), levels = selection$lv_names),
                   Membership = as.integer(counts))

  ggplot(df, aes(x = .data$LV, y = .data$Membership)) +
    geom_col() +
    theme_classic() +
    ylab("Membership") +
    theme(axis.title.x = element_blank(), text = element_text(size = 10))
}

#' ED6D / ED6I — one LV's top-affinity features across every sample.
#'
#' Columns are samples in the analysis/04_plier.R order (tissue, group,
#' timepoint, sex) and are annotated with all four; rows are the `n` features
#' with the highest affinity to `lv`, in affinity order.
ed6_marker_heatmap <- function(selection, lv, n, labels = NULL,
                               legend_margin_mm = 10,
                               colour_title = "Z-Score") {
  features <- plier_top_features(selection$fit, lv, n)
  mat <- selection$input$matrix[features, , drop = FALSE]
  annotation <- selection$input$metadata[, c("Sex", "Timepoint", "Modality", "Tissue")]

  row_labels <- if (is.null(labels)) rownames(mat) else labels[features]
  row_labels[is.na(row_labels)] <- features[is.na(row_labels)]

  function() {
    heatmap <- pheatmap::pheatmap(
      mat,
      breaks = seq(-5, 5, length.out = 100),
      color = gplots::colorpanel(101, "blue", "white", "red"),
      show_colnames = FALSE,
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      annotation_col = annotation,
      annotation_colors = plier_annotation_colors(),
      annotation_legend = FALSE,
      fontsize = 14,
      labels_row = row_labels,
      silent = TRUE
    )

    margin <- grid::unit(legend_margin_mm, "mm")
    legend <- ed6_annotation_legend(
      plier_annotation_colors(),
      rows = list(c("Tissue", "Modality", "Sex"), "Timepoint"),
      max_width_in = grid::convertWidth(
        grid::unit(1, "npc") - margin * 2, "in", valueOnly = TRUE)
    )
    # Split the page rather than appending a row to pheatmap's gtable: that
    # gtable's heights already sum to the full device, so an added row lands
    # below the page edge and the legend's second line is cut off.
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(
      2, 1, heights = grid::unit.c(grid::unit(1, "null"), legend$height))))

    # The heatmap, held short of the top so the colour bar has room for a title.
    grid::pushViewport(grid::viewport(layout.pos.row = 1))
    grid::pushViewport(grid::viewport(
      y = grid::unit(0, "npc"), just = "bottom",
      height = grid::unit(1, "npc") - grid::unit(5, "mm")))
    grid::grid.draw(heatmap$gtable)
    legend_left <- heatmap_legend_left(heatmap$gtable)
    grid::popViewport()
    draw_heatmap_legend_title(legend_left, colour_title, fontsize = 10)
    grid::popViewport()

    # The margin has to go on a viewport of its own: grid ignores x/width when
    # layout.pos.row places the viewport in a layout cell.
    grid::pushViewport(grid::viewport(layout.pos.row = 2))
    grid::pushViewport(grid::viewport(
      x = margin, just = "left",
      width = grid::unit(1, "npc") - margin * 2))
    grid::grid.draw(legend$grob)
    grid::popViewport(3)
  }
}

#' The annotation key ED6D and ED6I carry under the heatmap.
#'
#' pheatmap's own annotation legend is a vertical column on the right, one
#' variable above the next, which costs the width these two panels do not have —
#' 1,417 and 1,554 columns across 8 inches. This is the same key laid out along
#' the bottom instead: swatch, label, next, wrapping to a second line when the
#' row runs out of width.
#'
#' `colors` is plier_annotation_colors(); `rows` names which variables share a
#' line, in reading order. Timepoint has seven levels and gets a line of its own
#' so its scale reads as one sequence. A line that still overruns wraps.
#' Returns the grob and the height to reserve for it.
ed6_annotation_legend <- function(colors, rows = list(names(colors)),
                                  max_width_in = 7.5, fontsize = 9) {
  SWATCH_MM <- 3.4   # the key square
  PAD_MM <- 1.2      # square to its label
  ITEM_GAP_MM <- 4   # label to the next square
  TITLE_GAP_MM <- 1.8
  ROW_MM <- 5.4

  text_width <- function(label, bold = FALSE) {
    grid::convertWidth(
      grid::grobWidth(grid::textGrob(
        label, gp = grid::gpar(fontsize = fontsize,
                               fontface = if (bold) "bold" else "plain"))),
      "mm", valueOnly = TRUE
    )
  }

  entries_for <- function(variables) {
    entries <- list()
    for (variable in variables) {
      title <- paste0(variable, ":")
      entries[[length(entries) + 1L]] <- list(
        kind = "title", label = title,
        width = text_width(title, bold = TRUE) + TITLE_GAP_MM
      )
      keys <- colors[[variable]]
      for (level in names(keys)) {
        entries[[length(entries) + 1L]] <- list(
          kind = "key", label = level, fill = unname(keys[level]),
          width = SWATCH_MM + PAD_MM + text_width(level) + ITEM_GAP_MM
        )
      }
    }
    entries
  }

  # One line per declared row, wrapping only if a line does not fit.
  limit <- grid::convertWidth(grid::unit(max_width_in, "in"), "mm", valueOnly = TRUE)
  laid_out <- list()
  for (variables in rows) {
    line <- list()
    x <- 0
    for (entry in entries_for(variables)) {
      if (length(line) > 0L && x + entry$width > limit) {
        laid_out[[length(laid_out) + 1L]] <- line
        line <- list()
        x <- 0
      }
      entry$x <- x
      line[[length(line) + 1L]] <- entry
      x <- x + entry$width
    }
    if (length(line) > 0L) laid_out[[length(laid_out) + 1L]] <- line
  }
  rows <- laid_out

  children <- grid::gList()
  for (i in seq_along(rows)) {
    y <- grid::unit(1, "npc") - grid::unit((i - 0.5) * ROW_MM, "mm")
    for (entry in rows[[i]]) {
      children <- if (identical(entry$kind, "title")) {
        grid::gList(children, grid::textGrob(
          entry$label, x = grid::unit(entry$x, "mm"), y = y,
          just = c("left", "centre"),
          gp = grid::gpar(fontsize = fontsize, fontface = "bold")))
      } else {
        grid::gList(
          children,
          grid::rectGrob(
            x = grid::unit(entry$x, "mm"), y = y,
            width = grid::unit(SWATCH_MM, "mm"), height = grid::unit(SWATCH_MM, "mm"),
            just = c("left", "centre"),
            gp = grid::gpar(fill = entry$fill, col = "grey30", lwd = 0.5)),
          grid::textGrob(
            entry$label, x = grid::unit(entry$x + SWATCH_MM + PAD_MM, "mm"), y = y,
            just = c("left", "centre"), gp = grid::gpar(fontsize = fontsize)))
      }
    }
  }

  list(grob = grid::gTree(children = children),
       height = grid::unit(length(rows) * ROW_MM + 3, "mm"))
}

#' The gene symbol of each transcript, for ED6D's row labels.
ed6_transcript_symbols <- function() {
  f2g <- as.data.frame(MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE)
  f2g <- f2g[f2g$assay == "transcript-rna-seq", ]
  setNames(as.character(f2g$gene_symbol), as.character(f2g$feature_id))
}

# ---- ED6A — each selected RNA latent variable across every sample ----------

# One line per selected latent variable: its value in every sample the fit was
# run on, in the sample order analysis/04_plier.R wrote — adipose, then blood,
# then muscle, and within each tissue by group, timepoint and sex. What the panel
# shows is which tissue block a latent variable lives in and how flat it is
# elsewhere.
#
# The legacy drew these as eight separate 6 x 5 graphics, to be arranged by hand
# into the schematic that is panel A. They are one panel here, faceted in the
# same LV order the rest of the figure uses — the FIG5G precedent: what the
# manuscript places should be what the script wrote.

ed6a <- function() {
  skip <- skip_if_no_fit("ED6A")
  if (!is.null(skip)) return(skip)

  panel_init("ED6A")

  selection <- rna_plier()
  plier_report_selection(selection)

  B <- selection$fit$B
  rownames(B) <- plier_lv_names(nrow(B))

  df <- do.call(rbind, lapply(selection$lv_names, function(lv) {
    data.frame(
      LV = lv,
      Sample = seq_len(ncol(B)),
      Value = B[lv, ]
    )
  }))
  df$LV <- factor(df$LV, levels = selection$lv_names)

  p <- ggplot(df, aes(x = .data$Sample, y = .data$Value)) +
    geom_line(linewidth = 0.3) +
    geom_point(size = 0.3) +
    facet_wrap(~ LV, ncol = 2, scales = "free_y") +
    labs(x = "Sample", y = "Relative Expression") +
    theme_classic() +
    theme(text = element_text(size = 10), strip.background = element_blank())

  export_panel(p, "ED6A")
}

# ---- ED6Aii — the strongest pathways of the LVs ED6A highlights ------------

# One horizontal bar plot per LV in PLIER_RNA_PATHWAY_BARPLOT_LVS, bars being
# PLIER's AUC for that LV's top PLIER_PATHWAYS_PER_LV pathways. The legacy drew
# them as three separate 8 x 3 graphics to be placed beside the ED6A schematic,
# the same way it drew ED6A's line plots; they are one panel here for the same
# reason.
#
# Labels come off the pathway ids. The legacy hand-wrote five display strings per
# plot and matched them to the top-5 rows positionally, so a refit that moved one
# pathway relabelled the rest.

ed6aii <- function() {
  skip <- skip_if_no_fit("ED6Aii")
  if (!is.null(skip)) return(skip)

  panel_init("ED6Aii")

  LVS <- plier_env_ints("PLIER_RNA_PATHWAY_BARPLOT_LVS")
  PATHWAYS_PER_LV <- plier_env_int("PLIER_PATHWAYS_PER_LV")

  selection <- rna_plier()
  plier_check_lvs(selection$fit, LVS, "rna")

  auc <- selection$fit$Uauc

  one_lv <- function(lv) {
    ranked <- rownames(auc)[order(-auc[, lv])][seq_len(PATHWAYS_PER_LV)]
    ranked <- ranked[auc[ranked, lv] > 0]
    df <- data.frame(
      Pathway = factor(plier_pathway_label(ranked), levels = rev(plier_pathway_label(ranked))),
      Uauc = auc[ranked, lv]
    )
    message(sprintf("        LV%d: %d pathways", lv, nrow(df)))
    ggplot(df, aes(x = .data$Pathway, y = .data$Uauc, fill = .data$Pathway)) +
      geom_col() +
      coord_flip() +
      scale_fill_brewer(palette = "Reds", direction = -1) +
      labs(title = paste0("LV", lv), y = "Uauc") +
      theme_classic() +
      theme(legend.position = "none", axis.title.y = element_blank(),
            text = element_text(size = 12))
  }

  export_panel(Reduce(`/`, lapply(LVS, one_lv)), "ED6Aii")
}

# ---- ED6B — cross-tissue agreement of each selected RNA LV -----------------

# A latent variable is one vector over every sample of all three tissues. Split
# it by tissue, key the pieces on participant and timepoint, and correlate the
# pieces over what two tissues share: three correlations per LV. A latent
# variable that is one coordinated cross-tissue programme correlates positively
# in all three pairs; one that is a tissue's own signal does not.
#
# This is also what orders the LVs in FIG4G, FIG4H, ED6A and ED6C — ascending
# blood-vs-muscle correlation.

ed6b <- function() {
  skip <- skip_if_no_fit("ED6B")
  if (!is.null(skip)) return(skip)

  panel_init("ED6B")

  selection <- rna_plier()
  plier_report_selection(selection)

  export_panel(ed6_correlation_plot(selection), "ED6B")
}

# ---- ED6C — transcripts belonging to each selected RNA LV ------------------

# A transcript belongs to an LV when its affinity (PLIER's Z) is more than
# PLIER_RNA_MARKER_SD standard deviations above the mean affinity, taken down the
# Z column. It is the same membership FIG4H's pathways and ED6D's heatmap are
# selected through.

ed6c <- function() {
  skip <- skip_if_no_fit("ED6C")
  if (!is.null(skip)) return(skip)

  panel_init("ED6C")

  selection <- rna_plier()
  plier_report_selection(selection)

  export_panel(ed6_membership_plot(selection), "ED6C")
}

# ---- ED6D — the strongest transcripts of one RNA LV ------------------------

# Which LV is PLIER_RNA_MARKER_HEATMAP_LV; how many transcripts is
# PLIER_MARKER_HEATMAP_N. Columns are every sample the fit was run on, ordered
# tissue, group, timepoint, sex and annotated with all four, so what the panel
# shows is whether the LV's transcripts move together in one tissue or in three.
#
# The legacy drew all 100 of these as PNGs to choose between and re-drew the
# chosen one as the panel. Only the chosen one is a panel here.

ed6d <- function() {
  skip <- skip_if_no_fit("ED6D")
  if (!is.null(skip)) return(skip)

  panel_init("ED6D")

  LV <- plier_env_int("PLIER_RNA_MARKER_HEATMAP_LV")
  N <- plier_env_int("PLIER_MARKER_HEATMAP_N")

  selection <- rna_plier()
  plier_check_lvs(selection$fit, LV, "rna")
  message(sprintf("        LV%d, top %d transcripts", LV, N))

  export_panel(
    ed6_marker_heatmap(selection, LV, N, labels = ed6_transcript_symbols()),
    "ED6D"
  )
}

# ---- ED6E — metabolomics latent variables against the acute response -------

# FIG4G over the metabolomics fit: same grid, same construction, same
# analysis/04_plier.R computation. Two drawing details differ and are the
# legacy's rather than a decision — the cells are outlined in grey instead of
# black, and the bubbles are drawn at 1/1.4 of the radius, which is what keeps 7
# rows of larger circles inside the same 6 mm cell. Both are arguments to
# plier_bubble_heatmap().

ed6e <- function() {
  skip <- skip_if_no_fit("ED6E")
  if (!is.null(skip)) return(skip)

  panel_init("ED6E")

  selection <- metab_plier()
  plier_report_selection(selection)

  response <- selection$response
  rows <- selection$lv_names

  bubble <- plier_bubble_heatmap(
    effect = response$effect[rows, , drop = FALSE],
    size = response$neglog10_padj[rows, , drop = FALSE],
    comparisons = response$comparisons,
    cell_border = "grey60",
    radius_divisor = 1.4
  )

  draw <- function() {
    ComplexHeatmap::draw(bubble$heatmap,
                         heatmap_legend_list = bubble$legends,
                         merge_legend = FALSE)
  }

  export_panel(draw, "ED6E")
}

# ---- ED6F — what the selected metabolomics latent variables are made of ----

# FIG4H's graphic over the metabolomics fit, with one difference in what it
# selects: the RefMet prior is small enough after PLIER's minimum-set-size filter
# that every class with any association to a selected LV fits on the page, so the
# panel keeps all of them rather than the top few of each LV.
#
# Row labels come off the RefMet class ids. The legacy hand-wrote six ("Acyl
# Carnitines", "Ceramides", "O-PC", "Saturated FA", "SM", "Unsaturated FA") and
# matched them to the rows positionally.

ed6f <- function() {
  skip <- skip_if_no_fit("ED6F")
  if (!is.null(skip)) return(skip)

  panel_init("ED6F")

  selection <- metab_plier()
  plier_report_selection(selection)

  auc <- selection$fit$Uauc
  columns <- plier_lv_index(selection$lv_names)

  auc_mat <- auc[, columns, drop = FALSE]
  auc_mat <- auc_mat[apply(auc_mat, 1, max) > 0, , drop = FALSE]
  colnames(auc_mat) <- selection$lv_names
  message(sprintf("        %d RefMet classes associated with at least one of %d latent variables",
                  nrow(auc_mat), length(columns)))

  export_panel(
    plier_auc_heatmap(auc_mat, plier_pathway_label(rownames(auc_mat))),
    "ED6F"
  )
}

# ---- ED6G — cross-tissue agreement of each selected metabolomics LV --------

# ED6B over the metabolomics fit. Same three tissue pairs, same construction; see
# ed6_correlation_plot() above.

ed6g <- function() {
  skip <- skip_if_no_fit("ED6G")
  if (!is.null(skip)) return(skip)

  panel_init("ED6G")

  selection <- metab_plier()
  plier_report_selection(selection)

  export_panel(ed6_correlation_plot(selection), "ED6G")
}

# ---- ED6H — metabolites belonging to each selected metabolomics LV ---------

# ED6C over the metabolomics fit, and the one place the two arms genuinely
# differ: the cutoff is PLIER_METAB_MARKER_SD (2, against RNA's 3) and it is
# taken ACROSS a feature's affinities rather than down an LV's, because 456
# metabolites do not spread over 100 latent variables the way 14,000 transcripts
# do. Both are the legacy's.

ed6h <- function() {
  skip <- skip_if_no_fit("ED6H")
  if (!is.null(skip)) return(skip)

  panel_init("ED6H")

  selection <- metab_plier()
  plier_report_selection(selection)

  export_panel(ed6_membership_plot(selection), "ED6H")
}

# ---- ED6I — the strongest metabolites of one metabolomics LV ---------------

# ED6D over the metabolomics fit: PLIER_METAB_MARKER_HEATMAP_LV, the same
# PLIER_MARKER_HEATMAP_N features, the same sample order and annotation. Rows are
# RefMet names, so they need no lookup.

ed6i <- function() {
  skip <- skip_if_no_fit("ED6I")
  if (!is.null(skip)) return(skip)

  panel_init("ED6I")

  LV <- plier_env_int("PLIER_METAB_MARKER_HEATMAP_LV")
  N <- plier_env_int("PLIER_MARKER_HEATMAP_N")

  selection <- metab_plier()
  plier_check_lvs(selection$fit, LV, "metab")
  message(sprintf("        LV%d, top %d metabolites", LV, N))

  export_panel(ed6_marker_heatmap(selection, LV, N), "ED6I")
}

run_panels(list(ED6A = ed6a, ED6Aii = ed6aii, ED6B = ed6b, ED6C = ed6c,
                ED6D = ed6d, ED6E = ed6e, ED6F = ed6f, ED6G = ed6g,
                ED6H = ed6h, ED6I = ed6i))
