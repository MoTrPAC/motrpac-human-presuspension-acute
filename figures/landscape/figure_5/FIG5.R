#!/usr/bin/env Rscript
# Figure 5 — Transcription-factor motif enrichment and TF exercise responses
#
# Panels:  FIG5A  nominally significant transcripts, per tissue, arm and timepoint
#          FIG5B  motifs enriched at q < 0.05 over that same grid
#          FIG5C  the five strongest motifs of every comparison, across all 22
#          FIG5D  phosphosites on transcription factors responding in muscle
#          FIG5E  endurance against resistance transcript response at 3.5/4 hr
#          FIG5F  endurance against resistance phosphosite response at 3.5/4 hr
#          FIG5G  factors by exercise response and motif enrichment, three tissues
#          FIG5H  muscle motifs enriched for DEGs whose factor is phosphorylated
#
# Every panel but FIG5A reads the HOMER knownResults of the 22 comparisons and
# tfproanno.RDS, both vendored under figure_5/sources/. Nothing is refitted here.
# ST5a is the same enrichment over the same grid and is written by tables/ST5.R,
# through the helper these panels read.
#
#   Rscript figures/landscape/FIG5.R          every panel
#   Rscript figures/landscape/FIG5.R FIG5D    one panel

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(geomtextpath)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(here, "FIG5_ED7_helpers.R"))

# ---- shared ----------------------------------------------------------------

# The HOMER parse and the motif annotation sit behind seven of the eight panels
# and are memoised in helpers/FIG5_ED7.R — tf_enrichment_cached() and
# tf_annotation_cached() — because ED7 and tables/ST5.R read the same two calls
# with the same arguments.
#
# Nothing else is shared. FIG5A and FIG5E both touch the transcript DA tables,
# but those are lazy-loaded package objects and the two panels ask different
# questions of them, so there is no loader to memoise.

# The threshold FIG5B counts motifs against and FIG5C stars its cells at.
ENRICHMENT_ALPHA <- 0.05

# Column order in FIG5C, sub-plot order in FIG5G.
DISPLAY_TISSUE_ORDER <- c("Adipose", "Blood", "Muscle")

# FIG5A and FIG5B are one bar plot each over the comparison grid, faceted by
# tissue, and share their geometry.
#
# ggbarplot rather than geom_col: it carries theme_pubr(), which is the look the
# published panels have. The dodge width equals the bar width, so the EE and RE
# bars of a timepoint touch. preserve = "single" keeps a bar the same width
# whether or not its timepoint has a second arm - blood has no during-exercise
# resistance contrast.
#
# No outline (color = NA): a stroke is centred on the bar edge, so with abutting
# bars each one's outline would be drawn over its neighbour. The legacy's
# outlines came from ggplot2's default hue scale by accident.
BAR_WIDTH <- 0.7

# ---- FIG5A — nominally significant transcripts per comparison --------------

# One bar per HOMER comparison: the number of transcripts with p_value < 0.05 in
# that comparison's exercise_with_controls contrast, counted from the transcript
# DA tables. The x axis holds all six acute timepoints in every facet, so a
# tissue that was not sampled at a timepoint reads as a gap rather than as a
# shifted axis.
#
# The HOMER output headers carry a smaller "(of N)" count than these; the
# published panel plots the DA counts.

fig5a <- function() {
  panel_init("FIG5A")

  transcript_da <- list(
    muscle  = as.data.frame(MotrpacHumanPreSuspensionAnalysis::MUSCLE_TRNSCRPT_DA),
    adipose = as.data.frame(MotrpacHumanPreSuspensionAnalysis::ADIPOSE_TRNSCRPT_DA),
    blood   = as.data.frame(MotrpacHumanPreSuspensionAnalysis::BLOOD_TRNSCRPT_DA)
  )

  counts <- tf_comparisons()
  counts$n_genes <- vapply(seq_len(nrow(counts)), function(i) {
    da <- transcript_da[[counts$tissue[i]]]
    in_contrast <- da$contrast_type == "exercise_with_controls" &
      da$randomGroupCode == counts$randomGroupCode[i] &
      da$Timepoint == counts$Timepoint[i]
    if (!any(in_contrast)) {
      stop("no exercise_with_controls contrast in the ", counts$tissue[i],
           " transcript DA for ", counts$comparison[i], call. = FALSE)
    }
    sum(in_contrast & da$p_value < 0.05, na.rm = TRUE)
  }, integer(1))

  counts$Timepoint <- factor(counts$tp_code, levels = unname(TF_TIMEPOINT_CODES))
  counts$Tissue <- factor(counts$Tissue, levels = DISPLAY_TISSUE_ORDER)
  counts$Modality <- factor(counts$Modality, levels = c("EE", "RE"))

  modality_colors <- tf_annotation_colors()$Modality

  p <- ggbarplot(
    counts,
    x = "Timepoint", y = "n_genes",
    fill = "Modality", color = NA,
    width = BAR_WIDTH,
    position = position_dodge(width = BAR_WIDTH, preserve = "single")
  ) +
    facet_wrap(~ Tissue, ncol = 1) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = modality_colors) +
    labs(x = "Time Point", y = "Sig Count") +
    theme(axis.text = element_text(angle = 45, hjust = 1))

  export_panel(p, "FIG5A")
}

# ---- FIG5B — motifs enriched per comparison --------------------------------

# The same grid as FIG5A, counting enriched motifs rather than significant
# transcripts. Read together the two say whether a comparison found few factors
# because few genes changed.
#
# The legacy assembled these counts by writing out all 22 column names by hand,
# pasting the EE and RE counts of a timepoint into one "12;7" string, splitting
# it back apart, then blanking four cells by row number to undo the blood arm the
# paste had invented. The counts are taken off the q-value matrix here, and the
# comparison grid has no cell for a blood resistance during-exercise timepoint,
# so there is nothing to blank.

fig5b <- function() {
  panel_init("FIG5B")

  enrichment <- tf_enrichment_cached()

  counts <- enrichment$comparisons
  counts$n_enriched <- colSums(enrichment$q[, counts$comparison, drop = FALSE] <
                                 ENRICHMENT_ALPHA)
  counts$Timepoint <- factor(counts$tp_code, levels = unname(TF_TIMEPOINT_CODES))
  counts$Tissue <- factor(counts$Tissue, levels = DISPLAY_TISSUE_ORDER)
  counts$Modality <- factor(counts$Modality, levels = c("EE", "RE"))

  modality_colors <- tf_annotation_colors()$Modality

  p <- ggbarplot(
    counts,
    x = "Timepoint", y = "n_enriched",
    fill = "Modality", color = NA,
    width = BAR_WIDTH,
    position = position_dodge(width = BAR_WIDTH, preserve = "single")
  ) +
    facet_wrap(~ Tissue, ncol = 1) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = modality_colors) +
    labs(x = "Time point", y = "Significantly enriched motifs") +
    theme(axis.text = element_text(angle = 45, hjust = 1))

  export_panel(p, "FIG5B")
}

# ---- FIG5C — the strongest motifs of every comparison ----------------------

# Rows are the union of the five most significantly enriched motifs of each
# comparison. That union is far shorter than 22 x 5: the same factors top many of
# the comparisons, which is the point the panel makes.
#
# Columns run adipose, blood, muscle, endurance before resistance within each,
# timepoints in order. The legacy stated that order as a vector of 22 column
# numbers into a matrix built by cbind; it is derived from the comparison grid
# here, so a comparison appearing or disappearing cannot shift a tissue boundary
# without moving its own column.

fig5c <- function() {
  panel_init("FIG5C")

  TOP_MOTIFS_PER_COMPARISON <- 5

  enrichment <- tf_enrichment_cached()
  top_motifs <- tf_top_motifs(enrichment, n = TOP_MOTIFS_PER_COMPARISON)
  message(sprintf("        %d motifs, from the top %d of each of %d comparisons",
                  length(top_motifs), TOP_MOTIFS_PER_COMPARISON,
                  nrow(enrichment$comparisons)))

  comparisons <- enrichment$comparisons
  column_order <- order(
    match(comparisons$Tissue, DISPLAY_TISSUE_ORDER),
    match(comparisons$Modality, c("EE", "RE")),
    match(comparisons$Timepoint, names(TF_TIMEPOINT_CODES))
  )
  comparisons <- comparisons[column_order, ]

  enrich_mat <- -log10(enrichment$q[top_motifs, comparisons$comparison, drop = FALSE])

  enrich_colors <- gplots::colorpanel(101, "white", "firebrick")

  # An asterisk on every cell that clears the threshold FIG5B counts against, so
  # the panel says which cells are enriched rather than leaving it to the shade.
  stars <- ifelse(enrich_mat > -log10(ENRICHMENT_ALPHA), "*", "")
  message(sprintf("        %d of %d cells enriched at q < %g",
                  sum(stars == "*"), length(stars), ENRICHMENT_ALPHA))

  draw <- function() {
    heatmap <- pheatmap::pheatmap(
      enrich_mat,
      labels_row = tf_motif_label(rownames(enrich_mat)),
      cluster_cols = FALSE,
      color = enrich_colors,
      angle_col = 315,
      annotation_col = tf_column_annotation(comparisons),
      annotation_colors = tf_annotation_colors(),
      show_colnames = FALSE,
      display_numbers = stars,
      number_color = "black",
      fontsize_number = 10,
      silent = TRUE
    )$gtable
    # The cells are -log10 of HOMER's q, which FIG5B counts against q < 0.05. The
    # bar is labelled adj_p, as FIG4G and ED6E label theirs.
    draw_heatmap_with_legend_title(heatmap, "-log10(adj_p)")
  }

  export_panel(draw, "FIG5C")
}

# ---- FIG5D — transcription-factor phosphosites in muscle -------------------

# One row per phosphosite that sits on a gene one of the 471 HOMER motifs maps
# to, and that reaches FDR < 0.05 in at least one muscle exercise-with-controls
# contrast. Circle area is -log10(adj p), capped at 5; colour is log2FC.
#
# The site set is taken through the gene symbol against HUMAN_FEATURE_TO_GENE
# rather than through MUSCLE_PROT_PH_QC$feature_metadata. Both name the same
# sites; the DA table is the thing being plotted, so it is also what the row set
# is intersected against.

fig5d <- function() {
  panel_init("FIG5D")

  SITE_ALPHA <- 0.05

  enrichment <- tf_enrichment_cached()
  annotation <- tf_annotation_cached()
  tf_symbols <- sort(unique(stats::na.omit(annotation$gene_symbol)))

  feature_to_gene <- as.data.frame(
    MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE)
  phospho <- feature_to_gene[feature_to_gene$assay == "prot-ph" &
                               as.character(feature_to_gene$gene_symbol) %in% tf_symbols, ]
  phospho$feature_id <- as.character(phospho$feature_id)
  phospho$gene_symbol <- as.character(phospho$gene_symbol)

  muscle <- enrichment$comparisons[enrichment$comparisons$Tissue == "Muscle", ]
  stats <- tf_da_stats(MotrpacHumanPreSuspensionAnalysis::MUSCLE_PROT_PH_DA,
                       phospho$feature_id, muscle)

  keep <- apply(stats$adj_p, 1, min) < SITE_ALPHA
  if (!any(keep)) {
    stop("no transcription-factor phosphosite reaches FDR < ", SITE_ALPHA,
         " in muscle", call. = FALSE)
  }
  logFC <- stats$logFC[keep, , drop = FALSE]
  adj_p <- stats$adj_p[keep, , drop = FALSE]
  message(sprintf("        %d of %d transcription-factor phosphosites are significant",
                  nrow(logFC), nrow(stats$logFC)))

  labels <- tf_phosphosite_label(
    rownames(logFC),
    phospho$gene_symbol[match(rownames(logFC), phospho$feature_id)]
  )
  if (anyDuplicated(labels) > 0) {
    dupes <- unique(labels[duplicated(labels)])
    stop(length(dupes), " phosphosite label(s) are not unique, e.g. ",
         paste(utils::head(sort(dupes), 3), collapse = ", "),
         "\n  Two features would land on one row.", call. = FALSE)
  }
  rownames(logFC) <- labels
  rownames(adj_p) <- labels

  # Endurance before resistance, timepoints ascending. The legacy stated this as
  # c(4,5,6,1,2,3) into a cbind-built matrix.
  column_order <- order(match(muscle$Modality, c("EE", "RE")),
                        match(muscle$Timepoint, names(TF_TIMEPOINT_CODES)))
  muscle <- muscle[column_order, ]
  logFC <- logFC[order(rownames(logFC)), muscle$comparison, drop = FALSE]
  adj_p <- adj_p[rownames(logFC), muscle$comparison, drop = FALSE]

  annotation_df <- droplevels(tf_column_annotation(muscle)[, c("Modality", "Timepoint")])
  annotation_colors <- tf_annotation_colors()
  top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = annotation_df,
    col = list(
      Modality = annotation_colors$Modality[levels(annotation_df$Modality)],
      Timepoint = annotation_colors$Timepoint[levels(annotation_df$Timepoint)]
    ),
    which = "column",
    border = TRUE,
    gap = grid::unit(2, "pt"),
    annotation_name_gp = grid::gpar(fontsize = 12),
    annotation_legend_param = list(
      border = TRUE,
      title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 12)
    )
  )

  ht <- tf_bubble_heatmap(logFC, adj_p, top_annotation = top_annotation,
                          outline_alpha = SITE_ALPHA)
  message(sprintf("        %d of %d bubbles outlined at adj p < %g",
                  sum(adj_p < SITE_ALPHA, na.rm = TRUE), length(adj_p), SITE_ALPHA))

  export_panel(function() {
    ComplexHeatmap::draw(ht,
                         heatmap_legend_list = tf_bubble_legends(SITE_ALPHA),
                         merge_legend = TRUE)
  }, "FIG5D")
}

# ---- the EE-against-RE scatter, FIG5E and FIG5F ----------------------------

# One point per feature, x the resistance log2FC and y the endurance log2FC of
# the same 3.5/4 hr contrast, one fitted line per tissue labelled along itself.
# Read by these two panels only, so it lives here rather than in
# helpers/FIG5_ED7.R.

# Legacy tissue colours, kept. They DISAGREE with HUMAN_TISSUE_COLORS
# (muscle #abd9e9, adipose #ffffbf, blood #d7191c), and repainting a published
# panel is not a modernization. The package adipose is a pale yellow that would
# not carry a labelled regression line on white in any case.
TF_SCATTER_TISSUE_COLORS <- c(
  "Muscle"  = "darkblue",
  "Adipose" = "goldenrod",
  "Blood"   = "firebrick"
)

TF_SCATTER_TIMEPOINT <- "post_3.5_4_hr"

#' Endurance and resistance log2FC at one timepoint, for one tissue's features.
#'
#' @param da A per-tissue DA table.
#' @param feature_ids Features to keep.
#' @param labels Row labels, named by feature id.
#' @param tissue Display name, and the colour key.
#' @param comparisons The rows of tf_comparisons() for this tissue.
tf_scatter_frame <- function(da, feature_ids, labels, tissue, comparisons) {
  at_timepoint <- comparisons[comparisons$Timepoint == TF_SCATTER_TIMEPOINT, ]
  if (nrow(at_timepoint) != 2) {
    stop(tissue, " has ", nrow(at_timepoint), " arms at ", TF_SCATTER_TIMEPOINT,
         ", expected two", call. = FALSE)
  }
  stats <- tf_da_stats(da, feature_ids, at_timepoint)

  ee <- at_timepoint$comparison[at_timepoint$Modality == "EE"]
  re <- at_timepoint$comparison[at_timepoint$Modality == "RE"]
  data.frame(
    feature_id = rownames(stats$logFC),
    label = unname(labels[rownames(stats$logFC)]),
    RE = stats$logFC[, re],
    EE = stats$logFC[, ee],
    Tissue = tissue,
    stringsAsFactors = FALSE
  )
}

#' The scatter itself.
#'
#' @param points One row per feature: RE, EE, label, Tissue.
#' @param limits Shared x and y limits, as the legacy fixed them per panel.
#' @param title Panel title.
# Where each tissue's label sits on its own fitted line: hjust along it, vjust
# across it. Blood's and adipose's lines almost coincide near the origin, so
# sliding the labels apart is not enough on its own and one is lifted off its
# line.
TF_SCATTER_LABEL_HJUST <- c("Muscle" = 0.75, "Adipose" = 0.85, "Blood" = 0.1)
TF_SCATTER_LABEL_VJUST <- c("Muscle" = 0.5, "Adipose" = -0.9, "Blood" = 1.9)

tf_scatter_plot <- function(points, limits, title) {
  tissues <- unique(points$Tissue)
  unknown <- setdiff(tissues, names(TF_SCATTER_TISSUE_COLORS))
  if (length(unknown) > 0) {
    stop("tissue with no scatter colour: ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }

  points$label_hjust <- unname(TF_SCATTER_LABEL_HJUST[points$Tissue])
  points$label_vjust <- unname(TF_SCATTER_LABEL_VJUST[points$Tissue])

  ggplot2::ggplot(points, ggplot2::aes(x = .data$RE, y = .data$EE,
                                       colour = .data$Tissue,
                                       label = .data$label)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey50",
                         linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey40", linetype = "dashed") +
    ggplot2::geom_point() +
    ggrepel::geom_text_repel(max.overlaps = 15, show.legend = FALSE) +
    # hjust slides each tissue's label along its own fitted line. The lines
    # nearly coincide near the origin, which is where the default puts every
    # label; spacing them is what keeps all of them readable.
    geomtextpath::geom_labelsmooth(
      ggplot2::aes(label = .data$Tissue, hjust = .data$label_hjust,
                   vjust = .data$label_vjust),
      fill = "white", method = "lm", formula = y ~ x,
      size = 5, linewidth = 1, boxlinewidth = 0.8, show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = TF_SCATTER_TISSUE_COLORS[tissues]) +
    ggplot2::coord_cartesian(xlim = limits, ylim = limits) +
    ggplot2::guides(colour = "none") +
    ggplot2::labs(title = title,
                  x = "log2FC resistance", y = "log2FC endurance") +
    ggplot2::theme_classic()
}

# ---- FIG5E — transcript response, endurance against resistance -------------

# One point per transcript of a gene one of the 471 HOMER motifs maps to, in each
# of the three tissues, with a fitted line per tissue. A gene on the diagonal
# responds the same way to both arms.
#
# Every feature of the tissue's DA table whose gene symbol is a motif's gene is
# plotted, significant or not: the panel is about the correlation between the two
# arms, and thresholding either axis would bend it.
#
# The legacy selected these transcripts by grep()ing the tfproanno Ensembl id
# against the DA table's versioned feature ids, which matches by substring; the
# gene symbol is used here, so ENSG00000123456 cannot also select
# ENSG000001234567.

fig5e <- function() {
  panel_init("FIG5E")

  AXIS_LIMITS <- c(-2, 5.75)

  enrichment <- tf_enrichment_cached()
  annotation <- tf_annotation_cached()
  tf_symbols <- sort(unique(stats::na.omit(annotation$gene_symbol)))

  feature_to_gene <- as.data.frame(
    MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE)
  transcripts <- feature_to_gene[
    feature_to_gene$assay == "transcript-rna-seq" &
      as.character(feature_to_gene$gene_symbol) %in% tf_symbols, ]
  labels <- stats::setNames(as.character(transcripts$gene_symbol),
                            as.character(transcripts$feature_id))

  da_tables <- list(
    Muscle  = MotrpacHumanPreSuspensionAnalysis::MUSCLE_TRNSCRPT_DA,
    Adipose = MotrpacHumanPreSuspensionAnalysis::ADIPOSE_TRNSCRPT_DA,
    Blood   = MotrpacHumanPreSuspensionAnalysis::BLOOD_TRNSCRPT_DA
  )

  points <- do.call(rbind, lapply(names(da_tables), function(tissue) {
    tf_scatter_frame(
      da_tables[[tissue]], names(labels), labels, tissue,
      enrichment$comparisons[enrichment$comparisons$Tissue == tissue, ]
    )
  }))
  message(sprintf("        %d transcripts across %d tissues",
                  nrow(points), length(da_tables)))

  off_axis <- sum(points$RE < AXIS_LIMITS[1] | points$RE > AXIS_LIMITS[2] |
                    points$EE < AXIS_LIMITS[1] | points$EE > AXIS_LIMITS[2])
  if (off_axis > 0) {
    message(sprintf("        %d point(s) fall outside the fixed axes and are clipped",
                    off_axis))
  }

  p <- tf_scatter_plot(points, AXIS_LIMITS,
                       "Comparing EE vs RE in Transcriptomics P3.5/4H")

  export_panel(p, "FIG5E")
}

# ---- FIG5F — phosphosite response, endurance against resistance ------------

# FIG5E's question asked of the phosphoproteome. Blood is absent rather than
# omitted: there is no blood phosphoproteomics, so muscle and adipose are the
# whole comparison.

fig5f <- function() {
  panel_init("FIG5F")

  AXIS_LIMITS <- c(-0.8, 1.2)

  enrichment <- tf_enrichment_cached()
  annotation <- tf_annotation_cached()
  tf_symbols <- sort(unique(stats::na.omit(annotation$gene_symbol)))

  feature_to_gene <- as.data.frame(
    MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE)
  phospho <- feature_to_gene[
    feature_to_gene$assay == "prot-ph" &
      as.character(feature_to_gene$gene_symbol) %in% tf_symbols, ]
  labels <- stats::setNames(
    tf_phosphosite_label(as.character(phospho$feature_id),
                         as.character(phospho$gene_symbol)),
    as.character(phospho$feature_id))

  da_tables <- list(
    Muscle  = MotrpacHumanPreSuspensionAnalysis::MUSCLE_PROT_PH_DA,
    Adipose = MotrpacHumanPreSuspensionAnalysis::ADIPOSE_PROT_PH_DA
  )

  points <- do.call(rbind, lapply(names(da_tables), function(tissue) {
    tf_scatter_frame(
      da_tables[[tissue]], names(labels), labels, tissue,
      enrichment$comparisons[enrichment$comparisons$Tissue == tissue, ]
    )
  }))
  message(sprintf("        %d phosphosites across %d tissues",
                  nrow(points), length(da_tables)))

  off_axis <- sum(points$RE < AXIS_LIMITS[1] | points$RE > AXIS_LIMITS[2] |
                    points$EE < AXIS_LIMITS[1] | points$EE > AXIS_LIMITS[2])
  if (off_axis > 0) {
    message(sprintf("        %d point(s) fall outside the fixed axes and are clipped",
                    off_axis))
  }

  p <- tf_scatter_plot(points, AXIS_LIMITS,
                       "Comparing EE vs RE in Phosphoproteomics P3.5/4H")

  export_panel(p, "FIG5F")
}

# ---- the four-way classification FIG5G stacks ------------------------------

# Every motif that maps to a gene this study measured is placed in one of four
# boxes per tissue and ome: whether that gene responds to exercise in that ome,
# crossed with whether its motif is enriched for the tissue's differentially
# expressed genes. Read by FIG5G only, so it lives here.
#
# The legacy wrote each of the 32 counts as its own line of nested
# intersect(rownames(tfproanno[...]), rownames(tfproanno[...])) calls, and the
# three tissues did not agree on what they were counting: "measured" for
# phosphoproteomics was every prot-ph feature in HUMAN_FEATURE_TO_GENE rather
# than the ones in that tissue's DA table, so muscle and adipose shared a
# denominator. Every membership question here is asked of the tissue's own DA
# table.

TF_ACTIVITY_CLASSES <- c("Not Sig, Not Enriched", "Not Sig, Enriched",
                         "Sig, Not Enriched", "Sig, Enriched")

# No package equivalent; the legacy hexes, kept.
TF_ACTIVITY_COLORS <- stats::setNames(
  c("#97c976", "#dba24c", "#fd717a", "#ce69cd"),
  TF_ACTIVITY_CLASSES
)

#' The omes each tissue contributes, as display name -> (assay, DA object name).
TF_ACTIVITY_OMES <- list(
  Muscle = list(
    "Transcriptomics"    = c(assay = "transcript-rna-seq", da = "MUSCLE_TRNSCRPT_DA"),
    "Proteomics (MS)"    = c(assay = "prot-pr",            da = "MUSCLE_PROT_PR_DA"),
    "Phosphoproteomics"  = c(assay = "prot-ph",            da = "MUSCLE_PROT_PH_DA")
  ),
  Adipose = list(
    "Transcriptomics"    = c(assay = "transcript-rna-seq", da = "ADIPOSE_TRNSCRPT_DA"),
    "Proteomics (MS)"    = c(assay = "prot-pr",            da = "ADIPOSE_PROT_PR_DA"),
    "Phosphoproteomics"  = c(assay = "prot-ph",            da = "ADIPOSE_PROT_PH_DA")
  ),
  Blood = list(
    "Transcriptomics"     = c(assay = "transcript-rna-seq", da = "BLOOD_TRNSCRPT_DA"),
    "Proteomics (Olink)"  = c(assay = "prot-ol",            da = "BLOOD_PROT_OL_DA")
  )
)

#' Motifs whose gene has any feature in one tissue x ome, and which of those
#' have a significant one.
#'
#' @returns A list of two character vectors of motif names, `measured` and `sig`.
tf_ome_membership <- function(annotation, assay, da) {
  feature_to_gene <- as.data.frame(
    MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE)
  rows <- feature_to_gene[feature_to_gene$assay == assay, ]
  symbol_of <- stats::setNames(as.character(rows$gene_symbol),
                               as.character(rows$feature_id))

  measured_symbols <- unique(symbol_of[tf_measured_features(da)])
  sig_symbols <- unique(symbol_of[tf_significant_features(da)])

  motif_symbol <- annotation$gene_symbol
  list(
    measured = rownames(annotation)[!is.na(motif_symbol) &
                                      motif_symbol %in% measured_symbols],
    sig = rownames(annotation)[!is.na(motif_symbol) &
                                 motif_symbol %in% sig_symbols]
  )
}

#' The four counts, for every ome of one tissue.
#'
#' @param tissue "Muscle", "Adipose" or "Blood".
#' @param annotation tf_annotation() over the motif list.
#' @param enriched Motifs enriched anywhere in this tissue.
#' @returns A data.frame: Ome, Classification, Count.
tf_activity_counts <- function(tissue, annotation, enriched) {
  omes <- TF_ACTIVITY_OMES[[tissue]]
  if (is.null(omes)) {
    stop("no ome list for tissue '", tissue, "'", call. = FALSE)
  }

  rows <- lapply(names(omes), function(ome) {
    spec <- omes[[ome]]
    da <- get(spec[["da"]], envir = asNamespace("MotrpacHumanPreSuspensionAnalysis"))
    membership <- tf_ome_membership(annotation, spec[["assay"]], da)

    measured <- membership$measured
    sig <- intersect(membership$sig, measured)
    not_sig <- setdiff(measured, sig)

    data.frame(
      Ome = ome,
      Classification = TF_ACTIVITY_CLASSES,
      Count = c(
        length(setdiff(not_sig, enriched)),
        length(intersect(not_sig, enriched)),
        length(setdiff(sig, enriched)),
        length(intersect(sig, enriched))
      ),
      stringsAsFactors = FALSE
    )
  })

  counts <- do.call(rbind, rows)
  counts$Ome <- factor(counts$Ome, levels = names(omes))
  counts$Classification <- factor(counts$Classification,
                                  levels = TF_ACTIVITY_CLASSES)
  counts
}

# Axis labels, broken where the one-line names overlap.
TF_ACTIVITY_AXIS_LABELS <- c(
  "Transcriptomics"    = "Transcriptomics",
  "Proteomics (MS)"    = "Proteomics\n(MS)",
  "Proteomics (Olink)" = "Proteomics\n(Olink)",
  "Phosphoproteomics"  = "Phospho-\nproteomics"
)

#' The stacked bar plot itself.
tf_activity_plot <- function(counts, title) {
  ggplot2::ggplot(counts, ggplot2::aes(x = .data$Ome, y = .data$Count,
                                       fill = .data$Classification)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(values = TF_ACTIVITY_COLORS) +
    # Two lines rather than one: "Proteomics (MS)" and "Phosphoproteomics"
    # collide head-on across a 6 inch axis, and rotating them runs the longest
    # off the page.
    ggplot2::scale_x_discrete(labels = TF_ACTIVITY_AXIS_LABELS) +
    ggplot2::labs(title = title, x = NULL, y = "Transcription factors") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 12))
}

# ---- FIG5G — factors by exercise response and motif enrichment -------------

# Three stacked sub-plots, adipose over blood over muscle, sharing one legend.
# One bar per ome the tissue contributes. A transcription factor counts once per
# ome in which its gene is measured in that tissue; it is significant if any of
# that gene's features reaches FDR < 0.05 in any exercise-with-controls
# contrast, and enriched if its motif reaches q < 0.05 in any of that tissue's
# HOMER runs.
#
# The legacy drew adipose, blood and muscle as three PNGs (6G.I, 6G.II, 6G.III)
# to be stacked in Illustrator, each carrying its own copy of the legend. They
# are composited here, so what the manuscript places is what this panel wrote.
#
# The y scale is free: blood contributes two omes and about 280 factors, muscle
# three and about 430, and a shared scale would flatten adipose to nothing.

fig5g <- function() {
  panel_init("FIG5G")

  enrichment <- tf_enrichment_cached()
  annotation <- tf_annotation_cached()

  # Top to bottom in DISPLAY_TISSUE_ORDER, the order the legacy's three files
  # were numbered in.
  sub_plots <- lapply(DISPLAY_TISSUE_ORDER, function(tissue) {
    enriched <- tf_enriched_motifs(enrichment, tissue)
    counts <- tf_activity_counts(tissue, annotation, enriched)
    message(sprintf("        %-7s %d motifs enriched, %d factors over %d ome(s)",
                    tissue, length(enriched), sum(counts$Count),
                    nlevels(counts$Ome)))
    tf_activity_plot(counts, paste(tissue, "transcription factor activity"))
  })

  p <- patchwork::wrap_plots(sub_plots, ncol = 1) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "right")

  export_panel(p, "FIG5G")
}

# ---- FIG5H — muscle motifs enriched for DEGs of phosphorylated factors -----

# Twelve columns per row: the six muscle phosphorylation contrasts, then the six
# motif-enrichment ones. A row is a motif that is enriched somewhere (q < 0.05)
# AND carries a phosphosite that responds somewhere (FDR < 0.05), so every row
# is a factor for which both lines of evidence exist.
#
# The cell carries -log10 of the adjusted p value, as FIG5C does. The legacy used
# the natural log (threshold 2.995732, which is -log(0.05)) under a colour bar
# titled -log10. The scale runs 0 to 5, which on log10 is HOMER's q floor of
# 1e-5.
#
# Cells for which there is no measurement stay NA and are painted in pheatmap's
# na_col grey. The legacy carried them as -5, which falls outside the 0-5 breaks
# and lands on the same grey by a different route - and so did every value ABOVE
# 5, which is where the strongest enrichments are. Those are clipped to 5 here
# rather than left to render as missing data.

fig5h <- function() {
  panel_init("FIG5H")

  ALPHA <- 0.05
  SCALE_MAX <- 5
  TISSUE <- "Muscle"

  enrichment <- tf_enrichment_cached()
  annotation <- tf_annotation_cached()
  muscle <- enrichment$comparisons[enrichment$comparisons$Tissue == TISSUE, ]

  # Endurance before resistance, timepoints ascending: the order both blocks are
  # read in, so a column of one block sits above the matching column of the other.
  muscle <- muscle[order(match(muscle$Modality, c("EE", "RE")),
                         match(muscle$Timepoint, names(TF_TIMEPOINT_CODES))), ]

  motifs <- enrichment$motifs
  block <- function(values, prefix) {
    colnames(values) <- paste(prefix, muscle$Modality, muscle$tp_code)
    values
  }

  enrich_p <- block(enrichment$q[motifs, muscle$comparison, drop = FALSE], "Enrich")

  # The per-motif minimum is tf_min_adj_p_by_motif(), shared with ST5a: the table
  # and this panel report the same number for the same factor because they call
  # the same function, not because two copies of it agree.

  # The legacy also filled a six-column "TF gene" block, from the factor's own
  # transcript, and then drew columns 13:18 and 1:6 — phosphorylation and
  # enrichment. Nothing plotted the transcript block; it is not computed here.
  phos_p <- block(
    tf_min_adj_p_by_motif(annotation, "prot-ph",
                          MotrpacHumanPreSuspensionAnalysis::MUSCLE_PROT_PH_DA,
                          muscle),
    "TF phos")

  summary_mat <- cbind(enrich_p, phos_p)
  summary_log <- -log10(summary_mat)

  enrich_cols <- colnames(enrich_p)
  phos_cols <- colnames(phos_p)
  threshold <- -log10(ALPHA)

  has_max <- function(m) {
    out <- suppressWarnings(apply(m, 1, max, na.rm = TRUE))
    out[!is.finite(out)] <- -Inf
    out
  }
  keep <- has_max(summary_log[, enrich_cols, drop = FALSE]) > threshold &
    has_max(summary_log[, phos_cols, drop = FALSE]) >= threshold
  if (!any(keep)) {
    stop("no muscle motif is both enriched and phosphorylated at alpha ", ALPHA,
         call. = FALSE)
  }

  plot_mat <- summary_log[keep, c(phos_cols, enrich_cols), drop = FALSE]

  # Rows are clustered the way the published panel clustered them: on the
  # unclipped -log10 values, with unmeasured cells at -5. Clipping first changes the
  # distances and so the row order; the clip below is for the colour scale only.
  cluster_mat <- plot_mat
  cluster_mat[is.na(cluster_mat)] <- -5
  if (any(!is.finite(cluster_mat))) {
    stop("FIG5H: non-finite -log value in the clustering matrix", call. = FALSE)
  }
  row_tree <- stats::hclust(stats::dist(cluster_mat), method = "complete")

  clipped <- sum(plot_mat > SCALE_MAX, na.rm = TRUE)
  plot_mat[!is.na(plot_mat) & plot_mat > SCALE_MAX] <- SCALE_MAX
  message(sprintf(
    "        %d of %d motifs are both enriched and phosphorylated; %d cell(s) clipped at %d, %d not measured",
    nrow(plot_mat), nrow(summary_log), clipped, SCALE_MAX, sum(is.na(plot_mat))))

  # 101 colours need 102 break points. The legacy passed 101 of each, which leaves
  # the last colour unreachable; the ramp is otherwise the one it used.
  enrich_colors <- gplots::colorpanel(101, "white", "firebrick")

  draw <- function() {
    heatmap <- pheatmap::pheatmap(
      plot_mat,
      breaks = seq(0, SCALE_MAX, length.out = 102),
      color = enrich_colors,
      angle_col = 315,
      labels_row = toupper(tf_motif_label(rownames(plot_mat))),
      cluster_rows = row_tree,
      cluster_cols = FALSE,
      silent = TRUE
    )$gtable
    draw_heatmap_with_legend_title(heatmap, "-log10(adj_p)")
  }

  export_panel(draw, "FIG5H")
}

run_panels(list(
  FIG5A = fig5a,
  FIG5B = fig5b,
  FIG5C = fig5c,
  FIG5D = fig5d,
  FIG5E = fig5e,
  FIG5F = fig5f,
  FIG5G = fig5g,
  FIG5H = fig5h
))
