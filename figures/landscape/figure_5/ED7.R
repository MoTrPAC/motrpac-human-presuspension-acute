#!/usr/bin/env Rscript
# Extended Data 7 — transcription-factor supplement
#
# Panels:  ED7A  transcription-factor transcript responses, three tissues, first half
#          ED7B  the second half of the same heatmap
#
# One row per transcription-factor gene, 22 columns: every tissue x arm x
# timepoint the HOMER analysis covers. A gene not measured in a tissue carries
# log2FC 0 and adj p 1 there, which draws as no bubble - the same convention the
# legacy used. The row set does not fit on one page, so it is drawn as two panels
# rather than one.
#
# Built on helpers/FIG5_ED7.R rather than on a copy of it: ED7 is Figure 5's
# supplement and reads the same motif list, the same annotation and the same
# comparison grid, and two copies of those is how they would stop being the same.
#
#   Rscript figures/landscape/ED7.R          both panels
#   Rscript figures/landscape/ED7.R ED7B     one panel

suppressPackageStartupMessages({
  library(ComplexHeatmap)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(here, "FIG5_ED7_helpers.R"))

# ---- shared ----------------------------------------------------------------

ED7_ALPHA <- 1e-05
# Which bubbles are outlined. Separate from ED7_ALPHA, which decides which genes
# get a ROW at all: a gene earns its row by reaching 1e-5 somewhere, and each of
# its cells is marked on the ordinary 0.05, as FIG5D marks its own.
ED7_OUTLINE_ALPHA <- 0.05
ED7_DISPLAY_TISSUE_ORDER <- c("Adipose", "Blood", "Muscle")

#' Every tissue's transcription-factor transcript response, in one matrix pair.
#'
#' @param annotation tf_annotation() over the motif list.
#' @param comparisons tf_comparisons().
#' @returns A list of two gene x comparison matrices, `logFC` and `adj_p`.
tf_combined_transcript_stats <- function(annotation, comparisons) {
  tf_symbols <- sort(unique(stats::na.omit(annotation$gene_symbol)))

  feature_to_gene <- as.data.frame(
    MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE)
  transcripts <- feature_to_gene[
    feature_to_gene$assay == "transcript-rna-seq" &
      as.character(feature_to_gene$gene_symbol) %in% tf_symbols, ]
  symbol_of <- stats::setNames(as.character(transcripts$gene_symbol),
                               as.character(transcripts$feature_id))

  da_tables <- list(
    Muscle  = MotrpacHumanPreSuspensionAnalysis::MUSCLE_TRNSCRPT_DA,
    Adipose = MotrpacHumanPreSuspensionAnalysis::ADIPOSE_TRNSCRPT_DA,
    Blood   = MotrpacHumanPreSuspensionAnalysis::BLOOD_TRNSCRPT_DA
  )

  per_tissue <- lapply(names(da_tables), function(tissue) {
    stats <- tf_da_stats(da_tables[[tissue]], names(symbol_of),
                         comparisons[comparisons$Tissue == tissue, ])
    rownames(stats$logFC) <- symbol_of[rownames(stats$logFC)]
    rownames(stats$adj_p) <- symbol_of[rownames(stats$adj_p)]
    stats
  })
  names(per_tissue) <- names(da_tables)

  genes <- sort(Reduce(union, lapply(per_tissue, function(s) rownames(s$logFC))))

  combine <- function(field, absent) {
    out <- matrix(absent, nrow = length(genes), ncol = nrow(comparisons),
                  dimnames = list(genes, comparisons$comparison))
    for (stats in per_tissue) {
      m <- stats[[field]]
      out[rownames(m), colnames(m)] <- m
    }
    out
  }

  list(logFC = combine("logFC", 0), adj_p = combine("adj_p", 1))
}

#' The display column order: adipose, blood, muscle; endurance before resistance.
tf_combined_column_order <- function(comparisons) {
  comparisons[order(
    match(comparisons$Tissue, ED7_DISPLAY_TISSUE_ORDER),
    match(comparisons$Modality, c("EE", "RE")),
    match(comparisons$Timepoint, names(TF_TIMEPOINT_CODES))
  ), ]
}

# The matrix pair both halves are cut from, and the column order they are cut
# in. One pass, lazy and memoised: ED7A and ED7B ask for exactly the same thing,
# and a build of the whole figure would otherwise read every transcript DA table
# twice. The HOMER parse and the motif annotation behind it are memoised in
# helpers/FIG5_ED7.R, where FIG5 and tables/ST5.R read them too.
ed7_combined_stats <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    comparisons <- tf_combined_column_order(tf_enrichment_cached()$comparisons)
    stats <- tf_combined_transcript_stats(tf_annotation_cached(), comparisons)

    keep <- apply(stats$adj_p, 1, min) < ED7_ALPHA
    if (!any(keep)) {
      stop("no transcription-factor transcript reaches FDR < ", ED7_ALPHA,
           " in any tissue", call. = FALSE)
    }
    message(sprintf("        %d of %d transcription-factor genes reach FDR < %g",
                    sum(keep), length(keep), ED7_ALPHA))

    cache <<- list(
      comparisons = comparisons,
      logFC = stats$logFC[keep, comparisons$comparison, drop = FALSE],
      adj_p = stats$adj_p[keep, comparisons$comparison, drop = FALSE]
    )
    cache
  }
})

#' One half of the row set, as a two-matrix list.
#'
#' The heatmap is too tall for one page, so it is drawn as two. The legacy split
#' at `nrow / 2` and indexed with `1:half` and `(half + 1):n`; on an odd row count
#' those truncate to the same integer and the last row is dropped from both
#' halves. Split here so that the two halves are always the whole.
#'
#' @param half 1 or 2.
tf_combined_half <- function(stats, half) {
  stopifnot(half %in% c(1, 2))
  n <- nrow(stats$logFC)
  cut <- ceiling(n / 2)
  rows <- if (half == 1) seq_len(cut) else seq(cut + 1, n)
  list(logFC = stats$logFC[rows, , drop = FALSE],
       adj_p = stats$adj_p[rows, , drop = FALSE])
}

#' Tissue, modality and timepoint annotation for the combined columns.
tf_combined_annotation <- function(comparisons) {
  annotation_df <- droplevels(tf_column_annotation(comparisons)[
    , c("Tissue", "Modality", "Timepoint")])
  colors <- tf_annotation_colors()
  ComplexHeatmap::HeatmapAnnotation(
    df = annotation_df,
    col = list(
      Tissue = colors$Tissue[levels(annotation_df$Tissue)],
      Modality = colors$Modality[levels(annotation_df$Modality)],
      Timepoint = colors$Timepoint[levels(annotation_df$Timepoint)]
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
}

#' Everything ED7A and ED7B do apart from which half they draw.
tf_combined_panel <- function(panel, half) {
  stats <- ed7_combined_stats()

  this_half <- tf_combined_half(stats, half)
  message(sprintf("        half %d of 2: rows %d..%d of %d",
                  half,
                  match(rownames(this_half$logFC)[1], rownames(stats$logFC)),
                  match(rownames(this_half$logFC)[nrow(this_half$logFC)],
                        rownames(stats$logFC)),
                  nrow(stats$logFC)))

  ht <- tf_bubble_heatmap(this_half$logFC, this_half$adj_p,
                          top_annotation = tf_combined_annotation(stats$comparisons),
                          outline_alpha = ED7_OUTLINE_ALPHA)
  message(sprintf("        %d of %d bubbles outlined at adj p < %g",
                  sum(this_half$adj_p < ED7_OUTLINE_ALPHA, na.rm = TRUE),
                  length(this_half$adj_p), ED7_OUTLINE_ALPHA))

  export_panel(function() {
    ComplexHeatmap::draw(ht,
                         heatmap_legend_list = tf_bubble_legends(ED7_OUTLINE_ALPHA),
                         merge_legend = TRUE)
  }, panel)
}

# ---- ED7A — transcript responses, first half -------------------------------

# Every gene one of the 471 HOMER motifs maps to that reaches FDR < 1e-5 in at
# least one tissue, arm and timepoint. Circle area is -log10(adj p), capped at 5;
# colour is log2FC. Columns run adipose, blood, muscle, endurance before
# resistance within each. The split is on the alphabetical row order and is
# stated once, in tf_combined_half().

ed7a <- function() {
  panel_init("ED7A")
  tf_combined_panel("ED7A", half = 1)
}

# ---- ED7B — transcript responses, second half ------------------------------

# ED7A's second half: the same matrix, the same columns, the rows below the cut.

ed7b <- function() {
  panel_init("ED7B")
  tf_combined_panel("ED7B", half = 2)
}

run_panels(list(
  ED7A = ed7a,
  ED7B = ed7b
))
