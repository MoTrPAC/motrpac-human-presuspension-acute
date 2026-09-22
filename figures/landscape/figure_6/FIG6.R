#!/usr/bin/env Rscript
# Figure 6 — Multi-omic regulatory network analysis (SC-ION)
#
# Panels:  FIG6C  TFEB phosphosite log fold change in muscle, EE and RE
#          FIG6D  ORA of ChIP-bound TFEB targets up or down at 3.5-4 h, EE and
#                 RE                          (also writes ST6c through ST6g)
#          FIG6E  muscle transcript response of the shared SC-ION targets of
#                 TFEB and RXRA
#          FIG6F  lactic acid across the three tissues
# Tables:  ST6c   TFEB ChIP targets            (written by FIG6D)
#          ST6d   EE 3.5-4 h up                (written by FIG6D)
#          ST6e   EE 3.5-4 h down              (written by FIG6D)
#          ST6f   RE 3.5-4 h up                (written by FIG6D)
#          ST6g   RE 3.5-4 h down              (written by FIG6D)
#
# Not built here, see FIG6's not_built in config/panel_map.json:
#   FIG6A  SC-ION workflow schematic                  Illustrator drawing
#   FIG6B  merged EE+RE network and the TFEB          Cytoscape rendering; the
#          first-neighbour subnetwork                 edges and node scores it
#                                                     draws are ST6a and ST6b
#   FIG6G  TFEB signalling module schematic           BioRender illustration
#
# The network inference is not reproduced here. The merged network's edge and
# node tables are vendored under figure_6/sources/, FIG6E reads the edge table
# as given, and tables/ST6.R re-emits both; see docs/external_dependencies.md.
#
#   Rscript figures/landscape/FIG6.R          every panel
#   Rscript figures/landscape/FIG6.R FIG6D    one panel

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(ggplot2)
  library(patchwork)
  library(MotrpacHumanPreSuspensionAnalysis)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(here, "..", "lib", "table_export.R"))
source(file.path(here, "..", "lib", "highlights.R"))
source(file.path(here, "..", "single_feature_plots.R"))
source(file.path(here, "FIG6_ED8_helpers.R"))

# ---- shared ----------------------------------------------------------------

# There is no shared loader here. The two loads more than one panel makes —
# the muscle transcript DA (FIG6D, FIG6E) and the published TFEB ChIP sources
# (FIG6D) — are memoised in helpers/FIG6_ED8.R instead, because ED8 reads the
# same two with the same arguments. FIG6E's scion_edges() is this figure's
# alone and is read once.

# ---- FIG6C — TFEB phosphosite log fold change in muscle --------------------

# Every quantified TFEB phosphosite in the muscle phosphoproteome over the six
# exercise-with-controls contrasts, a star where adj p < 0.05. Rows are named
# TFEB.<site>; the colour ramp runs from the smallest logFC drawn through white
# at 0 to red at 1, as the legacy's did.

fig6c <- function() {
  panel_init("FIG6C")

  GENE <- "TFEB"

  da <- MUSCLE_PROT_PH_DA[
    MUSCLE_PROT_PH_DA$contrast_type == "exercise_with_controls" &
      MUSCLE_PROT_PH_DA$contrast_category %in% names(EXERCISE_GROUP_OF), ]
  symbols <- HUMAN_FEATURE_TO_GENE[HUMAN_FEATURE_TO_GENE$assay == "prot-ph",
                                   c("feature_id", "gene_symbol")]
  symbols <- symbols[!duplicated(symbols$feature_id), ]
  da <- merge(as.data.frame(da), symbols, by = "feature_id")
  da <- da[!is.na(da$gene_symbol) & da$gene_symbol == GENE, ]
  if (nrow(da) == 0) {
    stop("no ", GENE, " phosphosite in the muscle phosphoproteome DA", call. = FALSE)
  }
  # P19484_S122s -> TFEB.S122s
  da$site <- paste(GENE, sub("^[^_]+_", "", da$feature_id), sep = ".")

  matrices <- exercise_contrast_matrices(da, key = "site")
  message(sprintf("        %d %s phosphosites over %d contrasts",
                  nrow(matrices$fc), GENE, ncol(matrices$fc)))

  heatmap <- logfc_star_heatmap(matrices, breaks = c(min(matrices$fc, na.rm = TRUE), 0, 1))

  export_panel(draw_heatmap(heatmap), "FIG6C")
}

# ---- FIG6D — ORA of ChIP-bound TFEB targets at 3.5-4 h ---------------------

# Four inputs: the muscle transcripts with published TFEB ChIP evidence that
# are DA (adj p < 0.05) at 3.5-4 h after EE or RE, split by the sign of the
# fold change. Each is tested against every muscle transcript with a gene
# symbol over the full MOLECULAR_SIGNATURES collection through run_ORA(). The
# panel draws the curated selection in config/highlights.json: a dot per
# set x input sized by -log10(adj p), capped at 5, on a grey cell where the
# adjusted p clears 0.05.
#
# This panel also writes ST6c (the ChIP evidence over every muscle transcript,
# with the SC-ION prediction flag) and ST6d-ST6g (the four unfiltered ORA
# results).

ALPHA <- 0.05
TIMEPOINT <- "post_3.5_4_hr"
DOT_CAP <- 5

fig6d <- function() {
  panel_init("FIG6D")

  # ---- ST6c: ChIP evidence over every muscle transcript ----

  transcripts <- muscle_transcript_symbols()
  evidence <- tfeb_chip_evidence(transcripts$gene_symbol)
  scion_targets <- scion_tfeb_targets()
  chip <- data.frame(
    feature_id = transcripts$feature_id,
    gene_symbol = transcripts$gene_symbol,
    SCION.target = ifelse(transcripts$gene_symbol %in% scion_targets, "Y", "N"),
    evidence[, setdiff(names(evidence), "gene_symbol")],
    stringsAsFactors = FALSE
  )
  write_st6c(chip)

  # ---- the four inputs ----

  da <- muscle_transcript_exercise_da()
  background <- unique(as.character(da$gene_symbol[!is.na(da$gene_symbol)]))
  at_timepoint <- da[da$Timepoint == TIMEPOINT & !is.na(da$adj_p_value) &
                       da$adj_p_value < ALPHA &
                       da$feature_id %in% chip$feature_id[chip$any == "Y"], ]
  input_of <- function(category, direction) {
    rows <- at_timepoint[at_timepoint$contrast_category == category &
                           sign(at_timepoint$logFC) == direction, ]
    unique(as.character(rows$gene_symbol[!is.na(rows$gene_symbol)]))
  }
  inputs <- list(
    ST6d = list(name = "EE_3.5_4_hr_up", genes = input_of("EE-CON", 1)),
    ST6e = list(name = "EE_3.5_4_hr_down", genes = input_of("EE-CON", -1)),
    ST6f = list(name = "RE_3.5_4_hr_up", genes = input_of("RE-CON", 1)),
    ST6g = list(name = "RE_3.5_4_hr_down", genes = input_of("RE-CON", -1))
  )
  message(sprintf("        inputs: %s; background %d",
                  paste(vapply(inputs, function(x) sprintf("%s %d", x$name, length(x$genes)),
                               character(1)), collapse = ", "),
                  length(background)))

  # ---- the ORA, and ST6d-ST6g ----

  ora <- lapply(inputs, function(input) {
    as.data.frame(run_ORA(input$genes, background))
  })
  write_st6d(ora$ST6d)
  write_st6e(ora$ST6e)
  write_st6f(ora$ST6f)
  write_st6g(ora$ST6g)
  # The drawn columns are named for the input, not for the sub-table it was
  # written as.
  names(ora) <- vapply(inputs, function(x) x$name, character(1))

  # ---- the curated selection ----

  curated <- highlight_pathways("FIG6D")$set
  missing_sets <- setdiff(curated, unique(unlist(lapply(ora, function(x) x$set))))
  if (length(missing_sets) > 0) {
    message(sprintf("        FIG6D: %d curated set(s) in no ORA result and not drawn: %s",
                    length(missing_sets), paste(missing_sets, collapse = ", ")))
  }
  drawn <- setdiff(curated, missing_sets)
  if (length(drawn) == 0) {
    stop("none of the curated sets is in this run's ORA", call. = FALSE)
  }

  neg_log_p <- sapply(ora, function(result) {
    -log10(result$adj_p_value[match(drawn, result$set)])
  })
  rownames(neg_log_p) <- drawn
  message(sprintf("        drawing %d of %d curated sets over %d inputs",
                  nrow(neg_log_p), length(curated), ncol(neg_log_p)))

  # ---- plot ----

  cell <- grid::unit(5, "mm")
  scale_max <- max(neg_log_p, na.rm = TRUE)
  col_fun <- circlize::colorRamp2(breaks = c(0, scale_max), colors = c("white", "#57347d"))

  heatmap <- ComplexHeatmap::Heatmap(
    neg_log_p,
    col = col_fun,
    rect_gp = grid::gpar(type = "none"),
    heatmap_legend_param = list(title = "-log10(adjusted p-value)",
                                legend_direction = "vertical",
                                legend_width = grid::unit(50, "mm"),
                                at = c(0, round(scale_max, 2))),
    width = ncol(neg_log_p) * cell, height = nrow(neg_log_p) * cell,
    cluster_columns = FALSE, cluster_rows = TRUE,
    show_column_names = TRUE, column_names_side = "top", show_row_names = TRUE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      value <- neg_log_p[i, j]
      significant <- !is.na(value) && value > -log10(ALPHA)
      grid::grid.rect(x = x, y = y, width, height,
                      gp = grid::gpar(col = "#555555",
                                      fill = if (significant) "gray" else "white"))
      if (!is.na(value)) {
        grid::grid.circle(x = x, y = y, r = min(value, DOT_CAP) * 0.02,
                          gp = grid::gpar(fill = col_fun(value), col = "black"))
      }
    },
    column_title_rot = 90, row_title_rot = 0,
    column_names_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 8)
  )

  significance_legend <- ComplexHeatmap::Legend(
    labels = c("< 0.05", ">= 0.05"), title = "BH Adjusted p-value", type = "grid",
    legend_gp = grid::gpar(fill = c("gray", "white")), border = "black")

  export_panel(draw_heatmap(heatmap, legend_side = "left",
                            annotation_legend_list = list(significance_legend)),
               "FIG6D")
}

# ---- ST6c — TFEB ChIP targets ----------------------------------------------

# Every muscle transcript with a gene symbol, flagged Y/N for the SC-ION
# prediction and for each of the four ChIP sources. It is FIG6D's own input
# rather than a table of its own, which is why it is written here.
write_st6c <- function(chip) {
  export_table(chip[, table_spec("ST6c")$columns], "ST6c")
}

# ---- ST6d-ST6g — the four unfiltered ORA results ---------------------------

# One function per id so each is named where it is written, as ST6c is. The
# results are FIG6D's four inputs before the curated selection cuts them to the
# rows the panel draws.
write_st6d <- function(ora) {
  export_table(ora[, table_spec("ST6d")$columns], "ST6d")
}

write_st6e <- function(ora) {
  export_table(ora[, table_spec("ST6e")$columns], "ST6e")
}

write_st6f <- function(ora) {
  export_table(ora[, table_spec("ST6f")$columns], "ST6f")
}

write_st6g <- function(ora) {
  export_table(ora[, table_spec("ST6g")$columns], "ST6g")
}

# ---- FIG6E — shared SC-ION targets of TFEB and RXRA ------------------------

# The targets both TFEB and RXRA point at in the merged EE+RE network, over the
# six exercise-with-controls muscle transcript contrasts. Targets are the
# regulator's OUTGOING edges: a first-neighbour definition would also bring in
# TFEB's own upstream regulators (NFIC regulates TFEB and is not a target).
#
# The legacy's 3-way version with CHCHD3 is not reproduced: under the v2.0
# network CHCHD3 has no direct edge from TFEB.
#
# No guard on the edge table. It is vendored under figure_6/sources/ rather
# than produced by a fit, so panel_source() resolves it from the repo and
# FIG6_SCION_EDGES_CSV only points the panel at a newer copy.

fig6e <- function() {
  panel_init("FIG6E")

  REGULATORS <- c("TFEB", "RXRA")

  edges <- scion_edges()
  targets <- lapply(REGULATORS, scion_targets_of, edges = edges)
  names(targets) <- REGULATORS
  shared <- sort(Reduce(intersect, targets))
  message(sprintf("        %s; shared %d: %s",
                  paste(sprintf("%s targets %d", REGULATORS, lengths(targets)), collapse = ", "),
                  length(shared), paste(shared, collapse = ", ")))
  if (length(shared) == 0) {
    stop("TFEB and RXRA share no target in the edge table", call. = FALSE)
  }

  # The vendored target list and the edge table describe the same network.
  recorded <- scion_tfeb_targets()
  if (!setequal(recorded, targets$TFEB)) {
    message(sprintf(
      "        FIG6E: the edge table's TFEB targets (%d) differ from TFEB_targets_v2.xlsx (%d): %d only in the table, %d only in the file",
      length(targets$TFEB), length(recorded),
      length(setdiff(targets$TFEB, recorded)), length(setdiff(recorded, targets$TFEB))))
  }

  da <- muscle_transcript_exercise_da()
  da <- da[!is.na(da$gene_symbol) & da$gene_symbol %in% shared, ]
  da <- one_transcript_per_symbol(da, "FIG6E")
  absent <- setdiff(shared, unique(da$gene_symbol))
  if (length(absent) > 0) {
    message(sprintf("        %d shared target(s) with no muscle transcript DA, not drawn: %s",
                    length(absent), paste(absent, collapse = ", ")))
  }

  matrices <- exercise_contrast_matrices(da, key = "gene_symbol")

  heatmap <- logfc_star_heatmap(matrices, breaks = c(-1, 0, 1.5))

  export_panel(draw_heatmap(heatmap), "FIG6E")
}

# ---- FIG6F — lactic acid trajectories --------------------------------------

# The metabolite over the acute timepoints in adipose, blood and muscle:
# exercise-group means with 95% confidence intervals, a black point where the
# timepoint clears FDR 0.05.
#
# One entry of the single-feature catalog in single_feature_plots.R.

fig6f <- function() {
  panel_init("FIG6F")
  export_panel(single_feature_plot("FIG6F"), "FIG6F")
}

run_panels(list(
  FIG6C = fig6c,
  FIG6D = fig6d,
  FIG6E = fig6e,
  FIG6F = fig6f
))
