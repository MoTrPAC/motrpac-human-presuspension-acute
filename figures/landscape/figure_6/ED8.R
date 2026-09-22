#!/usr/bin/env Rscript
# Extended Data 8 — TFEB pathway regulation in response to EE and RE
#
# Panels:  ED8A  muscle transcript response of the ChIP-supported SC-ION
#                targets of TFEB
#          ED8B  TFEB transcript across the three tissues
#          ED8C  PP2A subunit transcripts quantified in both blood and muscle
#          ED8D  HSP90AA1 transcript across the three tissues
#          ED8E  ChIP-supported TFEB targets in GO-BP Regulation of Autophagy,
#                up at 3.5-4 h
#          ED8F  leucine and Leu-Ile in blood and muscle
#
# Figure 6's supplement: the TFEB ChIP evidence ED8A and ED8E select on is the
# same computation FIG6D writes ST6c from, and helpers/FIG6_ED8.R carries it.
# The published ChIP sources and the SC-ION target list are vendored under
# figure_6/sources/; nothing here refits anything.
#
#   Rscript figures/landscape/ED8.R          every panel
#   Rscript figures/landscape/ED8.R ED8C     one panel

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(ggplot2)
  library(patchwork)
  library(MotrpacHumanPreSuspensionAnalysis)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(here, "..", "lib", "highlights.R"))
source(file.path(here, "..", "lib", "single_feature_helpers.R"))
source(file.path(here, "..", "single_feature_plots.R"))
source(file.path(here, "FIG6_ED8_helpers.R"))

# ---- shared ----------------------------------------------------------------

# There is no shared loader here. The two loads more than one panel makes —
# the muscle transcript DA (ED8A, ED8E) and the published TFEB ChIP sources
# (ED8A, ED8E) — are memoised in helpers/FIG6_ED8.R instead, because FIG6 reads
# the same two with the same arguments. ED8C reads the blood and muscle
# transcript DA with different arguments again, so it shares nothing with them.

ALPHA <- 0.05
TIMEPOINT <- "post_3.5_4_hr"

# ---- ED8A — ChIP-supported SC-ION targets of TFEB --------------------------

# The SC-ION-predicted TFEB targets that at least one published ChIP source
# calls bound, over the six exercise-with-controls muscle transcript contrasts.
# Too many rows to label; the row order is the clustering's.
#
# The ChIP evidence is the same function FIG6D writes ST6c from, so the rows
# here are ST6c's SCION.target = Y and any = Y rows.

ed8a <- function() {
  panel_init("ED8A")

  targets <- scion_tfeb_targets()
  evidence <- tfeb_chip_evidence(targets)
  bound <- evidence$gene_symbol[evidence$any == "Y"]
  message(sprintf("        %d of %d SC-ION targets have ChIP evidence", length(bound), length(targets)))

  da <- muscle_transcript_exercise_da()
  da <- da[!is.na(da$gene_symbol) & da$gene_symbol %in% bound, ]
  da <- one_transcript_per_symbol(da, "ED8A")
  absent <- setdiff(bound, unique(da$gene_symbol))
  if (length(absent) > 0) {
    message(sprintf("        %d bound target(s) with no muscle transcript DA, not drawn: %s",
                    length(absent), paste(absent, collapse = ", ")))
  }

  matrices <- exercise_contrast_matrices(da, key = "gene_symbol")

  heatmap <- logfc_star_heatmap(matrices, breaks = c(-1, 0, 1.5), show_row_names = FALSE)

  export_panel(draw_heatmap(heatmap), "ED8A")
}

# ---- ED8B — TFEB transcript trajectories -----------------------------------

# Exercise-group means with 95% confidence intervals over the acute timepoints,
# a black point where the timepoint clears FDR 0.05.
#
# One entry of the single-feature catalog in single_feature_plots.R.

ed8b <- function() {
  panel_init("ED8B")
  export_panel(single_feature_plot("ED8B"), "ED8B")
}

# ---- ED8C — PP2A subunit transcripts ---------------------------------------

# The PP2A subunit genes in config/highlights.json over the exercise-with-
# controls transcript contrasts of blood (six EE, four RE) and muscle (three
# and three), side by side under a tissue annotation. A subunit quantified in
# one tissue only is not drawn, as in the legacy.

ed8c <- function() {
  panel_init("ED8C")

  genes <- highlight_gene_list("ED8C")

  tissue_da <- function(da, tissue) {
    da <- da[da$contrast_type == "exercise_with_controls" &
               da$contrast_category %in% names(EXERCISE_GROUP_OF), ]
    symbols <- HUMAN_FEATURE_TO_GENE[HUMAN_FEATURE_TO_GENE$assay == "transcript-rna-seq",
                                     c("feature_id", "gene_symbol")]
    symbols <- symbols[!duplicated(symbols$feature_id), ]
    da <- merge(as.data.frame(da), symbols, by = "feature_id")
    da <- da[!is.na(da$gene_symbol) & da$gene_symbol %in% genes &
               !grepl("PAR_Y", da$feature_id, fixed = TRUE), ]
    exercise_contrast_matrices(da, key = "gene_symbol", tissue = tissue)
  }
  blood <- tissue_da(BLOOD_TRNSCRPT_DA, "blood")
  muscle <- tissue_da(MUSCLE_TRNSCRPT_DA, "muscle")

  shared <- intersect(rownames(blood$fc), rownames(muscle$fc))
  message(sprintf("        %d of %d PP2A subunits quantified in both tissues (blood %d, muscle %d)",
                  length(shared), length(genes), nrow(blood$fc), nrow(muscle$fc)))
  if (length(shared) == 0) {
    stop("no PP2A subunit is quantified in both blood and muscle", call. = FALSE)
  }

  matrices <- list(
    fc = cbind(blood$fc[shared, , drop = FALSE], muscle$fc[shared, , drop = FALSE]),
    adj_p = cbind(blood$adj_p[shared, , drop = FALSE], muscle$adj_p[shared, , drop = FALSE]),
    annotation = rbind(blood$annotation, muscle$annotation)
  )

  heatmap <- logfc_star_heatmap(matrices, breaks = c(-1, 0, 1))

  export_panel(draw_heatmap(heatmap), "ED8C")
}

# ---- ED8D — HSP90AA1 transcript trajectories -------------------------------

# Exercise-group means with 95% confidence intervals over the acute timepoints,
# a black point where the timepoint clears FDR 0.05.

ed8d <- function() {
  panel_init("ED8D")
  export_panel(single_feature_plot("ED8D"), "ED8D")
}

# ---- ED8E — ChIP-supported TFEB targets in Regulation of Autophagy ---------

# The members of the set named in config/highlights.json among FIG6D's two
# "up" inputs: muscle transcripts with published TFEB ChIP evidence that are
# up (adj p < 0.05, logFC > 0) at 3.5-4 h after EE or RE. Drawn over all six
# exercise-with-controls contrasts.
#
# The PR read this gene list from an .rds it did not ship
# (tfeb_ora_leading_edge.rds); the definition above reproduces its 50 rows
# exactly.

ed8e <- function() {
  panel_init("ED8E")

  pathway <- highlight_pathways("ED8E")
  set_name <- pathway$set[1]
  members <- NULL
  for (collection in MOLECULAR_SIGNATURES) {
    if (set_name %in% names(collection)) members <- collection[[set_name]]
  }
  if (is.null(members)) {
    stop("MOLECULAR_SIGNATURES carries no set named ", set_name, call. = FALSE)
  }

  transcripts <- muscle_transcript_symbols()
  evidence <- tfeb_chip_evidence(transcripts$gene_symbol)
  bound_ids <- transcripts$feature_id[evidence$any == "Y"]

  da <- muscle_transcript_exercise_da()
  up <- da[da$Timepoint == TIMEPOINT & !is.na(da$adj_p_value) & da$adj_p_value < ALPHA &
             da$logFC > 0 & da$feature_id %in% bound_ids, ]
  genes <- sort(intersect(unique(up$gene_symbol), members))
  message(sprintf("        %d ChIP-supported targets up at 3.5-4 h are in %s (%d members)",
                  length(genes), set_name, length(members)))
  if (length(genes) == 0) {
    stop("no ChIP-supported TFEB target up at 3.5-4 h is in ", set_name, call. = FALSE)
  }

  da <- da[!is.na(da$gene_symbol) & da$gene_symbol %in% genes, ]
  da <- one_transcript_per_symbol(da, "ED8E")
  matrices <- exercise_contrast_matrices(da, key = "gene_symbol")

  heatmap <- logfc_star_heatmap(matrices, breaks = c(-1.5, 0, 1.5), row_fontsize = 7)

  export_panel(draw_heatmap(heatmap), "ED8E")
}

# ---- ED8F — leucine and Leu-Ile trajectories -------------------------------

# Exercise-group means with 95% confidence intervals over the acute timepoints,
# a black point where the timepoint clears FDR 0.05.

ed8f <- function() {
  panel_init("ED8F")
  export_panel(single_feature_plot("ED8F"), "ED8F")
}

run_panels(list(
  ED8A = ed8a,
  ED8B = ed8b,
  ED8C = ed8c,
  ED8D = ed8d,
  ED8E = ed8e,
  ED8F = ed8f
))
