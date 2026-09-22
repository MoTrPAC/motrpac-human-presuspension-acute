#!/usr/bin/env Rscript
# Extended Data 9 — Secretome supplement
#
# Panels:  ED9A  exerkine candidates shared across tissues, by exercise arm
#          ED9B  exerkine candidates not shown in FIG7A
#          ED9C  FLNA phosphosite S1533, muscle and adipose
#          ED9D  CCN1 transcript in adipose against muscle, per timepoint
#
# The COMPARTMENTS extracellular scores ED9A and ED9B read are vendored under
# figure_7/sources/; see docs/external_dependencies.md.
#
# Needs consortium data access. ED9D reads the adipose and muscle transcript
# qc matrices through load_qc().
#
#   Rscript figures/landscape/ED9.R          every panel
#   Rscript figures/landscape/ED9.R ED9B     one panel

suppressPackageStartupMessages({
  library(circlize)
  library(ComplexHeatmap)
  library(ComplexUpset)
  library(ggplot2)
  library(patchwork)
  # load_differential_analysis(), which plot_single_feature() calls, resolves its
  # lazy-loaded DA objects by name through the search path, so the Analysis
  # package has to be attached and not only namespace-qualified.
  library(MotrpacHumanPreSuspensionAnalysis)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(here, "..", "lib", "highlights.R"))
source(file.path(here, "..", "lib", "single_feature_helpers.R"))
source(file.path(here, "..", "single_feature_plots.R"))
source(file.path(here, "FIG7_ED9_helpers.R"))

# ---- shared ----------------------------------------------------------------

# There is no shared loader here. The one frame ED9A and ED9B both build — the
# exerkine candidates — is memoised in helpers/FIG7_ED9.R instead, because
# FIG7A reads it with the same arguments. ED9D's load_qc() is this figure's
# alone and is read once.

# ---- ED9A — exerkine candidates shared across tissues ----------------------

# The FIG7A candidate frame, one UpSet per arm (EE only, RE only, both), each
# gene counted once in the intersection of the tissues it is differentially
# abundant in. Bars are filled by the ome that carries the gene: within a
# tissue the highest-priority ome (protein, then transcript, then phosphosite);
# across tissues "Mixed" when they disagree.

OME_LABEL <- c("prot-pr" = "Protein", "transcript-rna-seq" = "Transcript",
               "prot-ph" = "Phospho")
OME_PRIORITY <- c(Protein = 1, Transcript = 2, Phospho = 3)
OME_FILL <- c(
  Protein    = unname(MotrpacHumanPreSuspensionAnalysis::HUMAN_OME_COLORS["prot-pr"]),
  Transcript = unname(MotrpacHumanPreSuspensionAnalysis::HUMAN_OME_COLORS["transcript-rna-seq"]),
  Phospho    = unname(MotrpacHumanPreSuspensionAnalysis::HUMAN_OME_COLORS["prot-ph"]),
  Mixed      = "grey50"
)
TISSUE_SETS <- c("Muscle", "Adipose", "Blood")
ARMS <- c("EE", "RE", "Both")

ed9a <- function() {
  panel_init("ED9A")

  candidates <- exerkine_candidate_frame()

  # One row per gene x tissue, carrying that tissue's best ome.
  per_tissue <- unique(data.frame(
    gene_symbol = candidates$gene_symbol,
    Tissue = unname(TISSUE_DISPLAY[candidates$tissue]),
    Omics = unname(OME_LABEL[candidates$assay]),
    Exercise2 = candidates$Exercise2,
    stringsAsFactors = FALSE
  ))
  per_tissue <- per_tissue[order(per_tissue$gene_symbol, per_tissue$Tissue,
                                 OME_PRIORITY[per_tissue$Omics]), ]
  per_tissue <- per_tissue[!duplicated(per_tissue[, c("gene_symbol", "Tissue")]), ]

  # One row per gene: presence in each tissue, one ome for the fill.
  genes <- unique(per_tissue$gene_symbol)
  presence <- data.frame(gene_symbol = genes, stringsAsFactors = FALSE)
  for (tissue in TISSUE_SETS) {
    presence[[tissue]] <- genes %in% per_tissue$gene_symbol[per_tissue$Tissue == tissue]
  }
  ome_by_gene <- tapply(per_tissue$Omics, per_tissue$gene_symbol, function(x) {
    if (length(unique(x)) == 1) unique(x) else "Mixed"
  })
  presence$Omics <- factor(unname(ome_by_gene[genes]), levels = names(OME_FILL))
  presence$Exercise2 <- factor(
    unname(tapply(per_tissue$Exercise2, per_tissue$gene_symbol, unique)[genes]),
    levels = ARMS
  )
  stopifnot(nrow(presence) == length(genes), !anyNA(presence$Omics))

  # ComplexUpset drops a set with no members and applies the stripes by
  # position, so each arm passes only the tissues it has, with their own colours.
  upset_plots <- lapply(ARMS, function(arm) {
    arm_presence <- presence[presence$Exercise2 == arm, c(TISSUE_SETS, "Omics")]
    sets <- TISSUE_SETS[colSums(arm_presence[, TISSUE_SETS, drop = FALSE]) > 0]
    ComplexUpset::upset(
      arm_presence,
      sets,
      sort_sets = FALSE,
      min_size = 1,
      base_annotations = list(
        "Intersection size" = intersection_size(
          counts = TRUE,
          mapping = aes(fill = Omics)
        ) +
          scale_fill_manual(values = OME_FILL,
                            guide = guide_legend(title = "Omics Layer"))
      ),
      set_sizes = upset_set_size(),
      width_ratio = 0.2,
      stripes = unname(TISSUE_COLORS[sets]),
      name = "Tissue Overlap"
    ) + ggtitle(paste("Exercise:", arm))
  })

  export_panel(function() print(wrap_plots(upset_plots, ncol = 1)), "ED9A")
}

# ---- ED9B — exerkine candidates not shown in FIG7A -------------------------

# The same frame as FIG7A, every gene that panel leaves out: an extracellular
# score under 4, or a plasma protein that falls after exercise. Cells carry
# both plasma directions, so all four direction pairs occur.

ed9b <- function() {
  panel_init("ED9B")

  candidates <- exerkine_candidate_frame()

  in_main <- candidates$score >= COMPARTMENTS_MAIN_SCORE & candidates$blood_direction == "Up"
  main_genes <- unique(candidates$gene_symbol[in_main])
  rest <- candidates[!candidates$gene_symbol %in% main_genes, , drop = FALSE]

  m <- candidate_heatmap_matrix(rest)

  heatmap <- candidate_heatmap(
    m, CANDIDATE_CELL_COLORS,
    candidate_row_annotation(rest, rownames(m), score_range = c(0, 5))
  )

  export_panel(heatmap, "ED9B")
}

# ---- ED9C — FLNA phosphosite S1533 -----------------------------------------

# One phosphosite on filamin A, in the two tissues that quantify it: muscle
# under the isoform-2 accession, adipose under the canonical one. Exercise-group
# means with 95% confidence intervals across the acute timepoints, under one
# collected legend. A point is filled black where that timepoint's adjusted p
# value is below 0.05.
#
# The plot is one entry of the single-feature catalog in single_feature_plots.R.

ed9c <- function() {
  panel_init("ED9C")
  export_panel(single_feature_plot("ED9C"), "ED9C")
}

# ---- ED9D — CCN1 transcript in adipose against muscle ----------------------

# Every participant with both biopsies at a timepoint, coloured by exercise
# group, with a least-squares line per group and its Pearson r and p in the
# facet. Pre-exercise pools all three groups.

ed9d <- function() {
  panel_init("ED9D")

  points <- ccn1_cross_tissue(MotrpacHumanPreSuspensionData::load_qc())
  stats <- ccn1_cross_tissue_stats(points)

  timepoints <- names(SECRETOME_TIMEPOINT_LABELS)[names(SECRETOME_TIMEPOINT_LABELS) %in% points$Timepoint]
  points$Timepoint <- factor(points$Timepoint, levels = timepoints,
                             labels = SECRETOME_TIMEPOINT_LABELS[timepoints])
  points$group <- factor(points$group, levels = c("CON", "EE", "RE"))

  captions <- aggregate(label ~ Timepoint, data = stats,
                        FUN = function(x) paste(x, collapse = "\n"))
  captions$Timepoint <- factor(captions$Timepoint, levels = timepoints,
                               labels = SECRETOME_TIMEPOINT_LABELS[timepoints])

  p <- ggplot(points, aes(x = Adipose, y = Muscle, color = group)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
    geom_text(data = captions, aes(label = label), color = "black",
              x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2, size = 3.2,
              inherit.aes = FALSE) +
    facet_wrap(~ Timepoint, ncol = 2, scales = "free") +
    scale_color_manual(name = "Group", values = EXERCISE_COLORS[c("CON", "EE", "RE")]) +
    labs(x = "Adipose CCN1", y = "Muscle CCN1") +
    theme_minimal(base_size = 15)

  export_panel(p, "ED9D")
}

run_panels(list(
  ED9A = ed9a,
  ED9B = ed9b,
  ED9C = ed9c,
  ED9D = ed9d
))
