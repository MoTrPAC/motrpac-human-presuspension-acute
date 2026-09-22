# plier_helpers.R — everything the cross-tissue PLIER figures share.
#
# Two figures read one pair of fits: FIG4G and FIG4H are the RNA arm, ED6A-ED6I
# are both arms. analysis/04_plier.R produces the fits and the per-LV response
# statistics; no figure script refits anything.
#
# This file sits in lib/ rather than with a figure because both figures read it
# AND the script that writes their input does — the same rule that puts
# da_overlap_helpers.R here.
#
# Sourced after panel_export.R, which defines landscape_root().
#
# What lives here:
#   plier_lvs(), plier_env_*()      the hand-set numbers out of figure_4/plier_lvs.env
#   plier_fit(), plier_input()      the fitted objects, read by name
#   plier_response()                the LV x comparison response grid
#   plier_markers()                 an LV's marker features, by affinity cutoff
#   plier_tissue_correlation()      an LV's cross-tissue agreement
#   plier_lv_display_order()        the order every panel draws its LVs in
#   plier_annotation_colors()       the package palettes under the display names
#   plier_bubble_heatmap()          the ComplexHeatmap FIG4G and ED6E both draw
#   plier_lv_response_stats()       analysis/04_plier.R's computation, here so the
#                                   grid the fit writes and the grid the panels
#                                   read are one definition
#
# Display names. The legacy renamed tissue, timepoint and exercise group to short
# figure labels by ASSIGNING OVER the package palettes positionally
# (`names(cols$Tissue) <- c("Blood","Muscle","Adipose")`). That is correct only
# while the package keeps that order. The maps below are keyed by the package's
# own names, so a reordered constant renames nothing and a renamed one fails.

suppressPackageStartupMessages({
  library(MotrpacHumanPreSuspensionAnalysis)
})

# ---- the hand-set numbers ---------------------------------------------------

# plier.env sits beside this file and is read on load, so a panel gets the
# latent-variable selections by sourcing this helper and nothing else. It is a
# shell file because it is a list of settings and reads like one, and because
# analysis/04_plier.R sources it directly; parsed here rather than run through a
# shell so an R session needs no subprocess.
#
# Values already exported in the environment WIN, which is what makes a one-off
# override possible:
#   PLIER_METAB_LVS=4,17 Rscript figures/landscape/ED6.R ED6E
# It is also how the fit script, which has sourced the file itself, is left
# alone.
.plier_load_env <- function() {
  path <- file.path(landscape_root(), "figure_4", "plier_lvs.env")
  if (!file.exists(path)) {
    stop("plier.env is missing from lib/", call. = FALSE)
  }
  for (line in readLines(path, warn = FALSE)) {
    line <- trimws(line)
    if (!startsWith(line, "export PLIER_")) next
    body <- sub("^export ", "", line)
    name <- sub("=.*$", "", body)
    if (nzchar(Sys.getenv(name, unset = ""))) next          # already set: leave it
    value <- sub("^[^=]*=", "", body)
    value <- sub('^"', "", sub('"$', "", value))
    # "${PLIER_X:-default}" -> default
    value <- sub("^\\$\\{[A-Za-z0-9_]+:-", "", value)
    value <- sub("\\}$", "", value)
    do.call(Sys.setenv, stats::setNames(list(value), name))
  }
  invisible(NULL)
}

.plier_load_env()

# ---- naming -----------------------------------------------------------------

PLIER_TISSUE_LABELS <- c(adipose = "Adipose", blood = "Blood", muscle = "Muscle")

PLIER_QC_TISSUE_LABELS <- c(
  "Human Adipose Powder" = "Adipose",
  "PaxGene RNA"          = "Blood",
  "Muscle"               = "Muscle"
)

PLIER_TIMEPOINT_LABELS <- c(
  pre_exercise      = "Pre",
  during_20_min     = "D20M",
  during_40_min     = "D40M",
  post_10_min       = "P10M",
  post_15_30_45_min = "P15-45M",
  post_3.5_4_hr     = "P3.5/4H",
  post_24_hr        = "P24H"
)

PLIER_MODALITY_LABELS <- c(ADUResist = "RE", ADUEndur = "EE", ADUControl = "CON")

PLIER_TISSUE_ORDER    <- c("Adipose", "Blood", "Muscle")
PLIER_MODALITY_ORDER  <- c("EE", "RE")
PLIER_TIMEPOINT_ORDER <- unname(PLIER_TIMEPOINT_LABELS)

#' Recode a vector of package-level values to the figures' display labels.
#'
#' Unmapped values are an error rather than an NA: an unrecognised timepoint in
#' the metadata means the study grew a timepoint, and a panel that silently drew
#' it as NA would be wrong in a way no check catches.
plier_relabel <- function(x, map, what) {
  x <- as.character(x)
  unknown <- setdiff(unique(x), names(map))
  if (length(unknown) > 0) {
    stop(what, " value with no display label: ", paste(unknown, collapse = ", "),
         "\n  Add it to PLIER_", toupper(what), "_LABELS in lib/plier_helpers.R.",
         call. = FALSE)
  }
  unname(map[x])
}

# ---- figure_4/plier_lvs.env ----------------------------------------------------------

# `file` names the env file that declares the variable, for the error only. The
# PLIER knobs are split in two: the hand-set latent-variable selections the
# panels draw are in figure_4/plier_lvs.env, and the fit's own parameters are in
# figure_4/plier.env, which analysis/04_plier.R reads. Sending a reader to the
# wrong one is worse than not naming a file at all.
plier_env <- function(name, file = "figure_4/plier_lvs.env") {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    stop(name, " is unset.\n  It is declared in ", file, ".", call. = FALSE)
  }
  value
}

plier_env_int <- function(name, file = "figure_4/plier_lvs.env") {
  as.integer(plier_env(name, file))
}
plier_env_num <- function(name, file = "figure_4/plier_lvs.env") {
  as.numeric(plier_env(name, file))
}

plier_env_ints <- function(name) {
  parts <- trimws(strsplit(plier_env(name), ",", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  out <- suppressWarnings(as.integer(parts))
  if (anyNA(out)) {
    stop(name, " is not a comma-separated list of integers: ", plier_env(name),
         call. = FALSE)
  }
  out
}

#' The latent variables an arm's panels are built on, in fit order.
#'
#' @param arm "rna" or "metab"
plier_lvs <- function(arm) {
  plier_env_ints(sprintf("PLIER_%s_LVS", toupper(plier_check_arm(arm))))
}

plier_marker_sd <- function(arm) {
  plier_env_num(sprintf("PLIER_%s_MARKER_SD", toupper(plier_check_arm(arm))))
}

plier_check_arm <- function(arm) {
  arm <- tolower(arm)
  if (!arm %in% c("rna", "metab")) {
    stop("unknown PLIER arm: ", arm, " (expected \"rna\" or \"metab\")", call. = FALSE)
  }
  arm
}

# ---- the fitted objects -----------------------------------------------------

#' Where analysis/04_plier.R writes its fits, and where the panels read them.
#'
#' PLIER_DIR overrides the default, which is under outputs/ so the fit and the
#' panels agree without anything being set.
plier_dir <- function() {
  dir <- Sys.getenv("PLIER_DIR", unset = "")
  if (!nzchar(dir)) {
    dir <- file.path(landscape_root(), "outputs", "fits", "plier")
  }
  dir
}

plier_path <- function(arm, what) {
  file.path(plier_dir(), sprintf("crosstissue_%s_%s.rds", plier_check_arm(arm), what))
}

plier_read <- function(arm, what) {
  path <- plier_path(arm, what)
  if (!file.exists(path)) {
    stop("the PLIER fit has not produced ", basename(path), ".\n  Expected it at ",
         path, "\n  Run Rscript figures/landscape/analysis/04_plier.R (hours for",
         " the RNA arm on a cold start).", call. = FALSE)
  }
  readRDS(path)
}

#' The PLIER result itself: $Z affinities, $B latent-variable values, $U and
#' $Uauc prior associations.
plier_fit <- function(arm) plier_read(arm, "plier")

#' What the fit was run on: $matrix (features x samples, z-scored per tissue) and
#' $metadata (samples x annotation, in matrix column order, display labels).
plier_input <- function(arm) plier_read(arm, "input")

#' The per-LV response grid: $effect (LV x comparison median difference from
#' control), $neglog10_padj (same shape), $comparisons (one row per column).
plier_response <- function(arm) plier_read(arm, "response")

# ---- selections off a fit ---------------------------------------------------

#' Name a fit's latent variables LV1..LVk.
plier_lv_names <- function(k) paste0("LV", seq_len(k))

#' Assert the hand-set LV numbers name something in this fit.
plier_check_lvs <- function(fit, lvs, arm) {
  k <- ncol(fit$Z)
  bad <- lvs[lvs < 1 | lvs > k]
  if (length(bad) > 0) {
    stop("figure_4/plier_lvs.env asks for LV", paste(bad, collapse = ", LV"),
         " but the ", arm, " fit has ", k, " latent variables.",
         "\n  A refit renumbers latent variables; re-check PLIER_",
         toupper(arm), "_LVS against it.", call. = FALSE)
  }
  invisible(lvs)
}

#' The features whose affinity to an LV is more than `sd_cutoff` standard
#' deviations above the mean, taken in affinity order.
#'
#' The legacy scaled the two arms' Z matrices differently — by column for RNA,
#' by row for metabolomics — and both are kept: which margin the cutoff is taken
#' over decides how many markers each LV gets, and those counts are what ED6C and
#' ED6H plot.
plier_markers <- function(fit, lvs, sd_cutoff, margin = c("column", "row")) {
  margin <- match.arg(margin)
  scored <- if (margin == "column") scale(fit$Z) else t(scale(t(fit$Z)))
  counts <- colSums(scored > sd_cutoff)
  out <- lapply(lvs, function(i) {
    n <- counts[i]
    if (is.na(n) || n < 1) {
      stop("LV", i, " has no feature above ", sd_cutoff,
           " SD of affinity — nothing for the panels to draw.", call. = FALSE)
    }
    rownames(fit$Z)[order(-fit$Z[, i])][seq_len(n)]
  })
  names(out) <- paste0("LV", lvs)
  out
}

#' The top-n features of an LV by affinity, whatever their score.
plier_top_features <- function(fit, lv, n) {
  rownames(fit$Z)[order(-fit$Z[, lv])][seq_len(n)]
}

# ---- cross-tissue agreement -------------------------------------------------

#' Correlate each LV's sample values between each pair of tissues.
#'
#' A latent variable is one vector over all samples of all three tissues; this
#' splits it by tissue, keys the pieces on participant and timepoint, and
#' correlates the pieces over the participant-timepoints two tissues share. It
#' is what ED6B and ED6G plot and what orders the LVs in every other panel.
plier_tissue_correlation <- function(B, metadata) {
  key <- paste(metadata$pid, metadata$Timepoint, sep = "_")
  by_tissue <- lapply(PLIER_TISSUE_ORDER, function(tissue) {
    keep <- metadata$Tissue == tissue
    m <- B[, keep, drop = FALSE]
    colnames(m) <- key[keep]
    m
  })
  names(by_tissue) <- PLIER_TISSUE_ORDER

  pairs <- list(
    Adipose_v_Blood  = c("Adipose", "Blood"),
    Adipose_v_Muscle = c("Adipose", "Muscle"),
    Blood_v_Muscle   = c("Blood", "Muscle")
  )
  out <- data.frame(row.names = rownames(B))
  for (nm in names(pairs)) {
    a <- by_tissue[[pairs[[nm]][1]]]
    b <- by_tissue[[pairs[[nm]][2]]]
    shared <- intersect(colnames(a), colnames(b))
    out[[nm]] <- vapply(seq_len(nrow(B)), function(i) {
      stats::cor(a[i, shared], b[i, shared])
    }, numeric(1))
  }
  out$LV <- rownames(B)
  out
}

#' The order every LV panel draws its latent variables in: ascending
#' blood-vs-muscle correlation, over ALL latent variables, then subset to the
#' selected ones.
#'
#' Taken over all of them rather than over the selection because that is what the
#' legacy did — the factor levels come from the full correlation frame — and
#' because it keeps the order stable when one LV is added to or dropped from the
#' selection.
plier_lv_display_order <- function(correlation, lv_names) {
  ordered <- correlation$LV[order(correlation$Blood_v_Muscle)]
  ordered[ordered %in% lv_names]
}

# ---- the response grid ------------------------------------------------------

#' The tissue x arm x timepoint cells a response is reported for.
#'
#' Derived from the metadata rather than stated as 22 column names: blood is the
#' only tissue with during-exercise samples and only its endurance arm has them,
#' which is a fact about the study the legacy restated as a hand-written vector
#' of 22 strings in four places.
plier_comparison_grid <- function(metadata) {
  cells <- unique(metadata[
    metadata$Timepoint != "Pre" & metadata$Modality %in% PLIER_MODALITY_ORDER,
    c("Tissue", "Modality", "Timepoint")
  ])
  cells <- cells[order(
    match(cells$Tissue, PLIER_TISSUE_ORDER),
    match(cells$Modality, PLIER_MODALITY_ORDER),
    match(cells$Timepoint, PLIER_TIMEPOINT_ORDER)
  ), ]
  cells$comparison <- paste(cells$Tissue, cells$Modality, cells$Timepoint, sep = "_")
  rownames(cells) <- cells$comparison
  cells
}

#' Participants with a pre-exercise sample in a tissue.
#'
#' Every response below is a within-participant difference from that
#' participant's own baseline, so a participant with no baseline in a tissue
#' contributes nothing there.
plier_baselined_subjects <- function(metadata, tissue) {
  in_tissue <- metadata[metadata$Tissue == tissue, ]
  pre <- table(in_tissue$pid[in_tissue$Timepoint == "Pre"])
  names(pre)[pre == 1]
}

#' The whole LV x comparison response grid, for every latent variable in a fit.
#'
#' For each LV, each tissue: subtract each participant's own pre-exercise value,
#' then compare each exercise arm to the controls at each timepoint. Reported as
#' the median difference between the arm and the controls (the colour of a
#' bubble) and as -log10 of the adjusted p (its size).
#'
#' Two departures from the legacy, both of them fixes §14.3 describes:
#'
#'   - It read `compare_means()`'s output BY ROW NUMBER — `adiposetestout[3, ]`
#'     is asserted to be the 3.5/4 hr endurance test — in 88 hand-written lines
#'     per matrix, four matrices. A tissue that gained or lost a timepoint would
#'     silently relabel every row after it. The rows are joined on their own
#'     Timepoint and group2 here.
#'   - The 22 columns were four hand-written vectors of names, and the metadata
#'     frame describing them a fifth. One grid, from plier_comparison_grid().
plier_lv_response_stats <- function(fit, metadata, progress_every = 10) {
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("ggpubr is needed for the response tests", call. = FALSE)
  }
  B <- fit$B
  rownames(B) <- plier_lv_names(nrow(B))
  stopifnot(identical(colnames(B), rownames(metadata)))

  comparisons <- plier_comparison_grid(metadata)
  tissues <- unique(comparisons$Tissue)

  # Per tissue: the rows to use, and each row's baseline row, resolved once
  # rather than per latent variable.
  frames <- lapply(tissues, function(tissue) {
    keep <- metadata$Tissue == tissue &
      metadata$pid %in% plier_baselined_subjects(metadata, tissue)
    df <- metadata[keep, c("pid", "Tissue", "Modality", "Timepoint"), drop = FALSE]
    df$sample <- rownames(metadata)[keep]
    baseline_row <- match(
      paste(df$pid, "Pre"),
      paste(df$pid, df$Timepoint)
    )
    list(df = df, baseline = baseline_row)
  })
  names(frames) <- tissues

  effect <- matrix(NA_real_, nrow = nrow(B), ncol = nrow(comparisons),
                   dimnames = list(rownames(B), comparisons$comparison))
  neglog10_padj <- effect

  for (i in seq_len(nrow(B))) {
    if (progress_every > 0 && i %% progress_every == 0) {
      message(sprintf("        LV %d / %d", i, nrow(B)))
    }
    for (tissue in tissues) {
      df <- frames[[tissue]]$df
      df$B <- B[i, df$sample]
      df$Diff <- df$B - df$B[frames[[tissue]]$baseline]
      post <- df[df$Timepoint != "Pre", ]

      tested <- ggpubr::compare_means(
        Diff ~ Modality, data = post, ref.group = "CON", group.by = "Timepoint"
      )
      tested$comparison <- paste(tissue, tested$group2, tested$Timepoint, sep = "_")

      cells <- comparisons[comparisons$Tissue == tissue, ]
      for (cell in cells$comparison) {
        arm <- cells[cell, "Modality"]
        tp  <- cells[cell, "Timepoint"]
        in_arm  <- post$Timepoint == tp & post$Modality == arm
        in_ctrl <- post$Timepoint == tp & post$Modality == "CON"
        effect[i, cell] <- stats::median(post$Diff[in_arm]) -
          stats::median(post$Diff[in_ctrl])
        hit <- match(cell, tested$comparison)
        neglog10_padj[i, cell] <- if (is.na(hit)) NA_real_ else -log10(tested$p.adj[hit])
      }
    }
  }

  list(effect = effect, neglog10_padj = neglog10_padj, comparisons = comparisons)
}

# ---- palettes and annotation ------------------------------------------------

#' The four annotation palettes, under the display labels, from the package
#' constants and keyed on the package's own names.
plier_annotation_colors <- function() {
  tissue <- MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS[names(PLIER_TISSUE_LABELS)]
  names(tissue) <- unname(PLIER_TISSUE_LABELS)

  timepoint <- MotrpacHumanPreSuspensionAnalysis::HUMAN_ACUTE_TIMEPOINT_COLORS[
    names(PLIER_TIMEPOINT_LABELS)
  ]
  names(timepoint) <- unname(PLIER_TIMEPOINT_LABELS)

  modality <- MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS[
    names(PLIER_MODALITY_LABELS)
  ]
  names(modality) <- unname(PLIER_MODALITY_LABELS)

  sex <- MotrpacHumanPreSuspensionAnalysis::HUMAN_SEX_COLORS[c("Male", "Female")]

  list(Tissue = tissue, Sex = sex, Timepoint = timepoint, Modality = modality)
}

#' The bubble heatmap FIG4G and ED6E both draw: one circle per LV x comparison,
#' coloured by the median difference from control and sized by -log10 adjusted p,
#' outlined where that p clears 0.05.
#'
#' The two legacy panels differ in two details and both are kept as arguments
#' rather than unified: the RNA panel outlines its cells in black and draws the
#' bubble at the full radius, the metabolomics panel outlines in grey and divides
#' the radius by 1.4. Unifying them is a decision about the figure, not a port.
plier_bubble_heatmap <- function(effect, size, comparisons,
                                 cell_border = "black",
                                 radius_divisor = 1,
                                 cell_size = grid::unit(6, "mm"),
                                 size_cap = 3,
                                 effect_limit = 0.5,
                                 base_fontsize = 14) {
  stopifnot(identical(dimnames(effect), dimnames(size)))

  annotation_df <- comparisons[colnames(effect), c("Tissue", "Modality", "Timepoint"), drop = FALSE]
  # As factors, in study order. The legacy handed ComplexHeatmap a frame of
  # plain character columns, so its annotation legends came out alphabetical:
  # P24H above P3.5/4H, in a figure whose whole subject is time.
  annotation_df$Tissue <- factor(annotation_df$Tissue, levels = PLIER_TISSUE_ORDER)
  annotation_df$Modality <- factor(annotation_df$Modality, levels = PLIER_MODALITY_ORDER)
  annotation_df$Timepoint <- factor(annotation_df$Timepoint, levels = PLIER_TIMEPOINT_ORDER)
  annotation_df <- droplevels(annotation_df)
  colors <- plier_annotation_colors()

  top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = annotation_df,
    col = colors[c("Tissue", "Modality", "Timepoint")],
    which = "column",
    border = TRUE,
    gap = grid::unit(2, "pt"),
    annotation_name_gp = grid::gpar(fontsize = 0.9 * base_fontsize),
    annotation_legend_param = list(
      border = TRUE,
      title_gp = grid::gpar(fontsize = 0.9 * base_fontsize, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 0.9 * base_fontsize)
    ),
    annotation_name_side = "left"
  )

  col_fun <- circlize::colorRamp2(
    c(-effect_limit, 0, effect_limit), c("blue", "white", "red")
  )

  heatmap <- ComplexHeatmap::Heatmap(
    matrix = matrix(NA_real_, nrow = nrow(effect), ncol = ncol(effect),
                    dimnames = dimnames(effect)),
    width  = cell_size * ncol(effect),
    height = cell_size * nrow(effect),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    rect_gp = grid::gpar(fill = NA, col = cell_border),
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    show_heatmap_legend = FALSE,
    top_annotation = top_annotation,
    cell_fun = function(j, i, x, y, width, height, fill) {
      fc <- effect[i, j]
      pv <- size[i, j]
      if (is.na(fc) || is.na(pv)) return(invisible(NULL))
      if (pv > size_cap) pv <- size_cap
      grid::grid.circle(
        x = x, y = y,
        r = grid::unit(pv / radius_divisor, "mm"),
        gp = grid::gpar(
          lwd = 2,
          fill = col_fun(fc),
          col = if (size[i, j] > -log10(0.05)) "black" else NA
        )
      )
    }
  )

  legends <- list(
    ComplexHeatmap::Legend(
      title = "log2FC", col_fun = col_fun,
      at = c(-effect_limit, 0, effect_limit),
      labels = as.character(c(-effect_limit, 0, effect_limit))
    ),
    ComplexHeatmap::Legend(
      title = "-log10(adj_p)", type = "points", pch = 16,
      size = grid::unit(c(0, 2, 4, 6), "mm"),
      labels = c("0", "1", "2", paste0(size_cap, "+"))
    )
  )

  list(heatmap = heatmap, legends = legends)
}

#' The AUC heatmap FIG4H and ED6F both draw: latent variables across, prior sets
#' down, a cell coloured by PLIER's AUC on a fixed 0-1 scale.
#'
#' pheatmap has no legend title, so the colour bar is drawn with none and "AUC"
#' is placed above it here: the heatmap goes in a viewport short of the top by
#' `title_mm`, and the label at the left edge of the legend column, which
#' heatmap_legend_left() reads off pheatmap's own widths.
plier_auc_heatmap <- function(mat, labels_row, title = "AUC",
                              title_mm = 5, fontsize = 10) {
  function() {
    heatmap <- pheatmap::pheatmap(
      mat,
      breaks = seq(0, 1, length.out = 101),
      color = gplots::colorpanel(101, "white", "firebrick"),
      cluster_cols = FALSE,
      angle_col = 0,
      labels_row = labels_row,
      silent = TRUE
    )$gtable

    draw_heatmap_with_legend_title(heatmap, title, title_mm = title_mm,
                                   fontsize = fontsize)
  }
}

# The display names the legacy wrote by hand, KEYED ON THE PATHWAY ID.
#
# The legacy carried them as bare vectors positionally matched to the rows — 23
# strings for FIG4H, 6 for ED6F, 5 per ED6Aii bar plot — so a refit that
# reordered one row relabelled every row after it, and ED6F's vector had already
# come adrift: it names Ceramides and SM, which its matrix no longer has, and
# omits Amino Acids and PC, which it does. Keyed, a name can only ever land on
# the pathway it was written for.
#
# The RNA strings are the FIG4H vector paired back to the ids of the union it was
# written against; the RefMet ones cover every class that survives the prior
# filter, spelled as the legacy spelled the ones it named. An id with no entry
# here falls through to the derived label, which is why this list only has to
# carry the names worth adjusting rather than all 14,200.
PLIER_PATHWAY_DISPLAY_NAMES <- c(
  "REACTOME_PROCESSING_OF_CAPPED_INTRON_CONTAINING_PRE_MRNA" = "Pre-mRNA containing Capped Intron Processing",
  "GOMF_TRANSCRIPTION_REGULATOR_ACTIVITY" = "Transcription Regulator Activity",
  "REACTOME_RNA_POLYMERASE_II_TRANSCRIPTION" = "RNA Polymerase II Transcription",
  "GOMF_SEQUENCE_SPECIFIC_DNA_BINDING" = "Sequence Specific DNA Binding",
  "GOBP_NUCLEAR_TRANSCRIBED_MRNA_CATABOLIC_PROCESS" = "Nuclear Transcribed mRNA Catabolic Process",
  "GOCC_PRERIBOSOME" = "Preribosome",
  "GOBP_RIBOSOMAL_LARGE_SUBUNIT_BIOGENESIS" = "Ribosomal Large Subunit Biogenesis",
  "MITOCARTA_Translation" = "Translation",
  "GOMF_RIBOSOME_BINDING" = "Ribosome Binding",
  "GOMF_SNORNA_BINDING" = "snoRNA Binding",
  "GOMF_G_PROTEIN_COUPLED_RECEPTOR_ACTIVITY" = "GPCR Activity",
  "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL" = "Lymphoid Cell Immunoregulatory Interactions",
  "GOCC_EXTERNAL_ENCAPSULATING_STRUCTURE" = "External Encapsulating Structure",
  "GOMF_DNA_BINDING_TRANSCRIPTION_FACTOR_ACTIVITY" = "DNA Binding TF Activity",
  "REACTOME_NGF_STIMULATED_TRANSCRIPTION" = "NGF Stimulated Transcription",
  "WP_GLUCOCORTICOID_RECEPTOR_PATHWAY" = "Glucocorticoid Receptor Pathway",
  "PID_ATF2_PATHWAY" = "ATF2 Pathway",
  "GOBP_RESPONSE_TO_FIBROBLAST_GROWTH_FACTOR" = "Fibroblast Growth Factor Response",
  "GOBP_REGULATION_OF_TRANSPORTER_ACTIVITY" = "Transporter Activity Regulation",
  "REACTOME_ATTENUATION_PHASE" = "Attenuation Phase",
  "GOMF_HEAT_SHOCK_PROTEIN_BINDING" = "Heat Shock Protein Binding",
  "GOCC_PHAGOPHORE_ASSEMBLY_SITE" = "Phagophore Assembly Site",
  "WP_GENES_RELATED_TO_PRIMARY_CILIUM_DEVELOPMENT_BASED_ON_CRISPR" = "Primary Cilium Development",
  "REFMET_Acyl carnitines" = "Acyl Carnitines",
  "REFMET_Amino acids" = "Amino Acids",
  "REFMET_Cer" = "Ceramides",
  "REFMET_LPC" = "LPC",
  "REFMET_O-PC" = "O-PC",
  "REFMET_O-PE" = "O-PE",
  "REFMET_PC" = "PC",
  "REFMET_PE" = "PE",
  "REFMET_Saturated FA" = "Saturated FA",
  "REFMET_SM" = "SM",
  "REFMET_TG" = "TG",
  "REFMET_Unsaturated FA" = "Unsaturated FA"
)

#' Turn a pathway or molecular-signature id into something readable.
#'
#' The hand-written name if PLIER_PATHWAY_DISPLAY_NAMES has one for this id,
#' otherwise derived from the id: the database prefix drops, underscores become
#' spaces, and the result is title case. Title case is what mangles an acronym —
#' "REFMET_PC" would read "Pc" — so every id whose name carries one is in the
#' list above.
plier_pathway_label <- function(ids) {
  out <- unname(PLIER_PATHWAY_DISPLAY_NAMES[ids])

  derived <- sub("^(BIOCARTA|KEGG_MEDICUS|PID|REACTOME|WP|GOBP|GOCC|GOMF|MITOCARTA|REFMET)_",
                 "", ids)
  derived <- gsub("_", " ", derived)
  derived <- tolower(derived)
  derived <- gsub("\\b([a-z])", "\\U\\1", derived, perl = TRUE)

  unnamed <- is.na(out)
  out[unnamed] <- derived[unnamed]
  out
}

# The order samples are laid out in for the marker heatmaps (ED6D, ED6I): tissue,
# then exercise group, then timepoint, then sex.
PLIER_SAMPLE_ORDER_MODALITY <- c("CON", "EE", "RE")

# ---- everything a panel needs off one arm -----------------------------------

#' One arm's fit, the selections figure_4/plier_lvs.env makes on it, and the order the
#' panels draw them in.
#'
#' Every panel of both figures opens with this call, so the LV set, the marker
#' cutoff and the display order are decided once. Panels then take what they
#' need: the bubble panels the response grid, the pathway panels $Uauc, the
#' marker heatmaps the input matrix.
#'
#' @return list with
#'   fit, input, response   the fitted objects
#'   lvs                    the selected LV numbers, in fit order
#'   lv_names               the same as "LV<n>", in DISPLAY order
#'   markers                named list of marker features per LV
#'   correlation            per-LV cross-tissue correlation, all LVs
plier_selection <- function(arm) {
  arm <- plier_check_arm(arm)
  fit <- plier_fit(arm)
  input <- plier_input(arm)
  response <- plier_response(arm)

  lvs <- plier_lvs(arm)
  plier_check_lvs(fit, lvs, arm)

  B <- fit$B
  rownames(B) <- plier_lv_names(nrow(B))
  correlation <- plier_tissue_correlation(B, input$metadata)

  markers <- plier_markers(
    fit, lvs, plier_marker_sd(arm),
    margin = if (arm == "rna") "column" else "row"
  )

  list(
    arm = arm,
    fit = fit,
    input = input,
    response = response,
    lvs = lvs,
    lv_names = plier_lv_display_order(correlation, names(markers)),
    markers = markers,
    correlation = correlation
  )
}

#' The LV number behind a display name ("LV33" -> 33).
plier_lv_index <- function(lv_names) as.integer(sub("^LV", "", lv_names))

#' Report what a panel drew.
#'
#' Every PLIER panel except ED6Aii, ED6D and ED6I prints this. It is the cheap
#' version of the check figure_4/plier_lvs.env's header asks for: after a refit,
#' LV33 is a different latent variable, and the marker counts are the first place
#' that shows.
plier_report_selection <- function(selection) {
  counts <- vapply(selection$markers[selection$lv_names], length, integer(1))
  message(sprintf("        %s: %s",
                  selection$arm,
                  paste(sprintf("%s(%d)", names(counts), counts), collapse = " ")))
  invisible(counts)
}
