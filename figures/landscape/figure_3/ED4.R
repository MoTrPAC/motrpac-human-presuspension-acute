#!/usr/bin/env Rscript
# Extended Data 4 — Resistance vs endurance exercise supplement
#
# Panels:  ED4A  differential features shared across tissue and exercise modality
#          ED4B  group-level pathway enrichment triangle, RE vs EE
#                                                        (also writes ST3d)
#          ED4C  REACTOME_HSF1_ACTIVATION feature heatmap, all three tissues
#          ED4D  HSPH1 muscle transcript trajectory
#          ED4E  diverging feature grid: FHOD1, MYOG, ZBTB4, NOL3, MYOT and an
#                lncRNA
#          ED4F  muscle features responding in opposite directions to EE and RE
#          ED4G  baseline features associated with each clinical trait
#                                                        (also writes ST3f)
# Tables:  ST3d  triangle selected pathways              (written by ED4B)
#          ST3f  clinical x omics all significant associations  (written by ED4G)
#
# Only ED4G needs consortium data access. ED4A reads one
# load_differential_analysis(epigen = TRUE) pass, which downloads the ATAC and
# methylCap DA tables (about 7.7 GB) from the public c2.0 CloudFront release on
# every run; nothing is cached.
#
# ED4G reads the baseline clinical x omics fit, built from
# MotrpacHumanPreSuspensionData and not in this repository: run
#
#   Rscript figures/landscape/figure_3/clinical_omics_fit.R
#
# once and it is cached. Without it ED4G reports SKIPPED and the other six
# panels still build.
#
#   Rscript figures/landscape/ED4.R          every panel
#   Rscript figures/landscape/ED4.R ED4F     one panel

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(ComplexUpset)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(tidyr)
  # Attached, not just namespace-qualified: load_differential_analysis() and
  # plot_feature_heatmap() resolve the per-tissue DA objects by name off the
  # search path.
  library(MotrpacHumanPreSuspensionAnalysis)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(here, "..", "lib", "table_export.R"))
source(file.path(here, "..", "lib", "highlights.R"))
source(file.path(here, "..", "lib", "da_overlap_helpers.R"))
source(file.path(here, "..", "lib", "single_feature_helpers.R"))
source(file.path(here, "..", "single_feature_plots.R"))
source(file.path(here, "FIG3_ED4_helpers.R"))

# ---- shared ----------------------------------------------------------------

# Nothing is shared here, deliberately. The seven panels read seven different
# things: ED4A the epigen-inclusive DA tables, ED4B the CAMERA-PR results, ED4C
# one signature's features, ED4D and ED4E single features, ED4F the muscle slice
# of the opposite-direction frame and ED4G the clinical x omics fit. The two
# frames ED4B and ED4F share with Figure 3 are built in helpers/FIG3_ED4.R, and
# ED4F takes selected_tissues = "muscle" where FIG3C takes "all", so even that
# is not one load.

# ---- ED4A — differential features shared across tissue and modality --------

# FIG2B's UpSet with the exercise arm added to the set definition: six sets,
# each a tissue crossed with endurance or resistance exercise, so an
# intersection says both where a feature responds and to which kind of
# exercise. Bars are stacked by ome and intersections smaller than 20 features
# are dropped, which is what keeps the panel readable at six sets.
#
# It sits in Extended Data 4 rather than with Figure 2 because it is an RE-vs-EE
# comparison, which is what Extended Data 4 supplements. The legacy filed it
# under figure_4/ and named it ED4a.

# Tissue-major, endurance before resistance within each. Stated so that the
# stripe colours below, which pair the two arms of a tissue, land on the rows
# they name.
MODALITY_SETS <- c(
  "Muscle_EE", "Muscle_RE",
  "Blood_EE", "Blood_RE",
  "Adipose_EE", "Adipose_RE"
)

ed4a <- function() {
  panel_init("ED4A")

  all_da <- cross_tissue_da(epigen = TRUE)
  sig <- significant_da(all_da)

  presence <- feature_presence_matrix(sig, by = "tissue_modality",
                                      sets = MODALITY_SETS)

  # Both arms of a tissue carry that tissue's colour, so the six stripes read as
  # three pairs. Derived from the set names rather than written out, which is
  # what keeps them attached to the rows if MODALITY_SETS is reordered.
  stripe_colors <- unname(
    MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS[
      tolower(sub("_.*$", "", MODALITY_SETS))]
  )
  if (anyNA(stripe_colors)) {
    stop("tissue not in HUMAN_TISSUE_COLORS: ",
         paste(MODALITY_SETS[is.na(stripe_colors)], collapse = ", "),
         call. = FALSE)
  }

  # From the package palette. The legacy hardcoded the six hex values, five of
  # which matched HUMAN_OME_COLORS and one of which did not: Transcriptomics was
  # "#4477AA" against the package's "#377EB8", so this panel drew transcriptomics
  # in a different blue from every other panel in the manuscript.
  ome_colors <- MotrpacHumanPreSuspensionAnalysis::HUMAN_OME_COLORS[OME_LEVELS] |>
    stats::setNames(OME_LEVELS)
  if (anyNA(ome_colors)) {
    stop("ome not in HUMAN_OME_COLORS: ",
         paste(OME_LEVELS[is.na(ome_colors)], collapse = ", "), call. = FALSE)
  }

  upset_plot <- ComplexUpset::upset(
    presence[, c(MODALITY_SETS, "Ome")],
    MODALITY_SETS,
    encode_sets = FALSE,
    width_ratio = 0.2,
    height_ratio = 1,
    set_sizes = FALSE,
    min_size = 20,
    base_annotations = list(
      "Number of \ndifferential features" = intersection_size(
        counts = TRUE,
        mapping = aes(fill = Ome),
        text = list(size = 4)
      ) +
        scale_fill_manual(values = ome_colors)
    ),
    themes = upset_modify_themes(
      list(intersections_matrix = theme(axis.title.x = element_blank()))
    ),
    sort_sets = FALSE,
    wrap = TRUE,
    stripes = stripe_colors
  )

  export_panel(function() print(upset_plot), "ED4A")
}

# ---- ED4B — group-level pathway enrichment triangle, RE vs EE --------------

# A curated set of molecular signatures, transcriptomics only, one row per
# signature and one column per post-exercise timepoint within each tissue. Each
# cell carries two triangles: the endurance arm nudged up, the resistance arm
# nudged down, pointing up or down with the sign of that arm's CAMERA-PR t
# statistic and sized by its magnitude. The grey background tile marks
# timepoints where the direct RE-vs-EE contrast is significant.
#
# This panel also writes ST3d, which is not a summary of the panel but the
# panel's own plotting frame: one row per triangle. Anything that changed the
# figure would change the table in the same edit, which is the point of building
# them together.

ed4b <- function() {
  panel_init("ED4B")

  # The curated signatures, from config/highlights.json. set_id is the five-digit
  # string CAMERA_RESULTS carries.
  interesting_pathways <- highlight_pathways("ED4B")$set_id

  enrich_wide <- re_vs_ee_enrichment_table()
  enrich_wide$set_short <- as.character(enrich_wide$set_short)

  selected_pathways_long <- enrich_wide |>
    # Transcriptomics only; overlaying several assays on one grid is unreadable.
    dplyr::filter(set_id %in% interesting_pathways,
                  assay %in% c("transcript-rna-seq")) |>
    # A separate column, not an overwrite: the 35-character cap exists so the y
    # axis stays readable, and ST3d has to carry the signature's real name.
    dplyr::mutate(set_label = truncate_set_label(set_short)) |>
    dplyr::rename(t_overall = t, adj_p_value_overall = adj_p_value) |>
    # One row per signature x timepoint x exercise arm.
    tidyr::pivot_longer(
      cols = c(t_ADUResist, t_ADUEndur,
               adj_p_value_ADUResist, adj_p_value_ADUEndur),
      names_to = c(".value", "condition"),
      names_pattern = "(t|adj_p_value)_(ADUResist|ADUEndur)"
    ) |>
    dplyr::mutate(
      # Direction and magnitude within the single arm.
      shape = ifelse(t > 0, "Upregulated", "Downregulated"),
      size = abs(t),
      # Significance of the direct RE-vs-EE contrast, shown as the tile.
      tile_filt = ifelse(adj_p_value_overall < 0.05,
                         "Adj P <= 0.05", "Adj P > 0.05")
    ) |>
    # Larger triangles are drawn first so smaller ones stay visible on top.
    dplyr::arrange(dplyr::desc(size)) |>
    dplyr::mutate(
      Timepoint_labels = dplyr::case_when(
        Timepoint == "post_10_min" ~ "P10M",
        Timepoint == "post_15_30_45_min" ~ "P15-45M",
        Timepoint == "post_3.5_4_hr" ~ "P3.5/4H",
        Timepoint == "post_24_hr" ~ "P24H",
        TRUE ~ Timepoint
      ),
      Timepoint_labels = factor(
        Timepoint_labels,
        levels = c("P10M", "P15-45M", "P3.5/4H", "P24H")
      )
    ) |>
    droplevels()

  write_st3d(selected_pathways_long)

  triangle_plots <- ggplot(selected_pathways_long,
                           aes(x = Timepoint_labels, y = set_label)) +
    geom_tile(aes(alpha = tile_filt),
              fill = "grey60",
              colour = NA,
              height = 0.95, width = 0.6) +
    geom_point(
      data = subset(selected_pathways_long, condition == "ADUEndur"),
      aes(fill = condition, shape = shape, size = size),
      stroke = 0.8,
      position = position_nudge(y = 0.2)
    ) +
    geom_point(
      data = subset(selected_pathways_long, condition == "ADUResist"),
      aes(fill = condition, shape = shape, size = size),
      stroke = 0.8,
      position = position_nudge(y = -0.2)
    ) +
    scale_alpha_manual(
      values = c("Adj P <= 0.05" = 0.3, "Adj P > 0.05" = 0),
      name = "RE-EE"
    ) +
    scale_fill_manual(
      values = MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS,
      # The palette keys on randomGroupCode; the arms are labelled EE and RE
      # everywhere else in the manuscript, so relabel the keys and leave the
      # values the colours are looked up by alone.
      labels = c("ADUEndur" = "EE", "ADUResist" = "RE"),
      name = "Group"
    ) +
    scale_shape_manual(
      values = c("Upregulated" = 24, "Downregulated" = 25),
      name = "Direction"
    ) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(),
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(linewidth = 1),

      axis.text.y.right = element_text(colour = "black", size = 11.5),
      axis.text.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(colour = "black", angle = 45,
                                 hjust = 1, vjust = 1, size = 12),
      strip.text = element_text(size = 11),
      plot.title = element_text(size = 15),

      # All four guides on one line: every guide below is nrow = 1, and the key
      # and text sizes are set to fill the page width rather than to fit within
      # it — the row runs from the left edge to the right one.
      legend.position = "bottom",
      legend.box = "horizontal",
      # Left-justified so the row starts flush with the panel rather than being
      # centred on it.
      legend.justification = "left",
      legend.box.just = "left",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.key.size = unit(3, "mm"),
      # The gap BETWEEN the four guides, so they read as four separate
      # categories rather than one run of keys...
      legend.spacing.x = unit(2, "mm"),
      # ...while the keys WITHIN a guide stay tight to their own title.
      legend.key.spacing.x = unit(1, "mm"),
      legend.box.spacing = unit(1, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      # The plot margin is set wide enough for the rotated timepoint labels,
      # which the legend row does not need; pull it back out to the page edge
      # so the row starts flush and the width goes to the keys.
      legend.box.margin = margin(l = -16),
      plot.margin = margin(t = 5, r = 8, b = 5, l = 24)
    ) +
    guides(
      fill = guide_legend(
        order = 1,
        nrow = 1,
        override.aes = list(shape = 22, size = 3, colour = "black")
      ),
      shape = guide_legend(order = 2, nrow = 1,
                           override.aes = list(size = 3)),
      size = guide_legend(
        title = "z-std",
        order = 3,
        nrow = 1,
        override.aes = list(shape = 24, fill = "grey60", color = "black")
      ),
      alpha = guide_legend(
        order = 4,
        nrow = 1,
        override.aes = list(fill = "grey60", color = "black",
                            shape = 22, size = 3)
      )
    ) +
    scale_y_discrete(position = "right") +
    facet_grid(~ tissue, scales = "free_x", space = "free_x") +
    ggtitle("Group-level pathway enrichment")

  export_panel(triangle_plots, "ED4B")
}

# ---- ST3d — one row per triangle -------------------------------------------

# Exported before the plot, from the frame the plot is about to be drawn from.
# The row order is the frame's own, largest triangle first, which is the order
# the panel draws them in so that smaller ones stay visible on top.
#
# set_short here is the full signature name. The axis uses set_label, which is
# the same name capped at 35 characters; a table carrying that cap would be
# shipping a property of the figure's y axis.
write_st3d <- function(selected_pathways_long) {
  export_table(
    as.data.frame(selected_pathways_long)[, table_spec("ST3d")$columns,
                                          drop = FALSE],
    "ST3d"
  )
}

# ---- ED4C — REACTOME_HSF1_ACTIVATION feature heatmap -----------------------

# Every transcript in the HSF1-activation signature, drawn as its z-standardised
# exercise effect against the time-matched control. Columns are tissue x
# exercise group x timepoint; a cell is starred where its adjusted p value is
# below 0.05. The heat-shock chaperones (HSPA1A/B, HSPH1, DNAJB1, HSP90AA1) are
# the block the panel exists to show.

ed4c <- function() {
  panel_init("ED4C")

  # The signature, from config/highlights.json, resolved by name through
  # SET_TO_ID so a collection update moves the id without editing this file.
  hsf1_set <- highlight_pathways("ED4C")$set
  hsf1_set_id <- MotrpacHumanPreSuspensionAnalysis::SET_TO_ID %>%
    filter(set == hsf1_set) %>%
    pull(set_id)

  if (length(hsf1_set_id) != 1) {
    stop(hsf1_set, " resolved to ", length(hsf1_set_id),
         " set ids, expected exactly one", call. = FALSE)
  }

  # return_drawing = TRUE draws with newpage = FALSE onto the device
  # export_panel() opens, as cmeans_feature_heatmap() does.
  hsf1_heatmap <- MotrpacHumanPreSuspensionAnalysis::plot_feature_heatmap(
    set_id = hsf1_set_id,
    selected_ome = "transcript-rna-seq",
    return_drawing = TRUE
  )

  export_panel(function() hsf1_heatmap$draw(), "ED4C")
}

# ---- ED4D — HSPH1 muscle transcript trajectory -----------------------------

# Exercise-group means with 95% confidence intervals across the acute
# timepoints. A point is filled black where that timepoint's adjusted p value is
# below 0.05. HSPH1 is a heat-shock chaperone and the RE-vs-EE contrast is what
# puts it in this figure.

ed4d <- function() {
  panel_init("ED4D")
  export_panel(single_feature_plot("ED4D"), "ED4D")
}

# ---- ED4E — diverging feature grid, muscle ---------------------------------

# Six further features whose responses separate between the resistance and
# endurance exercise groups, each drawn as the exercise-group mean with a 95%
# confidence interval across the acute timepoints. A point is filled black where
# that timepoint's adjusted p value is below 0.05. The six sub-plots sit two per
# row under one collected legend.

ed4e <- function() {
  panel_init("ED4E")
  export_panel(single_feature_plot("ED4E"), "ED4E")
}

# ---- ED4F — muscle features responding in opposite directions --------------

# FIG3C restricted to muscle. Same frame, same label rule, same fill scale; the
# tissue outline goes because there is one tissue left, and the freed aesthetic
# carries the timepoint instead, which FIG3C cannot show while it is spending
# shape on the point outline.
#
# Ported from Opposite_EE_vs_RE.Rmd, the muscle-only sibling of
# the per-tissue plot FIG3C draws. The frame is built by helpers/FIG3_ED4.R,
# shared because ED4 is Figure 3's supplement and the two draw the same
# comparison.
#
# The label rule runs before the muscle filter, so this panel labels the muscle
# points FIG3C labels rather than promoting new features into the space the other
# tissues left.

# Filled shapes only, so the ome fill is visible on every point.
TIMEPOINT_SHAPES = c(21, 22, 24, 23, 25, 21, 22)

ed4f = function() {
  panel_init("ED4F")

  opposite = opposite_ee_re_table(selected_tissues = "muscle")

  ome_scale = ome_fill_scale(opposite$assay)

  # Timepoints in protocol order, not alphabetical, and only the ones with a point.
  timepoint_levels = intersect(unname(TIMEPOINT_SHORT),
                              unique(opposite$timepoint_short))
  opposite$timepoint_short = factor(opposite$timepoint_short,
                                   levels = timepoint_levels)

  if (length(timepoint_levels) > length(TIMEPOINT_SHAPES)) {
    stop("more timepoints than shapes: ", length(timepoint_levels), call. = FALSE)
  }
  timepoint_shapes = stats::setNames(
  TIMEPOINT_SHAPES[seq_along(timepoint_levels)], timepoint_levels
  )

  message(sprintf("        ED4F: %d muscle feature(s), %d labelled, timepoints %s",
                  nrow(opposite), sum(nzchar(opposite$plot_label)),
                  paste(timepoint_levels, collapse = ", ")))

  # fill is mapped at the top level so the label layer can recolour itself with
  # after_scale(fill). The zero lines take inherit.aes = FALSE so they pick up
  # neither scale.
  p = ggplot(opposite, aes(x = .data$logFC_Endur, y = .data$logFC_Resist,
                          fill = .data$assay, shape = .data$timepoint_short)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey80",
               inherit.aes = FALSE) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey80",
               inherit.aes = FALSE) +
    geom_point(size = 3, color = "black", stroke = 0.6) +
    geom_text_repel(
      aes(label = .data$plot_label, color = after_scale(fill)),
      size = 3.5, box.padding = 0.3, point.padding = 0.2,
      max.overlaps = Inf, max.iter = 10000, force = 2, force_pull = 1,
      min.segment.length = 0, segment.size = 0.3, segment.alpha = 0.6,
      show.legend = FALSE
    ) +
    scale_fill_manual(name = "Ome", values = ome_scale$colors,
                      labels = ome_scale$labels) +
    scale_shape_manual(name = "Timepoint", values = timepoint_shapes) +
    # The fill legend draws its keys as the ome colour rather than inheriting the
    # first timepoint's shape.
    guides(fill = guide_legend(override.aes = list(shape = 21)),
           shape = guide_legend(override.aes = list(fill = "grey60"))) +
    labs(
      title = "Skeletal muscle EE vs RE",
      x = "log2FC endurance",
      y = "log2FC resistance"
    ) +
    theme_minimal() +
    theme(text = element_text(size = 14))

  export_panel(p, "ED4F")
}

# ---- ED4G — baseline features associated with each clinical trait ----------

# The clinical x omics fit, counted: how many features on each omics platform are
# associated with each baseline clinical trait, at FDR < 0.05, in each tissue.
# Two label columns carry the platform and the trait, three heat columns carry
# the counts. Each tissue column is coloured on its own scale, count divided by
# that column's largest count, over ColorBrewer Blues, as in the published
# table. A combination the fit did not cover is a grey square with no number,
# not a 0.
#
# This panel also writes ST3f, every significant association the counts are
# drawn from.

TISSUE_COLUMNS <- c("adipose", "blood", "muscle")

# Platform and trait display names, and the row order.
# Olink and mass spec are one "Proteomics" row. They are disjoint by tissue —
# Olink is assayed in blood only, mass spec in muscle and adipose only — so
# collapsing them loses nothing and each tissue's column still comes from the
# platform that tissue was actually measured on. FIG2A merges them for the same
# reason.
PLATFORM_LABELS <- c(
  "metab"              = "Metabolomics",
  "prot-ol"            = "Proteomics",
  "prot-pr"            = "Proteomics",
  "transcript-rna-seq" = "Transcriptomics"
)

# Row order within a platform, and it is THIS order rather than alphabetical.
# The metabolic and cardiorespiratory traits lead, strength follows: sorting the
# labels put leg press between Lactate and NEFA, which reads as though the
# strength measures belong with them.
#
# The strength traits are the rest-corrected 1RM columns of the curated screening
# table (visit ADU_SCP). lp_rest_corrected_rmat is the legacy's "Leg press
# strength" input: the legacy function run on it reproduces the fit exactly.
TRAIT_LABELS <- c(
  "HOMA_IR"                = "HOMA-IR",
  "Lactate"                = "Lactate",
  "NEFA"                   = "NEFA",
  "VO2max_L"               = "VO2peak\n(L/min)",
  "lp_rest_corrected_rmat" = "Leg press\nstrength",
  "peaktorq_iske"          = "Peak torque\n(ISKE)",
  "cp_rest_corrected_rmat" = "Chest press\nstrength",
  "le_rest_corrected_rmat" = "Leg extension\nstrength",
  "peak_handgrip"          = "Peak handgrip\nstrength"
)

# The three 1RM strength traits were measured at screening in the resistance
# exercise participants only, so their models run on that arm alone. Starred on
# the panel; the per-model n goes to the build log. Peak handgrip and isokinetic
# peak torque were measured in every arm and are not starred.
LOW_N_TRAITS <- c("cp_rest_corrected_rmat", "lp_rest_corrected_rmat",
                  "le_rest_corrected_rmat")

BLUES_9 <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6",
             "#4292C6", "#2171B5", "#08519C", "#08306B")
UNFITTED_GREY <- "#BDBDBD"
COLUMN_X <- c(omics_platform = 1, clinical_trait = 2, adipose = 3, blood = 4, muscle = 5)
COLUMN_TITLES <- c("Omics platform", "Clinical trait", "Adipose", "Blood", "Muscle")
STAR_X <- 0.21
STAR_Y <- -0.17

# Each platform block is its own box, separated from the next by white space,
# with the platform named outside it on the left. Inside a box, mid-grey lines
# separate every cell and a dark rule separates the traits from the counts.
GAP <- 0.45
GRID <- list(colour = "grey55", linewidth = 0.45)
BOX <- list(colour = "black", linewidth = 1.1)

heat_cells <- function(rows) {
  rows |>
    select("row", "y", all_of(TISSUE_COLUMNS)) |>
    pivot_longer(all_of(TISSUE_COLUMNS), names_to = "tissue", values_to = "count") |>
    group_by(.data$tissue) |>
    mutate(column_max = suppressWarnings(max(.data$count, na.rm = TRUE)),
           fill_value = dplyr::case_when(is.na(.data$count) ~ NA_real_,
                                         .data$column_max > 0 ~ .data$count / .data$column_max,
                                         TRUE ~ 0)) |>
    ungroup() |>
    mutate(x = unname(COLUMN_X[.data$tissue]),
           label = ifelse(is.na(.data$count), "", as.character(.data$count)),
           label_colour = ifelse(!is.na(.data$fill_value) & .data$fill_value > 0.6,
                                 "#F1F1F1", "black"))
}

block_extents <- function(rows) {
  rows |>
    group_by(.data$block, .data$omics_platform) |>
    summarise(top = max(.data$y) + 0.5, bottom = min(.data$y) - 0.5,
              y = mean(.data$y), .groups = "drop")
}

text_layers <- function(rows, heat) {
  list(
    geom_text(data = rows, aes(x = 2, y = .data$y, label = .data$clinical_trait), size = 5),
    geom_text(data = rows[rows$low_n, , drop = FALSE],
              aes(x = 2 + STAR_X, y = .data$y + STAR_Y),
              label = "*", size = 9, hjust = 0, vjust = 0.6),
    geom_text(data = block_extents(rows),
              aes(x = 1, y = .data$y, label = .data$omics_platform),
              size = 5, fontface = "bold"),
    geom_text(data = heat, aes(x = .data$x, y = .data$y, label = .data$label,
                               colour = .data$label_colour), size = 6)
  )
}

finish <- function(p) {
  p +
    scale_fill_gradientn(colours = BLUES_9, limits = c(0, 1), na.value = UNFITTED_GREY) +
    scale_colour_identity() +
    scale_x_continuous(breaks = COLUMN_X, labels = COLUMN_TITLES, position = "top",
                       limits = c(0.45, 5.55), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(add = 0.05)) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), axis.text.x = element_text(face = "bold"),
          legend.position = "none")
}

ed4g <- function() {
  # The manifest marks CLINICAL_OMICS_TABLE must_exist for this panel, so an
  # unrun fit is reported with the command that produces it rather than dying on
  # a missing path. A fit produced elsewhere therefore has to be named in
  # CLINICAL_OMICS_TABLE, not only in FIG3_CLINICAL_OMICS_TABLE, for this panel
  # to run.
  skip <- skip_if_no_fit("ED4G")
  if (!is.null(skip)) return(skip)

  panel_init("ED4G")

  associations <- clinical_omics_associations()

  # Every fitted trait is a candidate row. Which ones are DRAWN is decided below,
  # by whether they found anything — not by this list. Read from
  # config/landscape.env so the panel and the fit cannot disagree about what was
  # meant to be fitted.
  declared_traits <- trimws(strsplit(
    Sys.getenv("CLINICAL_OMICS_TRAITS", unset = ""), ",", fixed = TRUE)[[1]])
  declared_traits <- declared_traits[nzchar(declared_traits)]
  if (length(declared_traits) == 0) {
    stop("CLINICAL_OMICS_TRAITS is not set; it is defaulted in config/landscape.env",
         call. = FALSE)
  }
  unlabelled <- setdiff(declared_traits, names(TRAIT_LABELS))
  if (length(unlabelled) > 0) {
    stop("trait declared in CLINICAL_OMICS_TRAITS with no display label: ",
         paste(unlabelled, collapse = ", "),
         "\n  Add it to TRAIT_LABELS.", call. = FALSE)
  }

  associations <- associations[associations$clinical_trait %in% declared_traits, ,
                               drop = FALSE]
  if (nrow(associations) == 0) {
    stop("the fit carries none of the declared traits: ",
         paste(declared_traits, collapse = ", "), call. = FALSE)
  }

  significant <- associations |>
    filter(!is.na(.data$adj_p), .data$adj_p < 0.05)

  write_st3f(significant)

  # ---- the count grid ----

  # The full platform x trait grid, not only the combinations that produced a
  # hit, so a fitted cell with no hits is drawn as 0 and an unfitted one as grey.
  # The legacy tallied the significant rows and pivoted, so a tissue with no hits
  # anywhere would have dropped its whole column and broken the select() below.
  # Rows with no hit in any tissue are still dropped, further down.
  grid <- expand.grid(
    omics_platform = sort(unique(unname(PLATFORM_LABELS))),
    clinical_trait = declared_traits,
    stringsAsFactors = FALSE
  )

  # Collapsed to the display platform BEFORE counting, so the two proteomics
  # assays add into one row rather than producing two.
  significant$platform <- unname(PLATFORM_LABELS[significant$omics_platform])
  counts <- significant |>
    count(.data$platform, .data$clinical_trait, .data$tissue, name = "n") |>
    dplyr::rename(omics_platform = "platform")

  unexpected <- setdiff(unique(counts$tissue), TISSUE_COLUMNS)
  if (length(unexpected) > 0) {
    stop("significant associations in tissue(s) the panel does not declare: ",
         paste(sort(unexpected), collapse = ", "), call. = FALSE)
  }

  wide <- counts |>
    pivot_wider(names_from = "tissue", values_from = "n") |>
    right_join(grid, by = c("omics_platform", "clinical_trait"))

  # A cell was fitted when at least one feature in that tissue x platform x trait
  # has a p-value. A combination the fit skipped has no rows at all; one where no
  # feature's model could be fitted (no complete observations, or a constant
  # trait or sex) has rows and no p-values. Neither found zero associations, so
  # neither is a 0: the cell stays NA and is drawn as a grey square.
  fitted_cells <- associations |>
    filter(!is.na(.data$p_value)) |>
    mutate(omics_platform = unname(PLATFORM_LABELS[.data$omics_platform])) |>
    distinct(.data$omics_platform, .data$clinical_trait, .data$tissue)

  for (tissue in TISSUE_COLUMNS) {
    if (!tissue %in% names(wide)) wide[[tissue]] <- NA_integer_
    fitted_here <- fitted_cells[fitted_cells$tissue == tissue, , drop = FALSE]
    is_fitted <- paste(wide$omics_platform, wide$clinical_trait) %in%
      paste(fitted_here$omics_platform, fitted_here$clinical_trait)
    wide[[tissue]][is.na(wide[[tissue]]) & is_fitted] <- 0L
  }

  # A platform x trait row fitted in no tissue at all is not in the figure.
  never_fitted <- rowSums(!is.na(wide[, TISSUE_COLUMNS, drop = FALSE])) == 0
  wide <- wide[!never_fitted, , drop = FALSE]
  if (nrow(wide) == 0) {
    stop("no platform x trait combination was fitted", call. = FALSE)
  }

  # Only rows that found something. A platform x trait with no significant
  # association in ANY tissue is dropped rather than drawn as three zeros: at nine
  # traits that was two thirds of the grid, and an empty row says nothing a reader
  # can act on. The row set is therefore data-driven — a trait that starts
  # producing hits appears on its own, with no list to edit.
  found <- rowSums(wide[, TISSUE_COLUMNS, drop = FALSE], na.rm = TRUE) > 0
  dropped <- wide[!found, , drop = FALSE]
  if (nrow(dropped) > 0) {
    message(sprintf(
      "        ED4G: %d platform x trait row(s) drawn, %d dropped for having no ",
      sum(found), nrow(dropped)),
      "significant association in any tissue:")
    for (i in seq_len(nrow(dropped))) {
      message(sprintf("          %-16s %s", dropped$omics_platform[i],
                      gsub("\n", " ", unname(TRAIT_LABELS[dropped$clinical_trait[i]]))))
    }
  }
  wide <- wide[found, , drop = FALSE]
  if (nrow(wide) == 0) {
    stop("no platform x trait combination has a significant association in any ",
         "tissue", call. = FALSE)
  }

  unfitted <- which(is.na(wide[, TISSUE_COLUMNS, drop = FALSE]), arr.ind = TRUE)
  if (nrow(unfitted) > 0) {
    message(sprintf("        ED4G: %d cell(s) drawn grey — combination does not have ",
                    nrow(unfitted)), "sufficient samples to estimate:")
    for (i in seq_len(nrow(unfitted))) {
      message(sprintf("          %-16s %-24s %s",
                      wide$omics_platform[unfitted[i, "row"]],
                      gsub("\n", " ",
                           unname(TRAIT_LABELS[wide$clinical_trait[unfitted[i, "row"]]])),
                      TISSUE_COLUMNS[unfitted[i, "col"]]))
    }
  }

  # ---- plot ----

  # Recoded BEFORE sorting, so the panel's row order is the order of the labels a
  # reader sees. Sorting on the raw names put VO2max after the lowercase strength
  # ids for no visible reason.
  plot_df <- wide |>
    mutate(low_n = .data$clinical_trait %in% LOW_N_TRAITS,
           clinical_trait = unname(TRAIT_LABELS[.data$clinical_trait])) |>
    arrange(.data$omics_platform,
            match(.data$clinical_trait, unname(TRAIT_LABELS))) |>
    mutate(row_id = factor(dplyr::row_number(), levels = rev(dplyr::row_number())))

  # Per-model n for the figure legend, which is where the grey squares and the
  # asterisk are explained: the most complete feature's n in each fitted
  # combination, ranged over the starred traits and over the rest.
  model_n <- associations |>
    filter(!is.na(.data$p_value)) |>
    group_by(.data$clinical_trait, .data$tissue, .data$omics_platform) |>
    summarise(n = max(.data$n_obs), .groups = "drop")
  low_n_range <- range(model_n$n[model_n$clinical_trait %in% LOW_N_TRAITS])
  other_n_max <- max(model_n$n[!model_n$clinical_trait %in% LOW_N_TRAITS])

  if (any(plot_df$low_n)) {
    message(sprintf(
      "        ED4G: starred traits were measured in resistance exercise participants only: n = %d-%d per model, versus up to %d for the other traits",
      low_n_range[1], low_n_range[2], other_n_max))
  }

  rows <- plot_df |>
    mutate(row = dplyr::row_number(),
           block = cumsum(c(TRUE, .data$omics_platform[-1] != .data$omics_platform[-dplyr::n()])))
  n_rows <- nrow(rows)
  rows$last_in_block <- c(rows$block[-1] != rows$block[-n_rows], TRUE)

  rows$y <- -(rows$row - 1) - (rows$block - 1) * GAP
  heat <- heat_cells(rows)
  blocks <- block_extents(rows)
  inner <- rows[!rows$last_in_block, , drop = FALSE]

  p <- ggplot() +
    geom_tile(data = heat, aes(x = .data$x, y = .data$y, fill = .data$fill_value),
              colour = GRID$colour, linewidth = GRID$linewidth) +
    text_layers(rows, heat) +
    geom_segment(data = inner, aes(x = 1.5, xend = 2.5, y = .data$y - 0.5, yend = .data$y - 0.5),
                 colour = GRID$colour, linewidth = GRID$linewidth) +
    geom_segment(data = blocks, aes(x = 2.5, xend = 2.5, y = .data$bottom, yend = .data$top),
                 colour = BOX$colour, linewidth = 0.8) +
    geom_rect(data = blocks, aes(xmin = 1.5, xmax = 5.5, ymin = .data$bottom, ymax = .data$top),
              fill = NA, colour = BOX$colour, linewidth = BOX$linewidth)
  p <- finish(p)

  export_panel(p, "ED4G")
}

# ---- ST3f — every significant association ----------------------------------

# gene_symbol comes from the package feature-to-gene table, keyed on assay and
# feature_id. Metabolites map to no gene and get an empty cell.
write_st3f <- function(significant) {
  feature_to_gene <- as.data.frame(
    MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE
  ) |>
    filter(.data$assay %in% unique(significant$omics_platform)) |>
    transmute(omics_platform = .data$assay,
              feature_id = as.character(.data$feature_id),
              gene_symbol = as.character(.data$gene_symbol))

  st3f <- significant |>
    left_join(feature_to_gene, by = c("omics_platform", "feature_id"),
              relationship = "many-to-one") |>
    arrange(.data$clinical_trait, .data$tissue, .data$omics_platform,
            .data$adj_p) |>
    as.data.frame()

  export_table(st3f[, table_spec("ST3f")$columns, drop = FALSE], "ST3f")
}

run_panels(list(ED4A = ed4a, ED4B = ed4b, ED4C = ed4c, ED4D = ed4d,
                ED4E = ed4e, ED4F = ed4f, ED4G = ed4g))
