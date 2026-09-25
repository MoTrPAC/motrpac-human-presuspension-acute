# FIG3_ED4.R — the computations Figure 3 and Extended Data 4 share.
#
# Extended Data 4 is Figure 3's supplement, so the two draw the same three
# comparisons at different scopes. Each of these was its own helper file, and
# each was read by both figures:
#
#   re_vs_ee_*          FIG3A, FIG3B, ED4B
#   opposite_ee_re_*    FIG3C, ED4F
#   clinical_omics_*    FIG3EFG, ED4G, and tables/ST3.R and tables/EDT1.R
#
# One frame and one rule per comparison. Splitting them would let a panel and
# its supplement label different features, or key a join differently, as the
# thresholds moved.
#
# Everything is namespace-qualified so a figure script that sources this file
# inherits no attached packages from it. %>% comes from lib/panel_export.R,
# which every figure script sources first.

# =============================================================================
# RE vs EE — the joins FIG3A, FIG3B and ED4B are built on
# =============================================================================
#
# Two tables, one shape. Both put the endurance arm and the resistance arm of
# the same comparison side by side in one row, and attach the direct
# Endur_vs_Resist contrast for that same unit:
#
#   re_vs_ee_feature_table()      one row per tissue x assay x feature
#   re_vs_ee_enrichment_table()   one row per tissue x assay x molecular signature

#' Feature-level logFC and significance for RE, EE and the RE-vs-EE contrast.
#'
#' Widens the exercise_with_controls contrasts (each exercise group against its
#' own time-matched control) so that one row carries logFC_ADUEndur,
#' logFC_ADUResist and the two group-wise adjusted p-values, then inner-joins the
#' Endur_vs_Resist contrast for the same tissue, assay, feature and timepoint.
#' Each feature is reduced to the timepoint where the RE-vs-EE contrast is most
#' significant, and labelled by which arms it is significant in.
#'
#' load_differential_analysis() resolves the per-tissue DA objects by name off
#' the search path, so MotrpacHumanPreSuspensionAnalysis has to be attached with
#' library() by the calling script; a :: call alone is not enough.
#'
#' @returns A data.frame with logFC_ADUEndur, logFC_ADUResist,
#'   adj_p_value_ADUEndur, adj_p_value_ADUResist, the RE-vs-EE logFC and
#'   adj_p_value, the three "Sig"/"Non-Sig" flags, and Significant_in.
re_vs_ee_feature_table <- function() {
  if (!"package:MotrpacHumanPreSuspensionAnalysis" %in% search()) {
    stop(
      "re_vs_ee_feature_table() needs MotrpacHumanPreSuspensionAnalysis on the ",
      "search path: add library(MotrpacHumanPreSuspensionAnalysis) to the ",
      "figure script.",
      call. = FALSE
    )
  }

  da <- MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(
    selected_omes = "all",
    selected_tissues = "all",
    single_matrix = TRUE,
    verbose = FALSE
  )

  # The direct RE-vs-EE contrast, one row per tissue x assay x feature x timepoint.
  ee_re_sig <- da |>
    dplyr::filter(contrast_type == "Endur_vs_Resist") |>
    dplyr::select(tissue, assay, feature_id, Timepoint, adj_p_value, logFC)

  # Each exercise group against its time-matched control, one column per group.
  per_group <- da |>
    dplyr::filter(contrast_type == "exercise_with_controls") |>
    dplyr::select(tissue, assay, feature_id, Timepoint, randomGroupCode,
                  logFC, adj_p_value) |>
    tidyr::pivot_wider(
      names_from = c("randomGroupCode"),
      values_from = c("logFC", "adj_p_value")
    )

  # merge() rather than a dplyr join: it is an inner join that also sorts by the
  # key columns, and that sort is what fixes tissue order downstream.
  ee_re_wide <- merge(
    per_group, ee_re_sig,
    by = c("tissue", "assay", "feature_id", "Timepoint")
  )

  ee_re_wide <- ee_re_wide |>
    dplyr::mutate(
      sig_contrast  = ifelse(adj_p_value < 0.05, "Sig", "Non-Sig"),
      sig_ee_single = ifelse(adj_p_value_ADUEndur < 0.05, "Sig", "Non-Sig"),
      sig_re_single = ifelse(adj_p_value_ADUResist < 0.05, "Sig", "Non-Sig")
    )

  ee_re_wide |>
    dplyr::mutate(
      Significant_in = dplyr::case_when(
        sig_contrast == "Sig" & sig_ee_single == "Sig" &
          sig_re_single == "Non-Sig" ~ "Only EE",
        sig_contrast == "Sig" & sig_ee_single == "Non-Sig" &
          sig_re_single == "Sig" ~ "Only RE",
        sig_contrast == "Sig" & sig_ee_single == "Sig" &
          sig_re_single == "Sig" ~ "Both Groups",
        sig_contrast == "Sig" & sig_ee_single == "Non-Sig" &
          sig_re_single == "Non-Sig" ~ "Neither",
        TRUE ~ "Not RE vs EE Hit"
      )
    ) |>
    # Collapse timepoints: keep the one where the RE-vs-EE contrast is strongest.
    dplyr::group_by(tissue, feature_id) |>
    dplyr::slice_min(adj_p_value) |>
    dplyr::ungroup()
}


#' Set-level CAMERA-PR statistics for RE, EE and the RE-vs-EE contrast.
#'
#' Splits CAMERA_RESULTS$contrast_short into an exercise group and a timepoint,
#' widens the exercise_with_controls contrasts to one column per group, and
#' inner-joins the Endur_vs_Resist statistics for the same tissue, assay,
#' signature and timepoint. The join drops during_20_min and during_40_min,
#' which only the endurance arm has.
#'
#' @returns A data.frame with t_ADUEndur, t_ADUResist, adj_p_value_ADUEndur,
#'   adj_p_value_ADUResist, and the RE-vs-EE t and adj_p_value.
re_vs_ee_enrichment_table <- function() {
  camera <- MotrpacHumanPreSuspensionAnalysis::CAMERA_RESULTS

  # The left-hand side of "<group>.<timepoint> - <reference>" carries both the
  # exercise group and the timepoint, so it is split once and read twice.
  contrast_left <- stringr::str_split_fixed(camera$contrast_short, " - ", 2)[, 1]
  left_parts <- stringr::str_split_fixed(contrast_left, "\\.", 2)

  camera_tp_group <- camera |>
    dplyr::mutate(
      randomGroupCode = dplyr::case_when(
        contrast_type == "Endur_vs_Resist" ~ "ADUEndur - ADUResist",
        TRUE ~ paste0("ADU", left_parts[, 1])
      ),
      Timepoint = as.factor(left_parts[, 2])
    )

  enrich_results <- camera_tp_group |>
    dplyr::filter(contrast_type == "Endur_vs_Resist") |>
    dplyr::select(tissue, assay, set_id, set_short, adj_p_value,
                  contrast_short, t, randomGroupCode, Timepoint)

  per_group <- camera_tp_group |>
    dplyr::filter(contrast_type == "exercise_with_controls") |>
    dplyr::select(tissue, assay, set_id, set_short, Timepoint,
                  randomGroupCode, t, adj_p_value) |>
    tidyr::pivot_wider(
      names_from = c("randomGroupCode"),
      values_from = c("t", "adj_p_value")
    )

  merge(
    per_group, enrich_results,
    by = c("tissue", "assay", "set_id", "set_short", "Timepoint")
  )
}


#' Shorten a molecular-signature label for use as an axis label.
#'
#' @param x character; signature names.
#' @param max_chars integer; labels longer than this are cut and given a
#'   trailing ellipsis.
truncate_set_label <- function(x, max_chars = 35L) {
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1L, max_chars), "..."), x)
}


# =============================================================================
# Opposite-direction EE/RE — the frame FIG3C and ED4F are both drawn from
# =============================================================================
#
# The two panels are one computation at two scopes: FIG3C draws every tissue,
# ED4F draws muscle alone. Opposite_EE_vs_RE.Rmd drew them as
# `opposite_ee_vs_re_all_per_tissue.pdf` and its muscle-only sibling.

TIMEPOINT_SHORT = c(
"pre_exercise"      = "PRE",
"during_20_min"     = "D20M",
"during_40_min"     = "D40M",
"post_10_min"       = "P10M",
"post_15_30_45_min" = "P15-45M",
"post_3.5_4_hr"     = "P3.5/4H",
"post_24_hr"        = "P24H"
)

OME_LABELS = c(
"metab"              = "Metabolomics",
"transcript-rna-seq" = "Transcriptomics",
"prot-pr"            = "Proteomics",
"prot-ol"            = "Proteomics (Olink)",
"prot-ph"            = "Phosphoproteomics"
)

#' Features significantly up in one exercise arm and down in the other.
#'
#' One row per feature x tissue x timepoint, restricted to features the direct
#' Endur_vs_Resist contrast also calls significant at the same timepoint, and
#' carrying the hand-set label rule from the Rmd.
#'
#' @param selected_tissues Tissues to keep, or "all".
opposite_ee_re_table = function(selected_tissues = "all") {
  # epigen = FALSE: the panels label features by gene symbol or metabolite name,
  # and neither ATAC peaks nor methylation regions are labelled that way. The Rmd
  # loaded without them for the same reason.
  all_da = MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(
  selected_omes = "all",
  selected_tissues = "all",
  epigen = FALSE,
  single_matrix = TRUE,
  combine_with_featgene = TRUE,
  verbose = FALSE
)

  # The contrast label carries both sides of the comparison:
  #   "Endur.post_24_hr - Control.post_24_hr (delta-delta)"
  # Split on the first " - " to get the arm and timepoint being tested. The Rmd
  # did this with two nested tidyr::separate() calls; separate() is superseded,
  # and the timepoint contains dots ("post_3.5_4_hr") which is why the arm split
  # is anchored on the FIRST dot rather than on all of them.
  contrast_key = MotrpacHumanPreSuspensionAnalysis::CONTRAST_CONVERTER %>%
    dplyr::mutate(
      .lhs = sub(" - .*$", "", .data$contrast_short),
      key_group = sub("\\..*$", "", .data$.lhs),
      key_timepoint = sub("^[^.]*\\.", "", .data$.lhs)
    ) %>%
    dplyr::select("contrast", "key_group", "key_timepoint")

  # One row per feature x tissue x timepoint, the two arms side by side.
  per_arm = all_da %>%
    dplyr::filter(.data$contrast_type == "exercise_with_controls") %>%
    dplyr::inner_join(contrast_key, by = "contrast") %>%
    dplyr::mutate(ex_effect = dplyr::case_when(
      .data$z.std > 0 & .data$adj_p_value < 0.05 ~ "Up",
      .data$z.std < 0 & .data$adj_p_value < 0.05 ~ "Down",
      .default = "Neither"
    )) %>%
    dplyr::select("assay", "feature_id", "tissue", "key_timepoint", "key_group",
                  "ex_effect", "logFC") %>%
    tidyr::pivot_wider(
      id_cols = c("assay", "feature_id", "tissue", "key_timepoint"),
      names_from = "key_group",
      values_from = c("ex_effect", "logFC")
    )

  for (needed in c("ex_effect_Endur", "ex_effect_Resist",
                   "logFC_Endur", "logFC_Resist")) {
    if (!needed %in% names(per_arm)) {
      stop("the exercise-with-controls contrasts did not widen to ", needed,
           "\n  key_group carried: ",
           paste(sort(unique(contrast_key$key_group)), collapse = ", "),
           call. = FALSE)
    }
  }

  # The direct RE-vs-EE test for the same feature at the same timepoint. Joining
  # it restricts the panels to features the two arms are also significantly
  # different BETWEEN, not merely significant within each arm in opposite
  # directions.
  ee_re_hits = all_da %>%
    dplyr::filter(.data$contrast_type == "Endur_vs_Resist",
                  .data$adj_p_value < 0.05) %>%
    dplyr::inner_join(contrast_key, by = "contrast") %>%
    dplyr::select("assay", "feature_id", "tissue", "key_timepoint")

  # semi_join, not inner_join: this is a membership test and it takes no columns
  # from ee_re_hits, so a duplicated key there must not be able to duplicate a
  # point on the scatter.
  opposite = per_arm %>%
    dplyr::semi_join(ee_re_hits,
                     by = c("assay", "feature_id", "tissue", "key_timepoint")) %>%
    dplyr::filter(
      (.data$ex_effect_Endur == "Down" & .data$ex_effect_Resist == "Up") |
        (.data$ex_effect_Endur == "Up" & .data$ex_effect_Resist == "Down")
    )

  if (nrow(opposite) == 0) {
    stop("no feature responds in opposite directions to the two arms",
         call. = FALSE)
  }

  # gene_symbol comes off the DA table, which combine_with_featgene = TRUE has
  # already annotated. The Rmd re-joined HUMAN_FEATURE_TO_GENE here, which
  # duplicates rows wherever a feature maps to more than one gene.
  symbol_lookup = all_da %>%
    dplyr::select("assay", "feature_id", "tissue", "gene_symbol") %>%
    dplyr::distinct(.data$assay, .data$feature_id, .data$tissue,
                    .keep_all = TRUE)

  opposite = opposite %>%
    # many-to-one: a feature appears once per timepoint on the left and once in
    # the lookup, so one symbol row legitimately serves several points.
    dplyr::left_join(symbol_lookup, by = c("assay", "feature_id", "tissue"),
                     relationship = "many-to-one") %>%
    dplyr::mutate(
      timepoint_short = unname(TIMEPOINT_SHORT[.data$key_timepoint]),
      gene_symbol = as.character(.data$gene_symbol)
    )

  if (anyNA(opposite$timepoint_short)) {
    stop("timepoint with no short label: ",
         paste(sort(unique(opposite$key_timepoint[is.na(opposite$timepoint_short)])),
               collapse = ", "), call. = FALSE)
  }

  # Which points get a name. A hand-set distance rule from the Rmd, kept as
  # written: metabolites are labelled by feature_id because they have no gene
  # symbol, everything else by gene symbol, and the thresholds differ between the
  # two because the metabolite logFCs are on a wider scale.
  #
  # The rule is applied BEFORE any tissue filter, so ED4F labels the muscle
  # points FIG3C labels rather than promoting new ones into the gap.
  opposite = opposite %>%
    dplyr::mutate(plot_label = dplyr::case_when(
      .data$assay == "metab" &
        ((.data$logFC_Endur > 1 | .data$logFC_Resist < -1) |
           (.data$logFC_Endur < -0.5 & .data$logFC_Resist > 0.5)) ~
        as.character(.data$feature_id),
      .data$assay != "metab" &
        ((.data$logFC_Endur > 1 | .data$logFC_Resist < -1) |
           (.data$logFC_Endur < -0.3 & .data$logFC_Resist > 0.4)) ~
        .data$gene_symbol,
      .default = ""
    ))
  opposite$plot_label[is.na(opposite$plot_label)] = ""

  if (!identical(selected_tissues, "all")) {
    opposite = dplyr::filter(opposite,
                            as.character(.data$tissue) %in% selected_tissues)
    if (nrow(opposite) == 0) {
      stop("no opposite-direction feature in tissue(s): ",
           paste(selected_tissues, collapse = ", "), call. = FALSE)
    }
  }

  opposite
}

#' Look up palette entries by key, failing on one the palette does not carry.
lookup_colors = function(palette, keys, what) {
  hits = palette[keys]
  if (anyNA(hits)) {
    stop(what, " not in the package palette: ",
         paste(keys[is.na(hits)], collapse = ", "), call. = FALSE)
  }
  stats::setNames(unname(hits), keys)
}

#' The ome fill scale for one frame: display labels and their colours.
#'
#' HUMAN_OME_COLORS is not keyed by assay id: it carries the display names and
#' the SPECIFIC platform ids ("metab-u-hilicpos"), not the collapsed "metab" the
#' DA tables use. So the colour is looked up through the display name, which is
#' the same string the legend is labelled with.
ome_fill_scale = function(assays) {
  assays_present = sort(unique(as.character(assays)))
  unlabelled = setdiff(assays_present, names(OME_LABELS))
  if (length(unlabelled) > 0) {
    stop("assay with no display label: ", paste(unlabelled, collapse = ", "),
         "\n  Add it to OME_LABELS.", call. = FALSE)
  }
  labels = unname(OME_LABELS[assays_present])
  list(
    labels = labels,
    colors = stats::setNames(
    unname(lookup_colors(MotrpacHumanPreSuspensionAnalysis::HUMAN_OME_COLORS,
                         labels, "ome")),
    assays_present
  )
  )
}


# =============================================================================
# The clinical x omics fit, as FIG3EFG, ED4G, EDT1, ST3g and ST3h read it
# =============================================================================
#
# ED4G counts the significant associations; FIG3EFG plots the VO2peak column of
# them against the acute differential analysis; EDT1 reports selected features
# and does its own differential analysis join; ST3g and ST3h run ORA over them.
# All five read the table through clinical_omics_associations(); FIG3EFG also
# uses clinical_omics_with_da(), ST3g and ST3h clinical_omics_ora().
#
# The legacy duplicated the second of these verbatim across the two scripts
# (`vo2max_dumbbell_plot_combined.R:32-66` is a copy of
# `clinical_x_omics_landscape_fxn_final.R:365-399`), which is how a join key
# gets changed in one panel and not its neighbour.

#' The clinical x omics association table analysis/03_clinical_omics.R writes.
#'
#' @returns A data.frame: clinical_trait, tissue, omics_platform, feature_id,
#'   beta, se, t_stat, p_value, n_obs, adj_p.
clinical_omics_associations <- function() {
  override <- Sys.getenv("FIG3_CLINICAL_OMICS_TABLE", unset = "")
  path <- if (nzchar(override)) {
    override
  } else {
    Sys.getenv("CLINICAL_OMICS_TABLE", unset = "")
  }

  if (!nzchar(path)) {
    stop("neither FIG3_CLINICAL_OMICS_TABLE nor CLINICAL_OMICS_TABLE is set",
         "\n  config/landscape.env defaults the second.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("the clinical x omics fit is not on this machine: ", path,
         "\n  Build it with:  Rscript figures/landscape/analysis/03_clinical_omics.R",
         call. = FALSE)
  }

  associations <- utils::read.csv(path, sep = "\t", check.names = FALSE,
                                  stringsAsFactors = FALSE)

  required <- c("clinical_trait", "tissue", "omics_platform", "feature_id",
                "beta", "se", "t_stat", "p_value", "n_obs", "adj_p")
  missing <- setdiff(required, names(associations))
  if (length(missing) > 0) {
    stop("the clinical x omics fit is missing column(s): ",
         paste(missing, collapse = ", "),
         "\n  Refit with:  CLINICAL_OMICS_FORCE=TRUE Rscript ",
         "figures/landscape/analysis/03_clinical_omics.R", call. = FALSE)
  }
  associations
}

#' One trait x tissue x ome slice, joined to its acute differential analysis.
#'
#' One row per feature per timepoint, carrying the baseline association
#' statistics alongside the EE-CON, RE-CON and EE-RE contrasts for that feature
#' at that timepoint.
#'
#' Both sides are cut to one tissue and assay before the join on feature_id. The
#' legacy keyed it on tissue and feature_id alone, which is safe only for the one
#' slice it drew (muscle transcriptomics, where the ids are Ensembl gene ids that
#' collide with nothing) and would silently fan out on any other.
#'
#' @param trait,tissue,ome The slice to take.
clinical_omics_with_da <- function(trait, tissue, ome) {
  if (!"package:MotrpacHumanPreSuspensionAnalysis" %in% search()) {
    stop("clinical_omics_with_da() needs MotrpacHumanPreSuspensionAnalysis on ",
         "the search path: declare it in the panel's data_packages.",
         call. = FALSE)
  }

  slice <- clinical_omics_associations() |>
    dplyr::filter(clinical_trait == trait, tissue == !!tissue,
                  omics_platform == ome)
  if (nrow(slice) == 0) {
    stop("the fit has no rows for ", trait, " / ", tissue, " / ", ome,
         "\n  analysis/03_clinical_omics.R skipped or failed that combination; ",
         "its run reports which.", call. = FALSE)
  }

  da <- MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(
    selected_omes = ome,
    selected_tissues = tissue,
    single_matrix = TRUE,
    combine_with_featgene = TRUE,
    epigen = FALSE,
    verbose = FALSE
  ) |>
    dplyr::filter(contrast_type %in% c("exercise_with_controls",
                                       "Endur_vs_Resist")) |>
    dplyr::select("tissue", "assay", "feature_id", "gene_symbol",
                  "contrast_category", "Timepoint", "z.std", "adj_p_value")

  # One row per feature x timepoint, the three contrasts side by side.
  wide <- da |>
    dplyr::filter(contrast_category %in% c("EE-CON", "RE-CON", "EE-RE")) |>
    dplyr::mutate(is_sig = adj_p_value < 0.05) |>
    dplyr::select("feature_id", "gene_symbol", "Timepoint",
                  "contrast_category", "z.std", "is_sig") |>
    tidyr::pivot_wider(
      id_cols = c("feature_id", "gene_symbol", "Timepoint"),
      names_from = "contrast_category",
      values_from = c("is_sig", "z.std"),
      names_sep = "_"
    )

  dplyr::inner_join(wide, slice, by = "feature_id",
                    relationship = "many-to-one")
}

# ---- the ORA over the fit ---------------------------------------------------

# An ORA on a single gene is not a test. Combinations below this are skipped
# with a logged input size rather than reported; see ST3g.
CLINICAL_OMICS_ORA_MIN_INPUT <- 2L

#' Gene symbols for one tissue x ome slice of the fit.
#'
#' HUMAN_FEATURE_TO_GENE is keyed on assay AND feature_id. Joining on
#' feature_id alone fans a feature out across every assay carrying that id.
#'
#' @param slice One trait x tissue x ome slice of clinical_omics_associations().
clinical_omics_symbols <- function(slice) {
  feature_to_gene <- as.data.frame(
    MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE
  ) |>
    dplyr::filter(assay %in% unique(slice$omics_platform)) |>
    dplyr::transmute(omics_platform = as.character(assay),
                     feature_id = as.character(feature_id),
                     gene_symbol = as.character(gene_symbol))

  slice |>
    dplyr::mutate(feature_id = as.character(feature_id)) |>
    dplyr::left_join(feature_to_gene,
                     by = c("omics_platform", "feature_id"),
                     relationship = "many-to-one")
}

#' Over-representation analysis of one trait x tissue x ome slice of the fit.
#'
#' Two calls, one per sign of the association. The background is every gene the
#' slice fitted, so it is what was measured in that tissue and ome rather than
#' what responded.
#'
#' database is run_ORA()'s default, the full MOLECULAR_SIGNATURES roster. The
#' metabolite and phosphosite collections carry no gene symbols and fall out on
#' min_size, leaving C2, C5, CELLMARKER and MITOCARTA.
#'
#' @param trait,tissue,ome The slice to test.
#' @returns A data.frame: run_ORA()'s columns, prefixed by the three slice keys
#'   and `direction`. Zero rows for a direction with too few genes to test.
clinical_omics_ora <- function(trait, tissue, ome) {
  slice <- clinical_omics_associations() |>
    dplyr::filter(clinical_trait == trait, tissue == !!tissue,
                  omics_platform == ome)
  if (nrow(slice) == 0) {
    stop("the fit has no rows for ", trait, " / ", tissue, " / ", ome,
         "\n  analysis/03_clinical_omics.R skipped or failed that combination; ",
         "its run reports which.", call. = FALSE)
  }

  mapped <- clinical_omics_symbols(slice)
  fitted <- mapped |>
    dplyr::filter(!is.na(adj_p), !is.na(gene_symbol))
  background <- unique(fitted$gene_symbol)

  if (length(background) == 0) {
    stop("no feature of ", trait, " / ", tissue, " / ", ome,
         " maps to a gene symbol",
         "\n  HUMAN_FEATURE_TO_GENE has no rows for assay '", ome, "'.",
         call. = FALSE)
  }

  results <- lapply(c("positive", "negative"), function(direction) {
    sign_matches <- if (direction == "positive") {
      fitted$t_stat > 0
    } else {
      fitted$t_stat < 0
    }
    input <- unique(fitted$gene_symbol[fitted$adj_p < 0.05 & sign_matches])

    # Both sides come from one frame, so a stray input means the background was
    # built from a different set of rows than the foreground.
    stray <- setdiff(input, background)
    if (length(stray) > 0) {
      stop("ORA input outside its own background: ",
           paste(utils::head(sort(stray), 3), collapse = ", "), call. = FALSE)
    }

    if (length(input) < CLINICAL_OMICS_ORA_MIN_INPUT) {
      message(sprintf("  %s | %s | %s | %s: %d gene(s) — skipped",
                      trait, tissue, ome, direction, length(input)))
      return(NULL)
    }

    message(sprintf("  %s | %s | %s | %s: %d gene(s) against %d",
                    trait, tissue, ome, direction, length(input),
                    length(background)))

    enrichment <- MotrpacHumanPreSuspensionAnalysis::run_ORA(
      input = input,
      background = background
    ) |>
      as.data.frame()

    # run_ORA counts the overlap but does not name it. ST3g and ST3h report the
    # member genes, so they are recovered from the signature the set came from,
    # in the signature's own order.
    enrichment$geneID <- vapply(
      seq_len(nrow(enrichment)),
      function(i) {
        members <- MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES[[
          as.character(enrichment$database[i])
        ]][[enrichment$set[i]]]
        paste(intersect(members, input), collapse = "/")
      },
      character(1)
    )

    dplyr::mutate(enrichment, clinical_trait = trait, tissue = tissue,
                  omics_platform = ome, direction = direction, .before = 1)
  })

  dplyr::bind_rows(results)
}

# ---- the VO2peak dumbbell panels --------------------------------------------

# Row cap. The panels are drawn from every transcript in the fit, so this is
# what keeps them a readable size. Applied on |t|, so it is symmetric in the
# sign of the VO2peak association — see vo2max_panel_rows().
VO2MAX_MAX_ROWS <- 70

# Z-score axis, fixed rather than derived from each block's own range. The three
# response categories are read side by side, so a segment has to mean the same
# distance in all of them; a per-block range silently rescales that.
VO2MAX_X_LIMITS <- c(-11, 11)

# The three response categories, one column of FIG3EFG each. `sort_on` names the
# contrast a category's rows are ordered by, which has to be a contrast that
# category's features are actually significant in.
VO2MAX_PANELS <- list(
  ee_specific = list(title = "EE-specific", sort_on = "z.std_EE-CON",
               keep = function(df) df[["is_sig_EE-CON"]] &
                 !df[["is_sig_RE-CON"]] & df[["is_sig_EE-RE"]]),
  re_specific = list(title = "RE-specific", sort_on = "z.std_RE-CON",
               keep = function(df) !df[["is_sig_EE-CON"]] &
                 df[["is_sig_RE-CON"]] & df[["is_sig_EE-RE"]]),
  shared = list(title = "EE & RE shared", sort_on = "z.std_EE-CON",
               keep = function(df) df[["is_sig_EE-CON"]] &
                 df[["is_sig_RE-CON"]])
)

# The three muscle biopsy timepoints, top to bottom.
VO2MAX_TIMEPOINTS <- c(
  "post_15_30_45_min" = "post 15 min",
  "post_3.5_4_hr"     = "post 3.5 hr",
  "post_24_hr"        = "post 24 hr"
)

#' The rows one response category draws, split by timepoint.
#'
#' @param panel One of names(VO2MAX_PANELS).
#' @returns A named list of data.frames, one per timepoint in
#'   VO2MAX_TIMEPOINTS, each with a row_num column.
vo2max_panel_rows <- function(panel) {
  definition <- VO2MAX_PANELS[[panel]]
  if (is.null(definition)) {
    stop("no dumbbell definition for '", panel, "'", call. = FALSE)
  }

  combined <- clinical_omics_with_da("VO2max_L", "muscle", "transcript-rna-seq")

  # Significant baseline association, then the panel's own response pattern.
  # A missing significance flag is a contrast the feature has no row for, which
  # is not the same as a non-significant one; treated as not significant here,
  # as the legacy did, but stated rather than left to NA propagation.
  for (col in c("is_sig_EE-CON", "is_sig_RE-CON", "is_sig_EE-RE")) {
    if (!col %in% names(combined)) {
      stop("the differential analysis has no ", col, " column", call. = FALSE)
    }
    combined[[col]][is.na(combined[[col]])] <- FALSE
  }

  combined <- combined[!is.na(combined$adj_p) & combined$adj_p < 0.05, ,
                       drop = FALSE]
  combined <- combined[definition$keep(combined), , drop = FALSE]

  unknown <- setdiff(unique(combined$Timepoint), names(VO2MAX_TIMEPOINTS))
  if (length(unknown) > 0) {
    stop("timepoint not in the stated panel order: ",
         paste(sort(unknown), collapse = ", "), call. = FALSE)
  }

  # The 70 STRONGEST VO2peak associations in either direction, applied to the
  # whole category before the timepoint split — the position the legacy's cap
  # occupied, ordered differently.
  #
  # The legacy did `arrange(t_stat) %>% dplyr::slice(1:70)`
  # (precovid-analyses PR #103,
  # figure_3/vo2max_dumbbell_plot_combined.R:137-138, :156-157, :175-176).
  # arrange() sorts ASCENDING, so that keeps the most NEGATIVE associations
  # rather than the strongest: a one-sided cut on the very variable the VO2peak
  # tile column is coloured by. It never bound in the legacy, which drew from a
  # short hand-picked gene list; drawing from the whole fit it binds on every
  # category, and it took RE-specific from 24.9% positive associations to 0.0%
  # and shared from 64.3% to 12.9% — tile columns uniformly blue for a reason
  # nothing on the figure disclosed.
  #
  # Ordering on -abs(t_stat) keeps the same 70-row budget and spends it on
  # effect size instead of sign.
  combined <- combined[order(-abs(combined$t_stat)), , drop = FALSE]
  if (nrow(combined) > VO2MAX_MAX_ROWS) {
    combined <- combined[seq_len(VO2MAX_MAX_ROWS), , drop = FALSE]
  }

  # Rows with no gene symbol are labelled by feature_id. The curated list used
  # to remove them implicitly; without it they reach the axis.
  combined$gene_symbol <- as.character(combined$gene_symbol)
  no_symbol <- is.na(combined$gene_symbol) | !nzchar(combined$gene_symbol)
  combined$gene_symbol[no_symbol] <- as.character(combined$feature_id[no_symbol])

  by_timepoint <- lapply(names(VO2MAX_TIMEPOINTS), function(tp) {
    rows <- combined[combined$Timepoint == tp, , drop = FALSE]
    if (nrow(rows) == 0) return(rows)
    rows <- rows[order(rows[[definition$sort_on]]), , drop = FALSE]
    rows$Timepoint <- factor(unname(VO2MAX_TIMEPOINTS[tp]),
                             levels = unname(VO2MAX_TIMEPOINTS))
    rows$row_num <- seq_len(nrow(rows))
    rows
  })
  stats::setNames(by_timepoint, unname(VO2MAX_TIMEPOINTS))
}

# ---- drawing one dumbbell panel ---------------------------------------------

VO2MAX_CONTRAST_COLORS <- c("EE vs. CON" = "#d95f02", "RE vs. CON" = "#1b9e77")
VO2MAX_SIG_COLORS <- c("Significant" = "black", "Not significant" = "grey60")

#' One timepoint's dumbbell row block: the VO2peak tile column beside the
#' EE/RE segment plot.
#'
#' @param rows One element of vo2max_panel_rows().
#' @param x_limits,t_limits Shared across all nine FIG3EFG blocks, so the blocks
#'   are comparable down and across the grid.
vo2max_row_block <- function(rows, x_limits, t_limits,
                             title = NULL, show_x_title = FALSE) {
  y_limits <- range(rows$row_num) + c(-0.5, 0.5)

  vo2_tile <- ggplot2::ggplot(rows,
                              ggplot2::aes(x = 1, y = .data$row_num,
                                           fill = .data$t_stat)) +
    ggplot2::geom_point(shape = 22, size = 6, color = "grey30", stroke = 0.25) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = t_limits,
      name = expression(atop("VO"[2] * "peak", "association")),
      guide = ggplot2::guide_colorbar(order = 4)
    ) +
    ggplot2::scale_x_continuous(limits = c(0.5, 1.5), breaks = 1,
                                labels = expression("VO"[2] * "peak"),
                                position = "top") +
    ggplot2::scale_y_continuous(limits = y_limits, breaks = rows$row_num,
                                labels = rows$gene_symbol, expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 10),
      axis.text.x.top = ggplot2::element_text(angle = 90, hjust = 0.5,
                                              vjust = 0.5, face = "bold",
                                              size = 10),
      legend.position = "right",
      plot.margin = ggplot2::margin(5, 8, 5, 4),
      plot.title = ggplot2::element_text(size = 14, face = "italic")
    ) +
    ggplot2::ggtitle(title)

  points <- dplyr::bind_rows(
    data.frame(row_num = rows$row_num, contrast = "EE vs. CON",
               z = rows[["z.std_EE-CON"]],
               sig = ifelse(rows[["is_sig_EE-CON"]], "Significant",
                            "Not significant"),
               stringsAsFactors = FALSE),
    data.frame(row_num = rows$row_num, contrast = "RE vs. CON",
               z = rows[["z.std_RE-CON"]],
               sig = ifelse(rows[["is_sig_RE-CON"]], "Significant",
                            "Not significant"),
               stringsAsFactors = FALSE)
  )
  points$contrast <- factor(points$contrast,
                            levels = names(VO2MAX_CONTRAST_COLORS))
  points$sig <- factor(points$sig, levels = names(VO2MAX_SIG_COLORS))

  # Two zero-size, fully transparent layers carrying only the legend keys. The
  # real point layers are split by significance so the two can take different
  # borders, which means neither can produce a complete legend on its own.
  legend_anchor <- data.frame(x = 0, row_num = min(rows$row_num))

  dumbbell <- ggplot2::ggplot(rows, ggplot2::aes(y = .data$row_num)) +
    ggplot2::geom_hline(yintercept = rows$row_num, linetype = "dotted",
                        color = "grey85", linewidth = 0.3) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$`z.std_EE-CON`, xend = .data$`z.std_RE-CON`,
                   yend = .data$row_num),
      color = "grey45", linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data = points[points$sig == "Significant", ],
      ggplot2::aes(x = .data$z, fill = .data$contrast),
      shape = 21, color = "black", size = 6, stroke = 1.2, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = points[points$sig == "Not significant", ],
      ggplot2::aes(x = .data$z, fill = .data$contrast),
      shape = 21, color = "grey60", size = 6, stroke = 0.4, alpha = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = data.frame(contrast = factor(names(VO2MAX_CONTRAST_COLORS),
                                          levels = names(VO2MAX_CONTRAST_COLORS)),
                        x = legend_anchor$x, row_num = legend_anchor$row_num),
      ggplot2::aes(x = .data$x, y = .data$row_num, fill = .data$contrast),
      shape = 21, color = "black", size = 0, alpha = 0, show.legend = TRUE
    ) +
    ggplot2::geom_point(
      data = data.frame(sig = factor(names(VO2MAX_SIG_COLORS),
                                     levels = names(VO2MAX_SIG_COLORS)),
                        x = legend_anchor$x, row_num = legend_anchor$row_num),
      ggplot2::aes(x = .data$x, y = .data$row_num, color = .data$sig),
      shape = 21, fill = "white", size = 0, alpha = 0, show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = VO2MAX_CONTRAST_COLORS, limits = names(VO2MAX_CONTRAST_COLORS),
      drop = FALSE, name = "Contrast",
      guide = ggplot2::guide_legend(order = 2, override.aes = list(
        shape = 21, fill = unname(VO2MAX_CONTRAST_COLORS), color = "black",
        size = 8, alpha = 1))
    ) +
    ggplot2::scale_color_manual(
      values = VO2MAX_SIG_COLORS, limits = names(VO2MAX_SIG_COLORS),
      drop = FALSE, name = "Significance",
      guide = ggplot2::guide_legend(order = 3, override.aes = list(
        shape = 21, fill = "white", color = c("black", "grey60"),
        size = c(8, 8), alpha = c(1, 1), stroke = c(1.3, 0.4)))
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    ggplot2::scale_y_continuous(limits = y_limits, breaks = rows$row_num,
                                labels = NULL, expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = x_limits, clip = "off") +
    ggplot2::labs(x = if (show_x_title) "Z-score" else NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 12),
      axis.title.x = ggplot2::element_text(size = 12),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      legend.position = "right",
      legend.box = "vertical",
      plot.margin = ggplot2::margin(5, 5, 5, 4),
      legend.title = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 14)
    )

  vo2_tile + dumbbell + patchwork::plot_layout(ncol = 2, widths = c(0.05, 1))
}

#' FIG3EFG: the three response categories side by side, the three timepoints
#' down, with one legend, one z axis and one VO2peak colour scale for all nine
#' blocks.
vo2max_dumbbell_combined <- function() {
  categories <- names(VO2MAX_PANELS)
  timepoints <- unname(VO2MAX_TIMEPOINTS)
  rows <- stats::setNames(lapply(categories, vo2max_panel_rows), categories)

  n_rows <- vapply(categories, function(category) {
    vapply(timepoints, function(tp) nrow(rows[[category]][[tp]]), integer(1))
  }, integer(length(timepoints)))
  if (sum(n_rows) == 0) {
    stop("FIG3EFG: no feature meets any category's significance pattern",
         call. = FALSE)
  }

  blocks_with_rows <- Filter(function(df) nrow(df) > 0,
                             unlist(rows, recursive = FALSE))

  # A z outside VO2MAX_X_LIMITS would be drawn into the margin rather than
  # dropped, because the blocks set clip = "off", so it is named here.
  all_z <- unlist(lapply(blocks_with_rows, function(df) {
    c(df[["z.std_EE-CON"]], df[["z.std_RE-CON"]])
  }))
  outside <- which(abs(all_z) > max(abs(VO2MAX_X_LIMITS)))
  if (length(outside) > 0) {
    message(sprintf(
      "        FIG3EFG: %d point(s) outside the fixed z axis [%g, %g] — max |z| = %.2f",
      length(outside), VO2MAX_X_LIMITS[1], VO2MAX_X_LIMITS[2],
      max(abs(all_z[outside]))))
  }

  # One VO2peak scale across all nine blocks. Identical limits are also what
  # lets guides = "collect" merge the nine colour bars into one.
  all_t <- unlist(lapply(blocks_with_rows, function(df) df$t_stat))
  t_abs <- ceiling(max(abs(all_t), na.rm = TRUE))

  # A grid row is one timepoint across all three categories, as tall as the
  # category with the most rows there. A timepoint no category has rows at is
  # dropped; an empty block within a kept row is left blank under its title.
  heights <- apply(n_rows, 1, max)
  kept <- timepoints[heights > 0]
  if (length(kept) < length(timepoints)) {
    message("        FIG3EFG: no rows at ",
            paste(setdiff(timepoints, kept), collapse = ", "))
  }

  blocks <- list()
  for (tp in kept) {
    for (category in categories) {
      title <- sprintf("%s (%s)", VO2MAX_PANELS[[category]]$title, tp)
      df <- rows[[category]][[tp]]
      blocks[[length(blocks) + 1]] <- if (nrow(df) == 0) {
        message("        FIG3EFG: no rows for ", title)
        ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle(title)
      } else {
        vo2max_row_block(df, VO2MAX_X_LIMITS, c(-t_abs, t_abs), title = title)
      }
    }
  }

  patchwork::wrap_plots(blocks, ncol = length(categories), guides = "collect",
                        heights = heights[heights > 0]) +
    patchwork::plot_annotation(
      caption = "Z-score",
      theme = ggplot2::theme(
        plot.caption = ggplot2::element_text(hjust = 0.5, size = 18)
      )
    )
}
