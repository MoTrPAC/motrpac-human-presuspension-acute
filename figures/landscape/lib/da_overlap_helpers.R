# da_overlap_helpers.R — the significant-DA overlap sets that the cross-tissue
# UpSet and shared-feature panels are all built on.
#
# Lives in lib/ rather than with a figure because it belongs to no one figure:
# FIG2B, FIG2Bii, FIG2C, FIG2E and FIG2F read it out of Figure 2, ED2C, ED2E and
# ED2F out of Extended Data 2, and ED4A out of Extended Data 4. Same rule that
# puts single_feature_helpers.R here.
#
# Three steps, in order:
#
#   cross_tissue_da()          every exercise-with-controls DA call, epigen included
#   significant_da()           the FDR < 0.05 subset, with display names attached
#   feature_presence_matrix()  one row per feature, 0/1 per tissue (or per
#                              tissue x exercise modality)
#
# Nine panels across three figures begin from these three steps. One definition,
# read by all of them, is the point of this file.
#
# Everything is namespace-qualified so a figure script that sources this file
# inherits no attached packages from it.

# ---- ome display names -----------------------------------------------------

# The assay ids the differential analysis emits, mapped to the names the panels
# label with. Olink and mass-spec proteomics share a display name: the figures
# report "Proteomics" as one row, and each platform keeps its own denominator.
#
# Stated as a lookup rather than a case_when() with a .default. The legacy
# .default was "Metabolomics", so any assay not named above it — prot-clinical,
# or anything a future data release adds — was silently counted and coloured as
# a metabolite.
OME_DISPLAY_NAMES <- c(
  "transcript-rna-seq"   = "Transcriptomics",
  "prot-pr"              = "Proteomics",
  "prot-ol"              = "Proteomics",
  "prot-ph"              = "Phosphoproteomics",
  "epigen-atac-seq"      = "Chromatin Accessibility (ATAC)",
  "epigen-methylcap-seq" = "Methylation"
)

# Every metabolomics platform is one display row. DA metabolomics rows read
# assay = "metab", with the platform in a separate column.
ome_display_name <- function(assay) {
  assay <- as.character(assay)
  out <- unname(OME_DISPLAY_NAMES[assay])
  out[is.na(out) & grepl("^metab", assay)] <- "Metabolomics"
  if (anyNA(out)) {
    stop("assay with no display name: ",
         paste(sort(unique(assay[is.na(out)])), collapse = ", "),
         "\n  Add it to OME_DISPLAY_NAMES in lib/da_overlap_helpers.R.",
         call. = FALSE)
  }
  out
}

# Ome order wherever an ome is a factor: ATAC, metabolomics, methylation,
# phosphoproteomics, proteomics, transcriptomics. This is the order the UpSet
# panels stack their bars in and the order ED2C lays its sub-panels out, so the
# same ome sits in the same position in every panel that separates them.
#
# It is NOT FIG2A's row order, which is a different figure making a different
# point and is stated separately in that script.
OME_LEVELS <- c(
  "Chromatin Accessibility (ATAC)",
  "Metabolomics",
  "Methylation",
  "Phosphoproteomics",
  "Proteomics",
  "Transcriptomics"
)

# ---- the differential analysis ---------------------------------------------

#' Every exercise-with-controls DA call, epigenomics included.
#'
#' load_differential_analysis() resolves the per-tissue DA objects by name off
#' the search path, so MotrpacHumanPreSuspensionAnalysis has to be attached with
#' library() by the calling script; a :: call alone is not enough. The manifest's
#' data_packages field does that attach.
#'
#' @param epigen Include the ATAC and methylCap tables. TRUE for every panel
#'   that draws a cross-tissue overlap; the tables are ~7.7 GB, downloaded
#'   from the public c2.0 CloudFront release on every call.
#' @returns A data.frame, one row per tissue x assay x feature x contrast.
cross_tissue_da <- function(epigen = TRUE) {
  if (!"package:MotrpacHumanPreSuspensionAnalysis" %in% search()) {
    stop(
      "cross_tissue_da() needs MotrpacHumanPreSuspensionAnalysis on the search ",
      "path: declare it in the panel's data_packages in config/panel_map.json.",
      call. = FALSE
    )
  }

  da <- MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(
    selected_omes = "all",
    selected_tissues = "all",
    epigen = epigen,
    combine_with_featgene = TRUE,
    single_matrix = TRUE,
    verbose = FALSE
  )

  # Both epigen assays must come back; stop rather than draw from five omes.
  if (epigen) {
    got <- unique(as.character(da$assay))
    missing_epigen <- setdiff(c("epigen-atac-seq", "epigen-methylcap-seq"), got)
    if (length(missing_epigen) > 0) {
      stop("epigen = TRUE but the differential analysis came back without ",
           paste(missing_epigen, collapse = " and "),
           ".\n  The tables come from the public c2.0 CloudFront release; ",
           "check it is reachable and rebuild.",
           call. = FALSE)
    }
  }

  # contrast_short arrives as "<arm>.<timepoint> - <control>"; the panels label
  # with the left side only. The legacy dropped the right side with the same
  # gsub in every script.
  da$contrast_short <- sub(" - .*", "", da$contrast_short)

  # Bare column names, not .data$ pronouns: this file attaches no packages, and
  # .data is only in scope where dplyr is on the search path.
  da |>
    dplyr::filter(contrast_type == "exercise_with_controls")
}

#' Tissue display names.
#'
#' Explicit, not toTitleCase(): the three tissue names are a closed set and an
#' unrecognised one has to fail here rather than reach a set column or a
#' heatmap's column grid.
#'
#' Exposed separately from significant_da() because the panels that draw a
#' heatmap select their rows from the SIGNIFICANT frame and then read values off
#' the FULL differential analysis, which carries only the lowercase `tissue`.
tissue_display <- function(tissue) {
  known <- c("adipose" = "Adipose", "blood" = "Blood", "muscle" = "Muscle")
  out <- unname(known[as.character(tissue)])
  if (anyNA(out)) {
    stop("tissue with no display name: ",
         paste(sort(unique(as.character(tissue)[is.na(out)])), collapse = ", "),
         call. = FALSE)
  }
  out
}

#' The significant subset, with the display columns the overlap panels group on.
#'
#' Adds five columns and changes none:
#'
#'   Ome            the display name, Olink and MS proteomics collapsed
#'   Tissue         title case, for labelling
#'   ex_modality    "EE" or "RE"
#'   tissue_modality  "<Tissue>_<ex_modality>"
#'   overlap_id     the identity a feature is counted under across tissues
#'
#' @param da Output of cross_tissue_da().
#' @param fdr Adjusted-p cut. 0.05 in every panel that calls this.
significant_da <- function(da, fdr = 0.05) {
  sig <- da |>
    dplyr::filter(adj_p_value < fdr)

  sig$Ome <- ome_display_name(sig$assay)
  sig$Tissue <- tissue_display(sig$tissue)

  # The arm is the part of contrast_short before the first dot: "Endur.post_24_hr".
  arm <- sub("\\..*$", "", sig$contrast_short)
  known_arm <- c("Endur" = "EE", "Resist" = "RE")
  sig$ex_modality <- unname(known_arm[arm])
  if (anyNA(sig$ex_modality)) {
    stop("exercise arm not recognised in contrast_short: ",
         paste(sort(unique(sig$contrast_short[is.na(sig$ex_modality)])),
               collapse = ", "),
         call. = FALSE)
  }
  sig$tissue_modality <- paste0(sig$Tissue, "_", sig$ex_modality)

  sig$overlap_id <- overlap_identity(sig)
  sig
}

#' The identity a feature is counted under when tissues are compared.
#'
#' A protein measured on both Olink and mass spec is one protein, and the
#' cross-tissue overlaps have to count it once: blood proteomics is Olink,
#' muscle and adipose are mass spec, so without this every shared protein reads
#' as two tissue-specific ones. UniProt is the identity the two platforms share.
#'
#' Kept as a SEPARATE COLUMN rather than overwriting feature_id, which is what
#' the legacy did:
#'
#'   mutate(feature_id = case_when(Ome == "Proteomics" ~ uniprot, .default = feature_id))
#'
#' For prot-pr that is a no-op — HUMAN_FEATURE_TO_GENE has feature_id == uniprot
#' on all 8,548 mass-spec rows. For prot-ol it is not: feature_id is an Olink
#' panel id ("OID20051") and uniprot is an accession ("P61978"), so the rewrite
#' left every Olink row keyed by something no prot-ol row downstream matches.
#'
#' On the current data that costs nothing measurable — no protein reaches the
#' three-tissue set at all, because blood ships Olink only and muscle and
#' adipose ship mass spec only, so every protein overlap is mediated by this one
#' mapping and the three-way intersection of it is empty. Keeping the identity in
#' its own column is simply the version that stays correct if that changes:
#' grouping on overlap_id gives the intended cross-platform count, and every join
#' back to the DA tables and to HUMAN_FEATURE_TO_GENE still uses the key those
#' tables actually share.
overlap_identity <- function(sig) {
  id <- as.character(sig$feature_id)
  is_olink <- sig$assay == "prot-ol"
  if (any(is_olink)) {
    uniprot <- as.character(sig$uniprot)
    # An Olink target with no accession keeps its OID: it cannot be matched to a
    # mass-spec protein, and blanking it would merge every such target into one.
    usable <- is_olink & !is.na(uniprot) & nzchar(uniprot)
    id[usable] <- uniprot[usable]
  }
  id
}

# ---- the presence matrix ---------------------------------------------------

#' One row per feature, 0/1 per set, plus the ome it belongs to.
#'
#' The frame ComplexUpset::upset() takes: indicator columns for the sets, and
#' the Ome column its bars are filled by.
#'
#' @param sig Output of significant_da().
#' @param by "Tissue" for the three-tissue overlaps (FIG2B, FIG2C, FIG2E, ED2C),
#'   "tissue_modality" for the six tissue x arm sets (ED4A).
#' @param sets The set columns, in the order the caller wants them. Named
#'   explicitly so a set that happens to be empty in one data release still gets
#'   a column, rather than the panel losing a set and silently renumbering.
#' @returns A data.frame: overlap_id, one integer column per set, Ome.
feature_presence_matrix <- function(sig, by = c("Tissue", "tissue_modality"),
                                    sets) {
  by <- match.arg(by)

  present <- unique(data.frame(
    overlap_id = sig$overlap_id,
    set = as.character(sig[[by]]),
    stringsAsFactors = FALSE
  ))

  unexpected <- setdiff(unique(present$set), sets)
  if (length(unexpected) > 0) {
    stop("significant DA in sets the panel does not declare: ",
         paste(sort(unexpected), collapse = ", "),
         "\n  Add them to the `sets` argument, in the order they should be drawn.",
         call. = FALSE)
  }

  out <- data.frame(overlap_id = sort(unique(present$overlap_id)),
                    stringsAsFactors = FALSE)
  for (set in sets) {
    out[[set]] <- as.integer(out$overlap_id %in% present$overlap_id[present$set == set])
  }

  # One ome per feature. A feature can only carry one assay, and the display
  # names collapse Olink onto mass spec, so this is a lookup and not a choice —
  # except for a protein seen on both platforms, where either row gives
  # "Proteomics". Taking the first is therefore not arbitrary the way the
  # legacy's distinct(feature_id, .keep_all = TRUE) over the whole joined frame
  # was.
  ome_by_id <- unique(data.frame(
    overlap_id = sig$overlap_id,
    Ome = sig$Ome,
    stringsAsFactors = FALSE
  ))
  collisions <- ome_by_id$overlap_id[duplicated(ome_by_id$overlap_id)]
  if (length(collisions) > 0) {
    stop(length(collisions), " feature(s) carry two omes under one overlap id, ",
         "e.g. ", paste(utils::head(sort(unique(collisions)), 3), collapse = ", "),
         "\n  overlap_identity() has merged features it should not have.",
         call. = FALSE)
  }
  out$Ome <- factor(ome_by_id$Ome[match(out$overlap_id, ome_by_id$overlap_id)],
                    levels = OME_LEVELS)

  out
}

#' The gene symbols behind a set of overlap ids, for an ORA input or background.
#'
#' Joins back to HUMAN_FEATURE_TO_GENE on the key that table actually uses, so
#' Olink rows resolve through their OID and mass-spec rows through their
#' accession. Both are reachable because overlap_id is a separate column from
#' feature_id — see overlap_identity().
#'
#' @param sig Output of significant_da(), or cross_tissue_da() for a background.
#' @param ids overlap_id values to keep. NULL keeps everything in `sig`.
overlap_gene_symbols <- function(sig, ids = NULL) {
  if (!is.null(ids)) {
    sig <- sig[sig$overlap_id %in% ids, , drop = FALSE]
  }
  symbols <- as.character(sig$gene_symbol)
  sort(unique(symbols[!is.na(symbols) & nzchar(symbols)]))
}

# ---- the cross-tissue ORA --------------------------------------------------

# The six tissue-membership sets FIG2B tests for enrichment, and the background
# each is tested against. A set's background is the genes measured in every
# tissue that set requires: a feature can only be called shared between blood
# and muscle if it was assayed in both.
#
# Adipose_and_Blood_only is deliberately absent. The legacy tested these six and
# not the seventh, and the curated selections downstream name only these six.
CROSS_TISSUE_ORA_SETS <- list(
  Muscle_only             = list(present = "Muscle",
                                 absent  = c("Adipose", "Blood")),
  Blood_only              = list(present = "Blood",
                                 absent  = c("Adipose", "Muscle")),
  Blood_and_Muscle_only   = list(present = c("Blood", "Muscle"),
                                 absent  = "Adipose"),
  Adipose_only            = list(present = "Adipose",
                                 absent  = c("Blood", "Muscle")),
  Adipose_and_Muscle_only = list(present = c("Adipose", "Muscle"),
                                 absent  = "Blood"),
  Adipose_Blood_Muscle    = list(present = c("Adipose", "Blood", "Muscle"),
                                 absent  = character(0))
)

#' Over-representation analysis of each cross-tissue membership set.
#'
#' One run_ORA() call per set in CROSS_TISSUE_ORA_SETS, bound into one frame
#' with a `contrast` column naming the set. This is ST2e, and FIG2B, FIG2Bii and
#' FIG2F are three readings of it.
#'
#' The foreground and the background are both taken from the DA table's own
#' gene_symbol column, which combine_with_featgene = TRUE has already attached.
#' The legacy built them from two different sources — the foreground from a
#' re-join to HUMAN_FEATURE_TO_GENE on a rewritten feature_id, the background
#' straight off all_da — which is how a foreground can end up filtered more
#' aggressively than its own background. One source for both is what stops that.
#'
#' @param da Output of cross_tissue_da().
#' @param sig Output of significant_da(da).
#' @param presence Output of feature_presence_matrix(sig, by = "Tissue", ...).
#' @returns A data.frame, run_ORA()'s columns plus `contrast` and
#'   `statistic_column` (= -log10(p_value), as FIG2Bii and FIG2F draw).
cross_tissue_ora <- function(da, sig, presence) {
  # Genes measured in each tissue, from the full DA table rather than the
  # significant subset: the background is what was assayed, not what responded.
  background_by_tissue <- lapply(
    c(Adipose = "adipose", Blood = "blood", Muscle = "muscle"),
    function(tis) {
      rows <- da[as.character(da$tissue) == tis, , drop = FALSE]
      symbols <- as.character(rows$gene_symbol)
      unique(symbols[!is.na(symbols) & nzchar(symbols)])
    }
  )

  results <- list()
  for (set_name in names(CROSS_TISSUE_ORA_SETS)) {
    definition <- CROSS_TISSUE_ORA_SETS[[set_name]]

    keep <- rep(TRUE, nrow(presence))
    for (tis in definition$present) keep <- keep & presence[[tis]] == 1
    for (tis in definition$absent)  keep <- keep & presence[[tis]] == 0
    member_ids <- presence$overlap_id[keep]

    input <- overlap_gene_symbols(sig, member_ids)

    # Intersected across every tissue the set requires. A gene shared between
    # blood and muscle had to be measurable in both to be found there.
    background <- Reduce(intersect, background_by_tissue[definition$present])

    # run_ORA() requires input to be a subset of background. The two are built
    # from the same column of the same table, so a foreground gene outside its
    # own background means the membership set and the background disagree about
    # which tissues a feature was measured in — a real inconsistency, not a
    # filtering artefact, and it stops here rather than being silently dropped.
    stray <- setdiff(input, background)
    if (length(stray) > 0) {
      stop(set_name, ": ", length(stray), " foreground gene(s) outside their ",
           "own background, e.g. ", paste(utils::head(sort(stray), 3), collapse = ", "),
           call. = FALSE)
    }
    if (length(input) == 0) {
      stop(set_name, ": no gene symbols in the membership set", call. = FALSE)
    }

    results[[set_name]] <- MotrpacHumanPreSuspensionAnalysis::run_ORA(
      input = input,
      background = background,
      database = names(MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES)
    )
  }

  out <- dplyr::bind_rows(results, .id = "contrast")

  # -log10, as ED2E and the Figure 5 heatmaps use. The legacy took the natural
  # log here, which put this ORA on a different scale from every panel drawn
  # from it.
  out$statistic_column <- -log10(out$p_value)
  out$contrast <- factor(out$contrast, levels = names(CROSS_TISSUE_ORA_SETS))
  out
}
