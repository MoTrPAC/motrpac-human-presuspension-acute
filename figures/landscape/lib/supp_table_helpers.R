# supp_table_helpers.R — what the supplementary tables share across figures.
#
# Lives in lib/ rather than beside one supplementary table, because the
# things here belong to no single table: the assay roster and the order assays
# are written in are properties of the study, and a table that disagreed with
# another about either would be a defect rather than a choice.
#
# Sourced by every script that writes a per-assay sub-table, including the
# figure scripts that write one alongside a panel.

suppressPackageStartupMessages(library(MotrpacHumanPreSuspensionAnalysis))

#' The order every per-assay supplementary table is written in.
#'
#' Grouped by modality — epigenomics, transcriptomics, proteomics,
#' metabolomics — and alphabetical within each group, except proteomics, which
#' runs global, Olink, phospho: the order the platforms are introduced in the
#' manuscript rather than the order their codes happen to sort in.
#'
#' Written out rather than derived. There is no field in any package object that
#' encodes it, so deriving it would mean inferring modality from the assay code
#' and then hardcoding the modality order instead — the same list, one
#' indirection further from the reader. Being explicit also makes
#' order_by_assay() able to fail on an assay it has never seen, which a derived
#' order could not.
ASSAY_DISPLAY_ORDER <- c(
  "epigen-atac-seq", "epigen-methylcap-seq",
  "transcript-rna-seq",
  "prot-pr", "prot-ol", "prot-ph",
  # The differential-analysis objects collapse metabolomics to one assay and keep
  # the platform in a column of its own, so "metab" is an assay value a table can
  # legitimately carry. It sorts where the platforms start.
  "metab",
  "metab-t-acoa", "metab-t-amines", "metab-t-conv", "metab-t-imm-crt",
  "metab-t-ka", "metab-t-nuc", "metab-t-oxylipneg", "metab-t-tca",
  "metab-u-hilicpos", "metab-u-ionpneg", "metab-u-lrpneg", "metab-u-lrppos",
  "metab-u-rpneg", "metab-u-rppos"
)

#' Every assay a supplementary table has a row for.
#'
#' OME_TISSUE_CODE is the study's tissue x ome roster. Two kinds of entry are
#' dropped:
#'
#'   lab-*            clinical chemistry, reported as clinical measures rather
#'                    than as an omics assay.
#'   *-clinical       metab-t-clinical and prot-clinical, which load_qc() skips
#'                    unless load_clinical = TRUE. A table built from load_qc()
#'                    output would have no rows for them, so a roster carrying
#'                    them would promise rows nothing can fill.
#'
#' The roster is what makes a zero distinguishable from a blank: an assay with no
#' outliers is a row of empty cells, not an absent row.
#'
#' It is a SUBSET of ASSAY_DISPLAY_ORDER, not equal to it: the display order also
#' carries "metab", which is an assay in the differential-analysis objects but not
#' a row any roster-driven table has.
supp_table_assays <- function() {
  roster <- unique(OME_TISSUE_CODE$ome)
  roster <- roster[!grepl("^lab-", roster)]
  roster <- roster[!grepl("-clinical$", roster)]
  roster <- sort(roster)
  unplaced <- setdiff(roster, ASSAY_DISPLAY_ORDER)
  if (length(unplaced) > 0) {
    stop("assay(s) in OME_TISSUE_CODE with no place in ASSAY_DISPLAY_ORDER: ",
         paste(unplaced, collapse = ", "),
         "\n  Add them in the position the manuscript introduces them.",
         call. = FALSE)
  }
  roster
}

#' Order a data frame by ASSAY_DISPLAY_ORDER, refusing an assay it does not know.
#'
#' An unknown assay stops the build. Sorting it to the end or dropping it would
#' both be silent, and a new assay in the study is a thing to be placed
#' deliberately rather than discovered in a supplementary table.
#'
#' @param df        A data frame with an assay column.
#' @param assay_col Name of that column.
#' @param ...       Further columns to break ties on, in order.
order_by_assay <- function(df, assay_col = "assay", ...) {
  unknown <- setdiff(unique(df[[assay_col]]), ASSAY_DISPLAY_ORDER)
  if (length(unknown) > 0) {
    stop("assay(s) with no place in ASSAY_DISPLAY_ORDER: ",
         paste(sort(unknown), collapse = ", "),
         "\n  Add them to lib/supp_table_helpers.R, in the position the ",
         "manuscript introduces them.", call. = FALSE)
  }
  rank <- match(df[[assay_col]], ASSAY_DISPLAY_ORDER)
  df[order(rank, ...), , drop = FALSE]
}
