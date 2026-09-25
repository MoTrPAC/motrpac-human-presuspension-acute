#!/usr/bin/env Rscript
# Supplementary Table 3 — the standalone sub-tables
#
# Sub-tables:  ST3g  ORA of muscle transcripts associated with baseline VO2peak
#              ST3h  ORA of adipose transcripts associated with baseline HOMA-IR
#
# The rest of ST3 is a panel's own numbers and is written by that panel's
# figure script, which config/table_map.json records:
#   ST3a  EE muscle-only ORA features        FIG3.R  (FIG3B)
#   ST3b  RE muscle-only ORA features        FIG3.R  (FIG3B)
#   ST3c  all-muscle RE/EE-only ORA          FIG3.R  (FIG3B)
#   ST3d  triangle selected pathways         ED4.R   (ED4B)
#   ST3f  clinical x omics all significant   ED4.R   (ED4G)
#
# ST3e, DA features with opposite EE/RE direction, is not built in this
# repository: table_map.json records it as built elsewhere, with no script.
#
# Both sub-tables here read the baseline clinical x omics fit, which is not in
# this repository: run
#
#   Rscript figures/landscape/analysis/03_clinical_omics.R
#
# once and it is cached.
#
#   Rscript figures/landscape/tables/ST3.R          every sub-table here
#   Rscript figures/landscape/tables/ST3.R ST3g     one of them

suppressPackageStartupMessages({
  library(dplyr)
  library(GO.db)
  library(AnnotationDbi)
  library(MotrpacHumanPreSuspensionAnalysis)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
# panel_export.R for landscape_root() and the config/landscape.env read it
# performs on load, which is where CLINICAL_OMICS_TABLE comes from;
# table_export.R after it.
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(here, "..", "lib", "table_export.R"))
source(file.path(here, "FIG3_ED4_helpers.R"))

# ---- shared ----------------------------------------------------------------

# ---- run_ORA output in the shape ST3g and ST3h report ----------------------

# The columns are the ones the assembled workbook's two clinical x omics ORA
# sheets carry: ONTOLOGY, ID, Description, FoldEnrichment, pvalue, p.adjust,
# geneID, Count, and a direction column named after the trait.
#
# GO sets fill ONTOLOGY with BP, CC or MF and ID with the GO accession. The
# other collections fill ONTOLOGY with their database name and leave ID empty;
# MSigDB set names are the only identifier they have.
#
# GO.db is attached here rather than by a figure script. No panel needs it, and
# config/panel_packages.txt covers the panel arm only.

GO_DATABASES <- c(GOBP = "BP", GOCC = "CC", GOMF = "MF")

#' GO term name -> GO accession.
#'
#' MSigDB names a GO set for its term with the punctuation flattened
#' (`GOBP_CELLULAR_RESPIRATION`); SET_TO_ID's set_id is MSigDB's own index, not
#' the accession. Flattening GO.db's term names the same way inverts it for
#' 10,341 of the 10,461 GO sets. The rest are terms MSigDB carries under a name
#' GO.db no longer uses, and they get an empty ID.
go_accession_by_term <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      terms <- AnnotationDbi::Term(GO.db::GOTERM)
      cache <<- stats::setNames(names(terms),
                                toupper(gsub("[^A-Za-z0-9]+", "_", terms)))
    }
    cache
  }
})

#' One clinical x omics ORA, in the reported shape.
#'
#' @param enrichment clinical_omics_ora() output.
#' @param direction_column Name for the direction column, e.g. "vo2peak_dir".
#' @returns A data.frame of the significant sets, strongest first.
format_clinical_omics_ora <- function(enrichment, direction_column) {
  significant <- enrichment[!is.na(enrichment$adj_p_value) &
                              enrichment$adj_p_value < 0.05, , drop = FALSE]
  if (nrow(significant) == 0) {
    stop("no set reaches FDR < 0.05 in this enrichment", call. = FALSE)
  }

  database <- as.character(significant$database)
  is_go <- database %in% names(GO_DATABASES)

  # The set name with its database prefix removed, underscores as spaces. For a
  # GO set this is the term as the workbook prints it.
  description <- tolower(gsub("_", " ", sub("^[A-Z0-9]+_", "",
                                            as.character(significant$set))))

  accession <- unname(go_accession_by_term()[
    toupper(sub("^GO(BP|CC|MF)_", "", as.character(significant$set)))
  ])
  accession[!is_go | is.na(accession)] <- ""

  # Observed overlap over the overlap a set of this size would take by chance.
  fold_enrichment <- (significant$set_size_in_input / significant$input_size) /
    (significant$set_size / significant$background_size)

  formatted <- data.frame(
    ONTOLOGY = ifelse(is_go, unname(GO_DATABASES[database]), database),
    ID = accession,
    Description = description,
    FoldEnrichment = fold_enrichment,
    pvalue = significant$p_value,
    p.adjust = significant$adj_p_value,
    geneID = significant$geneID,
    Count = significant$set_size_in_input,
    stringsAsFactors = FALSE
  )
  formatted[[direction_column]] <- significant$direction

  # Strongest first, then the set name so the order is total. adj_p_value is
  # adjusted within collection, so p_value is what orders across them.
  formatted[order(formatted$p.adjust, formatted$pvalue,
                  formatted$Description), , drop = FALSE]
}

#' The schema check both sub-tables make before exporting.
check_ora_columns <- function(formatted, spec) {
  missing <- setdiff(spec$columns, names(formatted))
  if (length(missing) > 0) {
    stop(spec$table, ": the enrichment carries no column(s) ",
         paste(missing, collapse = ", "),
         "\n  It carries: ", paste(names(formatted), collapse = ", "),
         "\n  Reconcile `columns` for ", spec$table,
         " in config/table_map.json with format_clinical_omics_ora().",
         call. = FALSE)
  }
  invisible(TRUE)
}

# ---- ST3g — ORA of muscle transcripts associated with baseline VO2peak -----

# Standalone. ST3f reports every significant association in the fit and ED4G
# counts them; no panel draws this enrichment.
#
# The muscle transcriptomics x VO2peak slice of the fit, split on the sign of
# the association and tested against every gene the slice fitted. ST3h is the
# same table for adipose x HOMA-IR.
#
# Columns follow the assembled workbook's clinical_x_omics_SKM_transcript sheet.
# The rows are a superset of it: that sheet is Gene Ontology only, this one
# carries every collection run_ORA tests.
st3g <- function() {
  skip <- skip_if_no_fit("ST3g")
  if (!is.null(skip)) return(skip)
  spec <- table_init("ST3g")

  enrichment <- clinical_omics_ora("VO2max_L", "muscle", "transcript-rna-seq")
  formatted <- format_clinical_omics_ora(enrichment, "vo2peak_dir")
  check_ora_columns(formatted, spec)

  export_table(as.data.frame(formatted)[, spec$columns, drop = FALSE], "ST3g")
}

# ---- ST3h — ORA of adipose transcripts associated with baseline HOMA-IR ----

# The adipose transcriptomics x HOMA-IR slice of the fit, split on the sign of
# the association and tested against every gene the slice fitted. ST3g is the
# same table for muscle x VO2peak.
#
# Columns follow the assembled workbook's clinical_x_omics_AT_transcripts sheet.
st3h <- function() {
  skip <- skip_if_no_fit("ST3h")
  if (!is.null(skip)) return(skip)
  spec <- table_init("ST3h")

  enrichment <- clinical_omics_ora("HOMA_IR", "adipose", "transcript-rna-seq")
  formatted <- format_clinical_omics_ora(enrichment, "homa_ir_dir")
  check_ora_columns(formatted, spec)

  export_table(as.data.frame(formatted)[, spec$columns, drop = FALSE], "ST3h")
}

run_tables(list(ST3g = st3g, ST3h = st3h))
