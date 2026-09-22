#!/usr/bin/env Rscript
# Figure 1 — Study design and cohort overview
#
# Panels:  FIG1C  samples measured per ome, tissue, exercise group and timepoint
# Tables:  ST1c   participants per tissue and omics assay   (written by FIG1C)
#          ST2b   total features detected                   (written by FIG1C)
#
# Not built here, see FIG1's not_built in config/panel_map.json:
#   FIG1A  study design        BioRender illustration
#   FIG1B  sampling timeline   BioRender illustration
#
# Needs consortium data access. load_qc(epigen = TRUE) reads about 6 GB of QC
# matrices; the ATAC and methylCap objects are not in the data package and are
# downloaded into EPIGEN_QC_DIR on first use, then reused.
#
#   Rscript figures/landscape/FIG1.R          every panel
#   Rscript figures/landscape/FIG1.R FIG1C    one panel

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
  library(MotrpacHumanPreSuspensionData)
  library(MotrpacHumanPreSuspensionAnalysis)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(here, "..", "lib", "table_export.R"))
source(file.path(here, "..", "lib", "supp_table_helpers.R"))

# ---- shared ----------------------------------------------------------------

# Six omes, not seven: this panel folds the two proteomics platforms together,
# where ED1B keeps Proteomics (MS) and Proteomics (Olink) apart.
OME_LEVELS <- c(
  "Transcriptomics", "Metabolomics", "Proteomics", "Phosphoproteomics",
  "Chromatin Accessibility (ATAC)", "Methylation"
)

TIMEPOINT_LABELS <- c(
  pre_exercise      = "Pre-\nExercise",
  during_20_min     = "During\n20 Min",
  during_40_min     = "During\n40 Min",
  post_10_min       = "Post\n10 Min",
  post_15_30_45_min = "Post\n15/30/45 Min",
  post_3.5_4_hr     = "Post\n3.5/4 Hour",
  post_24_hr        = "Post\n24 Hour"
)

GROUP_LABELS <- c(ADUResist = "RE", ADUEndur = "EE", ADUControl = "CON")

# One load_qc(epigen = TRUE) pass behind the panel and both tables. Lazy and
# memoised rather than loaded at the top of the script: 6 GB should not be read
# to build a figure whose selected panels do not need it.
qc_data <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- MotrpacHumanPreSuspensionData::load_qc(
        epigen = TRUE,
        repo_local_dir = epigen_qc_dir()
      )
    }
    cache
  }
})

# One row per measured vial, with the participant, group and timepoint it
# belongs to. The panel rolls it up by ome, ST1c by assay.
acute_vials <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    qc_object <- qc_data()

    # Column names of each qc_norm matrix are the vial labels measured on that
    # tissue x ome.
    vial_blocks <- list()
    for (tiss in names(qc_object)) {
      for (assay1 in names(qc_object[[tiss]])) {
        vial_blocks[[length(vial_blocks) + 1L]] <-
          data.frame(vialLabel = names(qc_object[[tiss]][[assay1]][["qc_norm"]])) %>%
          mutate(tissue = tiss, assay = assay1)
      }
    }
    all_vials <- bind_rows(vial_blocks)

    # ADU_BAS is the acute visit. Restricting here rather than after counting
    # keeps a participant's post-training vials out of the acute bars, and out
    # of ST1c.
    cache <<- MotrpacHumanPreSuspensionData::load_pheno()$pheno_data %>%
      select("pid", "vialLabel", "visitcode", "randomGroupCode", "Timepoint") %>%
      right_join(all_vials, by = "vialLabel") %>%
      filter(.data$visitcode == "ADU_BAS")
    cache
  }
})

# ---- FIG1C — samples measured per ome, tissue, group and timepoint ---------

fig1c <- function() {
  panel_init("FIG1C")
  vials <- acute_vials()

  sample_counts <- vials %>%
    mutate(ome = case_when(
      str_detect(.data$assay, "metab") ~ "Metabolomics",
      str_detect(.data$assay, "prot-pr") ~ "Proteomics",
      str_detect(.data$assay, "prot-ol") ~ "Proteomics",
      str_detect(.data$assay, "prot-ph") ~ "Phosphoproteomics",
      str_detect(.data$assay, "rna") ~ "Transcriptomics",
      str_detect(.data$assay, "methyl") ~ "Methylation",
      str_detect(.data$assay, "atac") ~ "Chromatin Accessibility (ATAC)",
      .default = .data$assay
    )) %>%
    select("pid", "tissue", "Timepoint", "randomGroupCode", "ome") %>%
    # One row per participant x tissue x ome x timepoint. Without this a
    # participant measured on several metabolomics platforms would be counted
    # once per platform rather than once for Metabolomics.
    distinct() %>%
    mutate(
      ome = factor(.data$ome, levels = OME_LEVELS),
      tissue = factor(.data$tissue, levels = c("blood", "adipose", "muscle")),
      Group = factor(GROUP_LABELS[as.character(.data$randomGroupCode)],
                     levels = unname(GROUP_LABELS)),
      Timepoint = factor(TIMEPOINT_LABELS[as.character(.data$Timepoint)],
                         levels = unname(TIMEPOINT_LABELS))
    )

  write_st1c(vials)
  write_st2b(qc_data())

  p <- ggplot(sample_counts, aes(x = .data$tissue, fill = .data$ome)) +
    geom_bar(position = position_dodge(preserve = "single")) +
    facet_grid(Group ~ Timepoint, switch = "y") +
    scale_fill_manual(values = HUMAN_OME_COLORS[OME_LEVELS]) +
    labs(x = NULL, y = "Number of Samples", fill = "Ome") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "white"),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.key.size = unit(4, "mm"),
      legend.spacing.x = unit(2, "mm"),
      legend.spacing.y = unit(2, "mm"),
      legend.box.spacing = unit(2, "pt"),
      legend.box.margin = margin(2, 2, 2, 2, "pt")
    )

  export_panel(p, "FIG1C")
}

# ---- ST1c — participants per tissue and omics assay ------------------------

# Counted as distinct participants, not vials: the legend for this table reads
# "number of participants for which samples were analyzed", and a participant
# with two vials on one tissue x assay x timepoint is still one participant.
write_st1c <- function(vials) {
  st1c_long <- vials %>%
    distinct(.data$pid, .data$tissue, .data$assay,
             .data$randomGroupCode, .data$Timepoint) %>%
    count(.data$assay, .data$randomGroupCode, .data$tissue, .data$Timepoint,
          name = "n") %>%
    mutate(column = paste0(str_to_title(.data$tissue), "_", .data$Timepoint))

  # The manifest declares which tissue x timepoint columns this table has, and
  # they are not the full cross product: adipose and muscle are not sampled
  # during exercise, so those columns do not exist. A combination that appears
  # in the data and not in the schema is a real change in what was measured,
  # not a formatting detail, so it stops the build.
  st1c_declared <- table_spec("ST1c")$columns
  st1c_value_cols <- setdiff(st1c_declared, c("assay", "randomGroupCode"))
  st1c_unexpected <- setdiff(unique(st1c_long$column), st1c_value_cols)
  if (length(st1c_unexpected) > 0) {
    stop("ST1c: measured tissue x timepoint combinations the manifest does not ",
         "declare: ", paste(sort(st1c_unexpected), collapse = ", "),
         "\n  Add them to `columns` for ST1c in config/table_map.json if the ",
         "sampling really changed.", call. = FALSE)
  }

  st1c <- st1c_long %>%
    select("assay", "randomGroupCode", "column", "n") %>%
    pivot_wider(names_from = "column", values_from = "n")

  # Modality order, not alphabetical, and shared with every other per-assay
  # sub-table — see lib/supp_table_helpers.R.
  st1c <- order_by_assay(st1c, "assay", st1c[["randomGroupCode"]])

  # A declared column with no measurement anywhere is absent after pivoting
  # rather than present and empty. Adding it back keeps the header the schema
  # promises.
  for (missing_col in setdiff(st1c_value_cols, names(st1c))) {
    st1c[[missing_col]] <- NA_integer_
  }
  st1c <- st1c[, st1c_declared]

  export_table(st1c, "ST1c")
}

# ---- ST2b — total features per tissue and assay ----------------------------

# The row count of each qc_norm matrix, where ST1c counted the columns. The
# missingness filter the table's legend describes is already applied in the
# stored objects, so this is the filtered count and nothing is filtered again
# here. remove_redundant_metab is on by default in load_qc(), which is what puts
# a metabolite measured on several platforms into the one where its CV is lowest.
write_st2b <- function(qc_object) {
  feature_count_blocks <- list()
  for (tiss in names(qc_object)) {
    for (assay1 in names(qc_object[[tiss]])) {
      feature_count_blocks[[paste(tiss, assay1)]] <- data.frame(
        tissue = tiss, assay = assay1,
        n = nrow(qc_object[[tiss]][[assay1]][["qc_norm"]])
      )
    }
  }
  st2b_long <- bind_rows(feature_count_blocks) %>%
    filter(.data$n > 0) %>%
    mutate(column = str_to_title(.data$tissue))

  st2b_declared <- table_spec("ST2b")$columns
  st2b_tissue_cols <- setdiff(st2b_declared, "assay")
  st2b_unexpected <- setdiff(unique(st2b_long$column), st2b_tissue_cols)
  if (length(st2b_unexpected) > 0) {
    stop("ST2b: features measured in tissue(s) the manifest does not declare: ",
         paste(sort(st2b_unexpected), collapse = ", "),
         "\n  Add them to `columns` for ST2b in config/table_map.json.",
         call. = FALSE)
  }

  # Every assay in the study gets a row, as in ST1d: a blank cell means the
  # assay was not run on that tissue, which is a result rather than an omission.
  st2b <- data.frame(assay = supp_table_assays(), stringsAsFactors = FALSE) %>%
    left_join(select(st2b_long, "assay", "column", "n"), by = "assay") %>%
    filter(!is.na(.data$column) | !.data$assay %in% st2b_long$assay) %>%
    tidyr::pivot_wider(names_from = "column", values_from = "n") %>%
    as.data.frame()
  st2b[["NA"]] <- NULL
  for (missing_col in setdiff(st2b_tissue_cols, names(st2b))) {
    st2b[[missing_col]] <- NA_integer_
  }
  st2b <- order_by_assay(st2b, "assay")[, st2b_declared]

  export_table(st2b, "ST2b")
}

run_panels(list(FIG1C = fig1c))
