#!/usr/bin/env Rscript
# Extended Data 1 — Cohort and data overview
#
# Panels:  ED1A  genetic ancestry PCA, participants projected onto 1000 Genomes
#          ED1B  sample overlap across assays within muscle, blood and adipose
#          ED1C  sample counts for significant vs non-significant
#                (phospho)proteomic features
#          ED1D  covariate association with principal components, all
#                tissue x assay
# Tables:  ST1e   observations per named feature            (written by ED1C)
#
# ---------------------------------------------------------------------------
# ED1A cannot be run from this repository, and its input cannot be released.
#
# It reads a genotype PCA carrying 238 MoTrPAC participant IDs and 32 genotype
# principal components per individual. That is individual-level genetic data:
# it is not distributed with this repository and may not be published. Running
# ED1A requires an approved MoTrPAC consortium data-access request and a local
# copy of the file at figure_1/sources/, or ED1A_PCA_RDS pointing at one.
#
# Every other ED1 panel runs without it.
# ---------------------------------------------------------------------------
#
# Needs consortium data access. ED1B and ED1D read QC matrices through
# load_qc(epigen = TRUE); the ATAC and methylCap objects are not in the data
# package and are downloaded into EPIGEN_QC_DIR on first use, then reused.
#
#   Rscript figures/landscape/ED1.R          every panel
#   Rscript figures/landscape/ED1.R ED1C     one panel

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggh4x)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(tidyr)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(here, "..", "lib", "highlights.R"))
source(file.path(here, "..", "lib", "table_export.R"))
source(file.path(here, "..", "lib", "supp_table_helpers.R"))

# ---- shared ----------------------------------------------------------------

# The two data-dependently acquired omes. ED1C draws their count distribution
# and write_st1e() reports the counts feature by feature, so the constant is
# shared by the panel and the table rather than repeated in each.
PROT_OMES <- c("prot-pr", "prot-ph")

# The three load_qc() calls below are not shared: ED1B takes every ome with
# epigen, ED1C every ome without it, ED1D everything but methylCap and then
# rewrites the metabolomics entries. Each stays inside its panel.

# ---- ED1A — genetic ancestry PCA over 1000 Genomes -------------------------

# Plots gPC1 against gPC2 of a joint PCA over LD-pruned genotypes shared between
# the 1000 Genomes reference panel and MoTrPAC participants. Reference samples
# are coloured by 1000 Genomes population code; MoTrPAC samples are open
# triangles with marginal rugs. Legends sit under the plot.
#
# The genotype PCA it reads is individual-level data, is not in this repository
# and cannot be released; see the note at the top of this file.

ed1a <- function() {
  skip <- skip_if_gated(
    "ED1A",
    path = "figure_1/sources/1kg.motrpac.merged.maf_0_05.ld_pruned.pca.RDS",
    env_var = "ED1A_PCA_RDS",
    input = paste("a genotype PCA carrying 238 MoTrPAC participant IDs and 32",
                  "genotype principal components per individual"),
    access = "an approved MoTrPAC consortium data-access request"
  )
  if (!is.null(skip)) return(skip)

  spec <- panel_init("ED1A")

  # ---- inputs ----
  #
  # Nothing here computes a principal component. The three files below ship with
  # the repo, under this figure's sources/; each takes an environment-variable
  # override.
  #
  # The participants' coordinates come from the data hub by default — the released
  # resources/motrpac_human-precovid_1kg_pca.csv, downloaded and cached. That
  # resource carries no coordinate for any of the 2,504 reference samples, so the
  # 1000 Genomes cloud comes from the vendored snpgdsPCA object.
  # ED1A_PC_SOURCE=cached_pca takes the participants from it as well.

  pca_rds_path <- panel_source(spec$figure_dir,
                               "1kg.motrpac.merged.maf_0_05.ld_pruned.pca.RDS",
                               env_var = "ED1A_PCA_RDS")

  samples_path <- panel_source(spec$figure_dir, "1kg_samples.tsv",
                               env_var = "ED1A_1KG_SAMPLES")

  pop_colors_path <- panel_source(spec$figure_dir, "1kg_population_colors.csv",
                                  env_var = "ED1A_1KG_POP_COLORS")

  pc_source <- Sys.getenv("ED1A_PC_SOURCE", unset = "data_hub")
  if (!pc_source %in% c("data_hub", "cached_pca")) {
    stop("ED1A_PC_SOURCE must be data_hub or cached_pca, got '", pc_source, "'",
         call. = FALSE)
  }

  # The released participant PCs, pinned to a collection version: which
  # collection's coordinates the figure shows is not a thing to leave to whatever
  # is newest.
  ANCESTRY_PCA_GCS <- paste0(
    "gs://motrpac-data-hub/analysis/human-precovid-sed-adu/c1.3/resources/",
    "motrpac_human-precovid_1kg_pca.csv"
  )
  # Cached under staging/, NOT under this figure's sources/. The resource is one
  # row per participant vial label carrying 32 genotype principal components, so
  # a copy inside the tracked tree could be committed to a public repository.
  # staging/ is the download cache and is the only ignored location.
  ANCESTRY_PCA_CACHE <- file.path(landscape_root(), "staging", "ancestry_pca")

  # ---- data ----

  populations <- read.csv(samples_path, sep = "\t", header = TRUE,
                          quote = "\"", check.names = FALSE)
  population_colors <- read.csv(pop_colors_path, check.names = FALSE)

  pca <- readRDS(pca_rds_path)

  pca_colors <- stats::setNames(population_colors$Color, population_colors$Population)

  pca_coords <- data.frame(
    Sample.ID = pca$sample.id,
    PC1 = pca$eigenvect[, 1],
    PC2 = pca$eigenvect[, 2]
  )

  # The reference cloud, always from the cached object. merge() sorts on the join
  # key, which fixes the point draw order; the second merge also restricts the
  # panel to populations that have a colour.
  pca_1kg <- pca_coords %>%
    merge(populations, by.x = "Sample.ID", by.y = "Sample name") %>%
    dplyr::select("Sample.ID", "PC1", "PC2", Population = "Population code") %>%
    merge(population_colors, by = "Population") %>%
    dplyr::mutate(
      Consortium = "1000 Genomes",
      Population = factor(.data$Population, levels = population_colors$Population)
    )

  # ---- the participants ----------------------------------------------------

  #' The released ancestry PCs, downloaded from the data hub and cached.
  #'
  #' One row per participant, keyed on vialLabel, plus a leading
  #' `var_explained_pct` row that carries the percentage each component explains.
  #' That row is data about the columns rather than a sample, so it is split off
  #' here and supplies the axis labels.
  released_ancestry_pcs <- function() {
    dir.create(ANCESTRY_PCA_CACHE, recursive = TRUE, showWarnings = FALSE)

    # check_first = TRUE reuses the cached copy; the download happens once.
    released <- MotrpacBicQC::dl_read_gcp(ANCESTRY_PCA_GCS, sep = ",",
                                          tmpdir = ANCESTRY_PCA_CACHE,
                                          gsutil_path = Sys.getenv("GSUTIL", unset = "gsutil"),
                                          check_first = TRUE)
    released <- as.data.frame(released)

    variance_row <- released[released$vialLabel == "var_explained_pct", , drop = FALSE]
    if (nrow(variance_row) != 1) {
      stop("the released ancestry PCA has no var_explained_pct row: ",
           ANCESTRY_PCA_GCS, call. = FALSE)
    }
    samples <- released[released$vialLabel != "var_explained_pct", , drop = FALSE]

    list(
      coords = data.frame(
        Sample.ID = as.character(samples$vialLabel),
        PC1 = as.numeric(samples$PC1),
        PC2 = as.numeric(samples$PC2)
      ),
      var_explained = c(PC1 = as.numeric(variance_row$PC1),
                        PC2 = as.numeric(variance_row$PC2))
    )
  }

  if (identical(pc_source, "data_hub")) {
    released <- released_ancestry_pcs()
    pca_motrpac <- released$coords %>% dplyr::mutate(Consortium = "MoTrPAC")
    var_explained <- released$var_explained
    message(sprintf("[ED1A] %d participants from the released resource",
                    nrow(pca_motrpac)))
  } else {
    # Everything that is not a 1000 Genomes sample is a MoTrPAC one. This is the
    # 238-sample set the published panel drew, against the released cohort's 174.
    pca_motrpac <- pca_coords %>%
      dplyr::filter(!(.data$Sample.ID %in% populations[["Sample name"]])) %>%
      dplyr::mutate(Consortium = "MoTrPAC")
    var_explained <- c(PC1 = pca$varprop[1] * 100, PC2 = pca$varprop[2] * 100)
    message(sprintf("[ED1A] %d participants from the cached PCA object",
                    nrow(pca_motrpac)))
  }

  # ---- plot ----

  panel_theme <- theme_minimal(base_size = 18) +
    theme(
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.line.x.bottom = element_line(colour = "black", linewidth = 0.25),
      axis.line.y.left = element_line(colour = "black", linewidth = 0.25),
      legend.position = "bottom",
      legend.box = "vertical"
    )

  p <- ggplot() +
    geom_point(
      data = pca_1kg,
      aes(x = .data$PC1, y = .data$PC2, shape = .data$Consortium,
          color = .data$Population),
      size = 3, alpha = 0.5
    ) +
    geom_point(
      data = pca_motrpac,
      aes(x = .data$PC1, y = .data$PC2, shape = .data$Consortium),
      size = 3
    ) +
    geom_rug(
      data = pca_motrpac,
      aes(x = .data$PC1, y = .data$PC2),
      sides = "br", alpha = 0.25, length = grid::unit(2, "mm")
    ) +
    scale_shape_manual(values = c("1000 Genomes" = 20, "MoTrPAC" = 2)) +
    scale_color_manual(values = pca_colors) +
    # Percentages from the PCA the points came from.
    xlab(sprintf("gPC1 (%.2f%%)", var_explained[["PC1"]])) +
    ylab(sprintf("gPC2 (%.2f%%)", var_explained[["PC2"]])) +
    panel_theme

  export_panel(p, "ED1A")
}

# ---- ED1B — sample overlap across assays -----------------------------------

# One tile per participant x assay, coloured by ome where that participant has a
# measured sample and grey where they do not, faceted by exercise group and
# timepoint. Assay membership comes from the vial labels carried on each
# tissue x ome QC matrix; participant, group, sex and timepoint come from pheno.
# The three tissue plots are stacked into a single column sharing one legend.

ed1b <- function() {
  panel_init("ED1B")

  # ---- external inputs ----

  # The ATAC-seq and methylCap-seq QC objects are not shipped inside
  # MotrpacHumanPreSuspensionData; they are downloaded from the consortium bucket
  # into EPIGEN_QC_DIR, which config/landscape.env defaults to staging/raw-files/.
  # load_qc(epigen = TRUE) is what does the downloading, and reuses the cache on
  # every run after the first.
  epigen_cache <- epigen_qc_dir()

  # ---- data ----

  # Two orderings of the same seven omes: the first fixes factor level order
  # during reshaping, the second fixes the x-axis order of the tiles.
  ome_levels <- c(
    "Transcriptomics", "Metabolomics", "Proteomics (MS)", "Proteomics (Olink)",
    "Phosphoproteomics", "Chromatin Accessibility (ATAC)", "Methylation"
  )
  track_levels <- c(
    "Transcriptomics", "Metabolomics", "Proteomics (MS)", "Proteomics (Olink)",
    "Methylation", "Phosphoproteomics", "Chromatin Accessibility (ATAC)"
  )

  qc_object <- MotrpacHumanPreSuspensionData::load_qc(
    epigen = TRUE,
    gsutil = Sys.getenv("GSUTIL", unset = "gsutil"),
    repo_local_dir = epigen_cache
  )
  pheno_object <- MotrpacHumanPreSuspensionData::load_pheno()

  # Column names of each qc_norm matrix are the vial labels measured on that
  # tissue x ome. Blocks are stacked last-seen-first, which is what sets the
  # tie-break order of participants on the y-axis further down.
  vial_blocks <- list()
  for (tiss in names(qc_object)) {
    for (assay1 in names(qc_object[[tiss]])) {
      vial_blocks[[length(vial_blocks) + 1L]] <-
        data.frame(vialLabel = names(qc_object[[tiss]][[assay1]][["qc_norm"]])) %>%
        mutate(tissue = tiss, assay = assay1)
    }
  }
  all_vials <- bind_rows(rev(vial_blocks))

  pheno_vials_long_all_times <- pheno_object$pheno_data %>%
    select("pid", "vialLabel", "visitcode", "tempSampProfile",
           "randomGroupCode", "Sex", "Timepoint") %>%
    mutate(tempSampProfile = factor(
      .data$tempSampProfile,
      levels = c("All", "Early", "Middle", "Late")
    )) %>%
    right_join(all_vials, by = "vialLabel") %>%
    mutate(binary = 1) %>%
    filter(.data$visitcode == "ADU_BAS")

  pheno_vials_long <- pheno_vials_long_all_times %>%
    mutate(ome = case_when(
      str_detect(.data$assay, "metab") ~ "Metabolomics",
      str_detect(.data$assay, "prot-pr") ~ "Proteomics (MS)",
      str_detect(.data$assay, "prot-ol") ~ "Proteomics (Olink)",
      str_detect(.data$assay, "prot-ph") ~ "Phosphoproteomics",
      str_detect(.data$assay, "rna") ~ "Transcriptomics",
      str_detect(.data$assay, "methyl") ~ "Methylation",
      str_detect(.data$assay, "atac") ~ "Chromatin Accessibility (ATAC)",
      .default = .data$assay
    )) %>%
    mutate(ome = factor(.data$ome, levels = ome_levels)) %>%
    mutate(Group = factor(
      .data$randomGroupCode,
      levels = c("ADUEndur", "ADUResist", "ADUControl")
    )) %>%
    mutate(tissue = factor(
      .data$tissue,
      levels = c("muscle", "adipose", "blood")
    )) %>%
    mutate(Group = fct_recode(
      .data$Group,
      EE = "ADUEndur", RE = "ADUResist", CON = "ADUControl"
    )) %>%
    mutate(Sex = tolower(.data$Sex)) %>%
    mutate(Timepoint = fct_recode(
      .data$Timepoint,
      "Pre-Exercise" = "pre_exercise",
      "During 20 Min" = "during_20_min",
      "During 40 Min" = "during_40_min",
      "Post 10 Min" = "post_10_min",
      "Post 15/30/45 Min" = "post_15_30_45_min",
      "Post 3.5/4 Hour" = "post_3.5_4_hr",
      "Post 24 Hour" = "post_24_hr"
    )) %>%
    select(-"vialLabel", -"visitcode", -"assay") %>%
    distinct() %>%
    arrange(.data$Timepoint, .data$ome)

  # One row per participant x timepoint, one 0/1 column per ome.
  pheno_vials <- pheno_vials_long %>%
    pivot_wider(values_from = "binary", names_from = "ome",
                values_fill = 0, names_sep = ";") %>%
    arrange(.data$tissue, .data$randomGroupCode, .data$Sex)

  # ---- plot ----

  id_cols <- c("pid", "randomGroupCode", "Timepoint", "tissue", "Group",
               "tempSampProfile", "Sex")

  fill_values <- c(
    MotrpacHumanPreSuspensionAnalysis::HUMAN_OME_COLORS[track_levels],
    "No Data" = "lightgrey"
  )

  # ---- facet strip colours ----
  #
  # The tiles are faceted timepoint across and exercise group down, and those two
  # facets are the only place either variable is named — the axis text is blank.
  # Colouring the strips is what puts the canonical timepoint and group palettes
  # on this panel.
  #
  # The timepoint strips are relabelled per tissue further down — muscle's third
  # column reads "Post 3.5 Hour", adipose's "Post 4 Hour" — so the colour for a
  # column is looked up from the label that column ends up carrying. Every one of
  # those labels has its own key in HUMAN_ACUTE_TIMEPOINT_COLORS, and the keys
  # that name the same visit resolve to the same colour there (post_15_min,
  # post_30_min and post_45_min are one colour; post_3.5_4_hr and post_4_hr are
  # another), which is what keeps a timepoint the same colour across tissues.
  timepoint_label_to_key <- c(
    "Pre-Exercise"  = "pre_exercise",
    "During 20 Min" = "during_20_min",
    "During 40 Min" = "during_40_min",
    "Post 10 Min"   = "post_10_min",
    "Post 15 Min"   = "post_15_min",
    "Post 30 Min"   = "post_30_min",
    "Post 45 Min"   = "post_45_min",
    "Post 3.5 Hour" = "post_3.5_4_hr",
    "Post 4 Hour"   = "post_4_hr",
    "Post 24 Hour"  = "post_24_hr"
  )

  #' Strip fills for a set of timepoint labels, in the order the strips appear.
  #'
  #' @param labels The displayed strip labels, as drawn.
  timepoint_fills <- function(labels) {
    keys <- timepoint_label_to_key[labels]
    if (anyNA(keys)) {
      stop("no palette key for timepoint label(s): ",
           paste(labels[is.na(keys)], collapse = ", "),
           " — add them to timepoint_label_to_key", call. = FALSE)
    }
    strip_fills(MotrpacHumanPreSuspensionAnalysis::HUMAN_ACUTE_TIMEPOINT_COLORS,
                unname(keys), "timepoint(s)")
  }

  group_short_to_code <- c(EE = "ADUEndur", RE = "ADUResist", CON = "ADUControl")

  #' Palette colours for a set of keys, erroring on any the palette lacks.
  #'
  #' A missing key would otherwise leave a default grey strip, which reads as a
  #' deliberate colour rather than as a lookup that failed.
  strip_fills <- function(palette, keys, what) {
    hits <- palette[keys]
    if (anyNA(hits)) {
      stop(what, " not in the package palette: ",
           paste(keys[is.na(hits)], collapse = ", "), call. = FALSE)
    }
    unname(hits)
  }

  # Black or white strip text, whichever reads against the fill behind it.
  text_on <- function(fills) {
    rgb_mat <- grDevices::col2rgb(fills) / 255
    luminance <- 0.299 * rgb_mat[1, ] + 0.587 * rgb_mat[2, ] + 0.114 * rgb_mat[3, ]
    ifelse(luminance < 0.55, "white", "black")
  }

  tissue_plots <- list()

  for (tissue_select in c("Blood", "Muscle", "Adipose")) {
    tissue_select_lower <- tolower(tissue_select)

    df_long_building <- pheno_vials %>%
      filter(.data$tissue == tissue_select_lower) %>%
      mutate(across(all_of(c("Sex", track_levels, "pid")), as.character)) %>%
      arrange(.data$Transcriptomics) %>%
      pivot_longer(cols = -all_of(id_cols),
                   names_to = "Track", values_to = "Value") %>%
      mutate(Track = factor(.data$Track, levels = track_levels)) %>%
      mutate(Ome = factor(
        if_else(.data$Value == 1, as.character(.data$Track), "No Data"),
        levels = c(track_levels, "No Data")
      )) %>%
      arrange(desc(.data$tempSampProfile)) %>%
      mutate(pid = factor(.data$pid, levels = unique(.data$pid))) %>%
      mutate(Value = as.numeric(.data$Value))

    # Drop assays with no samples at all in this tissue rather than drawing a
    # column of grey.
    track_totals <- df_long_building %>%
      group_by(.data$tissue, .data$Track) %>%
      summarise(sum = sum(.data$Value), .groups = "drop")

    df_long <- df_long_building %>%
      left_join(track_totals, by = c("tissue", "Track")) %>%
      filter(.data$sum != 0)

    group_fills <- strip_fills(
      MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS,
      unname(group_short_to_code[levels(droplevels(df_long$Group))]),
      paste0(tissue_select, " exercise group(s)")
    )

    # Timepoint levels are relabelled positionally to the sampling schedule of
    # the tissue: muscle and adipose have four, blood has seven.
    if (tissue_select_lower == "muscle") {
      legend_pos <- "right"
      df_long <- df_long %>%
        mutate(Timepoint = factor(.data$Timepoint, labels = c(
          "Pre-Exercise", "Post 15 Min", "Post 3.5 Hour", "Post 24 Hour"
        )))
    } else if (tissue_select_lower == "adipose") {
      legend_pos <- "none"
      df_long <- df_long %>%
        mutate(Timepoint = factor(.data$Timepoint, labels = c(
          "Pre-Exercise", "Post 45 Min", "Post 4 Hour", "Post 24 Hour"
        )))
    } else {
      legend_pos <- "none"
      df_long <- df_long %>%
        mutate(Timepoint = factor(.data$Timepoint, labels = c(
          "Pre-Exercise", "During 20 Min", "During 40 Min", "Post 10 Min",
          "Post 30 Min", "Post 3.5 Hour", "Post 24 Hour"
        )))
    }

    # Column colours are keyed on the label each strip actually carries, after the
    # relabel above — not on the position of the level behind it. The two agree
    # today, but only because every tissue's schedule happens to be a prefix-like
    # subset of blood's; a tissue whose timepoints sat in a different order would
    # take blood's colours in blood's order and nothing would say so.
    tp_fills <- timepoint_fills(levels(droplevels(df_long$Timepoint)))

    tissue_plots[[tissue_select]] <-
      ggplot(df_long, aes(x = .data$Track, y = .data$pid, fill = .data$Ome)) +
      geom_tile(color = "#00000000", linewidth = 0.05) +
      scale_fill_manual(values = fill_values) +
      facet_grid2(
        Group ~ Timepoint, scales = "free_y", shrink = TRUE, space = "free_y",
        strip = strip_themed(
          background_x = lapply(tp_fills, function(f) {
            element_rect(fill = f, colour = "grey20")
          }),
          text_x = lapply(text_on(tp_fills), function(cl) {
            element_text(colour = cl, face = "bold")
          }),
          background_y = lapply(group_fills, function(f) {
            element_rect(fill = f, colour = "grey20")
          }),
          text_y = lapply(text_on(group_fills), function(cl) {
            element_text(colour = cl, face = "bold")
          })
        )
      ) +
      theme_classic() +
      ggtitle(paste(tissue_select, "Sample Overlap")) +
      ylab("Participant") +
      theme(
        legend.background = element_blank(),
        legend.position = legend_pos,
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank()
      )
  }

  p <- tissue_plots[["Muscle"]] +
    tissue_plots[["Blood"]] +
    tissue_plots[["Adipose"]] +
    plot_layout(guides = "collect", ncol = 1)

  export_panel(p, "ED1B")
}

# ---- ED1C — sample counts for (phospho)proteomic features ------------------

# Proteomics and phosphoproteomics are acquired data-dependently, so every
# feature has its own missingness pattern and its own effective n. For each
# feature x group x timepoint this counts the acute-visit samples with an
# observed value, splits those observations by whether the feature was
# differentially abundant anywhere in the exercise contrasts, and stacks the two
# distributions as a histogram per tissue x assay.
#
# This panel also writes ST1e, the per-feature observation counts for the
# features named in the manuscript. The panel is the distribution of those
# counts for prot-pr and prot-ph; the table reports them feature by feature and
# covers the other assays too. Both come out of the one load_qc() call below,
# and the proteomics half of the table IS all_group_stats — the object the panel
# is drawn from.

ed1c <- function() {
  panel_init("ED1C")

  # ---- data ----

  # Every non-epigenomics ome, not just the two the panel draws. ST1e reports
  # transcriptomics, Olink and the metabolomics platforms as well, and one load
  # answering both is the point of building them together. The panel's own
  # computation below is unchanged: it iterates PROT_OMES only.
  qc_data <- MotrpacHumanPreSuspensionData::load_qc(
    selected_omes = "all",
    epigen = FALSE
  )

  # Observed sample count per feature x group x timepoint, for every tissue that
  # has MS proteomics or phosphoproteomics.
  group_stats_list <- list()
  for (tissue_name in names(qc_data)) {
    for (ome_name in intersect(names(qc_data[[tissue_name]]), PROT_OMES)) {
      curr_data <- qc_data[[tissue_name]][[ome_name]][["qc_norm"]]
      if (nrow(curr_data) == 0) next

      # The acute bout is the ADU_BAS visit; training-visit samples do not
      # contribute to any of the contrasts plotted here.
      curr_meta <- qc_data[[tissue_name]][[ome_name]][["sample_metadata"]] %>%
        dplyr::filter(visitcode == "ADU_BAS")
      curr_data <- curr_data[, colnames(curr_data) %in% curr_meta$vialLabel,
                             drop = FALSE]
      curr_meta <- curr_meta[match(colnames(curr_data), curr_meta$vialLabel), ]

      vial_to_group <- setNames(curr_meta$randomGroupCode, curr_meta$vialLabel)
      vial_to_timepoint <- setNames(curr_meta$Timepoint, curr_meta$vialLabel)

      curr_data_long <- as.data.frame(t(curr_data))
      curr_data_long$Sample <- rownames(curr_data_long)
      curr_data_long <- tidyr::pivot_longer(
        curr_data_long, -Sample,
        names_to = "feature_id", values_to = "Value"
      )
      curr_data_long$randomGroupCode <- vial_to_group[curr_data_long$Sample]
      curr_data_long$Timepoint <- vial_to_timepoint[curr_data_long$Sample]

      group_stats <- curr_data_long %>%
        dplyr::group_by(randomGroupCode, feature_id, Timepoint) %>%
        dplyr::summarize(Count = sum(!is.na(Value)), .groups = "drop") %>%
        dplyr::mutate(tissue = tissue_name, ome = ome_name)

      group_stats_list[[paste(tissue_name, ome_name, sep = ".")]] <- group_stats
    }
  }
  all_group_stats <- dplyr::bind_rows(group_stats_list)

  # A feature x group x timepoint counts as "significant" when it appears as a
  # term of any delta-delta contrast the feature is significant in. Both the
  # exercise-vs-control and the endurance-vs-resistance families are included, and
  # every group.timepoint term of a significant contrast is marked, because the
  # baseline and control observations contribute to that estimate too.
  significant_prot_pr_ph_features <- MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(
    selected_omes = c("prot-pr", "prot-ph"),
    selected_tissues = "all",
    single_matrix = TRUE
  ) %>%
    dplyr::filter(contrast_type %in% c("exercise_with_controls", "Endur_vs_Resist"),
                  adj_p_value < 0.05) %>%
    dplyr::select(feature_id, tissue, assay, contrast) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      terms = stringr::str_extract_all(
        gsub("group_timepoint", "", contrast),
        "ADU[A-Za-z]+\\.[\\w\\.]+"
      )
    ) %>%
    tidyr::unnest(cols = "terms") %>%
    dplyr::mutate(
      randomGroupCode = stringr::str_split_fixed(terms, "\\.", 2)[, 1],
      Timepoint = stringr::str_split_fixed(terms, "\\.", 2)[, 2]
    ) %>%
    dplyr::rename(ome = assay) %>%
    dplyr::select(feature_id, tissue, ome, randomGroupCode, Timepoint) %>%
    dplyr::distinct()

  join_keys <- c("feature_id", "tissue", "ome", "randomGroupCode", "Timepoint")
  sig_prot_pr_ph_stats <- all_group_stats %>%
    dplyr::semi_join(significant_prot_pr_ph_features, by = join_keys)
  non_sig_prot_pr_ph_stats <- all_group_stats %>%
    dplyr::anti_join(significant_prot_pr_ph_features, by = join_keys)

  # Pre-aggregate rather than letting stat_bin do it. stat_bin emits a rectangle
  # for every bin x fill combination, including empty ones; those become
  # zero-height stroked rects in the PDF that Illustrator shows as stray lines.
  plot_data <- dplyr::bind_rows(
    dplyr::mutate(sig_prot_pr_ph_stats, Category = "Significant"),
    dplyr::mutate(non_sig_prot_pr_ph_stats, Category = "Non-significant")
  ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Category = factor(Category, levels = c("Significant", "Non-significant"))
    ) %>%
    dplyr::count(tissue, ome, Category, Count, name = "n_obs") %>%
    dplyr::filter(n_obs > 0)

  x_max <- max(plot_data$Count, na.rm = TRUE)

  write_st1e(qc_data, all_group_stats)

  # ---- plot ----

  strip_labels <- c(
    "prot-ph" = "Phosphoproteomics",
    "prot-pr" = "Proteomics",
    "adipose" = "Adipose",
    "muscle" = "Muscle",
    "blood" = "Blood"
  )
  relabel_strip <- function(x) {
    x <- as.character(x)
    out <- unname(strip_labels[x])
    ifelse(is.na(out), x, out)
  }

  p <- ggplot(plot_data, aes(x = Count, y = n_obs, fill = Category)) +
    # geom_col with no border: a stroked bar keeps its stroke width when the
    # figure is scaled down in Illustrator, so thin bars collapse into slivers.
    geom_col(width = 1, position = "stack") +
    scale_fill_manual(
      values = c("Significant" = "#E64B35", "Non-significant" = "#ADB5BD"),
      name = NULL
    ) +
    geom_vline(xintercept = 3, linetype = "dashed", color = "grey30",
               linewidth = 0.6) +
    scale_x_continuous(
      # bars are centred on integer counts, so pad by half a bin on each side
      limits = c(-0.5, x_max + 0.5),
      breaks = scales::breaks_width(5),
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(expand = c(0, 0, 0.05, 0)) +
    labs(
      title = "Observed sample counts for significant vs. non-significant (phospho)proteomic features",
      subtitle = "Each panel is one tissue × assay; each observation is one feature × group × timepoint",
      x = "Observed samples per group × timepoint for each feature",
      y = "Number of feature × group × timepoint observations"
    ) +
    facet_wrap(vars(tissue, ome), scales = "free_y",
               labeller = labeller(.default = relabel_strip)) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey40"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 13),
      strip.text = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 13),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
      plot.margin = margin(10, 16, 10, 10),
      legend.position = "bottom"
    )

  export_panel(p, "ED1C")
}

# ---- ED1D — covariate association with principal components ----------------

# For every tissue x assay matrix, the first five principal components are
# correlated against participant and design covariates with canonical
# correlation analysis. Each cell of a sub-heatmap is one covariate x PC. The
# metabolomics platforms are collapsed into a single matrix first, so
# metabolomics contributes one sub-heatmap per tissue rather than fourteen.
# The sub-heatmaps are stacked into one column per tissue.

ed1d <- function() {
  panel_init("ED1D")

  # ---- data ----

  # ATAC-seq is not shipped in MotrpacHumanPreSuspensionData; it is downloaded from
  # the consortium bucket into EPIGEN_QC_DIR, which config/landscape.env defaults to
  # staging/raw-files/. load_qc(epigen = TRUE) is what does the downloading.
  # methylCap is excluded here, so only the two ATAC matrices are fetched.
  epigen_local_dir <- epigen_qc_dir()

  non_methyl <- setdiff(
    MotrpacHumanPreSuspensionAnalysis::ome_available_list(),
    "epigen-methylcap-seq"
  )

  qc_norm_list <- MotrpacHumanPreSuspensionData::load_qc(
    selected_omes = non_methyl,
    epigen = TRUE,
    gsutil = Sys.getenv("GSUTIL", unset = "gsutil"),
    repo_local_dir = epigen_local_dir
  )

  # Replace the per-platform metabolomics entries with one combined matrix per
  # tissue. combine_qc_matrixes keys columns on pid..Timepoint, so the column
  # names are mapped back to vial labels to stay joinable with the sample
  # metadata.
  for (tissue in names(qc_norm_list)) {
    metab_platforms_only <- MotrpacHumanPreSuspensionData::load_qc(
      selected_tissues = tissue,
      selected_omes = MotrpacHumanPreSuspensionAnalysis::metab_only_list()
    )

    combined_metab_matrixes <- MotrpacHumanPreSuspensionData::combine_qc_matrixes(
      metab_platforms_only
    )
    # The combined matrix keeps only participants shared across platforms, so any
    # one platform's metadata covers every remaining column.
    metab_single_metadata <- metab_platforms_only[[tissue]][[1]][["sample_metadata"]]

    name_mapping <- metab_single_metadata %>%
      dplyr::mutate(pid_timepoint = interaction(pid, Timepoint, sep = "..")) %>%
      dplyr::filter(pid_timepoint %in% colnames(combined_metab_matrixes)) %>%
      dplyr::select(pid_timepoint, vialLabel) %>%
      tibble::deframe()
    colnames(combined_metab_matrixes) <-
      name_mapping[colnames(combined_metab_matrixes)]

    combined_metab_replacement <- list(
      "qc_norm" = combined_metab_matrixes,
      "sample_metadata" = metab_single_metadata
    )

    qc_norm_list[[tissue]][grep("metab", names(qc_norm_list[[tissue]]))] <- NULL
    qc_norm_list[[tissue]][["metab"]] <- combined_metab_replacement
  }

  # ---- plot ----

  # pheatmap draws as a side effect of being built, so the CCA loop runs against a
  # null device; only the returned gtables are kept.
  grDevices::pdf(NULL)

  cca_plot_list <- list()
  for (tissue in names(qc_norm_list)) {
    for (assay in names(qc_norm_list[[tissue]])) {
      data_to_pca <- qc_norm_list[[tissue]][[assay]][["qc_norm"]]
      if (assay == "prot-ph" || assay == "prot-pr") {
        data_to_pca <- data_to_pca %>% tidyr::drop_na()
      }
      if (nrow(data_to_pca) < 10) next

      sample_metadata <- qc_norm_list[[tissue]][[assay]][["sample_metadata"]]
      rownames(sample_metadata) <- sample_metadata$vialLabel

      cca_plot_list[[tissue]][[assay]] <-
        MotrpacHumanPreSuspensionAnalysis::plot_precovid_cca(
          data_to_pca,
          sample_metadata,
          custom_title = paste(tissue, assay)
        )
    }
  }

  grDevices::dev.off()

  grobs_by_tissue <- list()
  for (tissue in names(cca_plot_list)) {
    sublist <- list()
    for (assay in names(cca_plot_list[[tissue]])) {
      pheat <- cca_plot_list[[tissue]][[assay]]
      if (inherits(pheat, "pheatmap")) {
        sublist[[assay]] <- pheat$gtable
      }
    }
    if (length(sublist) > 0) {
      grobs_by_tissue[[tissue]] <- sublist
    }
  }

  # One column per tissue, that tissue's assays stacked vertically.
  tissue_panels <- lapply(grobs_by_tissue, function(grob_list) {
    gridExtra::arrangeGrob(grobs = grob_list, ncol = 1)
  })

  p <- gridExtra::arrangeGrob(grobs = tissue_panels, ncol = length(tissue_panels))

  export_panel(p, "ED1D")
}

# ---- ST1e — observations per named feature ---------------------------------

write_st1e <- function(qc_data, all_group_stats) {
  # The named features, from config/highlights.json.
  named_features <- highlight_named_features()

  # Gene is the symbol the manuscript refers to the feature by, and it is carried
  # in the free-text notes: "CTNND1 pS920" for a phosphosite, plain "CX3CL1" for a
  # protein. The symbol is the first token; the rest is the site.
  named_features$Gene <- sub("\\s.*$", "", trimws(named_features$notes))
  named_features$notes <- NULL

  # Three ways to count, because the assays differ in whether a feature can be
  # missing from a sample that was run:
  #
  #   no missingness   transcriptomics, Olink, metabolomics — every feature is
  #                    observed in every sample of that tissue x assay, so the
  #                    count is a property of the assay, not of the feature.
  #   missingness      prot-pr and prot-ph are acquired data-dependently, so each
  #                    feature has its own count. That is all_group_stats above.
  #
  # An assay in the list with no branch here would be dropped silently, so the
  # split is asserted rather than assumed.
  NO_MISSINGNESS_OMES <- c("transcript-rna-seq", "prot-ol")
  handled <- c(NO_MISSINGNESS_OMES, PROT_OMES, "metab")
  unhandled <- setdiff(unique(named_features$assay), handled)
  if (length(unhandled) > 0) {
    message("[ST1e] dropping ", sum(named_features$assay %in% unhandled),
            " row(s) for assay(s) with no counting rule: ",
            paste(sort(unhandled), collapse = ", "))
    named_features <- named_features[!named_features$assay %in% unhandled, ]
  }

  # Samples per tissue x assay x group x timepoint, for the assays where that is
  # the answer for every feature.
  sample_count_blocks <- list()
  for (tissue_name in names(qc_data)) {
    for (ome_name in names(qc_data[[tissue_name]])) {
      meta <- qc_data[[tissue_name]][[ome_name]][["sample_metadata"]]
      if (is.null(meta) || nrow(meta) == 0) next
      vials <- colnames(qc_data[[tissue_name]][[ome_name]][["qc_norm"]])
      sample_count_blocks[[paste(tissue_name, ome_name)]] <- meta %>%
        dplyr::filter(.data$visitcode == "ADU_BAS",
                      .data$vialLabel %in% vials) %>%
        dplyr::count(.data$randomGroupCode, .data$Timepoint, name = "n") %>%
        dplyr::mutate(tissue = tissue_name, assay = ome_name)
    }
  }
  assay_sample_counts <- dplyr::bind_rows(sample_count_blocks)

  # many-to-many by design: several named features share a tissue x assay, and
  # each fans out to that assay's group x timepoint counts. Declared so the fan-out
  # reads as intended rather than as dplyr's warning about an accidental one.
  st1e_flat <- named_features %>% dplyr::filter(.data$assay %in% NO_MISSINGNESS_OMES) %>%
    dplyr::left_join(assay_sample_counts, by = c("tissue", "assay"),
                     relationship = "many-to-many")

  # The list says "metab"; which platform a metabolite was quantified on is
  # recorded in the differential analysis, and the count follows the platform.
  metab_platforms <- MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(
    selected_omes = "metab", single_matrix = TRUE
  ) %>%
    dplyr::select("feature_id", "tissue", "assay", "platform") %>%
    dplyr::distinct()

  st1e_metab <- named_features %>% dplyr::filter(.data$assay == "metab") %>%
    dplyr::left_join(metab_platforms, by = c("feature_id", "tissue", "assay")) %>%
    dplyr::mutate(assay = .data$platform) %>%
    dplyr::select(-"platform") %>%
    dplyr::left_join(assay_sample_counts, by = c("tissue", "assay"),
                     relationship = "many-to-many")

  st1e_prot <- named_features %>% dplyr::filter(.data$assay %in% PROT_OMES) %>%
    dplyr::left_join(all_group_stats,
                     by = c("feature_id", "tissue", "assay" = "ome")) %>%
    dplyr::rename(n = "Count")

  st1e_long <- dplyr::bind_rows(st1e_flat, st1e_metab, st1e_prot) %>%
    dplyr::mutate(randomGroupCode = sub("^ADU", "", .data$randomGroupCode))

  st1e_declared <- table_spec("ST1e")$columns
  st1e_value_cols <- setdiff(st1e_declared,
                             c("feature_id", "tissue", "assay", "Gene"))
  st1e_long <- st1e_long %>%
    dplyr::mutate(column = paste(.data$randomGroupCode, .data$Timepoint, sep = "_"))

  unexpected <- setdiff(unique(stats::na.omit(st1e_long$column)), st1e_value_cols)
  if (length(unexpected) > 0) {
    stop("ST1e: group x timepoint combinations the manifest does not declare: ",
         paste(sort(unexpected), collapse = ", "),
         "\n  Add them to `columns` for ST1e in config/table_map.json.",
         call. = FALSE)
  }

  st1e <- st1e_long %>%
    dplyr::select("feature_id", "tissue", "assay", "Gene", "column", "n") %>%
    tidyr::pivot_wider(names_from = "column", values_from = "n") %>%
    as.data.frame()
  st1e[["NA"]] <- NULL
  for (missing_col in setdiff(st1e_value_cols, names(st1e))) {
    st1e[[missing_col]] <- NA_integer_
  }

  # Muscle first, then blood, then adipose - the order the manuscript walks the
  # tissues. Assay within tissue in the shared modality order, and feature within
  # that alphabetically, so the row order is fully determined.
  st1e <- st1e[order(match(st1e$tissue, c("muscle", "blood", "adipose")),
                     match(st1e$assay, ASSAY_DISPLAY_ORDER),
                     st1e$feature_id), st1e_declared]

  export_table(st1e, "ST1e")
}

run_panels(list(ED1A = ed1a, ED1B = ed1b, ED1C = ed1c, ED1D = ed1d))
