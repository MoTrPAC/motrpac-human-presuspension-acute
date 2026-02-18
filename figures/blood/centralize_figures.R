# Purpose: copy and rename individual figure panel files to a central location
# for easy access and modification for the blood manuscript.

# File paths should be relative to precovid-analyses/figures/blood/
files <- c(
  # Figure 2
  "2C" = "",
  "2D" = "feature_heatmaps/blood_transcript-rna-seq_exercise_with_controls_REACTOME_HSP90_CHAPERONE_CYCLE_FOR_STEROID_HORMONE.pdf",
  "2E" = "feature_heatmaps/muscle_transcript-rna-seq_exercise_with_controls_REACTOME_HSP90_CHAPERONE_CYCLE_FOR_STEROID_HORMONE.pdf",

  # Supplementary Figure 2 (SF2)
  "SF2A" = "",
  "SF2B" = "feature_heatmaps/blood_transcript-rna-seq_EE_vs_RE_custom_Figure_S2B.pdf",

  # Figure 3
  "3C" = "",
  "3E" = "feature_heatmaps/blood_transcript-rna-seq_EE_vs_RE_custom_Figure_3E.pdf",

  # Supplementary Figure 3 (SF3)
  "S3A" = "",
  "S3B" = "feature_heatmaps/blood_transcript-rna-seq_EE_vs_RE_naive_bcell_blood.pdf",
  "S3C" = "feature_heatmaps/blood_transcript-rna-seq_EE_vs_RE_B_cell_lymph_node.pdf",
  "S3D" = "feature_heatmaps/blood_transcript-rna-seq_EE_vs_RE_activated_memory_B_cell_blood.pdf",
  "S3E" = "feature_heatmaps/blood_transcript-rna-seq_exercise_with_controls_B_cell_lymph_node.pdf",

  # Figure 6
  "6C" = "",
  "6D_upper" = "",
  "6D_lower" = "",

  # Supplementary Figure 6 (SF6)
  "SF6A" = "",

  # Figure 7
  "7A" = "",

  # Supplementary Figure 7 (SF7)
  "SF7A" = "",
  "SF7B" = "",
  "SF7C" = "" # missing
)

# Remove missing files
files <- files[files != ""]

names(files) <- paste0(
  "figure_",
  names(files),
  sub(".*(\\.[^.]+$)", "\\1", files)
)

# Central location to save the files (relative to precovid-analyses/)
file_dir <- "figures/blood/final_figure_files"

if (!dir.exists(file_dir))
  dir.create(file_dir)

# Copy-rename figure files
file.copy(
  from = file.path("figures", "blood", files),
  to = file.path(file_dir, names(files))
)
