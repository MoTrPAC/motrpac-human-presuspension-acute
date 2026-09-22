#!/usr/bin/env Rscript
# Supplementary Table 6 — the merged EE+RE SC-ION network
#
# Sub-tables: ST6a  every edge: regulator (regulates) target, with the
#                   regulator's phosphosite for a TF, the target's feature id,
#                   modality, tissue, c-means cluster (blank for a hub edge)
#                   and weight. The source writes the literal string NA in site
#                   for a metabolite edge; that is missing.
#             ST6b  every node's Network Motif Score in the EE and RE networks
#                   and their mean, blank where the node is in one network only.
#
# Both are re-emitted from files vendored under figure_6/sources/, under the
# submitted workbook's headers. The inference, the permutation trim at edge
# weight 0.1 and the Cytoscape merge that produced them are not run here; see
# docs/external_dependencies.md.
#
# The rest of ST6 is a panel's own numbers and is written by that panel's
# figure script, which config/table_map.json records:
#   ST6c  TFEB ChIP targets     FIG6.R  (FIG6D)
#   ST6d  EE 3.5-4 h up         FIG6.R  (FIG6D)
#   ST6e  EE 3.5-4 h down       FIG6.R  (FIG6D)
#   ST6f  RE 3.5-4 h up         FIG6.R  (FIG6D)
#   ST6g  RE 3.5-4 h down       FIG6.R  (FIG6D)
#
# Neither sub-table here needs consortium data access: both read a vendored
# file and nothing else.
#
#   Rscript figures/landscape/tables/ST6.R          every sub-table here
#   Rscript figures/landscape/tables/ST6.R ST6b     one of them

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
# panel_export.R first, for panel_source() and landscape_root(); table_export.R
# last, so its manifest helpers are the ones in effect.
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(landscape_root(), "figure_6", "FIG6_ED8_helpers.R"))
source(file.path(here, "..", "lib", "table_export.R"))

# ---- shared ----------------------------------------------------------------

# There is no shared accessor. ST6a and ST6b are two sub-tables of one network,
# but they are read from two different vendored files and share no computation,
# so each entry does its own read and both ids stay selectable.

# ---- ST6a — the merged network's edges -------------------------------------

st6a <- function() {
  table_init("ST6a")

  edges <- scion_edges()
  edges$site[!is.na(edges$site) & edges$site == "NA"] <- NA
  message(sprintf("        %d edges: %s; %d hub edges without a cluster",
                  nrow(edges),
                  paste(sprintf("%s %d", names(table(edges$modality)), table(edges$modality)),
                        collapse = ", "),
                  sum(is.na(edges$cluster))))

  export_table(edges[, table_spec("ST6a")$columns], "ST6a")
}

# ---- ST6b — the merged network's node scores -------------------------------

st6b <- function() {
  table_init("ST6b")

  nms_path <- panel_source("figure_6", "Merged_Network_combined_NMS_v2.csv",
                           env_var = "FIG6_SCION_NMS_CSV")
  nms <- utils::read.csv(nms_path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("name", "nms.score.EE", "nms.score.RE", "nms.score.avg")
  missing_columns <- setdiff(required, names(nms))
  if (length(missing_columns) > 0) {
    stop("the NMS table is missing column(s): ", paste(missing_columns, collapse = ", "),
         call. = FALSE)
  }
  nodes <- data.frame(name = nms$name, nms.EE = nms$nms.score.EE, nms.RE = nms$nms.score.RE,
                      nms.score.avg = nms$nms.score.avg, stringsAsFactors = FALSE)
  message(sprintf("        %d nodes: %d in EE, %d in RE, %d in both",
                  nrow(nodes), sum(!is.na(nodes$nms.EE)), sum(!is.na(nodes$nms.RE)),
                  sum(!is.na(nodes$nms.score.avg))))

  export_table(nodes[, table_spec("ST6b")$columns], "ST6b")
}

run_tables(list(ST6a = st6a, ST6b = st6b))
