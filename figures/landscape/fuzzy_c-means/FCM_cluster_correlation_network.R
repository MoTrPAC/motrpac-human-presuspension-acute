## Author: Tyler Sagendorf
## Date: 2025-05-08
##
## Purpose: Construct a correlation network diagram to show connections between
## fuzzy c-means clusters within and between tissues separately by ome.

library(MotrpacHumanPreSuspensionData) # FCM_CAMERA
library(MotrpacHumanPreSuspension) # HUMAN_TISSUE_COLORS
library(dplyr)
library(ggplot2)
library(igraph)
library(ggnetwork)

# Matrix of p-values from the CAMERA-PR analysis of the fuzzy c-means
# clusters. Rows are molecular signatures, columns are clusters.
p_mat <- FCM_CAMERA %>%
  mutate(set = paste0(assay, ".", set),
         tissue_cluster = paste0(tissue, "_", cluster)) %>%
  tidyr::pivot_wider(id_cols = set,
                     names_from = tissue_cluster,
                     values_from = p_value) %>%
  tibble::column_to_rownames("set") %>%
  as.matrix()

# Generate matrix of Spearman correlations between clusters based on the
# enrichment analysis results.
pearson_cor <- function(mat, ome = NULL) {
  if (!is.null(ome)) {
    omes <- sub("(^[^\\.]+)\\..*", "\\1", rownames(mat))

    mat <- mat[omes == ome, ]

    keep_cols <- colMeans(!is.na(mat)) != 0
    mat <- mat[, keep_cols]
  }

  cor_mat <- cor(x = mat,
                 use = "pairwise.complete.obs",
                 method = "spearman")

  return(cor_mat)
}


## Spearman correlation network diagram ----------------------------------------
plot_network <- function(p_mat,
                         phospho = FALSE,
                         cor_cutoff = 0.5,
                         seed = NULL)
{
  omes <- c("transcript-rna-seq", "prot-pr", "prot-ph", "metab")

  plot_df <- lapply(omes, function(ome_i) {
    cor_mat <- pearson_cor(p_mat, ome = ome_i)
    cor_mat[cor_mat < cor_cutoff] <- 0
    rownames(cor_mat) <- colnames(cor_mat) <-
      paste0(ome_i, "_", rownames(cor_mat))

    graph_i <- igraph::graph_from_adjacency_matrix(
      adjmatrix = cor_mat,
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE
    )

    set.seed(seed)

    n <- ggnetwork(
      x = graph_i,
      layout = igraph::with_kk()
    ) %>%
      tidyr::separate_wider_delim(cols = name, delim = "_",
                                  names = c("ome", "tissue", "cluster"))

    return(n)
  }) %>%
    bind_rows() %>%
    mutate(ome = factor(ome,
                        levels = omes,
                        labels = c("Transcriptomics",
                                   "Proteomics",
                                   "Phosphoproteomics",
                                   "Metabolomics")))

  plot_rows <- 1L + phospho

  if (!phospho) {
    plot_df <- filter(plot_df, ome != "Phosphoproteomics") %>%
      droplevels.data.frame()
  }

  tissue_colors <- HUMAN_TISSUE_COLORS[c("adipose", "blood", "muscle")]

  ggplot(plot_df, aes(x = x, y = y,
                      xend = xend, yend = yend)) +
    geom_edges(aes(lwd = weight),
               color = "black",
               curvature = 0.2) +
    geom_nodes(aes(fill = tissue),
               size = 8,
               shape = 21) +
    geom_nodetext(aes(label = cluster)) +
    facet_wrap(~ ome, nrow = plot_rows) +
    scale_fill_manual(name = "Tissue",
                      values = tissue_colors,
                      breaks = names(tissue_colors)) +
    scale_linewidth(name = "Spearman\nCorrelation",
                    limits = c(cor_cutoff, 1),
                    breaks = scales::breaks_pretty(n = 4L)) +
    theme_blank(base_size = 10) +
    theme(panel.background = element_rect(color = "black"),
          strip.background = element_blank(),
          strip.text = element_text(size = rel(1.5)))
}

p <- plot_network(p_mat, cor_cutoff = 0.4, seed = 9001L)
p

folder <- file.path("figures", "landscape", "fuzzy_c-means")

ggsave(path = folder,
       filename = "FCM_cluster_correlation_network_diagram.pdf",
       plot = p, height = 6, width = 10, units = "in")


## Include phosphoproteomics
p2 <- plot_network(p_mat, phospho = TRUE, cor_cutoff = 0.4, seed = 9001)
p2

ggsave(path = folder,
       filename = "FCM_cluster_correlation_network_diagram_with_phospho.pdf",
       plot = p2, height = 7, width = 10, units = "in")
