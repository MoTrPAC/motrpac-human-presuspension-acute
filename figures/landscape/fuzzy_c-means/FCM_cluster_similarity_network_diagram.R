## NOTE: These plots are not used. Replaced by cluster correlation network
## diagrams.

## Author: Tyler Sagendorf
##
## Date: 2025-07-28
##
## Purpose: Construct a set similarity network diagram to show connections
## between fuzzy c-means clusters within and between tissues separately by ome.

library(MotrpacHumanPreSuspensionData) # FCM_CAMERA
library(MotrpacHumanPreSuspension) # HUMAN_TISSUE_COLORS
library(dplyr)
library(ggplot2)
library(igraph)
library(ggnetwork)
library(TMSig) # similarity()

# List of significant molecular signatures for each combination of tissue and
# cluster.
df <- FCM_CAMERA %>%
  filter(adj_p_value < 0.05 + 0.05 * (assay == "prot-ph")) %>%
  mutate(tissue_cluster = paste0(tissue, "_", cluster),
         set = as.character(set)) %>%
  select(tissue_cluster, assay, set)


# Input is a list of significant molecular signatures in each tissue/cluster
# combination.
calc_similarity <- function(df, ome, type = "jaccard") {
  ls <- df %>%
    filter(assay == ome) %>%
    summarise(.by = tissue_cluster,
              set_list = list(set)) %>%
    tibble::deframe()

  # Set similarity matrix
  s <- TMSig::similarity(x = ls,
                         type = type) %>%
    as.matrix()

  return(s)
}


## Set similarity network diagram ----
plot_network <- function(df,
                         type = "jaccard",
                         sim_cutoff = 0.5,
                         phospho = FALSE,
                         seed = NULL)
{
  omes <- c("transcript-rna-seq", "prot-pr", "prot-ph", "metab")

  plot_df <- lapply(omes, function(ome_i) {
    j <- calc_similarity(df, ome = ome_i, type = type)

    j[j < sim_cutoff] <- 0

    # Clusters must be fully unique for plotting
    rownames(j) <- colnames(j) <-
      paste0(ome_i, "_", rownames(j))

    graph_i <- igraph::graph_from_adjacency_matrix(
      adjmatrix = j,
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE
    )

    set.seed(seed)

    n <- ggnetwork(graph_i) %>%
      tidyr::separate_wider_delim(
        cols = name,
        delim = "_",
        names = c("ome", "tissue", "cluster")
      )

    return(n)
  }) %>%
    bind_rows() %>% # stack omes
    mutate(ome = factor(ome,
                        levels = omes,
                        labels = c("Transcriptomics",
                                   "Proteomics",
                                   "Phosphoproteomics",
                                   "Metabolomics")))

  plot_rows <- 1L + phospho

  # Remove phospho data
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
               size = 8, shape = 21) +
    geom_nodetext(aes(label = cluster)) +
    facet_wrap(~ ome, nrow = plot_rows) +
    scale_fill_manual(name = "Tissue",
                      values = tissue_colors,
                      breaks = names(tissue_colors)) +
    scale_linewidth(
      name = paste0(
        sub("(^.{1})(.*$)", "\\U\\1\\L\\2", perl = TRUE, type), # title case
        "\nSimilarity"
      ),
      limits = c(sim_cutoff, 1),
      breaks = scales::breaks_pretty(n = 4L)
    ) +
    theme_blank(base_size = 10) +
    theme(panel.background = element_rect(color = "black"),
          strip.background = element_blank(),
          strip.text = element_text(size = rel(1.5)))
}


folder <- file.path("figures", "landscape", "fuzzy_c-means")

## Overlap similarity ----
p_o <- plot_network(
  df = df,
  type = "overlap",
  sim_cutoff = 0.3,
  seed = 9001
)

p_o

ggsave(path = folder,
       filename = "FCM_cluster_overlap_similarity_network_diagram.pdf",
       plot = p_o,
       height = 5,
       width = 10,
       units = "in")

## Jaccard similarity ----
p_j <- plot_network(
  df = df,
  type = "jaccard",
  sim_cutoff = 0.3,
  seed = 9001
)

p_j

ggsave(path = folder,
       filename = "FCM_cluster_jaccard_similarity_network_diagram.pdf",
       plot = p_j,
       height = 5,
       width = 10,
       units = "in")
