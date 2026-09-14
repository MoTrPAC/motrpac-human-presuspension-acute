options(max.print=50)

ggplot2::theme_set(ggplot2::theme_classic(base_size = 6))
ggplot2::theme_update(strip.background = ggplot2::element_rect(fill = "snow2", linetype = "blank"),
                      strip.text = ggplot2::element_text(face = "bold"),
                      axis.text = ggplot2::element_text(color = "black"),
                      axis.title = ggplot2::element_text(face = "bold"),
                      axis.line = element_line(linewidth = rel(0.5), lineend = "square"),
                      axis.ticks = element_line(linewidth = rel(0.5), color = "black"),
                      legend.margin = margin(0, 0, 0, 0, "cm"))

aesthetically_pleasing_colors = list(
  colorBlindFriendly = rcartocolor::carto_pal(12, "Safe"),
  colorBlindFriendly2 = RColorBrewer::brewer.pal(8, "Set2")
)

attach(list(
  prettyLabels = function(x, short = FALSE) {
    x_og = as.character(x)
    x_new = x_og |>
      str_remove_all("pre_and_|(?<!^)ADU.*_(?=pre)") |>
      str_replace_all(c("15_30_45_min" = "15-45 min", "3.5_4_hr" = "3.5-4 hr", "24_hr" = "24 hr", "vs" = "/", "_" = " ", "ADUEndur" = "EE", "ADUResist" = "RE", "ADUControl" = "CON")) |>
      str_remove_all(" (?=hr|min)") |>
      str_replace_all(c("(?<=^)EE(?= )" = "EE:", "(?<=^)RE(?= )" = "RE:", "(?<=^)Ctrl(?= )" = "Ctrl:"))
    if (short) {
      x_new = str_replace_all(x_new, c("pre exercise" = "Pre", "post " = "P"))
    }
    return(x_new)
  },
  invertNames = function(x) {
    return(purrr::set_names(names(x), unname(x)))
  },
  pcaWrapper = function(counts, metadata, filter_by, covariates, protect) {
    pheno_df = metadata |>
      remove_rownames() |>
      column_to_rownames(var = "library") |>
      filter_by() |>
      select(any_of(covariates)) |>
      select(!where(~ n_distinct(.x) == 1)) |>
      droplevels()
    protect = intersect(protect, colnames(pheno_df))
    adj_counts = PeakGeneNet::adjustCovariateMatrix(counts[,rownames(pheno_df)] |> t(), pheno_df, protect) |>
      scale() |>
      prcomp(rank. = 10) |>
      pcaPostProcess(metadata) |>
      `attr<-`("in_model", colnames(pheno_df)) |>
      `attr<-`("protected", protect)
    return(adj_counts)
  },
  getPCAVariance = function(pca.out) {
    pca.out$variance = scales::percent(pca.out$sdev^2 / sum(pca.out$sdev^2), 0.01)[1:sum(grepl("^PC", colnames(pca.out$x)))] |>
      purrr::set_names(grep("^PC", colnames(pca.out$x), value = TRUE)) |>
      imap_chr(~ paste0(.y, " (", .x, ")"))
    return(pca.out)
  },
  pcaPostProcess = function(pca_obj, sample_data) {
    pca_obj = getPCAVariance(pca_obj)
    pca_obj$x = pca_obj$x |>
      as.data.frame() |>
      rownames_to_column(var = "library") |>
      left_join(sample_data, by = "library")
    return(pca_obj)
  },
  binarizeVar = function(fct, prefix = NULL) {
    # provide a factor to produce a binarized dataframe of that variable
    factor_options = levels(fct)
    factor_options = intersect(factor_options, unique(as.character(fct)))
    df = data.frame(og_fct = fct)
    for (x in factor_options) {
      if (!is.null(prefix)) {
        col_nm = paste0(prefix, "___", x)
      } else {
        col_nm = x
      }
      df[[col_nm]] = fct == x
    }
    df = df |>
      select(-1) |>
      mutate(across(everything(), as.numeric))
    return(df)
  }
), name = "misc")
