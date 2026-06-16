#' Create a Sankey Plot with Colored Nodes
#'
#' This function generates a Sankey plot from a given dataframe,
#'   visualizing relationships
#' between genes and terms. The node colors are determined
#'   by the provided palette.
#'
#' @param df A data frame containing gene-term associations with the
#'   following required columns:
#'   - `P.value`: Numeric, p-values for the terms.
#'   - `gene_count`: Numeric, count of genes associated with each term.
#'   - `Term`: Character, names of the terms.
#'   - `Odds.Ratio`: Numeric, odds ratio values.
#'   - `status`: Character, indicating if a node is "up", "down",
#'               or another category.
#' @param palette A named vector of colors.
#'   The colors should be named according to the
#'   categories in the `status` column, with special colors for "up" and "down".
#' @param gene_size : text size for genes
#' @param term_size : text size for terms
#'
#' @importFrom ggplot2 ggplot aes scale_x_discrete scale_size_manual
#' @importFrom ggplot2 theme_minimal theme element_blank scale_fill_manual
#' @importFrom ggsankeyfier geom_sankeynode geom_sankeyedge

#' @return A ggplot object representing the Sankey diagram.
sankey_plot_ <- function(
  df,
  v_space = "auto",
  term_size = 3,
  gene_size = 1.5,
  palette = NULL
) {

  if (is.null(palette)) {
    palette <- c(
      rlang::set_names(
        viridis::viridis(dplyr::n_distinct(df$Term)),
        setdiff(unique(df$Term), c("up", "down"))
      ),
      "up" = "#317AC1",
      "down" = "#E1A624"
    )
  }

  df <- df |>
    dplyr::group_by(.data$Term) |>
    tidyr::pivot_stages_longer(
      stages_from = c("Genes", "Term"),
      values_from = "gene_count",
      additional_aes_from = c("Term", "Odds.Ratio", "logP", "status")
    ) |>
    dplyr::mutate(
      status = ifelse(
        stage == "Genes",
        .data$status,
        as.character(.data$node)
      )
    )

  pos <- position_sankey(
    nudge_x = -0.06,
    v_space = v_space,
    order = "ascending"
  )

  sankey_plot <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$stage,
      y = .data$gene_count,
      group = .data$node,
      connector = .data$connector,
      edge_id = .data$edge_id
    )
  ) +
    ggsankeyfier::geom_sankeynode(
      v_space = v_space,
      aes(fill = .data$status)
    ) +
    ggsankeyfier::geom_sankeyedge(
      v_space = v_space,
      aes(fill = .data$Term, label = .data$status)
    ) +
    ggsankeyfier::geom_text(
      aes(label = .data$node, cex = .data$stage),
      position = pos,
      stat = "sankeynode",
      hjust = 1
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_discrete(
      expand = ggplot2::expansion(add = c(0.2, 0))
    ) +
    ggplot2::scale_size_manual(
      values = c(gene_size, term_size),
      guide = "none"
    ) +
    ggplot2::scale_fill_manual(
      values = palette,
      guide = "none"
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  sankey_plot
}

#' Create a Dot Plot for Gene-Term Associations
#'
#' This function generates a dot plot that visualizes gene-term associations
#' based on odds ratio, gene count, and statistical significance.
#'
#' @param df A data frame containing gene-term associations with the following
#'   required columns:
#'   - `P.value`: Numeric, p-values for the terms.
#'   - `gene_count`: Numeric, count of genes associated with each term.
#'   - `Term`: Character, names of the terms.
#'   - `Odds.Ratio`: Numeric, odds ratio values.
#'
#' @importFrom ggplot2 ggplot aes scale_color_viridis_c
#' @importFrom ggplot2 coord_cartesian theme_minimal theme
#' @importFrom ggplot2 element_blank element_text element_rect
#' @importFrom ggplot2 geom_point geom_rect
#' @return A ggplot object representing the dot plot.
#'
dot_plot <- function(df) {

  dot_plot <- df |>
    dplyr::group_by(
      .data$Term,
      .data$P.value,
      .data$yy,
      .data$status,
      .data$gene_count
    ) |>
    dplyr::mutate(Odds.Ratio = max(.data$Odds.Ratio)) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = .data$Odds.Ratio,
        y = .data$yy,
        size = .data$gene_count,
        color = .data$logP
      )
    ) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::geom_rect(ggplot2::aes(
      xmin = min(.data$Odds.Ratio) - 0.05 * max(.data$Odds.Ratio),
      ymin = 0,
      xmax = max(.data$Odds.Ratio) + 0.05 * max(.data$Odds.Ratio),
      ymax = max(.data$yy) + 0.2 * max(.data$yy)
    ),
    color = "black",
    fill = NA,
    linewidth = 0.1
    ) +
    scale_color_viridis_c(name = "-log10(P.value)") +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_text(size = 9),
      legend.key.size = ggplot2::unit(0.4, "cm"),
      legend.spacing.y = ggplot2::unit(0.2, "cm"),
      legend.background = ggplot2::element_rect(
        color = "black",
        fill = "white",
        linewidth = 0.5
      ),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  dot_plot
}


#' Create a Combined Sankey and Dot Plot
#'
#' This function generates a composite visualization that includes both a
#' Sankey plot and a dot plot, allowing for an intuitive representation of
#' gene-term associations.
#'
#' @param df A data frame containing gene-term associations with the following
#'   required columns:
#'   - `P.value`: Numeric, p-values for the terms.
#'   - `gene_count`: Numeric, count of genes associated with each term.
#'   - `Term`: Character, names of the terms.
#'   - `Odds.Ratio`: Numeric, odds ratio values.
#'   - `status`: Character, indicating node status
#'   (e.g., "up", "down", or other terms).
#' @param palette A named vector of colors for the terms and statuses.
#'   If `NULL`,
#'   the function generates a default palette using `viridis::viridis()`,
#'   with fixed colors for `"up"` (`#317AC1`) and `"down"` (`#E1A624`).
#'
#' @importFrom patchwork plot_layout
#' @importFrom ggplot2 scale_y_continuous ggplot_build
#' @importFrom dplyr arrange mutate
#' @importFrom forcats fct_reorder
#' @importFrom httr POST content status_code
#'
#' @return A ggplot object combining a Sankey plot and a dot plot.
#'
#' @export
#'
sankeydot_plot <- function(
  df,
  v_space = "auto",
  term_size = 3,
  gene_size = 1.5,
  palette = NULL
) {

  df <- df |>
    dplyr::mutate(
      logP = -log10(.data$P.value),
      Term = forcats::fct_reorder(.data$Term, .data$gene_count)
    ) |>
    dplyr::arrange(.data$gene_count, .data$Term)

  sp <- sankey_plot_(df, v_space, term_size, gene_size, palette)
  spb <- ggplot2::ggplot_build(sp)
  df <- df |> dplyr::arrange(.data$Term)
  df$yy <- spb$data[[2]]$yend_node
  ymax <- max(spb$data[[2]]$y, spb$data[[2]]$yend_node)
  ymax <- ymax + 0.1 * max(ymax)

  dp <- dot_plot(df)
  sdp <- sp + dp + patchwork::plot_layout(widths = c(6, 2)) &
    ggplot2::scale_y_continuous(limits = c(0, ymax))
  sdp
}


get_enrichr_url <- function(genes, description = "Ma liste de gènes") {
  # 1. On formate les gènes en mettant un gène par ligne (\n)
  gene_str <- paste(genes, collapse = "\n")

  # 2. On envoie la requête POST à l'API d'Enrichr
  response <- httr::POST(
    url = "https://maayanlab.cloud/Enrichr/addList",
    body = list(list = gene_str, description = description),
    encode = "multipart"
  )

  # 3. On extrait l'ID unique généré par Enrichr
  if (httr::status_code(response) == 200) {
    result <- httr::content(response, as = "parsed", type = "application/json")
    user_list_id <- result$userListId

    # 4. On génère l'URL finale pour l'utilisateur
    url <- paste0(
      "https://maayanlab.cloud/Enrichr/enrich?userListId=",
      user_list_id
    )
    return(url)
  }
  NULL # En cas d'erreur de l'API
}