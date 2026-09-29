
#' Shared CpG island categories used by all enrichment plots.
.cpg_island_levels <- c(
  "OpenSea", "Shelf", "S_Shelf", "S_Shore", "Island", "Shore",
  "N_Shore", "N_Shelf"
)

.cpg_island_colors <- c(
  "OpenSea" = "#A6CEE3", "Shelf" = "#1F78B4", "N_Shelf" = "#1F78B4",
  "S_Shelf" = "#1F78B4", "Shore" = "#B2DF8A", "N_Shore" = "#B2DF8A",
  "S_Shore" = "#B2DF8A", "Island" = "#33A02C"
)

.cpg_island_data <- function(mrs, index, type, fdr, tools) {
  if (type == "dmps") {
    return(get_dmps(mrs, index) |>
      dplyr::filter(.data$adj.P.Val < fdr) |>
      dplyr::mutate(
        tool = NA_character_,
        status = dplyr::if_else(
          .data$deltabetas > 0, "Hypermethylated", "Hypomethylated"
        )
      ))
  }

  get_dmrs(mrs, index, tools = tools) |>
    dplyr::filter(.data$HMFDR < fdr) |>
    tidyr::separate_longer_delim(.data$Relation_to_Island, delim = ";") |>
    dplyr::mutate(
      status = dplyr::if_else(
        .data$mean_deltabeta > 0, "Hypermethylated", "Hypomethylated"
      )
    )
}

.cpg_island_plot <- function(mrs, index, type, position, add_platform, fdr,
                             tools) {
  platform_name <- mrs@metadata$plateform
  dt <- .cpg_island_data(mrs, index, type, fdr, tools)

  if (nrow(dt) == 0) {
    warning("CpG island plot: input data frame is empty, nothing to plot.")
    return(NULL)
  }

  if (add_platform) {
    platform <- mrs@manifest |>
      as.data.frame() |>
      dplyr::select(Probe_ID, Relation_to_Island) |>
      dplyr::mutate(status = platform_name, tool = platform_name)
    dt <- dplyr::bind_rows(dt, platform)
  }

  dt <- dt |>
    dplyr::mutate(
      Relation_to_Island = dplyr::if_else(
        .data$Relation_to_Island %in% .cpg_island_levels,
        .data$Relation_to_Island, "OpenSea"
      ),
      Relation_to_Island = factor(
        .data$Relation_to_Island, levels = .cpg_island_levels
      ),
      status = factor(
        .data$status,
        levels = c(platform_name, "Hypomethylated", "Hypermethylated")
      )
    )

  if (type == "dmrs") {
    dt$tool <- factor(dt$tool, levels = c(platform_name, tools))
    dt$status_tool <- interaction(dt$tool, dt$status, sep = " - ", drop = TRUE)
    x_column <- "status_tool"
    x_label <- "DMR tool and methylation status"
  } else {
    x_column <- "status"
    x_label <- "Methylation status"
  }

  p_text <- ""
  if (add_platform) {
    contingency_table <- table(dt$status, dt$Relation_to_Island)
    p_values <- character()
    for (status in c("Hypermethylated", "Hypomethylated")) {
      if (all(c(platform_name, status) %in% rownames(contingency_table))) {
        mat <- contingency_table[c(platform_name, status), , drop = FALSE]
        mat <- mat[, colSums(mat) > 0, drop = FALSE]
        if (ncol(mat) > 1 && all(rowSums(mat) > 0)) {
          p_value <- stats::chisq.test(mat)$p.value
          if (is.finite(p_value)) {
            p_value <- if (p_value < 0.001) "p < 0.001" else
              paste0("p = ", format(p_value, digits = 3))
            p_values <- c(
              p_values,
              paste0(status, " vs ", platform_name, ": ", p_value)
            )
          }
        }
      }
    }
    p_text <- paste(p_values, collapse = "\n")
  }

  ggplot2::ggplot(dt, ggplot2::aes(x = .data[[x_column]],
                                  fill = .data$Relation_to_Island)) +
    ggplot2::geom_bar(position = position, width = 0.6) +
    ggplot2::scale_fill_manual(values = .cpg_island_colors) +
    ggplot2::labs(
      x = x_label,
      y = if (position == "stack") "Number of regions" else "Percentage (%)",
      fill = "CGI position",
      subtitle = if (p_text == "") NULL else
        paste("Chi-squared enrichment test vs background:\n", p_text)
    ) +
    ggplot2::theme_classic(base_size = 20) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8),
      axis.title = ggplot2::element_text(size = 12),
      plot.subtitle = ggplot2::element_text(size = 10, face = "italic",
                                            color = "darkgrey"),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 11)
    ) +
    ggplot2::coord_flip()
}

#' Plot DMP CpG island enrichment with platform background.
#' @param mrs A MethylResultSet object.
#' @param index Contrast index or name.
#' @param fdr Adjusted p-value threshold.
#' @return A ggplot2 object.
#' @export
cpgislands_plot_dmps_fill <- function(mrs, index, fdr = 0.05) {
  .cpg_island_plot(mrs, index, "dmps", "fill", TRUE, fdr, character())
}

#' Plot DMP CpG island counts without platform background.
#' @param mrs A MethylResultSet object.
#' @param index Contrast index or name.
#' @param fdr Adjusted p-value threshold.
#' @return A ggplot2 object.
#' @export
cpgislands_plot_dmps_stack <- function(mrs, index, fdr = 0.05) {
  .cpg_island_plot(mrs, index, "dmps", "stack", FALSE, fdr, character())
}

#' Plot DMR CpG island enrichment by tool with platform background.
#' @param mrs A MethylResultSet object.
#' @param index Contrast index or name.
#' @param fdr DMR FDR threshold.
#' @param tools DMR tools to include.
#' @return A ggplot2 object.
#' @export
cpgislands_plot_dmrs_fill <- function(
    mrs, index, fdr = 0.05,
    tools = c("dmrcate", "ipdmr", "combp", "dmrff")) {
  .cpg_island_plot(mrs, index, "dmrs", "fill", FALSE, fdr, tools)
}

#' Plot DMR CpG island counts by tool without platform background.
#' @param mrs A MethylResultSet object.
#' @param index Contrast index or name.
#' @param fdr DMR FDR threshold.
#' @param tools DMR tools to include.
#' @return A ggplot2 object.
#' @export
cpgislands_plot_dmrs_stack <- function(
    mrs, index, fdr = 0.05,
    tools = c("dmrcate", "ipdmr", "combp", "dmrff")) {
  .cpg_island_plot(mrs, index, "dmrs", "stack", FALSE, fdr, tools)
}

#' Plot CpG island enrichment (legacy interface).
#' @param mrs A MethylResultSet object.
#' @param index Contrast index or name.
#' @param what Either "dmps" or "dmrs".
#' @param position Either "fill" or "stack".
#' @param fdr Significance threshold.
#' @param tools DMR tools to include.
#' @return A ggplot2 object.
#' @export
cpgislands_plot <- function(mrs, index, what = "dmps", position = "fill",
                            fdr = 0.05,
                            tools = c("dmrcate", "ipdmr", "combp", "dmrff")) {
  if (!what %in% c("dmps", "dmrs") || !position %in% c("fill", "stack")) {
    stop("what must be 'dmps' or 'dmrs' and position must be 'fill' or 'stack'")
  }
  .cpg_island_plot(mrs, index, what, position, position == "fill", fdr, tools)
}

# Historical singular name retained for callers using the old API.
cpg_island <- cpgislands_plot





