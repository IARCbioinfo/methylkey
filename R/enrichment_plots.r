
#' Shared CpG island categories used by all enrichment plots.
.cpg_island_levels <- c(
  "OpenSea", "Shelf", "S_Shelf", "S_Shore", "Island", "Shore",
  "N_Shore", "N_Shelf"
)

.cpg_island_priority <- c(
  "Island", "Shore", "N_Shore", "S_Shore", "Shelf", "N_Shelf",
  "S_Shelf", "OpenSea"
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

  get_dmrs(mrs, index, tools = tools, max_fdr = fdr) |>
    dplyr::mutate(
      Relation_to_Island = vapply(
        strsplit(as.character(.data$Relation_to_Island), ";", fixed = TRUE),
        function(categories) {
          priority_matches <- .cpg_island_priority[
            .cpg_island_priority %in% trimws(categories)
          ]
          if (length(priority_matches) == 0L) "OpenSea" else
            priority_matches[[1L]]
        },
        character(1)
      ),
      status = dplyr::if_else(
        .data$mean_deltabeta > 0,
        "Hypermethylated",
        "Hypomethylated",
        missing = "Unknown"
      )
    )
}

.cpg_island_plot <- function(mrs, index, type, position, add_platform, fdr,
                             tools, dn = 1) {
  if (length(dn) != 1L || is.na(dn) || !dn %in% 0:3) {
    stop("dn must be one of 0, 1, 2, or 3")
  }

  add_platform <- isTRUE(add_platform)
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

  status_levels <- c("Hypomethylated", "Hypermethylated", "Unknown")
  #if (type == "dmps") {
    status_levels <- c(platform_name, status_levels)
  #}

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
        levels = status_levels
      )
    )

  if (type == "dmrs") {
    tool_levels <- if (add_platform) c(platform_name, tools) else tools
    dt$tool <- factor(dt$tool, levels = tool_levels)
    dt$status_tool <- interaction(dt$tool, dt$status, sep = " - ", drop = TRUE)
    x_column <- "status_tool"
    x_label <- "DMR tool and methylation status"
    y_label <- ifelse(
      position == "stack",
      "Number of probes",
      "Percentage (%)"
    )
    if (add_platform) {
      lv <- levels(dt$status_tool)
      dup <- paste(platform_name, platform_name, sep = " - ")
      levels(dt$status_tool)[lv == dup] <- platform_name
    }
    
  } else {
    x_column <- "status"
    x_label <- "Methylation status"
    y_label <- ifelse(
      position == "stack",
      "Number of probes",
      "Percentage (%)"
    )
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

  plot <- ggplot2::ggplot(dt, ggplot2::aes(x = .data[[x_column]],
                                         fill = .data$Relation_to_Island)) +
    ggplot2::geom_bar(position = position, width = 0.6)

  label_data <- dt
  if (dn %in% c(1, 3)) {
    label_position <- if (position == "fill") {
      ggplot2::position_fill(vjust = 0.5)
    } else {
      ggplot2::position_stack(vjust = 0.5)
    }
    plot <- plot + ggplot2::geom_text(
      data = label_data,
      ggplot2::aes(label = ggplot2::after_stat(count)),
      stat = "count",
      position = label_position,
      size = 3
    )
  }

  if (dn %in% c(2, 3)) {
    total_data <- label_data |>
      dplyr::group_by(.data[[x_column]]) |>
      dplyr::summarise(total = dplyr::n(), .groups = "drop") |>
      dplyr::mutate(
        label_y = if (position == "fill") 1.04 else
          .data$total + pmax(.data$total * 0.05, 0.5)
      )
    plot <- plot + ggplot2::geom_text(
      data = total_data,
      ggplot2::aes(
        x = .data[[x_column]], y = .data$label_y, label = .data$total
      ),
      inherit.aes = FALSE,
      hjust = 0,
      size = 3
    )
    plot <- plot + ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.2))
    )
  }

  plot +
    ggplot2::scale_fill_manual(values = .cpg_island_colors) +
    ggplot2::labs(
      x = x_label,
      y = y_label,
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
#' @param dn Number display mode: 0 none, 1 counts within segments,
#'   2 bar totals, or 3 both.
#' @return A ggplot2 object.
#' @export
cpgislands_plot_dmps_fill <- function(mrs, index, fdr = 0.05, dn = 1) {
  .cpg_island_plot(mrs, index, "dmps", "fill", TRUE, fdr, character(), dn)
}

#' Plot DMP CpG island counts without platform background.
#' @param mrs A MethylResultSet object.
#' @param index Contrast index or name.
#' @param fdr Adjusted p-value threshold.
#' @param dn Number display mode: 0 none, 1 counts within segments,
#'   2 bar totals, or 3 both.
#' @return A ggplot2 object.
#' @export
cpgislands_plot_dmps_stack <- function(mrs, index, fdr = 0.05, dn = 1) {
  .cpg_island_plot(mrs, index, "dmps", "stack", FALSE, fdr, character(), dn)
}

#' Plot DMR CpG island enrichment by tool with platform background.
#' @param mrs A MethylResultSet object.
#' @param index Contrast index or name.
#' @param fdr DMR FDR threshold.
#' @param tools DMR tools to include.
#' @param dn Number display mode: 0 none, 1 counts within segments,
#'   2 bar totals, or 3 both.
#' @return A ggplot2 object.
#' @export
cpgislands_plot_dmrs_fill <- function(
    mrs, index, fdr = 0.05,
    tools = c("dmrcate", "ipdmr", "combp", "dmrff"), dn = 1) {
  .cpg_island_plot(mrs, index, "dmrs", "fill", TRUE, fdr, tools, dn)
}

#' Plot DMR CpG island counts by tool without platform background.
#' @param mrs A MethylResultSet object.
#' @param index Contrast index or name.
#' @param fdr DMR FDR threshold.
#' @param tools DMR tools to include.
#' @param dn Number display mode: 0 none, 1 counts within segments,
#'   2 bar totals, or 3 both.
#' @return A ggplot2 object.
#' @export
cpgislands_plot_dmrs_stack <- function(
    mrs, index, fdr = 0.05,
    tools = c("dmrcate", "ipdmr", "combp", "dmrff"), dn = 1) {
  .cpg_island_plot(mrs, index, "dmrs", "stack", FALSE, fdr, tools, dn)
}

#' Plot CpG island enrichment (legacy interface).
#' @param mrs A MethylResultSet object.
#' @param index Contrast index or name.
#' @param what Either "dmps" or "dmrs".
#' @param position Either "fill" or "stack".
#' @param fdr Significance threshold.
#' @param tools DMR tools to include.
#' @param dn Number display mode: 0 none, 1 counts within segments,
#'   2 bar totals, or 3 both.
#' @return A ggplot2 object.
#' @export
cpgislands_plot <- function(mrs, index, what = "dmps", position = "fill",
                            fdr = 0.05,
                            tools = c("dmrcate", "ipdmr", "combp", "dmrff"),
                            dn = 1) {
  if (!what %in% c("dmps", "dmrs") || !position %in% c("fill", "stack")) {
    stop("what must be 'dmps' or 'dmrs' and position must be 'fill' or 'stack'")
  }
  .cpg_island_plot(mrs, index, what, position, position == "fill", fdr, tools,
                   dn)
}

# Historical singular name retained for callers using the old API.
cpg_island <- cpgislands_plot
