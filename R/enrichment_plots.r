
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

.cpg_island_data <- function(
  mrs,
  index = 1,
  type = "dmps",
  q = 0.05,
  method = c("fdr", "BH", "bonferroni"),
  tools = c("dmrcate", "ipdmr", "combp", "dmrff"),
  combine = NULL,
  min_size = 1,
  max_size = Inf,
  min_cpgs = 2,
  max_cpgs = Inf
) {
  if (type == "dmps") {
    return(get_dmps(mrs, index, q, method) |>
      dplyr::mutate(
        tool = NA_character_,
        status = dplyr::if_else(
          .data$deltabetas > 0, "Hypermethylated", "Hypomethylated"
        )
      ))
  }

  get_dmrs(
    mrs,
    index,
    q = q,
    method = method,
    tools = tools,
    combine = combine,
    min_size = min_size,
    max_size = max_size,
    min_cpgs = min_cpgs,
    max_cpgs = max_cpgs
  ) |>
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

.cpg_island_plot <- function(mrs, index, type, position, add_platform,
                             q = 0.05, method = c("fdr", "BH", "bonferroni"),
                             tools = c("dmrcate", "ipdmr", "combp", "dmrff"),
                             combine = NULL, min_size = 1, max_size = Inf,
                             min_cpgs = 1, max_cpgs = Inf, dn = 1,
                             separate_by_tool = TRUE) {
  method <- match.arg(method, c("fdr", "BH", "bonferroni"))
  if (length(dn) != 1L || is.na(dn) || !dn %in% 0:3) {
    stop("dn must be one of 0, 1, 2, or 3")
  }
  if (length(separate_by_tool) != 1L || is.na(separate_by_tool) ||
      !is.logical(separate_by_tool)) {
    stop("separate_by_tool must be TRUE or FALSE")
  }
  if (type == "dmrs" && !is.null(combine) && separate_by_tool) {
    warning(
      "`combine` is incompatible with `separate_by_tool = TRUE`; ",
      "setting `separate_by_tool = FALSE`.",
      call. = FALSE
    )
    separate_by_tool <- FALSE
  }
  if (length(index) == 0L || anyNA(index)) {
    stop("index must contain at least one valid contrast index or name")
  }
  if (!is.numeric(index) && !is.character(index)) {
    stop("index must contain numeric indices or character names")
  }
  if (anyDuplicated(index)) {
    stop("index values must be unique")
  }
  if (is.numeric(index) &&
      any(!is.finite(index) | index < 1 | index > length(mrs@dmps) |
            index != as.integer(index))) {
    stop("numeric index values must identify existing contrasts")
  }
  if (is.character(index) && any(!nzchar(index))) {
    stop("index names must not be empty")
  }
  if (is.character(index) &&
      any(!index %in% names(mrs@dmps))) {
    stop("all index names must identify existing contrasts")
  }

  add_platform <- isTRUE(add_platform)
  platform_name <- mrs@metadata$plateform
  index_labels <- if (is.numeric(index)) {
    contrast_names <- names(mrs@dmps)[index]
    ifelse(is.na(contrast_names), as.character(index), contrast_names)
  } else {
    as.character(index)
  }
  dt <- dplyr::bind_rows(lapply(seq_along(index), function(i) {
    .cpg_island_data(mrs, index[[i]], type, q, method, tools, combine,
  min_size, max_size, min_cpgs, max_cpgs) |>
      dplyr::mutate(index = index_labels[[i]])
  }))

  if (nrow(dt) == 0) {
    warning("CpG island plot: input data frame is empty, nothing to plot.")
    return(NULL)
  }

  if (add_platform) {
    platform_data <- mrs@manifest |>
      as.data.frame() |>
      dplyr::select(Probe_ID, Relation_to_Island) |>
      dplyr::mutate(status = platform_name, tool = platform_name)
    platform <- dplyr::bind_rows(lapply(index_labels, function(label) {
      platform_data |>
        dplyr::mutate(index = label)
    }))
    dt <- dplyr::bind_rows(dt, platform)
  }

  status_levels <- c("Hypomethylated", "Hypermethylated", "Unknown")
  status_levels <- c(platform_name, status_levels)

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

  y_label <- ifelse(
    position == "stack",
    "Number of probes",
    "Percentage (%)"
  )

  plot_data <- dt
  if (add_platform && length(index_labels) > 1L) {
    plot_data <- dplyr::filter(
      plot_data,
      .data$status != platform_name |
        .data$index == index_labels[[1L]]
    )
  }
  plot_data$plot_status <- ifelse(
    plot_data$status == platform_name,
    "All probes",
    as.character(plot_data$status)
  )

  plot_status_levels <- c("Hypermethylated", "Hypomethylated", "Unknown")
  plot_colors <- c(
    "Hypermethylated" = "#2E86AB",
    "Hypomethylated" = "#A23B72",
    "All probes" = "#000000",
    "Unknown" = "#000000"
  )
  if (type == "dmrs" && separate_by_tool) {
    tool_column <- if ("tool" %in% names(plot_data)) {
      "tool"
    } else if ("tools" %in% names(plot_data)) {
      "tools"
    } else {
      stop("DMR plot data must contain a 'tool' or 'tools' column.")
    }
    tool_values <- as.character(plot_data[[tool_column]])
    present_tools <- unique(tool_values[
      plot_data$plot_status != "All probes" &
        !is.na(tool_values) & nzchar(tool_values)
    ])
    used_tools <- if (tool_column == "tool") {
      tools[tools %in% present_tools]
    } else {
      sort(present_tools)
    }
    plot_category <- ifelse(
      plot_data$plot_status == "All probes",
      "All probes",
      paste(plot_data$plot_status, plot_data$index,
          tool_values, sep = "\034")
    )
    contrast_categories <- unlist(lapply(plot_status_levels, function(status) {
      unlist(lapply(index_labels, function(label) {
        paste(status, label, used_tools, sep = "\034")
      }), use.names = FALSE)
    }), use.names = FALSE)
    category_labels <- unlist(lapply(plot_status_levels, function(status) {
      unlist(lapply(index_labels, function(label) {
        paste(label, used_tools, sep = ": ")
      }), use.names = FALSE)
    }), use.names = FALSE)
    x_label <- if (length(index_labels) > 1L) {
      "Contrast and DMR tool"
    } else {
      "DMR tool"
    }
  } else {
    plot_category <- ifelse(
      plot_data$plot_status == "All probes",
      "All probes",
      paste(plot_data$plot_status, plot_data$index, sep = "\034")
    )
    contrast_categories <- unlist(lapply(plot_status_levels, function(status) {
      paste(status, index_labels, sep = "\034")
    }), use.names = FALSE)
    category_labels <- rep(index_labels, times = length(plot_status_levels))
    x_label <- "Contrast"
  }
  plot_levels <- contrast_categories
  category_colors <- rep(plot_colors[plot_status_levels], each = if (
    type == "dmrs" && separate_by_tool
  ) {
    length(index_labels) * length(used_tools)
  } else {
    length(index_labels)
  })
  if (add_platform && any(plot_data$plot_status == "All probes")) {
    plot_levels <- c(plot_levels, "All probes")
    category_labels <- c(category_labels, "All probes")
    category_colors <- c(category_colors, plot_colors[["All probes"]])
  }
  plot_levels <- rev(plot_levels)
  category_labels <- rev(category_labels)
  category_colors <- rev(category_colors)
  keep_categories <- plot_levels %in% plot_category
  plot_levels <- plot_levels[keep_categories]
  category_labels <- category_labels[keep_categories]
  category_colors <- category_colors[keep_categories]
  plot_data$plot_category <- factor(plot_category, levels = plot_levels)
  axis_labels <- paste0(
    "<span style='color:", category_colors, "'>",
    gsub("&", "&amp;", gsub("<", "&lt;", category_labels, fixed = TRUE),
         fixed = TRUE),
    "</span>"
  ) |>
    stats::setNames(plot_levels)
  x_column <- "plot_category"

  #stat <- data.frame(
  #  index = character(),
  #  status = character(),
  #  background = character(),
  #  Relation_to_Island = character(),
  #  observed = integer(),
  #  expected = numeric(),
  #  residual = numeric(),
  #  stdres = numeric(),
  #  statistic = numeric(),
  #  df = numeric(),
  #  p.value = numeric(),
  #  method = character(),
  #  data.name = character(),
  #  stringsAsFactors = FALSE
  #)
  # if (add_platform) {
  #   for (label in index_labels) {
  #     index_data <- dplyr::filter(dt, .data$index == label)
  #     contingency_table <- table(index_data$status,
  #                                index_data$Relation_to_Island)
  #     for (status in c("Hypermethylated", "Hypomethylated")) {
  #       if (all(c(platform_name, status) %in% rownames(contingency_table))) {
  #         mat <- contingency_table[c(platform_name, status), , drop = FALSE]
  #         mat <- mat[, colSums(mat) > 0, drop = FALSE]
  #         if (ncol(mat) > 1 && all(rowSums(mat) > 0)) {
  #           test <- stats::chisq.test(mat)
  #           stat <- dplyr::bind_rows(
  #             stat,
  #             data.frame(
  #               index = label,
  #               status = status,
  #               background = platform_name,
  #               Relation_to_Island = rep(
  #                 colnames(mat), each = nrow(mat)
  #               ),
  #               observed = as.vector(test$observed),
  #               expected = as.vector(test$expected),
  #               residual = as.vector(test$residuals),
  #               stdres = as.vector(test$stdres),
  #               statistic = as.numeric(test$statistic),
  #               df = as.numeric(test$parameter),
  #               p.value = test$p.value,
  #               method = test$method,
  #               data.name = test$data.name,
  #               stringsAsFactors = FALSE
  #             )
  #           )
  #         }
  #       }
  #     }
  #   }
  # }

  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = .data[[x_column]], fill = .data$Relation_to_Island
  )) +
    ggplot2::geom_bar(position = position, width = 0.6)

  legend_status <- c(
    "Hypermethylated", "Hypomethylated",
    if (add_platform) "EPIC v2"
  )
  legend_colors <- c(
    "Hypermethylated" = plot_colors[["Hypermethylated"]],
    "Hypomethylated" = plot_colors[["Hypomethylated"]],
    "EPIC v2" = plot_colors[["All probes"]]
  )
  legend_data <- data.frame(
    plot_category = factor(
      rep(plot_levels[[1L]], length(legend_status)),
      levels = plot_levels
    ),
    legend_status = factor(legend_status, levels = legend_status),
    y = 0
  )
  plot <- plot + ggplot2::geom_point(
    data = legend_data,
    ggplot2::aes(
      x = .data$plot_category, y = .data$y, colour = .data$legend_status
    ),
    inherit.aes = FALSE,
    alpha = 0,
    show.legend = TRUE
  )

  label_data <- plot_data
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
      dplyr::group_by(.data[[x_column]], .data$plot_status) |>
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

  plot <- plot +
    ggplot2::scale_x_discrete(
      limits = plot_levels,
      labels = axis_labels
    ) +
    ggplot2::scale_colour_manual(
      name = "Methylation status",
      values = legend_colors[legend_status],
      breaks = legend_status
    ) +
    ggplot2::scale_fill_manual(values = .cpg_island_colors) +
    ggplot2::labs(
      x = x_label,
      y = y_label,
      fill = "CGI position"
    ) +
    ggplot2::theme_classic(base_size = 20) +
    ggplot2::theme(
      axis.text.y = ggtext::element_markdown(size = 8 * 2 / 3),
      axis.title = ggplot2::element_text(size = 12),
      plot.subtitle = ggplot2::element_text(size = 10, face = "italic",
                                            color = "darkgrey"),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 11)
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(alpha = 1, shape = 15, size = 5)
      )
    ) +
    ggplot2::coord_flip()
  #plot$stat <- stat
  plot
}

#' Plot DMP CpG island enrichment with platform background.
#' @param mrs A MethylResultSet object.
#' @param index One or more contrast indices or names, in display order.
#' @param max_fdr Adjusted p-value threshold.
#' @param dn Number display mode: 0 none, 1 counts within segments,
#'   2 bar totals, or 3 both.
#' @return A ggplot2 object with a \code{stat} data frame containing
#'   chi-squared test details. Contrast labels are colored by methylation
#'   status with a matching legend; the platform background appears once as
#'   \code{"All probes"}.
#' @export
cpgislands_plot_dmps_fill <- function(mrs, index, q = 0.05, method = c("fdr", "BH", "bonferroni"), dn = 1) {
  .cpg_island_plot(mrs, index, "dmps", "fill", TRUE, q, method = method, dn = dn)
}

#' Plot DMP CpG island counts without platform background.
#' @param mrs A MethylResultSet object.
#' @param index One or more contrast indices or names, in display order.
#' @param max_fdr Adjusted p-value threshold.
#' @param dn Number display mode: 0 none, 1 counts within segments,
#'   2 bar totals, or 3 both.
#' @return A ggplot2 object with an empty \code{stat} data frame. Contrast
#'   labels are colored by methylation status with a matching legend.
#' @export
cpgislands_plot_dmps_stack <- function(mrs, index, q = 0.05, method = c("fdr", "BH", "bonferroni"), dn = 1) {
  .cpg_island_plot(mrs, index, "dmps", "stack", FALSE, q, method = method, dn = dn)
}

#' Plot DMR CpG island enrichment by tool with platform background.
#' @param mrs A MethylResultSet object.
#' @param index One or more contrast indices or names, in display order.
#' @param q DMR FDR threshold.
#' @param tools DMR tools to include.
#' @param combine Optional expression describing which tools should contribute
#'   to each combined DMR region, for example
#'   \code{"only(dmrcate,dmrff)"}.
#' @param min_size Minimum DMR size in base pairs.
#' @param max_size Maximum DMR size in base pairs.
#' @param min_cpgs Minimum number of probes in each input DMR.
#' @param max_cpgs Maximum number of probes in the resulting region.
#' @param separate_by_tool Logical; whether to plot each DMR tool separately
#'   (default \code{TRUE}). If \code{combine} is not \code{NULL}, this is
#'   forced to \code{FALSE} with a warning.
#' @param dn Number display mode: 0 none, 1 counts within segments,
#'   2 bar totals, or 3 both.
#' @return A ggplot2 object with a \code{stat} data frame containing
#'   chi-squared test details. Contrast labels are colored by methylation
#'   status with a matching legend; the platform background appears once as
#'   \code{"All probes"}.
#' @export
cpgislands_plot_dmrs_fill <- function(
  mrs,
  index,
  q = 0.05,
  method = c("fdr", "BH", "bonferroni"),
  tools = c("dmrcate", "ipdmr", "combp", "dmrff"),
  combine = NULL,
  min_size = 1,
  max_size = Inf,
  min_cpgs = 2,
  max_cpgs = Inf,
  dn = 1,
  separate_by_tool = TRUE) {
  .cpg_island_plot(mrs, index, "dmrs", "fill", TRUE, q, method = method, tools,
                   combine, min_size, max_size, min_cpgs, max_cpgs,
                   dn = dn, separate_by_tool = separate_by_tool)
}

#' Plot DMR CpG island counts by tool without platform background.
#' @param mrs A MethylResultSet object.
#' @param index One or more contrast indices or names, in display order.
#' @param q DMR FDR threshold.
#' @param tools DMR tools to include.
#' @param combine Optional expression describing which tools should contribute
#'   to each combined DMR region, for example
#'   \code{"only(dmrcate,dmrff)"}.
#' @param min_size Minimum DMR size in base pairs.
#' @param max_size Maximum DMR size in base pairs.
#' @param min_cpgs Minimum number of probes in each input DMR.
#' @param max_cpgs Maximum number of probes in the resulting region.
#' @param separate_by_tool Logical; whether to plot each DMR tool separately
#'   (default \code{TRUE}). If \code{combine} is not \code{NULL}, this is
#'   forced to \code{FALSE} with a warning.
#' @param dn Number display mode: 0 none, 1 counts within segments,
#'   2 bar totals, or 3 both.
#' @return A ggplot2 object with an empty \code{stat} data frame. Contrast
#'   labels are colored by methylation status with a matching legend.
#' @export
cpgislands_plot_dmrs_stack <- function(
  mrs,
  index,
  q = 0.05,
  method = c("fdr", "BH", "bonferroni"),
  tools = c("dmrcate", "ipdmr", "combp", "dmrff"),
  combine = NULL,
  min_size = 1,
  max_size = Inf,
  min_cpgs = 2,
  max_cpgs = Inf,
  dn = 1,
  separate_by_tool = TRUE) {
  .cpg_island_plot(mrs, index, "dmrs", "stack", FALSE, q, method = method, tools,
                   combine, min_size, max_size, min_cpgs, max_cpgs,
                   dn = dn, separate_by_tool = separate_by_tool)
}









#' Plot CpG island enrichment (legacy interface).
#' @param mrs A MethylResultSet object.
#' @param index One or more contrast indices or names, in display order.
#' @param what Either "dmps" or "dmrs".
#' @param position Either "fill" or "stack".
#' @param max_fdr DMR FDR threshold.
#' @param tools DMR tools to include.
#' @param separate_by_tool Logical; whether to plot each DMR tool separately
#'   when \code{what = "dmrs"} (default \code{TRUE}).
#' @param dn Number display mode: 0 none, 1 counts within segments,
#'   2 bar totals, or 3 both.
#' @return A ggplot2 object with a \code{stat} data frame containing
#'   chi-squared test details when a platform background is included.
#'   Contrast labels are colored by methylation status with a matching legend,
#'   and the platform background appears once as \code{"All probes"}.
#' @export
cpgislands_plot <- function(mrs, index, what = "dmps", position = "fill",
                            max_fdr = 0.05, method = c("fdr", "BH", "bonferroni"),
                            tools = c("dmrcate", "ipdmr", "combp", "dmrff"),
                            dn = 1, separate_by_tool = TRUE) {
  if (!what %in% c("dmps", "dmrs") || !position %in% c("fill", "stack")) {
    stop("what must be 'dmps' or 'dmrs' and position must be 'fill' or 'stack'")
  }
  .cpg_island_plot(mrs, index, what, position, position == "fill", max_fdr, method = method, tools,
                   dn = dn, separate_by_tool = separate_by_tool)
}

# Historical singular name retained for callers using the old API.
cpg_island <- cpgislands_plot
