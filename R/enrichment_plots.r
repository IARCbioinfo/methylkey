
#' Plot CpG island enrichment for DMPs
#'
#' This function plots the distribution of differentially methylated probes
#' across CpG island annotation categories for a given methylation result set.
#' It compares the probe categories of the selected DMPs against the platform
#' background and returns a bar plot with chi-squared enrichment statistics.
#'
#' @param mrs A methylation result set object
#'   containing metadata, manifest and DMPs.
#' @param index Index or name of the DMP set
#'   to extract from \code{mrs} via \code{getDMPs}.
#' @param what Either "dmps" or "dmrs" (définit le type d'analyse)
#' @param position Position argument passed to \code{geom_bar}
#'   ("stack" or "fill").
#' @return A ggplot2 object representing CpG island enrichment.
#'
#' @importFrom dplyr select mutate bind_rows
#' @importFrom tidyr separate_longer_delim
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual
#' @importFrom ggplot2 labs theme_classic theme coord_flip element_text
#' @importFrom stats chisq.test
#'
#' @export
cpgislands_plot <- function(mrs, index, what = "dmps", position = "fill") {

  platform_name <- mrs@metadata$plateform
  cgi_levels <- c("OpenSea", "Shelf", "Shore", "Island")

  manifest <- mrs@manifest |>
    as.data.frame() |>
    dplyr::select(.data$Probe_ID, .data$Relation_to_Island) |>
    dplyr::mutate(status = platform_name)

  if (what == "dmps") {
    dt <- get_dmps(mrs, index) |>
      dplyr::mutate(ID = .data$Probe_ID) |>
      dplyr::mutate(status = ifelse(
        .data$deltabetas > 0,
        "Hypermethylated",
        "Hypomethylated"
      ))
  } else if (what == "dmrs") {
    dt <- get_dmrs(mrs, index) |>
      dplyr::mutate(deltabetas = .data$mean_deltabeta) |>
      tidyr::separate_longer_delim(.data$Relation_to_Island, delim = ";") |>
      dplyr::mutate(status = ifelse(
        .data$deltabetas > 0,
        "Hypermethylated",
        "Hypomethylated"
      ))
  } else {
    stop("what must be either 'dmps' or 'dmrs'")
  }

  dt <- dt |> dplyr::bind_rows(manifest)

  # 1. Prepare data for the test and the plot
  dt_clean <- dt |>
    dplyr::select(.data$Probe_ID, .data$status, .data$Relation_to_Island) |>
    dplyr::mutate(Relation_to_Island = ifelse(
      .data$Relation_to_Island %in% cgi_levels,
      .data$Relation_to_Island, "OpenSea"
    )) |>
    dplyr::mutate(Relation_to_Island = factor(
      .data$Relation_to_Island,
      levels = cgi_levels
    ))

  dt_clean <- dt_clean |>
    dplyr::mutate(status = factor(
      .data$status,
      levels = c(platform_name, "Hypomethylated", "Hypermethylated")
    ))

  # 2. Chi-squared test for enrichment of DMPs in CGI
  #   categories compared to platform background
  contingency_table <- table(dt_clean$status, dt_clean$Relation_to_Island)

  p_text <- ""

  # Test for hyper (only if the group exists in the data)
  if ("Hypermethylated" %in% rownames(contingency_table) &&
      platform_name %in% rownames(contingency_table)
  ) {
    mat_hyper <- contingency_table[c(platform_name, "Hypermethylated"), ]
    p_hyper <- stats::chisq.test(mat_hyper)$p.value
    p_hyper_txt <- ifelse(p_hyper < 0.001, "p < 0.001",
      paste0("p = ", format(p_hyper, digits = 3))
    )
    p_text <- paste0("Hyper vs ", platform_name, ": ", p_hyper_txt)
  }

  # Test for hypo (only if the group exists in the data)
  if ("Hypomethylated" %in% rownames(contingency_table) &&
      platform_name %in% rownames(contingency_table)
  ) {
    mat_hypo <- contingency_table[c(platform_name, "Hypomethylated"), ]
    p_hypo <- stats::chisq.test(mat_hypo)$p.value
    p_hypo_txt <- ifelse(p_hypo < 0.001, "p < 0.001",
      paste0("p = ", format(p_hypo, digits = 3))
    )

    # Concatenation of p-values for the plot subtitle
    sep <- ifelse(p_text == "", "", "\n")
    p_text <- paste0(p_text, sep, "Hypo vs ", platform_name, ": ", p_hypo_txt)
  }

  # 3. GGeneration of the plot
  dt_clean |>
    ggplot2::ggplot(
      ggplot2::aes(x = .data$status, fill = .data$Relation_to_Island)
    ) +
    ggplot2::geom_bar(position = position, width = 0.6) +
    ggplot2::scale_fill_manual(
      values = c("OpenSea" = "#A6CEE3",
                 "Shelf"   = "#1F78B4",
                 "Shore"   = "#B2DF8A",
                 "Island"  = "#33A02C")
    ) +
    ggplot2::labs(x = "Methylation status",
      y = ifelse(position == "stack", "Number of DMPs", "Percentage (%)"),
      fill = "CGI position", subtitle = paste(
        "Chi-squared enrichment test vs background:\n", p_text
      )
    ) +
    ggplot2::theme_classic(base_size = 20) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8, face = "bold"),
      axis.title  = ggplot2::element_text(size = 12),
      plot.subtitle = ggplot2::element_text(
        size = 10,
        face = "italic",
        color = "darkgrey"
      ),
      legend.title = ggplot2::element_text(size = 12),
      legend.text  = ggplot2::element_text(size = 11)
    ) +
    ggplot2::coord_flip()
}
