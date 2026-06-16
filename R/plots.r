#' barplots_dist_to_tss
#'
#' Display barplot of deltabetas distribution around TSS
#'
#' @param dmps data.frame of dmps
#' @param bin size of bin
#'
#' @return plot
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_hline labs
#' @importFrom ggplot2 theme_minimal scale_x_discrete theme theme_minimal element_text
#' @importFrom dplyr mutate group_by count
#'
#' @export
#'
barplots_dist_to_tss <- function(dmps, bin = 50) {

  foo <- dmps |>
    dplyr::mutate(status = ifelse(.data$deltabetas > 0, "hyper", "hypo")) |>
    dplyr::mutate(status = ifelse(.data$adj.P.Val < 0.05,
                                  paste0(.data$status, "*"),
                                  .data$status)
    ) |>
    tidyr::separate_longer_delim("distToTSS", delim = ";") |>
    dplyr::mutate(distToTSS = as.numeric(.data$distToTSS)) |>
    dplyr::mutate(distToTSS = ggplot2::cut_number(.data$distToTSS, bin))

  bounds <- stringr::str_match(levels(foo$distToTSS), "\\[?\\(*(.*),(.*)\\]")
  bounds <- data.frame(
    name = bounds[, 1],
    lower = as.numeric(bounds[, 2]),
    upper = as.numeric(bounds[, 3])
  ) |>
    dplyr::mutate(labels = dplyr::case_when(
      lower < 0 & upper > 0 ~ "TSS",
      lower < -200 & upper > -200 ~ "-200",
      lower < -500 & upper > -500 ~ "-500",
      lower < -1000 & upper > -1000 ~ "-1000",
      lower < 200 & upper > 200 ~ "200",
      lower < 500 & upper > 500 ~ "500",
      lower < 1000 & upper > 1000 ~ "1000",
      TRUE ~ ""
    ))

  moy <- table(foo$distToTSS) |> mean() / 2

  foo |>
    mutate(status = factor(
      .data$status, levels = c("hyper*", "hyper", "hypo*", "hypo")
    )) |>
    dplyr::group_by(.data$status, .data$distToTSS) |>
    dplyr::count() |>
    dplyr::mutate(n = ifelse(grepl(.data$status, "hypo"), -n, n)) |>
    ggplot2::ggplot(aes(fill = .data$status, x = .data$distToTSS, y = n)) +
    ggplot2::geom_bar(position = "stack", stat = "identity") +
    ggplot2::geom_hline(yintercept = moy, color = "red") +
    ggplot2::geom_hline(yintercept = -moy, color = "red") +
    ggplot2::labs(
      x = "distance To Tss",
      y = "Count",
      title  = "Methylation change by regions"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_discrete(labels  = bounds$labels) +
    ggplot2::theme(
      legend.text = ggplot2::element_text(hjust  = 1, size = 12),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size  = 12)
    )
}

#' Manhattan
#'
#' Manhattan plot of DMPs
#'
#' The table must contain a column 'chr' for chromosome
#' and a column pos for the chromosomique position.
#'
#' @param tab dmps table (dataframe)
#' @param sig cutoff significance
#'
#' @return plot
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous scale_color_manual labs
#' @importFrom ggplot2 theme_minimal theme element_text scale_size_continuous
#' @importFrom dplyr rename_with mutate inner_join group_by summarize
#' @importFrom forcats as_factor
#'
#' @export
manhattan <- function(df, sig = 5e-8) {

  df <- df  |>
    rename_with(
      ~ case_when(
        . == "HMFDR" ~ "adj.P.Val",
        . == "seqnames" ~ "chr",
        TRUE ~ .
      )
    ) |>
    mutate(
      pos = if ("pos" %in% names(dplyr::pick(dplyr::everything()))) {
        .data$pos
      } else {
        .data$start + (.data$width / 2)
      }
    ) |>
    mutate(
      chr = factor(
        .data$chr,
        levels = c(paste0("chr", seq(1:22)), "chrX", "chrY")
      )
    )

  data_cum <- df |>
    dplyr::filter(!is.na(.data$chr)) |>
    dplyr::group_by(.data$chr) |>
    dplyr::summarise(max_bp  = max(.data$pos)) |>
    dplyr::mutate(
      bp_add = lag(cumsum(as.numeric(.data$max_bp)), default = 0)
    ) |>
    dplyr::select(.data$chr, .data$bp_add)

  gwas_data <- df  |>
    dplyr::inner_join(data_cum, by  = "chr") |>
    dplyr::mutate(bp_cum  = sum(.data$pos, .data$bp_add))

  axis_set <- gwas_data |>
    dplyr::group_by(.data$chr) |>
    dplyr::summarize(center  = mean(.data$bp_cum))

  ylim <- min(gwas_data$adj.P.Val)
  sig <- abs(floor(log10(ylim))) + 2

  ggplot2::ggplot(
    gwas_data,
    ggplot2::aes(
      x  = .data$bp_cum,
      y  = -log10(.data$adj.P.Val),
      color = forcats::as_factor(.data$chr),
      size = -log10(.data$adj.P.Val)
    )
  ) +
    ggplot2::geom_hline(
      yintercept  = -log10(sig),
      color  = "grey40",
      linetype  = "dashed"
    ) +
    ggplot2::geom_point(alpha  = 0.75) +
    ggplot2::scale_x_continuous(
      label = axis_set$chr,
      breaks = axis_set$center
    ) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, ylim)
    ) +
    ggplot2::scale_color_manual(
      values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))
    ) +
    ggplot2::scale_size_continuous(range = c(0.5, 3)) +
    ggplot2::labs(x  = NULL, y  = "-log<sub>10</sub>(p)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position  = "none",
      panel.border  = ggplot2::element_blank(),
      panel.grid.major.x  = ggplot2::element_blank(),
      panel.grid.minor.x  = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_text(angle  = 60, size  = 8, vjust  = 0.5)
    )
}


# TODO: debug this function ----
#' Create a genome track plot with CpG, DMRs, and probe annotations
#'
#' This function generates a genome track plot that includes
#'   various annotations, such as CpG islands,
#'   DMRs (Differentially Methylated Regions),
#'   and probe annotations. It takes a data frame 'dt' containing
#'   the probe information, a matrix 'betas' containing beta values,
#'   and additional parameters for specifying the genome region to display.
#'
#' @param dt A data frame containing the probe information.
#' @param betas A matrix containing beta values.
#' @param genome The genome assembly version (e.g., "hg19").
#' @param chromosome The chromosome to display (e.g., "Chr1").
#' @param start The start position of the genome region to display.
#' @param end The end position of the genome region to display.
#'
#' @return A genome track plot displaying CpG islands, DMRs,
#'   and probe annotations.
#'
#' @importFrom Gviz IdeogramTrack AnnotationTrack DataTrack plotTracks
#' @importFrom dplyr pull group_by summarise inner_join arrange mutate select
#' @importFrom tibble as_tibble
my_track <- function(
  dt,
  sample_sheet,
  betas,
  genome = "hg19",
  chromosome = "Chr1",
  start = 0,
  end = 5000
) {

  cpgs <- dt |> dplyr::pull("Probe_ID") |> unlist()

  betas_tbl <- betas[cpgs, , drop = FALSE] |>
    tibble::as_tibble(rownames = "Probe_ID")

  foo <- dt |>
    dplyr::group_by(.data$Probe_ID, .data$chr, .data$pos) |>
    dplyr::summarise(.groups = "drop") |>
    dplyr::inner_join(betas_tbl, by = "Probe_ID") |>
    dplyr::arrange(.data$pos) |>
    dplyr::mutate(start = .data$pos, end = .data$pos) |>
    dplyr::select(-c("Probe_ID", "pos"))

  if (nrow(foo) < 1) return(NULL)

  dmrcate <- dt |>
    dplyr::group_by(.data$dmrtool, .data$chr, .data$start, .data$end) |>
    dplyr::summarise(no.cpgs = n(), .groups = "drop") |>
    dplyr::filter(.data$dmrtool == "dmrcate")
  dmrff <- dt |>
    dplyr::group_by(.data$dmrtool, .data$chr, .data$start, .data$end) |>
    dplyr::summarise(no.cpgs = n(), .groups = "drop") |>
    dplyr::filter(.data$dmrtool == "dmrff")
  combp <- dt |>
    dplyr::group_by(.data$dmrtool, .data$chr, .data$start, .data$end) |>
    dplyr::summarise(no.cpgs = n(), .groups = "drop") |>
    dplyr::filter(.data$dmrtool == "combp")
  ipdmr <- dt |>
    dplyr::group_by(.data$dmrtool, .data$chr, .data$start, .data$end) |>
    dplyr::summarise(no.cpgs = n(), .groups = "drop") |>
    dplyr::filter(.data$dmrtool == "ipdmr")

  #Ideogram track
  itrack <- Gviz::IdeogramTrack(genome  = genome, chromosome  = chromosome)

  #Annotation track, title  = "CpG"
  a_track0 <- Gviz::AnnotationTrack(
    cpgIslands,
    genome  = genome,
    name  = "CpG",
    chr  = chromosome,
    from  = start,
    to  = end
  )

  #Annotation tracks
  a_track1 <- Gviz::AnnotationTrack(
    range  = dmrcate,
    name  = "DMRcate",
    chr  = chromosome,
    from  = start,
    to  = end
  )
  a_track2 <- Gviz::AnnotationTrack(
    range  = dmrff,
    name  = "DMRff",
    chr  = chromosome,
    from  = start,
    to  = end
  )
  a_track3 <- Gviz::AnnotationTrack(
    range  = combp,
    name  = "combp",
    chr  = chromosome,
    from  = start,
    to  = end
  )
  a_track4 <- Gviz::AnnotationTrack(
    range  = ipdmr,
    name  = "ipdmr",
    chr  = chromosome,
    from  = start,
    to  = end
  )

  foo <- GenomicRanges::makeGRangesFromDataFrame(
    foo,
    keep.extra.columns = TRUE
  )

  dtrack <- Gviz::DataTrack(
    foo,
    name = "probes",
    groups = sample_sheet$Group,
    type = "confint",
    showSampleNames = TRUE,
    cex.sampleNames = 0.6,
    genome = genome,
    chr = chromosome
  )

  Gviz::plotTracks(dtrack)

  Gviz::plotTracks(
    list(
      itrack,
      gtrack,
      dtrack,
      a_track0,
      a_track1,
      a_track2,
      a_track3,
      a_track4
    ),
    chr = chromosome,
    from = start,
    to  = end
  )

}

#' Create a DMR plot for visualizing beta values with confidence intervals
#'
#' This function generates a DMR (Differentially Methylated Region)
#'   plot for visualizing beta values with confidence intervals.
#'   It takes a data frame 'dt' containing the required data,
#'   a grouping variable 'group', and a title for the plot.
#'
#' @param dt A data frame containing the required data.
#' @param group A grouping variable for distinguishing different groups
#'   in the plot.
#' @param title A title for the plot.
#'
#' @return A ggplot2 plot displaying beta values with confidence intervals
#'   for different groups.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon
#' @importFrom ggplot2 theme element_text labs ggtitle
my_dmrplot <- function(dt, group, title) {

  ggplot2::ggplot(
    dt,
    ggplot2::aes(
      x = .data$Probe_ID,
      y = .data$betas,
      ymin = (.data$betas - .data$sd / 2),
      ymax = (.data$betas + .data$sd / 2),
      group = get(group),
      fill = get(group)
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::geom_ribbon(alpha = 0.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 50, hjust = 1, size = 10)
    ) +
    ggplot2::labs(fill = group) +
    ggplot2::ggtitle(title)
}


#' TODO: debug this function ----
#' Create a Manhattan Plot for GWAS Data
#'
#' This function generates a Manhattan plot for visualizing
#' genome-wide association study (GWAS) data. It takes a data frame
#'   containing chromosome (chr), position (pos), and p-value (P.Value)
#'   columns as input.
#'
#' @param df A data frame containing GWAS data with columns chr (chromosome),
#'   pos (position), and P.Value (p-value).
#' @param sig The significance threshold for highlighting points on the plot
#'   (default is 5e-8).
#'
#' @return A Manhattan plot visualizing GWAS data.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous scale_color_manual labs
#' @importFrom ggplot2 theme_minimal theme element_text scale_size_continuous
#' @importFrom dplyr rename_with mutate inner_join group_by summarize
#' @importFrom forcats as_factor
#' @importFrom stats lag
manhattan <- function(df, sig = 5e-8) {

  df <- df  |>
    data.frame() |>
    dplyr::mutate(
      chr = factor(
        .data$chr,
        levels = c(paste0("chr", seq(1:22)), "chrX", "chrY")
      )
    )

  data_cum <- df |>
    dplyr::group_by(.data$chr) |>
    dplyr::summarise(max_bp  = max(.data$pos)) |>
    dplyr::mutate(
      bp_add = stats::lag(cumsum(as.numeric(.data$max_bp)), default = 0)
    ) |>
    dplyr::select("chr", "bp_add")

  gwas_data <- df  |>
    dplyr::inner_join(data_cum, by  = "chr") |>
    dplyr::mutate(bp_cum  = .data$pos + .data$bp_add)

  axis_set <- gwas_data |>
    dplyr::group_by(.data$chr) |>
    dplyr::summarize(center  = mean(.data$bp_cum))

  ylim <- gwas_data |>
    dplyr::filter(.data$P.Value == min(.data$P.Value)) |>
    dplyr::mutate(ylim = abs(floor(log10(.data$P.Value))) + 2) |>
    dplyr::pull(ylim)

  ggplot2::ggplot(
    gwas_data,
    ggplot2::aes(
      x = .data$bp_cum,
      y = -log10(.data$P.Value),
      color = forcats::as_factor(.data$chr), size = -log10(.data$P.Value)
    )
  ) +
    ggplot2::geom_hline(
      yintercept = -log10(sig),
      color = "grey40",
      linetype = "dashed"
    ) +
    ggplot2::geom_point(alpha  = 0.75) +
    ggplot2::scale_x_continuous(
      label = axis_set$chr,
      breaks = axis_set$center
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits  = c(0, ylim)) +
    ggplot2::scale_color_manual(
      values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))
    ) +
    ggplot2::scale_size_continuous(range  = c(0.5, 3)) +
    ggplot2::labs(
      x = NULL,
      y = "-log<sub>10</sub>(p)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      panel.border = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle  = 60, size  = 8, vjust  = 0.5)
    )
}
