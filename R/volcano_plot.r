#' Volcano Plot for DMPs/DMRs
#'
#' Publication-ready volcano plot showing log2(fold-change) vs -log10(p-value)
#'
#' @param df Data frame with DMPs/DMRs
#'   (columns: deltabetas or maxdiff, adj.P.Val or HMFDR)
#' @param fdr_threshold FDR significance threshold (default 0.05)
#' @param deltabeta_threshold Delta beta threshold for highlighting
#'   (default 0, no highlighting)
#' @param title Plot title (default: "Volcano Plot")
#' @param subtitle Plot subtitle (optional)
#' @param label_probes Optional vector of probe IDs to annotate.
#'   If NULL, labels are chosen automatically.
#' @param max_labels Maximum number of labels to display
#'   when `label_probes` is NULL.
#' @param label_column Column used for point labels.
#'   Default is `genesUniq`, use `Probe_ID` to display probe IDs.
#'
#' @return ggplot object
#'
#' @importFrom dplyr rename_with case_when mutate
#' @importFrom ggplot2 ggplot geom_point geom_vline geom_hline
#' @importFrom ggplot2 scale_color_manual labs theme element_text
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggpubr theme_pubr
#'
#' @export
volcano <- function(df, fdr_threshold = 0.05, deltabeta_threshold = 0,
                    title = "Volcano Plot", subtitle = NULL,
                    label_probes = NULL, max_labels = 10,
                    label_column = "genesUniq") {

  if (!any(c("mean_deltabeta", "deltabetas") %in% names(df))) {
    stop("df must contain either 'deltabetas' or 'mean_deltabeta'.")
  }
  if (!any(c("adj.P.Val", "HMFDR") %in% names(df))) {
    stop("df must contain either 'adj.P.Val' or 'HMFDR'.")
  }
  if (!any(c("Probe_ID", "ID") %in% names(df))) {
    stop("df must contain either 'Probe_ID' or 'ID'.")
  }

  # Standardize column names to work with both DMPs and DMRs
  df <- df |>
    dplyr::rename_with(
      ~ dplyr::case_when(
        . == "HMFDR" ~ "adj.P.Val",
        . == "mean_deltabeta" ~ "deltabetas",
        . == "ID" ~ "Probe_ID",
        TRUE ~ .
      )
    ) |>
    dplyr::mutate(
      log10_pval = -log10(.data$adj.P.Val),
      # Categorize points for coloring
      significance = dplyr::case_when(
        .data$adj.P.Val < fdr_threshold &
          abs(deltabetas) > abs(deltabeta_threshold) ~
          ifelse(deltabetas > 0, "Up-regulated", "Down-regulated"),
        TRUE ~ "Not significant"
      )
    )

  if (!is.null(label_probes) && !"Probe_ID" %in% names(df)) {
    stop("Data frame must contain a 'Probe_ID' column to label points.")
  }

  if (!label_column %in% names(df)) {
    stop(sprintf("Label column '%s' not found in data frame.", label_column))
  }

  if (is.null(label_probes) && max_labels > 0) {
    label_probes <- volcano_label_probes(
      df,
      max_fdr = fdr_threshold,
      min_delta_beta = deltabeta_threshold,
      max_labels = max_labels
    )
  }

  # Create volcano plot with ggpubr styling
  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data$deltabetas,
    y = .data$log10_pval
  )) +
    # Background points (non-significant)
    ggplot2::geom_point(
      data = dplyr::filter(df, .data$significance == "Not significant"),
      color = "#CCCCCC",
      size = 2.5,
      alpha = 0.6
    ) +
    # Significant points
    ggplot2::geom_point(
      data = dplyr::filter(
        df,
        .data$significance %in% c("Up-regulated", "Down-regulated")
      ),
      ggplot2::aes(color = .data$significance), size = 2.5, alpha = 0.8
    ) +
    # Threshold lines
    ggplot2::geom_vline(
      xintercept = c(-deltabeta_threshold, deltabeta_threshold),
      linetype = "dashed",
      color = "#888888",
      linewidth = 0.5,
      alpha = 0.7
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(fdr_threshold),
      linetype = "dashed",
      color = "#888888",
      linewidth = 0.5,
      alpha = 0.7
    ) +
    # Colors
    ggplot2::scale_color_manual(
      values = c("Up-regulated" = "#2E86AB", "Down-regulated" = "#A23B72"),
      name = "Status"
    )

  if (!is.null(label_probes) && length(label_probes) > 0) {
    label_data <- df |>
      dplyr::filter(.data$Probe_ID %in% label_probes) |>
      dplyr::mutate(label_text = as.character(.data[[label_column]])) |>
      dplyr::filter(
        !is.na(.data$deltabetas),
        !is.na(.data$log10_pval),
        !is.na(.data$label_text),
        .data$label_text != ""
      )

    if (nrow(label_data) > 0) {
      p <- p +
        ggrepel::geom_text_repel(
          data = label_data,
          ggplot2::aes(label = .data$label_text),
          size = 3.5,
          box.padding = 0.35,
          point.padding = 0.3,
          segment.color = "#4D4D4D",
          segment.size = 0.3,
          max.overlaps = Inf
        )
    }
  }

  p +
    # Labels and theme
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Delta Beta (Effect Size)",
      y = bquote("-log"[10]*"(FDR)")
    ) +
    ggpubr::theme_pubr(base_size = 11, border = TRUE) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(hjust = 0, size = 10),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold")
    )
}

#' Select probes to label on a volcano plot
#'
#' Chooses probes eligible for annotation based on
#'   FDR and delta beta thresholds,
#' en alternant entre les meilleures probes par FDR
#'   et les plus fortes par effet.
#'
#' @param df Data frame of DMPs/DMRs. Must contain columns `Probe_ID`,
#' `adj.P.Val` or `HMFDR`, and `deltabetas` or `maxdiff`.
#' @param max_fdr Maximum FDR threshold for eligible probes.
#' @param min_delta_beta Minimum absolute delta beta for eligible probes.
#' @param max_labels Maximum number of probes to select.
#'
#' @return Character vector of selected `Probe_ID`s,
#'   ordered by alternating FDR / delta beta selection.
#'
#' @importFrom dplyr rename_with case_when arrange filter pull desc
#' @export
volcano_label_probes <- function(
    df,
    max_fdr = 0.05,
    min_delta_beta = 0.2,
    max_labels = 10) {

  df <- df |>
    dplyr::rename_with(
      ~ dplyr::case_when(
        . == "HMFDR" ~ "adj.P.Val",
        . == "maxdiff" ~ "deltabetas",
        TRUE ~ .
      )
    )

  if (!"Probe_ID" %in% names(df)) {
    stop("Data frame must contain a 'Probe_ID' column.")
  }

  eligible <- df |>
    dplyr::filter(
      !is.na(.data$Probe_ID),
      !is.na(.data$adj.P.Val),
      !is.na(.data$deltabetas),
      .data$adj.P.Val <= max_fdr,
      abs(.data$deltabetas) >= min_delta_beta
    )

  if (nrow(eligible) == 0) {
    return(character(0))
  }

  fdr_list <- eligible |>
    dplyr::arrange(
      .data$adj.P.Val,
      dplyr::desc(abs(.data$deltabetas)),
      .data$Probe_ID
    ) |>
    dplyr::pull("Probe_ID")

  delta_list <- eligible |>
    dplyr::arrange(
      dplyr::desc(abs(.data$deltabetas)),
      .data$adj.P.Val,
      .data$Probe_ID
    ) |>
    dplyr::pull("Probe_ID")

  selected <- character(0)
  i <- 1
  j <- 1
  while (length(selected) < max_labels &&
           (i <= length(fdr_list) || j <= length(delta_list))) {
    if (i <= length(fdr_list)) {
      candidate <- fdr_list[i]
      if (!candidate %in% selected) {
        selected <- c(selected, candidate)
      }
      i <- i + 1
    }
    if (length(selected) >= max_labels) break

    if (j <= length(delta_list)) {
      candidate <- delta_list[j]
      if (!candidate %in% selected) {
        selected <- c(selected, candidate)
      }
      j <- j + 1
    }
  }

  selected <- unique(selected)
  selected[seq_len(min(length(selected), max_labels))]
}

#' Select one representative probe per gene
#'
#' Pour chaque gène fourni, choisit la probe la plus significative
#'   et la plus forte en effet.
#'
#' @param df Data frame of DMPs/DMRs containing probe and gene annotations.
#' @param genes Character vector of gene symbols to select.
#' @param gene_column Optional gene column name. If NULL,
#'   tries common gene annotation fields.
#'
#' @return Named character vector of selected probes,
#'   names correspond to input genes.
#'
#' @importFrom dplyr rename_with case_when select mutate filter
#' @importFrom dplyr arrange desc slice_head
#' @importFrom tidyr separate_longer_delim
#' @export
#'
best_probe_by_gene <- function(df, genes, gene_column = "genesUniq") {
  df <- df |>
    dplyr::rename_with(
      ~ dplyr::case_when(
        . == "HMFDR" ~ "adj.P.Val",
        . == "maxdiff" ~ "deltabetas",
        TRUE ~ .
      )
    )

  if (!"Probe_ID" %in% names(df)) {
    stop("Data frame must contain a 'Probe_ID' column.")
  }

  candidates <- c(
    "genesUniq", "geneNames", "gene", "Gene", "SYMBOL", "symbol",
    "Hugo_Symbol", "Gene_Symbol", "geneName"
  )
  if (is.null(gene_column)) {
    gene_column <- intersect(candidates, names(df))
    if (length(gene_column) == 0) {
      stop("Unable to find a gene annotation column in the data frame.")
    }
    gene_column <- gene_column[1]
  }

  if (!gene_column %in% names(df)) {
    stop(sprintf("Gene column '%s' not found in data frame.", gene_column))
  }

  gene_map <- df |>
    dplyr::select(
      .data$Probe_ID,
      .data$adj.P.Val,
      .data$deltabetas,
      dplyr::all_of(gene_column)
    ) |>
    dplyr::mutate(
      gene_value = as.character(.data[[gene_column]]),
      gene_value = ifelse(is.na(.data$gene_value), "", .data$gene_value)
    ) |>
    tidyr::separate_longer_delim(.data$gene_value, delim = ";") |>
    dplyr::mutate(
      gene_value = trimws(.data$gene_value),
      gene_value_lower = tolower(.data$gene_value)
    ) |>
    dplyr::filter(.data$gene_value != "")

  if (nrow(gene_map) == 0) {
    return(setNames(rep(NA_character_, length(genes)), genes))
  }

  gene_list_lower <- tolower(genes)

  probe_per_gene <- lapply(seq_along(genes), function(idx) {
    gene_lower <- gene_list_lower[[idx]]
    subset <- gene_map |>
      dplyr::filter(.data$gene_value_lower == gene_lower)

    if (nrow(subset) == 0) {
      return(NA_character_)
    }

    best_fdr <- subset |>
      dplyr::arrange(
        .data$adj.P.Val,
        dplyr::desc(abs(.data$deltabetas)),
        .data$Probe_ID
      ) |>
      dplyr::slice_head(n = 1)

    best_deltabeta <- subset |>
      dplyr::arrange(
        dplyr::desc(abs(.data$deltabetas)),
        .data$adj.P.Val,
        .data$Probe_ID
      ) |>
      dplyr::slice_head(n = 1)

    if (best_fdr$Probe_ID == best_deltabeta$Probe_ID) {
      best_fdr$Probe_ID
    } else {
      combined <- subset |>
        dplyr::mutate(
          rank_fdr = dplyr::min_rank(.data$adj.P.Val),
          rank_deltabeta = dplyr::min_rank(dplyr::desc(abs(.data$deltabetas))),
          combined_rank = .data$rank_fdr + .data$rank_deltabeta
        ) |>
        dplyr::arrange(
          .data$combined_rank,
          .data$adj.P.Val,
          dplyr::desc(abs(.data$deltabetas)),
          .data$Probe_ID
        ) |>
        dplyr::slice_head(n = 1)
      combined$Probe_ID
    }
  })

  selected <- unlist(probe_per_gene, use.names = FALSE)
  names(selected) <- genes
  selected
}