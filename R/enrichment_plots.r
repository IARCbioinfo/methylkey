
#' Plot CpG island enrichment for DMPs
#'
#' This function plots the distribution of differentially methylated probes
#' across CpG island annotation categories for a given methylation result set.
#' It compares the probe categories of the selected DMPs against the platform
#' background and returns a bar plot with chi-squared enrichment statistics.
#'
#' @param mrs A methylation result set object containing metadata,
#'   manifest and DMPs.
#' @param index Index or name of the DMP set to extract from \code{mrs}
#'   via \code{getDMPs}.
#' @param position Position argument passed to \code{geom_bar}
#'  ("stack" or "fill").
#' @return A ggplot2 object representing CpG island enrichment.
#'
#' @export
cpgislands_plot <- function(mrs, index, what = "dmps", position = "fill") {

  platform_name <- mrs@metadata$plateform
  cgi_levels <- c("OpenSea", "Shelf", "Shore", "Island")

  manifest <- mrs@manifest |>
    as.data.frame() |>
    dplyr::select(Probe_ID, CGIposition) |>
    dplyr::mutate(status = platform_name)

  if (what == "dmps") {
    dt <- get_dmps(mrs, index) |>
      mutate(ID = Probe_ID) |>
      mutate(status = ifelse(
        deltabetas > 0,
        "Hypermethylated",
        "Hypomethylated"
      ))
  } else if (what == "dmrs") {
    dt <- get_dmrs(mrs, index) |>
      mutate(deltabetas = mean_deltabeta) |>
      tidyr::separate_rows(CGIposition, sep = ";") |>
      mutate(status = ifelse(
        deltabetas > 0,
        "Hypermethylated",
        "Hypomethylated"
      ))
  } else {
    stop("what must be either 'dmps' or 'dmrs'")
  }

  dt <- dt |> dplyr::bind_rows(manifest)

  # 1. Préparation des données pour le test et le plot
  dt_clean <- dt |>
    dplyr::select(Probe_ID, status, CGIposition) |>
    mutate(CGIposition = ifelse(CGIposition %in% cgi_levels, CGIposition, "OpenSea")) |>
    mutate(CGIposition = factor(
      CGIposition,
      levels = cgi_levels
    ))

  # S'assurer que l'ordre sur le graphique met la Puce en premier ou en dernier
  dt_clean <- dt_clean |>
    mutate(status = factor(
      status, levels = c(platform_name, "Hypomethylated", "Hypermethylated")
    ))

  # 2. Calculs statistiques (Comparaisons ciblées face au Background)
  contingency_table <- table(dt_clean$status, dt_clean$CGIposition)

  p_text <- ""

  # Test pour les Hyper (uniquement si le groupe existe dans les données)
  if ("Hypermethylated" %in% rownames(contingency_table) &&
      platform_name %in% rownames(contingency_table)
  ) {
    mat_hyper <- contingency_table[c(platform_name, "Hypermethylated"), ]
    p_hyper <- chisq.test(mat_hyper)$p.value
    p_hyper_txt <- ifelse(p_hyper < 0.001, "p < 0.001",
      paste0("p = ", format(p_hyper, digits = 3))
    )
    p_text <- paste0("Hyper vs ", platform_name, ": ", p_hyper_txt)
  }

  # Test pour les Hypo (uniquement si le groupe existe dans les données)
  if ("Hypomethylated" %in% rownames(contingency_table) &&
      platform_name %in% rownames(contingency_table)
  ) {
    mat_hypo <- contingency_table[c(platform_name, "Hypomethylated"), ]
    p_hypo <- chisq.test(mat_hypo)$p.value
    p_hypo_txt <- ifelse(p_hypo < 0.001, "p < 0.001",
      paste0("p = ", format(p_hypo, digits = 3))
    )

    # Concatenation avec le test précédent
    sep <- ifelse(p_text == "", "", "\n")
    p_text <- paste0(p_text, sep, "Hypo vs ", platform_name, ": ", p_hypo_txt)
  }

  # 3. Génération du graphique
  dt_clean |>
    ggplot(aes(x = status, fill = CGIposition)) +
    geom_bar(position = position, width = 0.6) +
    scale_fill_manual(
      values = c("OpenSea" = "#A6CEE3",
                 "Shelf"   = "#1F78B4",
                 "Shore"   = "#B2DF8A",
                 "Island"  = "#33A02C")
    ) +
    labs(x = "Methylation status",
      y = ifelse(position == "stack", "Number of DMPs", "Percentage (%)"),
      fill = "CGI position", subtitle = paste(
        "Chi-squared enrichment test vs background:\n", p_text
      )
    ) +
    theme_classic(base_size = 20) +
    theme(
      axis.text.y = element_text(size = 8, face = "bold"),
      axis.title  = element_text(size = 12),
      plot.subtitle = element_text(
        size = 10,
        face = "italic",
        color = "darkgrey"
      ),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 11)
    ) +
    coord_flip()
}


