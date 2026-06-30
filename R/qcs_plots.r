#' Generic function for plotting channels
#'
#' This generic function provides a common interface for plotting channels. The
#' behavior of this function depends on the class of the input data.
#'
#' @param x Input data. Depending on its class, behavior varies.
#'
#' @return A plot representing the channels.
#'
#' @export
setGeneric("plot_channels", function(x, num_positions  = 1000) {
  standardGeneric("plot_channels")
})


#' Plot channels from a list
#'
#' This method of \code{plot_channels} is specifically for plotting channels
#' when the input data is a list of channel sets. It extracts and organizes
#' the channels, creates a boxplot, and facets the plot by sentrix.
#'
#' @param x A list where each element represents a set of channels.
#'
#' @return A ggplot2 boxplot representing the channels.
#'
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_boxplot theme scale_fill_manual
#' @importFrom ggplot2 scale_y_continuous ylab facet_grid theme element_text
#'
#' @export
setMethod(
  "plot_channels",
  signature("list"),
  definition <- function(x, num_positions  = 1000) {

    samples_names <- names(x)
    if (identical(samples_names, names(x))) {
      samples_names = gsub(".*_R", "R", samples_names)
    }

    sampled_data <- map(x, ~ {
      if (num_positions < nrow(.)) {
        dplyr::sample_n(., num_positions)
      } else {
        .
      }
    })

    purrr::map_df(sampled_data, ~ tibble::tibble(
      UG = .x$UG,
      UR = .x$UR
    ), .id = "name") |>
      dplyr::mutate(
        samp = rep(samples_names, each = num_positions)
      ) |>
      tidyr::pivot_longer(
        cols = -c("name", "samp"),
        names_to = "channel",
        values_to = "Intensity"
      ) |>
      dplyr::mutate(
        sentrix = stringr::str_remove(.data$name, "_.*")
      ) |>
      ggplot2::ggplot() +
      ggplot2::geom_boxplot(aes(x = samp, y = Intensity, fill = channel)) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 90,
          vjust = 0.5, hjust = 1
        )
      ) +
      ggplot2::scale_fill_manual(values = c("#2ecc71", "#c0392b")) +
      ggplot2::scale_y_continuous(trans = "log2") +
      ggplot2::ylab("Log2 Intensity") +
      ggplot2::facet_grid(. ~ sentrix, scales  = "free_x")
  }
)


#' Plot channels from an RGChannelSet object
#'
#' This method of \code{plot_channels} is specifically for plotting channels
#' when the input data is an RGChannelSet object. It creates boxplots
#' for red and green channels.
#'
#' @param x An RGChannelSet object containing microarray channel data.
#'
#' @return A ggplot2 plot with separate boxplots for red and green channels.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom graphics par axis
#'
#' @export
setMethod(
  "plot_channels",
  signature("RGChannelSet"),
  definition = function(x) {

    if (!requireNamespace("minfi", quietly = TRUE)) {
      stop("Package 'minfi' is required for this function.",
        "Please install it.",
        call. = FALSE
      )
    }
    sample_sheet = SummarizedExperiment::colData(x)
    nsamp = ncol(x)
    ylab <- "log2 intensity of both green and red channel"
    par(xaxt = "n")

    boxplot(log2(minfi::getRed(x) + 1),
      col = "red",
      boxwex  = 0.25,
      at = 1:nsamp - 0.175,
      ylab = ylab,
      labels = sample_sheet$samples, cex = 0.5
    )

    boxplot(log2(minfi::getGreen(x) + 1),
      col = "green",
      boxwex  = 0.25,
      at = 1:nsamp + 0.175,
      axis = FALSE,
      add = TRUE,
      cex = 0.5
    )

    graphics::par(xaxt = "s")
    graphics::axis(1, at = 1:nsamp,
      labels = SummarizedExperiment::colData(x)$Basename,
      tick = TRUE,
      las = 2,
      cex.axis = 0.8
    )
    recordPlot()
  }
)



#' Generic function for plotting channels
#'
#' This generic function provides a common interface for plotting channels. The
#' behavior of this function depends on the class of the input data.
#'
#' @param x Input data. Depending on its class, behavior varies.
#'
#' @return A plot representing the channels.
#'
#' @export
setGeneric(
  "plot_channels2",
  function(x) standardGeneric("plot_channels2")
)

#' Call Plot channels from an RGChannelSet object for each sentrix
#'
#' This method of \code{plot_channels} is specifically for plotting channels
#' when the input data is an RGChannelSet object. It creates separate boxplots
#' for each sentrix
#'
#' @param x An RGChannelSet object containing microarray channel data.
#'
#' @return A list of boxplots.
#'
#' @export
setMethod(
  "plot_channels2",
  signature("RGChannelSet"),
  definition  = function(x) {

    nsubplots <- length(sentrix)
    stx <- gsub("_.*", "", colnames(x)) |> unique()

    lapply(1:nsubplots, function(i) {
      subsel <- grepl(stx[i], colnames(x))
      substx <- x[, subsel]
      plot_channels(substx)
    })

  }
)

#' Create a scatter plot with custom annotations
#'
#' This function generates a scatter plot with custom
#' annotations using ggplot2. The plot displays data from a data frame,
#' where 'mean_oob_grn' and 'mean_oob_red' are used as x and y coordinates,
#' and 'name' is used for annotations.
#'
#' @param df A data frame containing the data to be plotted,
#' with columns 'mean_oob_grn', 'mean_oob_red', and 'name'.
#'
#' @return A ggplot2 scatter plot with custom annotations.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_text geom_abline xlab ylab
#'
#' @export
plot_background <- function(df) {
  ggplot2::ggplot(df,
    ggplot2::aes(
      x = .data$mean_oob_grn,
      y = .data$mean_oob_red,
      label = .data$name
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_text(hjust = -0.1, vjust = 0.1) +
    ggplot2::geom_abline(intercept = 0, slope  = 1, linetype = "dotted") +
    ggplot2::xlab("Green Background") +
    ggplot2::ylab("Red Background")
}

#' Create a bar plot to visualize detection fractions
#'
#' This function generates a bar plot to visualize detection
#' fractions using ggplot2. It plots data from a data frame,
#' where 'name' is used on the x-axis and 'frac_dt' on the y-axis.
#'
#' @param df A data frame containing the data to be plotted,
#' with columns 'name' and 'frac_dt'.
#'
#' @return A ggplot2 bar plot for visualizing detection fractions.
#'
#' @importFrom ggplot2 ggplot geom_bar aes theme element_text
#'
#' @export
plot_detection <- function(df) {

  ggplot2::ggplot(df) +
    ggplot2::geom_bar(
      ggplot2::aes(
        x = .data$name,
        y = .data$frac_dt
      ),
      stat = "identity"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
}

#' Create a bar plot to visualize missing values (NA)
#'
#' This function generates a bar plot to visualize the number of
#' missing values (NA) for different categories using ggplot2.
#' It plots data from a data frame, where 'name' is used on the x-axis
#' and 'num_na_cg' on the y-axis.
#'
#' @param df A data frame containing the data to be plotted,
#' with columns 'name' and 'num_na_cg'.
#'
#' @return A ggplot2 bar plot for visualizing missing values (NA).
#'
#' @importFrom ggplot2 ggplot geom_bar aes theme element_text
#'
#' @export
plot_na <- function(df) {
  ggplot2::ggplot(df) +
    ggplot2::geom_bar(
      ggplot2::aes(x = .data$name, y = .data$num_na_cg),
      stat = "identity"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
}

#' densityPlot of betas values
#'
#' This function generates a ggplot2 density plot for matrix columns.
#' Each column is represented by a density line and colored by
#'   \code{sampGroups}.
#'
#' @param dat A matrix of betas.
#' @param sampGroups A vector of group labels for each column in \code{dat}.
#' @param main Plot title.
#' @param xlab Label for the x-axis.
#' @param pal Color palette vector.
#' @param xlim x-axis limits.
#' @param ylim y-axis limits.
#' @param legend Whether to draw a legend.
#' @return A ggplot2 object representing density estimates for each sample.
#'
#' @importFrom ggplot2 ggplot aes geom_line labs scale_color_manual
#' @importFrom dplyr bind_rows
#' @importFrom stats density
#'
#' @export
density_plot <- function(
    meth, group = "sentrix_id", main  = "",
    subset = 10000,
    xlab  = "Beta values distribution",
    pal  = RColorBrewer::brewer.pal(8, "Dark2"), ...) {

  if (!is(meth, "Betas")) {
    stop("argument 'meth' must be a Betas object'")
  }

  if (group %in% colnames(colData(meth))) {
    group <- colData(meth)[[group]]
  } else {
    warning(paste0("Group variable '", group, "' not found in colData.",
      "Using default group labels."
    ))
    group <- colData(meth)["sentrix_id"]
  }

  betas <- get_betas(meth)
  sset <- sample(nrow(betas), size = subset)
  betas <- betas[sset, ]

  d <- apply(betas, 2, function(x) stats::density(as.vector(x), na.rm  = TRUE))

  df <- dplyr::bind_rows(lapply(seq_along(d), function(i) {
    data.frame(
      x = d[[i]]$x,
      y = d[[i]]$y,
      sample_group = group[i],
      sample = if (!is.null(colnames(betas)))
        colnames(betas)[i] else as.character(i),
      stringsAsFactors = FALSE
    )
  }))

  ggplot2::ggplot(df, ggplot2::aes(
    x = .data$x,
    y = .data$y,
    color = .data$sample_group,
    group = interaction(.data$sample, .data$sample_group)
  )) +
    ggplot2::geom_line() +
    ggplot2::labs(x = xlab, y = "Density", title = main, color = "Groups") +
    ggplot2::scale_color_manual(
      values = rep(pal, length.out = length(levels(group)))
    )
}

#' violin_plot1
#'
#' Display betas in a violin plot by group
#'
#' @param betas matrix of betas
#' @param group variable group
#' @param numPositions number of random CpG used
#'
#' @return plot
#'
#' @importFrom ggplot2 ggplot aes geom_violin theme element_text guides ylab
#' @importFrom stats runif
#' @importFrom tidyr pivot_longer
#'
#' @export
#'
violin_plot <- function(betas, group, num_positions = 1000) {

  cpg <- betas[stats::runif(num_positions, 1, nrow(betas)), ]
  cpg <- data.frame(group, t(cpg)) |>
    tidyr::pivot_longer(!group, names_to  = "Probe_ID", values_to  = "betas")

  ggplot2::ggplot(cpg, ggplot2::aes(x = group, y = betas, fill = group)) +
    ggplot2::geom_violin() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 50, hjust = 1, size = 14)
    ) +
    ggplot2::guides(fill = "none") +
    ggplot2::ylab("Methylation %") + xlab("")

}

#' Create a scatter plot to visualize channel switching
#'
#' This function generates a scatter plot to visualize channel
#'   switching using ggplot2. It plots data from a data frame,
#'   where 'InfI_switch_G2R' and 'InfI_switch_R2G' are used as
#'   x and y coordinates.
#'
#' @param df A data frame containing the data to be plotted,
#'   with columns 'InfI_switch_G2R' and 'InfI_switch_R2G'.
#'
#' @return A ggplot2 scatter plot for visualizing channel switching.
#'
#'
#' @export
plot_channel_switch <- function(df) {
  ggplot(df) +
    geom_point(aes(.data$InfI_switch_G2R, .data$InfI_switch_R2G))
}