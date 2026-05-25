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
#setGeneric("plot_Channels", function(x) standardGeneric("plot_Channels") )
setGeneric("plot_Channels", function(x, numPositions  = 1000) {
  standardGeneric("plot_Channels")
})


#' Plot channels from a list
#'
#' This method of \code{plot_Channels} is specifically for plotting channels
#' when the input data is a list of channel sets. It extracts and organizes
#' the channels, creates a boxplot, and facets the plot by sentrix.
#'
#' @param x A list where each element represents a set of channels.
#'
#' @return A ggplot2 boxplot representing the channels.
#'
#' @export
setMethod(
  "plot_Channels",
  signature("list"),
  definition = function(x, numPositions  = 1000) {

    sampNames = names(x)
    if (identical(sampNames, names(x))) sampNames = gsub(".*_R", "R", sampNames)

    sampled_data <- map(x, ~ {
      if (numPositions < nrow(.)) {
        sample_n(., numPositions)
      } else {
        .
      }
  })

  bind_rows(
    map(sampled_data, ~ tibble(UG  = .$UG, UR  = .$UR)), .id  = "name") |>
    mutate(samp  = rep(sampNames, each  = numPositions)) |>
    tidyr::pivot_longer(-c(name, samp),
      names_to  = "channel",
      values_to  = "Intensity") |>
    mutate(sentrix  = gsub("_.*", "", name)) |>
    ggplot() +
    geom_boxplot(aes(x  = samp, y  = Intensity, fill  = channel)) +
    theme(axis.text.x  = element_text(angle  = 90, vjust  = 0.5, hjust  = 1)) +
    scale_fill_manual(values  = c("#2ecc71", "#c0392b")) +
    scale_y_continuous(trans  = 'log2') +
    ylab("Log2 Intensity") +
    facet_grid(. ~ sentrix, scales  = "free_x")
})


#' Plot channels from an RGChannelSet object
#'
#' This method of \code{plot_Channels} is specifically for plotting channels
#' when the input data is an RGChannelSet object. It creates boxplots
#' for red and green channels.
#'
#' @param x An RGChannelSet object containing microarray channel data.
#'
#' @return A ggplot2 plot with separate boxplots for red and green channels.
#'
#' @export
setMethod(
  "plot_Channels",
  signature("RGChannelSet"),
  definition  = function(x) {
    nsamp = ncol(x)
    ylab <- "log2 intensity of both green and red channel"
    par(xaxt = 'n')

    boxplot(log2(getRed(x)+1),
      col  = "red",
      boxwex  = 0.25,
      at = 1:nsamp - 0.175,
      ylab = ylab,
      labels = sampleSheet$samples, cex = 0.5)

    boxplot(log2(getGreen(x)+1),
      col  = "green",
      boxwex  = 0.25,
      at = 1:nsamp + 0.175,
      axis = F,
      add = T,
      cex = 0.5)

    par(xaxt = 's')
    axis(1, at = 1:nsamp,
      labels = SummarizedExperiment::colData(x)$Basename,
      tick = TRUE,
      las = 2,
      cex.axis = 0.8)
    recordPlot()
})



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
  "plot_Channels2",
  function(x) standardGeneric("plot_Channels2"))

#' Call Plot channels from an RGChannelSet object for each sentrix
#'
#' This method of \code{plot_Channels} is specifically for plotting channels
#' when the input data is an RGChannelSet object. It creates separate boxplots
#' for each sentrix
#'
#' @param x An RGChannelSet object containing microarray channel data.
#'
#' @return A list of boxplots.
#'
#' @export
setMethod(
  "plot_Channels2",
  signature("RGChannelSet"),
  definition  = function(x) {

  nsubplots <- length(sentrix)
  stx <- gsub("_.*", "", colnames(x)) |> unique()

  lapply(1:nsubplots, function(i) {
    subsel <- grepl(stx[i], colnames(x))
    substx <- x[, subsel]
    plot_Channels(substx)
  })

})

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
#' @export
plot_background <- function(df) {
  ggplot(df,
         aes(x  = mean_oob_grn, y = mean_oob_red, label  = barcode)) +
    geom_point() + geom_text(hjust  = -0.1, vjust  = 0.1) +
    geom_abline(intercept  = 0, slope  = 1, linetype  = 'dotted') +
    xlab('Green Background') + ylab('Red Background')
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
#' @export
plot_Detection <- function(df) {

  ggplot(df) + geom_bar(aes(x = barcode,y = frac_dt), stat = "identity") +
    theme(axis.text.x  = element_text(angle  = 90, vjust  = 0.5, hjust = 1))
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
#' @export
plot_NA <- function(df){
  ggplot(df) + geom_bar(aes(x = barcode,y = num_na_cg), stat = "identity") +
    theme(axis.text.x  = element_text(angle  = 90, vjust  = 0.5, hjust = 1))
}


#' densityPlot of betas values
#'
#' This function generates a ggplot2 density plot for matrix columns.
#' Each column is represented by a density line and colored by \code{sampGroups}.
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
#' @export
#'
density_plot <- function(meth, group = "sentrix_id", main  = "",
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

  betas <- getBetas(meth)
  sset <- sample(nrow(betas), size = subset)
  betas <- betas[sset, ]

  d <- apply(betas, 2, function(x) density(as.vector(x), na.rm  = TRUE))
  #if (missing(ylim)) ylim <- range(sapply(d, function(i) range(i$y)))
  #if (missing(xlim)) xlim <- range(sapply(d, function(i) range(i$x)))

  df <- bind_rows(lapply(seq_along(d), function(i) {
    data.frame(
      x = d[[i]]$x,
      y = d[[i]]$y,
      sampGroup = group[i],
      sample = if (!is.null(colnames(betas)))
        colnames(betas)[i] else as.character(i),
      stringsAsFactors = FALSE
    )
  }))

  p <- ggplot(df, aes(x = x,
                      y = y,
                      color = sampGroup,
                      group = interaction(sample, sampGroup))) +
    geom_line() +
    labs(x = xlab, y = "Density", title = main, color = "Groups") +
    scale_color_manual(values = rep(pal, length.out = length(levels(group))))

  p
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
#' @export
#'
violin_plot <- function(betas, group, numPositions = 1000) {

  cpg <- betas[runif(numPositions,1,nrow(betas)),]
  #cpg <- cbind(group = as.character(group),data.table(t(cpg)) )
  #cpg <- melt( cpg, id.vars = "group" )
  #group <- as.character(group)
  cpg <- data.frame(group, t(cpg)) |>
    pivot_longer(!group, names_to  = "Probe_ID", values_to  = "betas")

  #colnames(cpg)<-c("group","Probe_ID","betas")

  p1 <- ggplot(cpg, aes(x = group, y = betas, fill = group)) + geom_violin() +
    theme(axis.text.x = element_text(angle = 50,hjust = 1, size = 14)) +
    guides(fill  = FALSE) +
    ylab("Methylation %") + xlab("")

  return(p1)
}




#' Create a scatter plot to visualize channel switching
#'
#' This function generates a scatter plot to visualize channel switching using ggplot2. It plots data from a data frame, where 'InfI_switch_G2R' and 'InfI_switch_R2G' are used as x and y coordinates.
#'
#' @param df A data frame containing the data to be plotted, with columns 'InfI_switch_G2R' and 'InfI_switch_R2G'.
#'
#' @return A ggplot2 scatter plot for visualizing channel switching.
#'
#' @export
plot_channel_switch <- function(df){
  ggplot(df) + geom_point(aes(InfI_switch_G2R, InfI_switch_R2G))
}