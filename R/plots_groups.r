#' Create a Scatter Plot for Visualizing Beta Values
#'
#' This function generates a scatter plot for visualizing beta values using multidimensional scaling (MDS). It takes beta values, group information, sample names, color palette, point size, and other optional parameters for customization.
#'
#' @param betas A matrix of beta values (rows for CpG sites, columns for samples).
#' @param group A vector specifying group labels for samples.
#' @param sample_names A vector of sample names to label the points (optional).
#' @param palette A color palette for group labels.
#' @param size The size of points in the scatter plot.
#' @param showlabels Logical, whether to show sample labels or not.
#' @param ellipse Logical, whether to add ellipses around groups (optional).
#'
#' @return A scatter plot of beta values using multidimensional scaling.
#'
#' @export
my_scatterPlot <- function(
    mvals,
    group,
    sample_names,
    palette,
    size,
    showlabels,
    ellipse) {

  if (is.numeric(group)) return()
  #get lines in betas without na
  subset <- which(!is.na(rowSums(mvals)))
  if (length(subset) > 100000){subset <- sample(subset,size = 100000)}
  mds <- t(mvals[subset, ]) |> dist() |> cmdscale() |> as_tibble()
  colnames(mds) <- c("Dim.1", "Dim.2")
  if (!showlabels) sample_names <- NULL

  #K-means clustering
  #clust <- kmeans(mds, 3)$cluster |> as.factor()
  #mds <- mds |> mutate(cluster  = clust)
  mds <- mds |> mutate(grp = group)

  # Plot MDS
  ggscatter(mds, x  = "Dim.1", y  = "Dim.2",
            label  = sample_names,
            color  = "grp",
            palette  = palette,
            size  = size,
            ellipse  = ellipse,
            #ellipse.type  = "convex",
            repel  = TRUE,
            max.overlaps = 50)
}
