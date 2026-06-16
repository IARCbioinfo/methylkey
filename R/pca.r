#' PCA analysis
#'
#' Estimate correlation covariates vs. PCA axes
#'
#' @param mvals matrix of m values
#' @param pdata sampleSheet of data
#'
#' @importFrom stats prcomp
#'
#' @return betas matrix
#'
#' @export
#'
makepca <- function(mvals, pdata) {

  tmvals <- t(mvals) # n x p required for prcomp
  sel <- which(apply(tmvals, 2, var) == 0)

  if (length(sel)) {
    pca <- stats::prcomp(tmvals[, -sel], scale = TRUE, center = TRUE)
  } else {
    pca <- stats::prcomp(tmvals, scale = TRUE, center = TRUE)
  }
  rownames(pca$x) <- pdata$samples

  pca
}



#' Plot PCA Contribution
#'
#' Visualizes the contribution of principal components (PC) in a PCA analysis.
#'
#' @param pca An object of class 'prcomp' containing PCA results.
#' @param nPC Number of principal components to include in the plot.
#' @param title Title for the plot.
#'
#' @importFrom ggplot2 ggplot geom_col labs theme element_text
#'
#' @return A ggplot object representing the PCA contribution plot.
#'
#' @export
plot_pca_contribution <- function(pca, npc = 4, title = "") {
  v_var <- summary(pca)$importance[2, 1:npc]
  ggplot(data.frame(pca_axis = names(v_var), val = v_var)) +
    geom_col(aes(x = reorder(.data$pca_axis, 1:npc), y = .data$val)) +
    labs(
      title = title,
      subtitle = "PCA contribution",
      x = "PC",
      y = "Proportion expl. var."
    ) +
    theme(plot.title = element_text(lineheight = .8, face = "bold"))
}



#' Estimate PCA Correlation
#'
#' Computes correlation statistics for the contribution of variables
#'   to principal components (PC) in a PCA analysis.
#'
#' @param pca An object of class 'prcomp' containing PCA results.
#' @param pdata A data frame containing variables for which
#'   the correlation with PCs will be estimated.
#' @param nPC Number of principal components to include in the analysis.
#'
#' @return A matrix of p-values representing the correlation
#'   significance for each variable with each PC.
#'
#' @importFrom dplyr group_by pull
#' @importFrom rstatix wilcox_test kruskal_test cor_test
#' @importFrom tidyr gather
#' @importFrom tibble rownames_to_column
#' @importFrom purrr map
#' @importFrom stats reformulate
#'
#' @export
estimate_pca_corr <- function(pca, pdata, npc) {

  npc <- min(npc, nrow(pdata))

  tab <- data.frame(pca$x[, 1:npc]) |>
    tibble::rownames_to_column("samples") |>
    tidyr::gather("PC", "contrib", -"samples") |>
    merge(pdata, by = "samples")
  nb_lvl <- tab |>
    purrr::map(levels) |>
    sapply(length)
  is_num <- tab |>
    dplyr::select(-"contrib") |>
    purrr::map(is.numeric) |>
    unlist()
  # Wilcox_test for categorical variables with 2 levels
  wilcox <- sapply(names(which(nb_lvl == 2)), function(variable) {
    tab |>
      dplyr::group_by(.data$PC) |>
      rstatix::wilcox_test(stats::reformulate(variable, "contrib")) |>
      dplyr::pull("p")
  })
  # Kruskal test for categorical variables with more than 2 levels
  kruskal <- sapply(names(which(nb_lvl > 2)), function(variable) {
    tab |>
      dplyr::group_by(.data$PC) |>
      rstatix::kruskal_test(stats::reformulate(variable, "contrib")) |>
      dplyr::pull("p")
  })
  # cor.test (kendall) for numeric variables
  kendall <- sapply(names(which(is_num == TRUE)), function(variable) {
    tab |>
      dplyr::group_by(.data$PC) |>
      rstatix::cor_test(!!variable, "contrib", method = "kendall") |>
      dplyr::pull("p")
  })
  dt_pval <- cbind(wilcox, kruskal, kendall) |> t()
  colnames(dt_pval) <- paste0("PC", 1:npc)
  dt_pval
}
