#' CpG with too many NA
#'
#' select probes with percentage of missing values superior to CpGlimit
#'
#' @param betas matrix of betas
#' @param nalimit maximum proportion of NA accepted
#'
#' @importFrom dplyr filter
#'
#' @return List of probes to exclude
#'
cpg_na_excl <- function(betas, nalimit = 0.2) {

  na_row <- apply(betas, 1, function(i) sum(is.na(i)))
  rownames(betas)[na_row > nalimit * ncol(betas)]
}

#' impute missing values
#'
#' impute missing values with pamr for probes with less than 'nalimit'
#' proportion of NA values
#'
#' @param betas matrix of betas
#' @param nalimit maximum proportion of NA accepted
#'
#' @import pamr
#'
#' @return betas matrix
impute_na <- function(betas, nalimit = 0.2) {

  if (!requireNamespace("pamr", quietly = TRUE)) {
    stop(
      "Package 'pamr' is required for imputeNA(). ",
      "Please install it with install.packages('pamr').",
      call. = FALSE
    )
  }

  foo <- list(x = as.matrix(betas), y = colnames(betas))
  betas <- pamr::pamr.knnimpute(foo, rowmax = nalimit, colmax = 0.95)$x
  betas
}

#' replace_by_mean
#'
#' Replace missing values by mean for each group
#'
#' @param betas matrix of betas
#' @param group group name
#'
#' @return betas matrix
replace_by_mean <- function(betas, groups) {

  # select rows containing missing values
  prob_na_count <- which(apply(betas, 1, function(i) sum(is.na(i))) > 0)

  # calculate means by group
  mean_betas <- function(betas) {
    rmeans <- matrix(
      rowMeans(betas, na.rm = TRUE),
      ncol = ncol(betas),
      nrow = nrow(betas)
    )
    betas[is.na(betas)] <- rmeans[is.na(betas)]
    betas
  }

  # replace missing values by means
  for (group in levels(groups)) {
    sel <- groups == group
    if (!sum(sel) > 1) {
      stop("You should have at least 2 samples by group !")
    }
    betas[prob_na_count, sel] <- mean_betas(betas[prob_na_count, sel])
  }

  # if there is remaining na, because there is no value for a group
  prob_na_count <- which(apply(betas, 1, function(i) sum(is.na(i))) > 0)
  if (length(prob_na_count)) {
    betas <- betas[-prob_na_count, ]
  }

  betas
}