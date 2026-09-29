#' Calculate Delta betas between two groups
#'
#' @param betas array of betas values
#' @param design design matrix
#' @param cmtx contrast matrix
#' @param contrast_name column
#'
#' @importFrom MatrixGenerics rowMeans
#'
#' @return vector
get_delta_betas <- function(betas, design, cmtx, contrast_name) {

  score <- as.vector(design %*% cmtx[, contrast_name])

  group1 <- round(score, 6) == 1
  group0 <- round(score, 6) == -1

  group1_means <- MatrixGenerics::rowMeans(
    betas[, group1, drop = FALSE],
    na.rm = TRUE
  )
  group0_means <- MatrixGenerics::rowMeans(
    betas[, group0, drop = FALSE],
    na.rm = TRUE
  )

  (group1_means - group0_means) * 100
}


##' Calculate Delta betas between two groups
##'
##' @param betas array of betas values
##' @param design design matrix
##' @param factor_prefix prefix for factor columns
##' @param level_col column for the level of interest
##'
##' @importFrom MatrixGenerics rowMeans
##'
##' @return vector
# get_delta_betas <- function(betas, design, factor_prefix, level_col) {
#
#  factor_cols <- grep(
#    paste0("^", factor_prefix),
#     colnames(design),
#     value = TRUE
#   )

#   group1 <- design[, level_col] == 1
#   group0 <- rowSums(design[, factor_cols, drop = FALSE]) == 0

#   group1_means <- MatrixGenerics::rowMeans(
#     betas[, group1, drop = FALSE],
#     na.rm = TRUE
#   )

#   group0_means <- MatrixGenerics::rowMeans(
#     betas[, group0, drop = FALSE],
#     na.rm = TRUE
#   )

#   (group1_means - group0_means) * 100
# }

#' Calculate mvalues
#'
#' @param betas array of betas values
#'
#' @return matrix
beta2m <- function(betas) {
  log2(betas / (1 - betas))
}

#' Calculate mvalues
#'
#' @param mvalues array of M values
#'
#' @return matrix
#'
m2beta <- function(mvalues) {
  (2^mvalues) / (1 + 2^mvalues)
}


#' Filter probes from list of probes
#'
#' @param filters file containing probe list
#' @param plateform plateform
#' @return vector
#'
#' @importFrom RCurl url.exists
#' @importFrom readr read_tsv
cpg_excl <- function(
  filters = NULL,
  plateform = "IlluminaHumanMethylationEPIC",
  snp = TRUE,
  cross = TRUE,
  xy = TRUE
) {

  if (is.null(filters)) {
    # from list
    return(get_default_probeList(plateform, snp, cross, xy))
  }
  # from files
  probes <- c()
  for (file in filters) {
    if (!file.exists(file) && !RCurl::url.exists(file)) {
      stop(paste0(file, " not found"))
    }
    probes <- unique(c(probes, readr::read_tsv(file)[[1]]))
  }
  probes
}

#' Generic function to extract column data from an object
#'
#' @param x An object RGset or sdfs list.
#' @param outfile file output name
#'
#' @export
#'
setGeneric("to_geo_submission",
  function(
    x,
    outfile,
    digits = 4
  ) {
    standardGeneric("to_geo_submission")
  }
)


#' generate files to submit data to Geo database
#'
#' @param RGset an RGset object
#' @param outfile file output name
#'
#' @importFrom readr write_tsv
#' @importFrom dplyr as_tibble
#'
#' @export
setMethod(
  "to_geo_submission",
  signature(x = "RGChannelSet", outfile = "character", digits = "numeric"),
  definition = function(x, outfile = "rawdata2Geo.tsv", digits = 4) {

    if (!requireNamespace("minfi", quietly = TRUE)) {
      stop("Package 'minfi' is required for this function.",
        "Please install it.",
        call. = FALSE
      )
    }
    mset <- minfi::preprocessRaw(x)
    meth <- minfi::getMeth(mset)
    colnames(meth) <-
      paste(minfi::pData(x)$Basename, "Methylated signal")
    unmeth <- minfi::getUnmeth(mset)
    colnames(unmeth) <-
      paste(minfi::pData(x)$Basename, "Unmethylated signal")
    pval <- minfi::detectionP(x)
    colnames(pval) <- paste(minfi::pData(x)$Basename, "Detection Pval")

    msi <- cbind(unmeth, meth, pval)
    n <- ncol(pval) + 1
    m <- ncol(pval) * 2 + 1
    cols <- c(1, n, m) + rep(0:(n - 2), each = 3)
    msi <- msi[, cols]
    readr::write_tsv(dplyr::as_tibble(msi), quote = FALSE, file = outfile)

  }
)

#' generate files to submit data to Geo database
#'
#' @param Betas a Betas object
#' @param outfile file output name
#' @param digits number of decimal places to round to
#'
#' importFrom assertthat asserthat
#' importFrom readr write_tsv
#'
#' @export
setMethod(
  "to_geo_submission",
  signature(x = "list", outfile = "character", digits = "numeric"),
  definition = function(x, outfile = "betas2Geo.tsv", digits = 4) {

    if (is.null(names(x)) || any(names(x) == "")) {
      stop("`sdfs` must be a named list")
    }

    # --- betas and detection p-values, extract from sdf ---
    betas_list <- lapply(x, getBetas)
    pvals_list <- lapply(x,
      function(sdf) sesame::pOOBAH(sdf, return.pval = TRUE)
    )

    # same set of probes for everyone
    # (union, not intersection: we keep everything)
    all_probes <- Reduce(union, lapply(betas_list, names))

    betas_matrix <- sapply(betas_list, function(x) x[all_probes])
    pvals_matrix <- sapply(pvals_list, function(x) x[all_probes])
    rownames(betas_matrix) <- all_probes
    rownames(pvals_matrix) <- all_probes

    # security : same columns, same order
    pvals_matrix <- pvals_matrix[, colnames(betas_matrix), drop = FALSE]

    # --- construction dof combined format ---
    combined <- data.frame(ID_REF = all_probes, check.names = FALSE)
    for (s in colnames(betas_matrix)) {
      combined[[s]] <- round(betas_matrix[, s], digits)
      combined[[paste0(s, ".Detection Pval")]] <-
        round(pvals_matrix[, s], digits)
    }

    write.table(
      combined,
      file = outfile,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    invisible(combined)
  }
)

#' h_mean
#' Calculate the harmonic mean of a numeric vector
#' @param x A numeric vector
#' @param na.rm Logical, whether to remove NA values
#'   before calculation (default: TRUE)
#' @return The harmonic mean of the input vector
h_mean <- function(x, na.rm = TRUE) {
  # Suppress NA si asked
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)

  # To avoid division by zero
  if (any(x == 0)) return(0)

  # Calcul classique
  (1 / mean(1 / x))
}
