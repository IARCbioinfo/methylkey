#' Calculate Delta betas between two groups
#'
#' @param betas array of betas values
#' @param group group
#' @param case case
#' @param control control
#'
#' @importFrom MatrixGenerics rowMeans
#'
#' @return vector
get_delta_betas <- function(betas, group) {

  group1_means <- MatrixGenerics::rowMeans(betas[, group == 1, drop = FALSE], na.rm = TRUE)
  group0_means <- MatrixGenerics::rowMeans(betas[, group == 0, drop = FALSE], na.rm = TRUE)

  (group1_means - group0_means) * 100
}


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
    outfile
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
  signature(x = "RGChannelSet", outfile = "character"),
  definition = function(x, outfile = "rawdata2Geo.tsv") {

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
#'
#' importFrom assertthat asserthat
#' importFrom readr write_tsv
#'
#' @export
setMethod(
  "to_geo_submission",
  signature(x = "Betas", outfile = "character"),
  definition = function(x, outfile = "betas2Geo.tsv") {

    assertthat::assert_that(
      class(x) == "Betas", msg = "x must be a Betas object"
    )
    assertthat::assert_that(
      is.character(outfile), msg = "outfile must be a character string"
    )

    message("Extracting beta values and metadata...")
    betas <- methylkey::get_betas(x, masked = FALSE, na = FALSE, sex = TRUE)
    pdata <- colData(x)

    # Create output matrix with sample names as column headers
    colnames(betas) <- pdata$samples

    # Convert to tibble and add probe ID as first column
    output_df <- tibble::as_tibble(betas, rownames = "Probe_ID")

    readr::write_tsv(output_df, file = outfile)
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
