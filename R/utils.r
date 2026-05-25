#' getDeltaBetas
#'
#' calculate Delta betas between two groups
#'
#' @param betas array of betas values
#' @param group group
#' @param case case
#' @param control control
#'
#' @return vector
#'
getDeltaBetas <- function(betas, group) {
  (rowMeans(betas[, group == 1, drop = FALSE]) -
    rowMeans(betas[, group == 0, drop = FALSE])) * 100
}


#' beta2m
#'
#' calculate mvalues
#'
#' @param betas array of betas values
#'
#' @return matrix
#'
beta2m <- function(betas) {
  return(log2(betas / (1 - betas)))
}

#' beta2m
#'
#' calculate mvalues
#'
#' @param mvalues array of M values
#'
#' @return matrix
#'
m2beta <- function(mvalues) {
  return(2^mvalues / (1 + 2^mvalues))
}


#' CpGexcl
#'
#' Filter probes from list of probes
#'
#' @param filters file containing probe list
#' @param plateform plateform
#'
#' @return vector
#' 
CpGexcl <- function(
  filters = NULL,
  plateform = "IlluminaHumanMethylationEPIC",
  SNP = TRUE,
  CROSS = TRUE,
  XY = TRUE
) {

  probes <- c()
  if (is.null(filters)) {
    # from list
    probes <- getDefaultProbeList(plateform, SNP, CROSS, XY)
  } else {
    # from files
    for (file in filters) {
      if (!file.exists(file) & !RCurl::url.exists(file)) {
        stop(paste0(file, " not found"))
      }
      probes <- unique(c(probes, readr::read_tsv(file)[[1]]))
    }
  }

  return(probes)
}

#' Generic function to extract column data from an object
#'
#' @param x An object RGset or sdfs list.
#' @param outfile file output name
#'
#' @export
#' 
setGeneric("toGeoSubmission", function(x, outfile) standardGeneric("toGeoSubmission"))


#' generate files to submit data to Geo database
#'
#' @param RGset an RGset object
#' @param outfile file output name
#'
#' @export
#'
setMethod(
  "toGeoSubmission",
  signature(x = "RGChannelSet", outfile = "character"),
  definition = function(x, outfile = "rawdata2Geo.tsv") {

    Mset <- minfi::preprocessRaw(RGset)
    meth <- minfi::getMeth(Mset)
    colnames(meth) <- paste(minfi::pData(RGset)$Basename, "Methylated signal")
    unmeth <- minfi::getUnmeth(Mset)
    colnames(unmeth) <- paste(minfi::pData(RGset)$Basename, "Unmethylated signal")
    pval <- minfi::detectionP(RGset)
    colnames(pval) <- paste(minfi::pData(RGset)$Basename, "Detection Pval")

    MSI <- cbind(unmeth, meth, pval)
    n <- ncol(pval) + 1
    m <- ncol(pval) * 2 + 1
    cols <- c(1, n, m) + rep(0:(n - 2), each = 3)
    MSI <- MSI[, cols]
    readr::write_tsv(as_tibble(MSI), quote = FALSE, file = outfile)

  }
)


#' generate files to submit data to Geo database
#'
#' @param Betas a Betas object
#' @param outfile file output name
#'
#' @export
#'
setMethod(
  "toGeoSubmission",
  signature(x = "Betas", outfile = "character"),
  definition = function(x, outfile = "betas2Geo.tsv") {

    assertthat::assert_that(class(x) == "Betas", msg = "x must be a Betas object")
    assertthat::assert_that(is.character(outfile), msg = "outfile must be a character string")

    message("Extracting beta values and metadata...")
    betas <- methylkey::getBetas(x, masked = FALSE, na = FALSE, sex = TRUE)
    pdata <- colData(x)

    # Create output matrix with sample names as column headers
    colnames(betas) <- pdata$samples

    # Convert to tibble and add probe ID as first column
    output_df <- tibble::as_tibble(betas, rownames = "Probe_ID")

    message(sprintf(
      "Writing %d probes x %d samples to %s",
      nrow(betas),
      ncol(betas),
      outfile
    ))
    readr::write_tsv(output_df, file = outfile)

    message("GEO submission file created successfully.")
  }
)

#' h_mean
#' Calculate the harmonic mean of a numeric vector
#' @param x A numeric vector
#' @param na.rm Logical, whether to remove NA values
#'   before calculation (default: TRUE)
#' @return The harmonic mean of the input vector
#'
#' @export
h_mean <- function(x, na.rm = TRUE) {
  # Suppress NA si asked
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)

  # To avoid division by zero
  if (any(x == 0)) return(0)

  # Calcul classique
  return(1 / mean(1 / x))
}
