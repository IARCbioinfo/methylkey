#' methylkey: DNA Methylation Analysis from Illumina Arrays
#'
#' \code{methylkey} provides a comprehensive Bioconductor package for analysis of DNA methylation
#' data from Illumina methylation arrays (27k, 450k, EPIC, and mouse arrays).
#'
#' The package offers:
#' \itemize{
#'   \item Preprocessing pipelines using \code{sesame} and \code{minfi} packages
#'   \item Quality control and visualization functions
#'   \item Differential methylation analysis (DMP and DMR detection)
#'   \item Support for multiple array platforms and normalization methods
#'   \item Integration with SummarizedExperiment for data management
#' }
#'
#' Main workflow:
#' \enumerate{
#'   \item Prepare a sample sheet with sample metadata and barcode information
#'   \item Load IDAT files using \code{\link{sesame2Betas}} or \code{\link{minfi2Betas}}
#'   \item Perform quality control and visualization
#'   \item Run differential methylation analysis with \code{\link{methyldiff}}
#'   \item Generate reports with \code{\link{methylkey_report}}
#' }
#'
#' @docType package
#' @name methylkey-package
#' @aliases methylkey
#' @import assertthat
#' @import BiocGenerics
#' @import dplyr
#' @import readr
#' @import tidyr
#' @import purrr
#' @import ggplot2
#' @import stringr
#' @import rstatix
#' @import ggpubr
#' @import ggvenn
#' @import randomForest
#' @import SummarizedExperiment
#' @import sesame
#' @import sesameData
#' @import Gviz
#' @import ENmix
#' @import dmrff
#' @importFrom BiocParallel MulticoreParam
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom S4Vectors metadata metadata<- DataFrame
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom wateRmelon agep
#' @importFrom minfi read.metharray.exp mapToGenome getSex preprocessRaw fixMethOutliers getQC getBeta preprocessFunnorm detectionP estimateCellCounts densityPlot
#' @importFrom FlowSorted.Blood.EPIC estimateCellCounts2
#' @export SummarizedExperiment
#' @export DataFrame
#' @export metadata
#' @export metadata<-
NULL

#> NULL
