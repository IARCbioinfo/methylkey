#' methylkey: DNA Methylation Analysis from Illumina Arrays
#'
#' \code{methylkey} provides a comprehensive Bioconductor package for analysis
#'   of DNA methylation
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
#'   \item Load IDAT files using \code{\link{sesame2Betas}}
#'     or \code{\link{minfi2Betas}}
#'   \item Perform quality control and visualization
#'   \item Run differential methylation analysis with \code{\link{methyldiff}}
#'   \item Generate reports with \code{\link{methylkey_report}}
#' }
#'
#' @name methylkey-package
#' @aliases methylkey
#' @aliases MethylResultSet
#' @aliases add_combp add_dmrcate add_dmrff add_ipdmr
#' @aliases dmrtools_dmrcate dmrtools_dmrff dmrtools_ipdmr dmrtools_combp
#' @aliases as_dataframe
#' @aliases barplots_dist_to_tss best_probe_by_gene circosplot colData
#' @aliases cpgislands_plot density_plot my_scatter_plot my_violin manhattan
#' @aliases estimate_pca_corr format_sample_sheet
#' @aliases get_betas get_betas_ranges get_dmps get_dmps_ranges
#' @aliases get_dmrs get_dmrs_ranges get_genes get_manifest get_mask
#' @aliases get_model_name get_mvals get_number_of_dmrs get_plateform
#' @aliases get_results makepca minfi2betas sesame2betas
#' @aliases nbprobes nbsamples nbsentrix new_betas plot_background plot_channel
#' @aliases switch plot_channels plot_channels2 plot_detection plot_na plot
#' @aliases pca_contribution probes samples sankeydot_plot sentrix setColData
#' @aliases subset_mval to_geo_submission tools_density tools_venn variables
#' @aliases violin_plot volcano volcano_label_probes
#' @import methods
#' @importFrom dplyr .data
#' @importFrom S4Vectors metadata `metadata<-`
#' @export metadata
#' @export metadata<-
"_PACKAGE"
