% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmrtools.r
\name{searchDMR_dmrcate}
\alias{searchDMR_dmrcate}
\title{Search for Differentially Methylated Regions (DMRs) using DMRcate}
\usage{
searchDMR_dmrcate(
  dmps,
  fdr = 0.05,
  maxgap = 1000,
  pcutoff = 0.2,
  genome = "hg38"
)
}
\arguments{
\item{dmps}{A data frame containing CpG site information, including chromosome, position, strand, statistical test statistic, delta-betas, adjusted p-values, and probe IDs.}

\item{fdr}{The False Discovery Rate (FDR) threshold for identifying significant DMRs. Default is 0.05.}

\item{maxgap}{The maximum gap between neighboring CpGs to consider when defining DMRs. Default is 1000.}

\item{pcutoff}{The p-value threshold for filtering CpG sites before DMR analysis. Default is 0.2.}

\item{genome}{The genome build to use for genomic annotation. Default is "hg19."}
}
\value{
A data frame with information about identified DMRs, including chromosome, start and end positions, and associated CpG probe IDs.
}
\description{
This function searches for Differentially Methylated Regions (DMRs) in DNA methylation data using the DMRcate package. It takes a set of CpG sites with associated statistical information and annotates them for DMR analysis.
}
