% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MethylKey.r
\name{getPlateform}
\alias{getPlateform}
\title{Determine the methylation array platform from the number of rows in a matrix}
\usage{
getPlateform(mat)
}
\arguments{
\item{mat}{A matrix for which to determine the platform.}
}
\value{
The identified platform as a character string, or "unknown" if the platform cannot be determined.
}
\description{
This function attempts to determine the methylation array platform (e.g., IlluminaHumanMethylationEPIC, IlluminaHumanMethylation450k, IlluminaHumanMethylation27k, IlluminaMouseMethylation285k) based on the number of rows in a matrix.
}
