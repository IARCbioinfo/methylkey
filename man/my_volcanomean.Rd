% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.r
\name{my_volcanomean}
\alias{my_volcanomean}
\title{Create a volcano plot to visualize differential methylation analysis results}
\usage{
my_volcanomean(df)
}
\arguments{
\item{df}{A data frame containing the data for differential methylation analysis.}
}
\value{
A ggplot2 volcano plot showing differential methylation analysis results.
}
\description{
This function generates a volcano plot to visualize differential methylation analysis results using ggplot2. It takes a data frame 'dt' as input, groups the data by 'dmrdmrtool', 'ID', and 'fdr', calculates summary statistics, and creates a plot.
}
