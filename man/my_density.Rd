% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.r
\name{my_density}
\alias{my_density}
\title{Create a density plot of DMR lengths grouped by dmrtool}
\usage{
my_density(dt)
}
\arguments{
\item{dt}{A data frame containing the data.}
}
\value{
A ggplot2 density plot showing DMR length distributions grouped by 'dmrtool'.
}
\description{
This function generates a density plot of DMR (Differentially Methylated Region) lengths grouped by 'dmrtool'. It takes a data frame 'dt' as input, calculates the length of each DMR, and then creates the density plot.
}
