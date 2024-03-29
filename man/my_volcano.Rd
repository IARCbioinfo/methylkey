% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.r
\name{my_volcano}
\alias{my_volcano}
\title{Create a volcano plot to visualize statistical results}
\usage{
my_volcano(df)
}
\arguments{
\item{df}{A data frame containing the data for statistical results with columns 'deltabetas' and 'P.Value'.}
}
\value{
A ggplot2 volcano plot showing statistical results.
}
\description{
This function generates a volcano plot to visualize statistical results using ggplot2. It takes a data frame 'df' as input, where 'deltabetas' and 'P.Value' are used as x and y coordinates. The points in the plot are colored based on significance and effect size.
}
