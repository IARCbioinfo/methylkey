% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.r
\name{plot_NA}
\alias{plot_NA}
\title{Create a bar plot to visualize missing values (NA)}
\usage{
plot_NA(df)
}
\arguments{
\item{df}{A data frame containing the data to be plotted, with columns 'name' and 'num_na_cg'.}
}
\value{
A ggplot2 bar plot for visualizing missing values (NA).
}
\description{
This function generates a bar plot to visualize the number of missing values (NA) for different categories using ggplot2. It plots data from a data frame, where 'name' is used on the x-axis and 'num_na_cg' on the y-axis.
}
