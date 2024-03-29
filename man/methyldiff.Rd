% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regression.r
\name{methyldiff}
\alias{methyldiff}
\title{methyldiff}
\usage{
methyldiff(
  se = NULL,
  model = "~group",
  case = NULL,
  control = NULL,
  method = "ls",
  niter = 50,
  ncore = 2,
  qval = 0.05
)
}
\arguments{
\item{model}{model to apply with sva (string) eg: "~group+gender"}

\item{case}{case group}

\item{control}{control group}

\item{method}{\if{html}{\out{<limma>}} ls or robust}

\item{niter}{\if{html}{\out{<limma>}} number of iteration for robust method}

\item{ncore}{\if{html}{\out{<limma>}} number of core to use for limma}

\item{qval}{cutoff for significant probes}

\item{betas}{matrix of betas}

\item{pdata}{sample Sheet (dataframe)}

\item{sva}{run sva correction before regression}

\item{modelSVA}{model to use for batch correction}
}
\value{
toptable for all probes
}
\description{
differential analysis,
}
