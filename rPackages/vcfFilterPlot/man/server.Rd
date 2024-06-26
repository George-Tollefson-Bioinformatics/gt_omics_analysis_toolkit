% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcfFilterPlotFunctions.R
\name{server}
\alias{server}
\title{Shiny Server Function for VCF Missingness Analysis}
\usage{
server(input, output, session)
}
\arguments{
\item{input, }{output, session Reactive elements from Shiny UI.}
}
\description{
Handles the server-side logic for the VCF Missingness Analysis Shiny application.
It processes user inputs, performs filtering based on set thresholds, 
and generates plots for sample and site missingness as well as QUAL score distributions.
}
\keyword{internal}
