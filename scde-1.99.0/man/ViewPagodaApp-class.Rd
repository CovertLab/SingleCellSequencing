% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/functions.R
\docType{class}
\name{ViewPagodaApp-class}
\alias{ViewPagodaApp}
\alias{ViewPagodaApp-class}
\title{A Reference Class to represent the PAGODA application}
\description{
This ROOK application class enables communication with the client-side ExtJS framework and Inchlib HTML5 canvas libraries to create the graphical user interface for PAGODA
Refer to the code in \code{\link{make.pagoda.app}} for usage example
}
\section{Fields}{

\describe{
\item{\code{results}}{Output of the pathway clustering and redundancy reduction}

\item{\code{genes}}{List of genes to display in the Detailed clustering panel}

\item{\code{pathways}}{}

\item{\code{mat}}{Matrix of posterior mode count estimates}

\item{\code{matw}}{Matrix of weights associated with each estimate in \code{mat}}

\item{\code{goenv}}{Gene set list as an environment}

\item{\code{renv}}{Global environment}

\item{\code{name}}{Name of the application page; for display as the page title}

\item{\code{trim}}{Trim quantity used for Winsorization for visualization}

\item{\code{batch}}{Any batch or other known confounders to be included in the visualization as a column color track}
}}

