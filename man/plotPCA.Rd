% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plotPCA}
\alias{plotPCA}
\title{Principle component analysis (PCA) or multidimensional scaling plot between samples}
\usage{
plotPCA(exp_data, groups, do.legend=F, plot_time_points=F, log=T, do.MDS=F, cols=NULL, epsilon = 1, return_plotMatrix = FALSE, ...)
}
\arguments{
\item{exp_data}{expression matrix, columns represent the samples and rows the genes (or transcripts)}

\item{groups}{numeric vector describing the relationships of columns to certain groups (like biological replicates)}

\item{do.legend}{boolean value defining if a color legend should be ploted}

\item{plot_time_points}{boolean value indicating if the grafting time points shall be plotted}

\item{log}{boolean value defining if the expression data should be log transformed before the calculation of PCA or MDS}

\item{do.MDS}{boolean value defines if PCA or MDS should be calculated. Default: TRUE}

\item{cols}{character vector defining the colors to use for the different samples}

\item{epsilon}{numeric value indicates pseudocount used if \code{log=TRUE}}

\item{return_plotMatrix}{boolean value indicates if PCA or MDS results should be returned}

\item{...}{optional plot parameters \link{plot}}
}
\value{
a \code{matrix}
}
\description{
Plot samples on a two-dimensional scatterplot in which the distances on the plot indicate the expression differences between the samples
}
\author{
Alexander Gabel
}
