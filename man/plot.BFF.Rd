% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.BFF}
\alias{plot.BFF}
\title{Plot Bayes Factor Function}
\usage{
\method{plot}{BFF}(x, plot = TRUE, ...)
}
\arguments{
\item{x}{a BFF object}

\item{plot}{whether plot should be generated.
Default to \code{TRUE}. If \code{FALSE} a data
frame with the support, prior ordinates, and
posterior ordinates is returned instead.}

\item{...}{additional arguments to the plotting
function. These include: \describe{
 \item{"title"}{tile of the figure}
 \item{"xlab"}{x-axis label of the figure}
 \item{"ylab"}{y-axis label of the figure}
 \item{"add_segments"}{whether effect size
 segments should be added to the figure. Available only
 for standardized effect sizes. Defaults to \code{TRUE}}
}}
}
\value{
either a ggplot2 object if \code{plot = TRUE} or a data.frame
with a Bayes factor function if \code{plot = FALSE}
}
\description{
Creates a Bayes factor function plot
of an BFF object. The BFF object needs to be specified
with \code{omega = NULL} or .
}
\seealso{
[z_test_BFF()], [t_test_BFF()], [chi2_test_BFF()], [f_test_BFF()]
}
