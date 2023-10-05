\name{summary.intPGOcc}

\alias{summary.intPGOcc}
\alias{print.intPGOcc}
\alias{plot.intPGOcc}

\title{Methods for intPGOcc Object}

\description{
  Methods for extracting information from fitted single species integrated occupancy (\code{intPGOcc}) model.
}

\usage{
\method{summary}{intPGOcc}(object, quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{intPGOcc}(x, \dots)
\method{plot}{intPGOcc}(x, param, density = TRUE, \dots)
}

\arguments{
  \item{object, x}{object of class \code{intPGOcc}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{param}{parameter name for which to generate a traceplot. Valid names are 
               \code{"beta", "alpha"}.}
  \item{density}{logical value indicating whether to also generate a density plot for each parameter
                 in addition to the MCMC traceplot.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{intPGOcc}, including methods to the generic functions \code{\link[base]{print}}, \code{\link[base]{summary}}, and \code{\link[base]{plot}}. 
}

\value{
  No return value, called to display summary information of a \code{intPGOcc} object.
}

\keyword{model}
