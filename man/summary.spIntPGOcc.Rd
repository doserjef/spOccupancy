\name{summary.spIntPGOcc}

\alias{summary.spIntPGOcc}
\alias{print.spIntPGOcc}
\alias{plot.spIntPGOcc}

\title{Methods for spIntPGOcc Object}

\description{
  Methods for extracting information from fitted single-species spatial integrated occupancy (\code{spIntPGOcc}) model.
}

\usage{
\method{summary}{spIntPGOcc}(object, quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{spIntPGOcc}(x, \dots)
\method{plot}{spIntPGOcc}(x, param, density = TRUE, \dots)
}

\arguments{
  \item{object, x}{object of class \code{spIntPGOcc}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{param}{parameter name for which to generate a traceplot. Valid names are 
               \code{"beta", "alpha", "theta"}.}
  \item{density}{logical value indicating whether to also generate a density plot for each parameter
                 in addition to the MCMC traceplot.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{spIntPGOcc}, including methods to the generic functions \code{\link[base]{print}}, \code{\link[base]{summary}}, and \code{\link[base]{plot}}. 
}

\value{
  No return value, called to display summary information of a \code{spIntPGOcc} object.
}

\keyword{model}
