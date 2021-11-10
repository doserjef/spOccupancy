\name{summary.spPGOcc}

\alias{summary.spPGOcc}
\alias{print.spPGOcc}

\title{Methods for spPGOcc Object}

\description{
  Methods for extracting information from fitted single species spatial 
  occupancy (\code{spPGOcc}) model.
}

\usage{
\method{summary}{spPGOcc}(object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{spPGOcc}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{spPGOcc}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{spPGOcc}, including methods to the generic functions 
  \code{\link[base]{print}} and \code{\link[base]{summary}}. 
}

\value{
  No return value, called to display summary information of a \code{spPGOcc} object.
}

\keyword{model}
