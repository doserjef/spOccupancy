\name{summary.intPGOcc}

\alias{summary.intPGOcc}
\alias{print.intPGOcc}

\title{Methods for intPGOcc Object}

\description{
  Methods for extracting information from fitted single species integrated occupancy (\code{intPGOcc}) model.
}

\usage{
\method{summary}{intPGOcc}(object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), 
	digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{intPGOcc}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{intPGOcc}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{intPGOcc}, including methods to the generic functions \code{\link[base]{print}} and \code{\link[base]{summary}}. 
}

\value{
  No return value, called to display summary information of a \code{intPGOcc} object.
}

\keyword{model}
