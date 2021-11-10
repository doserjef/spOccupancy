\name{summary.PGOcc}

\alias{summary.PGOcc}
\alias{print.PGOcc}

\title{Methods for PGOcc Object}

\description{
  Methods for extracting information from fitted single species occupancy (\code{PGOcc}) model.
}

\usage{
\method{summary}{PGOcc}(object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), 
	digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{PGOcc}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{PGOcc}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{PGOcc}, including methods to the generic functions \code{\link[base]{print}} and \code{\link[base]{summary}}. 
}

\value{
  No return value, called to display summary information of a \code{PGOcc} object.
}


\keyword{model}
