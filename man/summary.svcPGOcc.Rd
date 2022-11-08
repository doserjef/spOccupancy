\name{summary.svcPGOcc}

\alias{summary.svcPGOcc}
\alias{print.svcPGOcc}

\title{Methods for svcPGOcc Object}

\description{
  Methods for extracting information from fitted single-species 
  spatially-varying coefficient occupancy (\code{svcPGOcc}) model.
}

\usage{
\method{summary}{svcPGOcc}(object, quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{svcPGOcc}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{svcPGOcc}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{svcPGOcc}, including methods to the generic functions 
  \code{\link[base]{print}} and \code{\link[base]{summary}}. 
}

\value{
  No return value, called to display summary information of a \code{svcPGOcc} object.
}

\keyword{model}
