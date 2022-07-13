\name{summary.stPGOcc}

\alias{summary.stPGOcc}
\alias{print.stPGOcc}

\title{Methods for stPGOcc Object}

\description{
  Methods for extracting information from fitted multi-season single-species spatial occupancy (\code{stPGOcc}) model.
}

\usage{
\method{summary}{stPGOcc}(object, quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{stPGOcc}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{stPGOcc}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{stPGOcc}, including methods to the generic functions \code{\link[base]{print}} and \code{\link[base]{summary}}. 
}

\value{
  No return value, called to display summary information of a \code{stPGOcc} object.
}


\keyword{model}
