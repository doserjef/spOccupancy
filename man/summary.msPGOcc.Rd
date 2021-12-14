\name{summary.msPGOcc}

\alias{summary.msPGOcc}
\alias{print.msPGOcc}

\title{Methods for msPGOcc Object}

\description{
  Methods for extracting information from fitted multi-species occupancy (\code{msPGOcc}) model.
}

\usage{
\method{summary}{msPGOcc}(object, level = 'both', quantiles = c(0.025, 0.5, 0.975), 
	digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{msPGOcc}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{msPGOcc}.}
  \item{level}{a quoted keyword that indicates the level to summarize the 
    posterior predictive check. Valid key words are: \code{"community"}, 
    \code{"species"}, or \code{"both"}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{msPGOcc}, including methods to the generic functions \code{\link[base]{print}} and \code{\link[base]{summary}}. 
}

\value{
  No return value, called to display summary information of a \code{msPGOcc} object.
}

\keyword{model}
