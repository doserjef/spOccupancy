\name{summary.spMsPGOcc}

\alias{summary.spMsPGOcc}
\alias{print.spMsPGOcc}

\title{Methods for spMsPGOcc Object}

\description{
  Methods for extracting information from fitted multi-species spatial occupancy (\code{spMsPGOcc}) model.
}

\usage{
\method{summary}{spMsPGOcc}(object, level, quantiles = c(0.025, 0.5, 0.975), 
	 digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{spMsPGOcc}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{spMsPGOcc}.}
  \item{level}{a quoted keyword that indicates the level to summarize the 
    posterior predictive check. Valid key words are: \code{"community"}, 
    \code{"species"}, or \code{"both"}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{spMsPGOcc}, including methods to the generic functions \code{\link[base]{print}} and \code{\link[base]{summary}}. 
}

\value{
  No return value, called to display summary information of a \code{spMsPGOcc} object.
}

\keyword{model}
