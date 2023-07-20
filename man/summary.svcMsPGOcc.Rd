\name{summary.svcMsPGOcc}

\alias{summary.svcMsPGOcc}
\alias{print.svcMsPGOcc}

\title{Methods for svcMsPGOcc Object}

\description{
  Methods for extracting information from fitted multi-species spatially-varying coefficient occupancy model.
}

\usage{
\method{summary}{svcMsPGOcc}(object, level, quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{svcMsPGOcc}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{svcMsPGOcc}.}
  \item{level}{a quoted keyword that indicates the level to summarize the 
    model results. Valid key words are: \code{"community"}, 
    \code{"species"}, or \code{"both"}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{svcMsPGOcc}, including methods to the generic functions \code{\link[base]{print}} and \code{\link[base]{summary}}. 
}

\value{
  No return value, called to display summary information of a \code{svcMsPGOcc} object.
}

\keyword{model}
