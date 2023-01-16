\name{summary.intMsPGOcc}

\alias{summary.intMsPGOcc}
\alias{print.intMsPGOcc}

\title{Methods for intMsPGOcc Object}

\description{
  Methods for extracting information from fitted integrated multi-species occupancy (\code{intMsPGOcc}) models.
}

\usage{
\method{summary}{intMsPGOcc}(object, level = 'both', quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{intMsPGOcc}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{intMsPGOcc}.}
  \item{level}{a quoted keyword that indicates the level to summarize the 
    model results. Valid key words are: \code{"community"}, 
    \code{"species"}, or \code{"both"}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{intMsPGOcc}, including methods to the generic functions \code{\link[base]{print}} and \code{\link[base]{summary}}. 
}

\value{
  No return value, called to display summary information of a \code{intMsPGOcc} object.
}

\keyword{model}
