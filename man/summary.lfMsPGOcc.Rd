\name{summary.lfMsPGOcc}

\alias{summary.lfMsPGOcc}
\alias{print.lfMsPGOcc}

\title{Methods for lfMsPGOcc Object}

\description{
  Methods for extracting information from a fitted latent factor multi-species occupancy model (\code{lfMsPGOcc}).
}

\usage{
\method{summary}{lfMsPGOcc}(object, level = 'both', quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{lfMsPGOcc}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{lfMsPGOcc}.}
  \item{level}{a quoted keyword that indicates the level to summarize the 
    model results. Valid key words are: \code{"community"}, 
    \code{"species"}, or \code{"both"}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{lfMsPGOcc}, including methods to the generic functions \code{\link[base]{print}} and \code{\link[base]{summary}}. 
}

\value{
  No return value, called to display summary information of a \code{lfMsPGOcc} object.
}

\keyword{model}
