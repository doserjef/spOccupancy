\name{summary.ppcOcc}

\alias{summary.ppcOcc}

\title{Methods for ppcOcc Object}

\description{
  Methods for extracting information from posterior predictive check objects of class \code{ppcOcc}.
}

\usage{
\method{summary}{ppcOcc}(object, level, digits = max(3L, getOption("digits") - 3L), \dots)
}

\arguments{
  \item{object}{object of class \code{ppcOcc}.}
  \item{level}{a quoted keyword for multi-species models that indicates 
    the level to summarize the posterior predictive check. Valid key words
    are: \code{"community"}, \code{"species"}, or \code{"both"}.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted posterior predictive
  check objects of class \code{ppcOcc}, including methods to the generic function
  \code{\link[base]{summary}}. 
}

\value{
  No return value, called to display summary information of a \code{ppcOcc} object.
}

