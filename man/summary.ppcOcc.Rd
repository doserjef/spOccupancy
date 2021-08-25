\name{summary.ppcOcc}

\alias{summary.ppcOcc}

\title{Methods for ppcOcc Object}

\description{
  Methods for extracting information from posterior predictive check objects of class \code{ppcOcc}.
}

\usage{
\method{summary}{ppcOcc}(object, digits = max(3L, getOption("digits") - 3L), \dots)
}

\arguments{
  \item{object}{object of class \code{ppcOcc}.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted posterior predictive
  check objects of class \code{ppcOcc}, including methods to the generic function
  \code{\link[base]{summary}}. 
}
