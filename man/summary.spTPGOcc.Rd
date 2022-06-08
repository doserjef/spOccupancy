\name{summary.spTPGOcc}

\alias{summary.spTPGOcc}
\alias{print.spTPGOcc}

\title{Methods for spTPGOcc Object}

\description{
  Methods for extracting information from fitted multi-season single-species spatial occupancy (\code{spTPGOcc}) model.
}

\usage{
\method{summary}{spTPGOcc}(object, quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{spTPGOcc}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{spTPGOcc}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{spTPGOcc}, including methods to the generic functions \code{\link[base]{print}} and \code{\link[base]{summary}}. 
}

\value{
  No return value, called to display summary information of a \code{spTPGOcc} object.
}


\keyword{model}
