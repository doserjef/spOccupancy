\name{summary.svcPGBinom}

\alias{summary.svcPGBinom}
\alias{print.svcPGBinom}

\title{Methods for svcPGBinom Object}

\description{
  Methods for extracting information from fitted single-species spatially-varying 
  coefficient binomial model (\code{svcPGBinom}).
}

\usage{
\method{summary}{svcPGBinom}(object, quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{svcPGBinom}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{svcPGBinom}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{svcPGBinom}, including methods to the generic functions 
  \code{\link[base]{print}} and \code{\link[base]{summary}}. 
}

\value{
  No return value, called to display summary information of a \code{svcPGBinom} object.
}

\keyword{model}
