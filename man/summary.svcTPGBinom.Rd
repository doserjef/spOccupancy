\name{summary.svcTPGBinom}

\alias{summary.svcTPGBinom}
\alias{print.svcTPGBinom}

\title{Methods for svcTPGBinom Object}

\description{
  Methods for extracting information from fitted multi-season single-species 
  spatially-varying coefficient binomial model (\code{svcTPGBinom}).
}

\usage{
\method{summary}{svcTPGBinom}(object, quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{svcTPGBinom}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{svcTPGBinom}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{svcTPGBinom}, including methods to the generic functions 
  \code{\link[base]{print}} and \code{\link[base]{summary}}. 
}

\value{
  No return value, called to display summary information of a \code{svcTPGBinom} object.
}

\keyword{model}
