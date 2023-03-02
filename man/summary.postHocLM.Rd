\name{summary.postHocLM}

\alias{summary.postHocLM}
\alias{print.postHocLM}

\title{Methods for postHocLM Object}

\description{
  Methods for extracting information from fitted posthoc linear models (\code{postHocLM}).
}

\usage{
\method{summary}{postHocLM}(object, quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{postHocLM}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{postHocLM}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{postHocLM}, including methods to the generic functions \code{\link[base]{print}} and \code{\link[base]{summary}}. 
}

\value{
  No return value, called to display summary information of a \code{postHocLM} object.
}


\keyword{model}
