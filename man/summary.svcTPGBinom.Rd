\name{summary.svcTPGBinom}

\alias{summary.svcTPGBinom}
\alias{print.svcTPGBinom}
\alias{plot.svcTPGBinom}

\title{Methods for svcTPGBinom Object}

\description{
  Methods for extracting information from fitted multi-season single-species 
  spatially-varying coefficient binomial model (\code{svcTPGBinom}).
}

\usage{
\method{summary}{svcTPGBinom}(object, quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{svcTPGBinom}(x, \dots)
\method{plot}{svcTPGBinom}(x, param, density = TRUE, \dots)
}

\arguments{
  \item{object, x}{object of class \code{svcTPGBinom}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{param}{parameter name for which to generate a traceplot. Valid names are 
               \code{"beta", "beta.star", "sigma.sq.psi", "theta"}.}
  \item{density}{logical value indicating whether to also generate a density plot for each parameter
                 in addition to the MCMC traceplot.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{svcTPGBinom}, including methods to the generic functions 
  \code{\link[base]{print}}, \code{\link[base]{summary}}, \code{\link[base]{plot}}. 
}

\value{
  No return value, called to display summary information of a \code{svcTPGBinom} object.
}

\keyword{model}
