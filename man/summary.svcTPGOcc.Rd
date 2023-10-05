\name{summary.svcTPGOcc}

\alias{summary.svcTPGOcc}
\alias{print.svcTPGOcc}
\alias{plot.svcTPGOcc}

\title{Methods for svcTPGOcc Object}

\description{
  Methods for extracting information from fitted multi-season single-species spatially-varying coefficient occupancy (\code{svcTPGOcc}) model.
}

\usage{
\method{summary}{svcTPGOcc}(object, quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{svcTPGOcc}(x, \dots)
\method{plot}{svcTPGOcc}(x, param, density = TRUE, \dots)
}

\arguments{
  \item{object, x}{object of class \code{svcTPGOcc}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{param}{parameter name for which to generate a traceplot. Valid names are 
               \code{"beta", "beta.star", "sigma.sq.psi", "alpha", "alpha.star", "sigma.sq.p", "theta"}.}
  \item{density}{logical value indicating whether to also generate a density plot for each parameter
                 in addition to the MCMC traceplot.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{svcTPGOcc}, including methods to the generic functions \code{\link[base]{print}}, \code{\link[base]{summary}}, and \code{\link[base]{plot}}. 
}

\value{
  No return value, called to display summary information of a \code{svcTPGOcc} object.
}


\keyword{model}
