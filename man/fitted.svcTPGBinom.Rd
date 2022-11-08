\name{fitted.svcTPGBinom}

\alias{fitted.svcTPGBinom}

\title{Extract Model Fitted Values for svcTPGBinom Object}

\description{
  Method for extracting model fitted values from a fitted multi-season single-species spatially-varying coefficients binomial model (\code{svcTPGBinom}). 
}

\usage{
\method{fitted}{svcTPGBinom}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{svcTPGBinom}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values for fitted model objects of class \code{svcTPGBinom}. 
}

\value{
  A three-dimensional matrix of fitted values for use in Goodness of Fit assessments. Dimensions correspond to MCMC samples, sites, and primary time periods.
}

\keyword{model}


