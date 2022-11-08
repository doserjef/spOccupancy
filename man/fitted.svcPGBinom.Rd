\name{fitted.svcPGBinom}

\alias{fitted.svcPGBinom}

\title{Extract Model Fitted Values for svcPGBinom Object}

\description{
  Method for extracting model fitted values from a fitted single-species spatially-varying coefficients binomial model (\code{svcPGBinom}). 
}

\usage{
\method{fitted}{svcPGBinom}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{svcPGBinom}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values for fitted model objects of class \code{svcPGBinom}. 
}

\value{
  A two-dimensional matrix of fitted values for use in Goodness of Fit assessments. Dimensions correspond to MCMC samples and sites.
}

\keyword{model}


