\name{fitted.sfJSDM}

\alias{fitted.sfJSDM}

\title{Extract Model Fitted Values for sfJSDM Object}

\description{
  Method for extracting model fitted values and probability values from a fitted spatial factor joint species distribution model (\code{sfJSDM}). 
}

\usage{
\method{fitted}{sfJSDM}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{sfJSDM}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values and probability values for fitted model objects of class \code{sfJSDM}. 
}

\value{
  A list comprised of: 

  \item{z.samples}{A three-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, species, and sites.}
  \item{psi.samples}{A three-dimensional numeric array of probability values. Array dimensions correspond to MCMC samples, species, and sites.}
}

\keyword{model}

