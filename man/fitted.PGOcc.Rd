\name{fitted.PGOcc}

\alias{fitted.PGOcc}

\title{Extract Model Fitted Values for PGOcc Object}

\description{
  Method for extracting model fitted values and detection probabilities from a fitted single-species occupancy (\code{PGOcc}) model. 
}

\usage{
\method{fitted}{PGOcc}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{PGOcc}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values and detection probabilities for fitted model objects of class \code{PGOcc}. 
}

\value{
  A list comprised of: 

  \item{y.rep.samples}{A three-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, sites, and replicates.}

  \item{p.samples}{A three-dimensional numeric array of detection probability values. Array dimensions correspond to MCMC samples, sites, and replicates.}
}

\keyword{model}

