\name{fitted.spPGOcc}

\alias{fitted.spPGOcc}

\title{Extract Model Fitted Values for spPGOcc Object}

\description{
  Method for extracting model fitted values and detection probabilities from a fitted single-species spatial occupancy (\code{spPGOcc}) model. 
}

\usage{
\method{fitted}{spPGOcc}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{spPGOcc}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values and detection probabilities for fitted model objects of class \code{spPGOcc}. 
}

\value{
  A list comprised of: 

  \item{y.rep.samples}{A three-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, sites, and replicates.}
  \item{p.samples}{A three-dimensional numeric array of detection probability values. Array dimensions correspond to MCMC samples, sites, and replicates.}
}

\keyword{model}

