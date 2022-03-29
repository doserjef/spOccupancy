\name{fitted.msPGOcc}

\alias{fitted.msPGOcc}

\title{Extract Model Fitted Values for msPGOcc Object}

\description{
  Method for extracting model fitted values and detection probability values from a fitted multi-species occupancy (\code{msPGOcc}) model. 
}

\usage{
\method{fitted}{msPGOcc}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{msPGOcc}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values and detection probability values for fitted model objects of class \code{msPGOcc}. 
}

\value{
  A list comprised of: 

  \item{y.rep.samples}{A four-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, species, sites, and replicates.}
  \item{p.samples}{A four-dimensional numeric array of detection probability values. Array dimensions correspond to MCMC samples, species, sites, and replicates.}
}

\keyword{model}

