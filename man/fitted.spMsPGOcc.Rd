\name{fitted.spMsPGOcc}

\alias{fitted.spMsPGOcc}

\title{Extract Model Fitted Values for spMsPGOcc Object}

\description{
  Method for extracting model fitted values and detection probability values from a fitted multi-species spatial occupancy (\code{spMsPGOcc}) model. 
}

\usage{
\method{fitted}{spMsPGOcc}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{spMsPGOcc}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values and detection probability values for fitted model objects of class \code{spMsPGOcc}. 
}

\value{
  A list comprised of: 

  \item{y.rep.samples}{A four-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, species, sites, and replicates.}
  \item{p.samples}{A four-dimensional numeric array of detection probability values. Array dimensions correspond to MCMC samples, species, sites, and replicates.}
}

\keyword{model}

