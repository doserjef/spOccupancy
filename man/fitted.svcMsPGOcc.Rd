\name{fitted.svcMsPGOcc}

\alias{fitted.svcMsPGOcc}

\title{Extract Model Fitted Values for svcMsPGOcc Object}

\description{
  Method for extracting model fitted values and detection probability values from a fitted multi-species spatially varying coefficient occupancy (\code{svcMsPGOcc}) model. 
}

\usage{
\method{fitted}{svcMsPGOcc}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{svcMsPGOcc}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values and detection probability values for fitted model objects of class \code{svcMsPGOcc}. 
}

\value{
  A list comprised of: 

  \item{y.rep.samples}{A four-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, species, sites, and replicates.}

  \item{p.samples}{A four-dimensional numeric array of detection probability values. Array dimensions correspond to MCMC samples, species, sites, and replicates.}
}

\keyword{model}

