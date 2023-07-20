\name{fitted.tMsPGOcc}

\alias{fitted.tMsPGOcc}

\title{Extract Model Fitted Values for tMsPGOcc Object}

\description{
  Method for extracting model fitted values and detection probability values from a fitted multi-species multi-season occupancy (\code{tMsPGOcc}) model. 
}

\usage{
\method{fitted}{tMsPGOcc}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{tMsPGOcc}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values and detection probability values for fitted model objects of class \code{tMsPGOcc}. 
}

\value{
  A list comprised of: 

  \item{y.rep.samples}{A five-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, species, sites, primary time period, and replicates.}
  \item{p.samples}{A five-dimensional numeric array of detection probability values. Array dimensions correspond to MCMC samples, species, sites, primary time period, and replicates.}
}

\keyword{model}

