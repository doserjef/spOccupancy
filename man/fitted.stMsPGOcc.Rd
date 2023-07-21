\name{fitted.stMsPGOcc}

\alias{fitted.stMsPGOcc}

\title{Extract Model Fitted Values for stMsPGOcc Object}

\description{
  Method for extracting model fitted values and detection probability values from a fitted multi-species multi-season spatial occupancy (\code{stMsPGOcc}) model. 
}

\usage{
\method{fitted}{stMsPGOcc}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{stMsPGOcc}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values and detection probability values for fitted model objects of class \code{stMsPGOcc}. 
}

\value{
  A list comprised of: 

  \item{y.rep.samples}{A five-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, species, sites, primary time period, and replicates.}
  \item{p.samples}{A five-dimensional numeric array of detection probability values. Array dimensions correspond to MCMC samples, species, sites, primary time period, and replicates.}
}

\keyword{model}

