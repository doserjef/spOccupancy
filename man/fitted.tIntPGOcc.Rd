\name{fitted.tIntPGOcc}

\alias{fitted.tIntPGOcc}

\title{Extract Model Fitted Values for tIntPGOcc Object}

\description{
  Method for extracting model fitted values and detection probabilities from a fitted multi-season single-species integrated occupancy (\code{tIntPGOcc}) model. 
}

\usage{
\method{fitted}{tIntPGOcc}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{tIntPGOcc}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values and detection probabilities for fitted model objects of class \code{tIntPGOcc}. 
}

\value{
  A list comprised of: 

  \item{y.rep.samples}{a list of four-dimensional numeric arrays of fitted values for each data set for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, sites, primary time periods, and replicates.}

  \item{p.samples}{a list of four-dimensional numeric arrays of detection probability values for each data set. Array dimensions correspond to MCMC samples, sites, primary time periods, and replicates.}
}

\keyword{model}

