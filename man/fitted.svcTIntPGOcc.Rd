\name{fitted.svcTIntPGOcc}

\alias{fitted.svcTIntPGOcc}

\title{Extract Model Fitted Values for svcTIntPGOcc Object}

\description{
  Method for extracting model fitted values and detection probabilities from a fitted multi-season single-species spatially-varying coefficient integrated occupancy (\code{svcTIntPGOcc}) model. 
}

\usage{
\method{fitted}{svcTIntPGOcc}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{svcTIntPGOcc}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values and detection probabilities for fitted model objects of class \code{svcTIntPGOcc}. 
}

\value{
  A list comprised of: 
  
  \item{y.rep.samples}{a list of four-dimensional numeric arrays of fitted values for each data set for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, sites, primary time periods, and replicates.}

  \item{p.samples}{a list of four-dimensional numeric arrays of detection probability values for each data set. Array dimensions correspond to MCMC samples, sites, primary time periods, and replicates.}
}

\keyword{model}

