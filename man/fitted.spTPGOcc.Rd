\name{fitted.spTPGOcc}

\alias{fitted.spTPGOcc}

\title{Extract Model Fitted Values for spTPGOcc Object}

\description{
  Method for extracting model fitted values and detection probabilities from a fitted multi-season single-species spatial occupancy (\code{spTPGOcc}) model. 
}

\usage{
\method{fitted}{spTPGOcc}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{spTPGOcc}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values and detection probabilities for fitted model objects of class \code{spTPGOcc}. 
}

\value{
  A list comprised of: 

  \item{y.rep.samples}{A four-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, sites, primary time periods, and replicates.}

  \item{p.samples}{A four-dimensional numeric array of detection probability values. Array dimensions correspond to MCMC samples, sites, primary time periods, and replicates.}
}

\keyword{model}

