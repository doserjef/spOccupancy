\name{fitted.spPGOcc}

\alias{fitted.spPGOcc}

\title{Extract Model Fitted Values for spPGOcc Object}

\description{
  Method for extracting model fitted values from a fitted single species spatial occupancy (\code{spPGOcc}) model. 
}

\usage{
\method{fitted}{spPGOcc}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{spPGOcc}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values for fitted model objects of class \code{spPGOcc}. 
}

\value{
  A three-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, sites, and replicates.
}

\keyword{model}

