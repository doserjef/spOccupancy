\name{fitted.lfMsPGOcc}

\alias{fitted.lfMsPGOcc}

\title{Extract Model Fitted Values for lfMsPGOcc Object}

\description{
  Method for extracting model fitted values and detection probability values from a fitted latent factor multi-species occupancy (\code{lfMsPGOcc}) model. 
}

\usage{
\method{fitted}{lfMsPGOcc}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{lfMsPGOcc}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values and detection probability values for fitted model objects of class \code{lfMsPGOcc}. 
}

\value{
  A list comprised of: 
  \item{y.rep.samples}{A four-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, species, sites, and replicates.}
  \item{p.samples}{A four-dimensional numeric array of detection probability values. Array dimensions correspond to MCMC samples, species, sites, and replicates.}
}

\keyword{model}

