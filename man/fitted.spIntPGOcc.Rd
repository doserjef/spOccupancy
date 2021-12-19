\name{fitted.spIntPGOcc}

\alias{fitted.spIntPGOcc}

\title{Extract Model Fitted Values for spIntPGOcc Object}

\description{
  Method for extracting model fitted values from a fitted single-species integrated spatial occupancy (\code{spIntPGOcc}) model. 
}

\usage{
\method{fitted}{spIntPGOcc}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{spIntPGOcc}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values for fitted model objects of class \code{spIntPGOcc}. 
}

\value{
  A list of three-dimensional numeric arrays of fitted values for each individual data source for use in Goodness of Fit assessments.
}

\keyword{model}

