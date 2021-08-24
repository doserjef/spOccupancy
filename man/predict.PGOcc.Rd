\name{predict.PGOcc}
\alias{predict.PGOcc}
\title{Function for prediction at new locations for single species occupancy models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new locations given an object of class `PGOcc`. 
}

\usage{
\method{predict}{PGOcc}(object, X.0, sub.sample, ...)
}

\arguments{

  \item{object}{an object of class PGOcc}

  \item{X.0}{the design matrix for prediction locations.}

  \item{sub.sample}{an optional list that specifies the samples to include in the composition sampling algorithm. Valid tags are \code{start}, \code{end}, and \code{thin}. Given the values associated with tags, the sample subset is selected using \code{seq(as.integer(start), as.integer(end), by = as.integer(thin))}}

  \item{...}{currently no additional arguments}
}

\value{
An object of class \code{predict.PGOcc}...
}

\examples{
}
