\name{residuals.PGOcc}

\alias{residuals.PGOcc}

\title{Occupancy and detection residuals for \code{PGOcc} models}

\description{
  Method for calculating occupancy and detection residuals for single-species occupancy models (\code{PGOcc}) following the approach of Wright et al. (2019). 
}

\usage{
\method{residuals}{PGOcc}(object, n.post.samples = 100, ...)
}

\arguments{
  \item{object}{object of class \code{PGOcc}.}
  \item{n.post.samples}{the number of posterior MCMC samples to calculate the residuals for. By default this is set to 100. If set to a value less than the total number of MCMC samples saved for the model, residuals will be calculated for a random subset of the total MCMC samples. Maximum value is the total number of MCMC samples saved.}
  \item{\dots}{currently no additional arguments}
}

\references{
  Wright, W. J., Irvine, K. M., & Higgs, M. D. (2019). Identifying occupancy model inadequacies: can residuals separately assess detection and presence?. Ecology, 100(6), e02703.
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}
}


\value{
  A list comprised of: 

  \item{occ.resids}{a matrix of occupancy residuals with first dimension equal to \code{n.post.samples} and second dimension equal to the number of sites in the data set.}

  \item{det.resids}{a three-dimensional array of detection residuals with first dimension equal to \code{n.post.samples}, second dimension equal to the number of sites in the data set, and third dimension equal to the maximum number of repeat visits. Note detection residuals are only calculated for a given site and MCMC iteration when the species is present.} 
}
