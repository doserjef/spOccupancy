\name{simIntOcc}
\alias{simIntOcc}
\title{Simulate Single-Species Occupancy Data from Multiple Data Sources}

\description{
  The function \code{simIntOcc} simulates single-species occupancy data from multiple data sources for simulation studies, power assessments, or function testing of integrated occupancy models. Data can be optionally be simulated with a spatial Gaussian Process on the occupancy process following an exponential correlation function.
}

\usage{
simIntOcc(n.data, J.x, J.y, J.obs, K.obs, beta, alpha, sigma.sq = 2, 
	  phi = 3/0.5, sp = FALSE)
}
\arguments{

\item{n.data}{Number of detection-nondetection data sources}

\item{J.x}{Number of cells on x-axis}

\item{J.y}{Number of cells on y-axis}

\item{J.obs}{Number of unique observed cells}

\item{K.obs}{Number of repeat visits}

\item{beta}{Occupancy covariate values}

\item{alpha}{Detection covariate values}

\item{sigma.sq}{Spatial variance parameter}

\item{phi}{Spatial range parameter}

\item{sp}{Logical indicating if simulate spatial model or not.}
}

\value{
stuff
}

\examples{
# Update
}
