\name{simOcc}
\alias{simOcc}
\title{Simulate Single-Species Occupancy Data}

\description{
  The function \code{simOcc} simulates single-species occupancy data for simulation studies, power assessments, or function testing. Data can be optionally be simulated with a spatial Gaussian Process following an exponential correlation function.
}

\usage{
simOcc(J.x, J.y, n.rep, beta, alpha, sigma.sq = 2, phi = 3/0.5, sp = FALSE)
}
\arguments{
\item{J.x}{Number of cells on x-axis}

\item{J.y}{Number of cells on y-axis}

\item{n.rep}{Number of repeat visits}

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
J.x <- 10
J.y <- 10
n.rep <- rep(4, J.x * J.y)
beta <- c(0.5, -0.15)
alpha <- c(0.7, 0.4)
phi <- 3 / .6
sigma.sq <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
             sigma.sq = sigma.sq, phi = phi, sp = TRUE)
}
