\name{simOcc}
\alias{simOcc}
\title{Simulate Single-Species Occupancy Data}
\usage{
simOcc(J.x, J.y, K, beta, alpha, sigma.sq = 2, phi = 3/0.5, sp = FALSE)
}
\arguments{
\item{J.x}{Number of cells on x-axis}

\item{J.y}{Number of cells on y-axis}

\item{K}{Number of repeat visits}

\item{beta}{Occupancy covariate values}

\item{alpha}{Detection covariate values}

\item{sigma.sq}{Spatial variance parameter}

\item{phi}{Spatial range parameter}

\item{sp}{Logical indicating if simulate spatial model or not.}
}
\value{
stuff
}
\description{
Simulate Single-Species Occupancy Data
}
\examples{
J.x <- 10
J.y <- 10
K <- rep(4, J.x * J.y)
beta <- c(0.5, -0.15)
alpha <- c(0.7, 0.4)
phi <- 3 / .6
sigma.sq <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, K = K, beta = beta, alpha = alpha,
             sigma.sq = sigma.sq, phi = phi, sp = TRUE)
}
