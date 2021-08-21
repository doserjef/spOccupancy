\name{simMsOcc}
\alias{simMsOcc}
\title{Simulate Multi-Species Occupancy Data}

\description{
  The function \code{simMsOcc} simulates multi-species occupancy data for simulation studies, power assessments, or function testing. Data can be optionally be simulated with a spatial Gaussian Process following an exponential correlation function.
}

\usage{
simMsOcc(J.x, J.y, K, N, beta, alpha, sigma.sq = 2, phi = 3/0.5, sp = FALSE)
}

\arguments{
\item{J.x}{Number of cells on x-axis}

\item{J.y}{Number of cells on y-axis}

\item{K}{Number of repeat visits}

\item{N}{Number of species}

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
J.x <- 8
J.y <- 8
J <- J.x * J.y
K <- sample(2:4, size = J, replace = TRUE)
N <- 10
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, -0.15)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 0.3)
# Detection
alpha.mean <- c(0.5, 0.2)
tau.sq.alpha <- c(0.2, 0.3)
p.det <- length(alpha.mean)
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
# Spatial parameters if desired
phi <- runif(N, 3/1, 3/.1)
sigma.sq <- runif(N, 0.3, 3)
sp = FALSE

dat <- simMsOcc(J.x = J.x, J.y = J.y, K = K, N = N, beta = beta, alpha = alpha,
		phi = phi, sigma.sq = sigma.sq, sp = sp)
}
