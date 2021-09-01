\name{simMsOcc}
\alias{simMsOcc}
\title{Simulate Multi-Species Occupancy Data}

\description{
  The function \code{simMsOcc} simulates multi-species occupancy data for simulation studies, power assessments, or function testing. Data can be optionally simulated with a spatial Gaussian Process following an exponential correlation function.
}

\usage{
simMsOcc(J.x, J.y, n.rep, N, beta, alpha, sigma.sq = 2, phi = 3/0.5, sp = FALSE)
}

\arguments{

\item{J.x}{a single numeric value indicating the number of sites to simulate detection-nondetection data along the horizontal axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{J.y}{a single numeric value indicating the number of sites to simulate detection-nondetection data along the vertical axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{n.rep}{a numeric vector of length \eqn{J = J.x \times J.y}{J = J.x x J.y} indicating the number of repeat visits at each of the \eqn{J}{J} sites.}

\item{N}{a single numeric value indicating the number of species to simulate detection-nondetection data.}

\item{beta}{a numeric matrix with \eqn{N}{N} rows containing the intercept and regression coefficient parameters for the occupancy portion of the multispecies occupancy model. Each row corresponds to the regression coefficients for a given species.} 

\item{alpha}{a numeric matrix with \eqn{N}{N} rows containing the intercept and regression coefficient parameters for the detection portion of the multispecies occupancy model. Each row corresponds to the regression coefficients for a given species.}

\item{sigma.sq}{a numeric vector of length \eqn{N}{N} containing the spatial variance parameter for each species. Ignored when \code{sp = FALSE}. }

\item{phi}{a numeric vector of length \eqn{N}{N} containing the spatial range parameter for each species. Ignored when \code{sp = FALSE}.}

\item{sp}{a logical value indicating whether to simulate a spatially-explicit occupancy model with a Gaussian process and exponential correlation function. By default set to \code{FALSE}.}
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  A list comprised of: 

  \item{X}{a \eqn{J \times p.occ}{J x p.occ} numeric design matrix for the occupancy portion of the model.}
  \item{X.p}{a three-dimensional numeric array with dimensions corresponding to sites, repeat visits, and number of detection regression coefficients. This is the design matrix used for the detection portion of the occupancy model.}
  \item{coords}{a \eqn{J \times J}{J x J} numeric matrix of coordinates of each occupancy site. Required for spatial models.}
  \item{D}{a \eqn{J \times J}{J x J} numeric matrix consisting of the distances between all occupancy sites.}
  \item{w}{a \eqn{N \times J}{N x J} matrix of the spatial random effects for each species. Only used to simulate data when \code{sp = TRUE}.}
  \item{psi}{a \eqn{N \times J}{N x J} matrix of the occupancy probabilities for each species at each site.}
  \item{z}{a \eqn{N \times J}{N x J} matrix of the latent occupancy states for each species at each site.}
  \item{p}{a \code{N x J x max(n.rep)} array of the detection probabilities for each species at each site site and replicate combination. Sites with fewer than \code{max(n.rep)} replicates will contain \code{NA} values.}
  \item{y}{a \code{N x J x max(n.rep)} array of the raw detection-nondetection data for each species at each site and replicate combination. Sites with fewer than \code{max(n.rep)} replicates will contain \code{NA} values.}


}

\examples{
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, size = J, replace = TRUE)
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
sp <- TRUE

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
		phi = phi, sigma.sq = sigma.sq, sp = sp)
}
