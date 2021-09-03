\name{simOcc}
\alias{simOcc}
\title{Simulate Single-Species Occupancy Data}

\description{
  The function \code{simOcc} simulates single-species occupancy data for simulation studies, power assessments, or function testing. Data can be optionally simulated with a spatial Gaussian Process following an exponential correlation function.
}

\usage{
simOcc(J.x, J.y, n.rep, beta, alpha, sigma.sq = 2, phi = 3/0.5, sp = FALSE)
}

\arguments{
\item{J.x}{a single numeric value indicating the number of sites to simulate detection-nondetection data along the horizontal axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{J.y}{a single numeric value indicating the number of sites to simulate detection-nondetection data along the vertical axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{n.rep}{a numeric vector of length \eqn{J = J.x \times J.y}{J = J.x x J.y} indicating the number of repeat visits at each of the \eqn{J}{J} sites.}

\item{beta}{a numeric vector containing the intercept and regression coefficient parameters for the occupancy portion of the single-species occupancy model.}

\item{alpha}{a numeric vector containing the intercept and regression coefficient parameters for the detection portion of the single-species occupancy model.}

\item{sigma.sq}{a numeric value indicating the spatial variance parameter. Ignored when \code{sp = FALSE}.}

\item{phi}{a numeric value indicating the spatial range parameter. Ignored when \code{sp = FALSE}. }

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
  \item{w}{a \eqn{J \times 1}{J x 1} matrix of the spatial random effects. Only used to simulate data when \code{sp = TRUE}.}
  \item{psi}{a \eqn{J \times 1}{J x 1} matrix of the occupancy probabilities for each site.}
  \item{z}{a length \eqn{J}{J} vector of the latent occupancy states at each site.}
  \item{p}{a \code{J x max(n.rep)} matrix of the detection probabilities for each site and replicate combination. Sites with fewer than \code{max(n.rep)} replicates will contain \code{NA} values.}
  \item{y}{a \code{J x max(n.rep)} matrix of the raw detection-nondetection data for each site and replicate combination.}
}

\examples{
set.seed(400)
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
