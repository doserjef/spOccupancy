\name{simMsOcc}
\alias{simMsOcc}
\title{Simulate Multi-Species Detection-Nondetection Data}

\description{
  The function \code{simMsOcc} simulates multi-species detection-nondetection data for simulation studies, power assessments, or function testing. Data can be optionally simulated with a spatial Gaussian Process in the occurrence portion of the model, as well as an option to allow for species correlations using a factor modeling approach. Non-spatial random intercepts can also be included in the detection or occurrence portions of the occupancy model.
}

\usage{
simMsOcc(J.x, J.y, n.rep, n.rep.max, N, beta, alpha, psi.RE = list(), 
         p.RE = list(), sp = FALSE, svc.cols = 1, cov.model, 
	 sigma.sq, phi, nu, factor.model = FALSE, n.factors, 
         range.probs, shared.spatial = FALSE, ...)
}

\arguments{

\item{J.x}{a single numeric value indicating the number of sites to simulate detection-nondetection data along the horizontal axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{J.y}{a single numeric value indicating the number of sites to simulate detection-nondetection data along the vertical axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{n.rep}{a numeric vector of length \eqn{J = J.x \times J.y}{J = J.x x J.y} indicating the number of repeat visits at each of the \eqn{J}{J} sites.}

\item{n.rep.max}{a single numeric value indicating the maximum number of replicate surveys. This is an optional argument, with its default value set to \code{max(n.rep)}. This can be used to generate data sets with different types of missingness (e.g., simulate data across 20 days (replicate surveys) but sites are only sampled a maximum of ten times each).}

\item{N}{a single numeric value indicating the number of species to simulate detection-nondetection data.}

\item{beta}{a numeric matrix with \eqn{N}{N} rows containing the intercept and regression coefficient parameters for the occurrence portion of the multi-species occupancy model. Each row corresponds to the regression coefficients for a given species.} 

\item{alpha}{a numeric matrix with \eqn{N}{N} rows containing the intercept and regression coefficient parameters for the detection portion of the multi-species occupancy model. Each row corresponds to the regression coefficients for a given species.}

\item{psi.RE}{a list used to specify the non-spatial random intercepts included in the occurrence portion of the model. The list must have two tags: \code{levels} and \code{sigma.sq.psi}. \code{levels} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the number of levels there are in each intercept. \code{sigma.sq.psi} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the variances for each random effect. If not specified, no random effects are included in the occurrence portion of the model.}

\item{p.RE}{a list used to specify the non-spatial random intercepts included in the detection portion of the model. The list must have two tags: \code{levels} and \code{sigma.sq.p}. \code{levels} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the number of levels there are in each intercept. \code{sigma.sq.p} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the variances for each random effect. If not specified, no random effects are included in the detection portion of the model.}

\item{sp}{a logical value indicating whether to simulate a spatially-explicit occupancy model with a Gaussian process. By default set to \code{FALSE}.}

\item{svc.cols}{a vector indicating the variables whose effects will be
  estimated as spatially-varying coefficients. \code{svc.cols} is an
  integer vector with values indicating the order of covariates specified
  in the model formula (with 1 being the intercept if specified).}

\item{cov.model}{a quoted keyword that specifies the covariance function used to model the spatial dependence structure among the latent occurrence values.  Supported covariance model key words are: \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and \code{"gaussian"}.}

\item{sigma.sq}{a numeric vector of length \eqn{N}{N} containing the spatial variance parameter for each species. Ignored when \code{sp = FALSE} or when \code{factor.model = TRUE}.}

\item{phi}{a numeric vector of length \eqn{N}{N} containing the spatial decay parameter for each species. Ignored when \code{sp = FALSE}. If \code{factor.model = TRUE}, this should be of length \code{n.factors}.}

\item{nu}{a numeric vector of length \eqn{N}{N} containing the spatial smoothness parameter for each species. Only used when \code{sp = TRUE} and \code{cov.model = 'matern'}. If \code{factor.model = TRUE}, this should be of length \code{n.factors}.} 

\item{factor.model}{a logical value indicating whether to simulate data following a factor modeling approach that explicitly incoporates species correlations. If \code{sp = TRUE}, the latent factors are simulated from independent spatial processes. If \code{sp = FALSE}, the latent factors are simulated from standard normal distributions.}

\item{n.factors}{a single numeric value specifying the number of latent factors to use to simulate the data if \code{factor.model = TRUE}.}

\item{range.probs}{a numeric vector of length \code{N} where each value should fall between 0 and 1, and indicates the probability that one of the \code{J} spatial locations simulated is within the simulated range of the given species. If set to 1, every species has the potential of being present at each location.}

\item{shared.spatial}{a logical value indicating used to specify whether a common spatial process should be estimated for all species instead of the factor modeling approach.}

\item{...}{currently no additional arguments}
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  A list comprised of: 

  \item{X}{a \eqn{J \times p.occ}{J x p.occ} numeric design matrix for the occurrence portion of the model.}
  \item{X.p}{a three-dimensional numeric array with dimensions corresponding to sites, repeat visits, and number of detection regression coefficients. This is the design matrix used for the detection portion of the occupancy model.}
  \item{coords}{a \eqn{J \times 2}{J x 2} numeric matrix of coordinates of each occupancy site. Required for spatial models.}
  \item{w}{a \eqn{N \times J}{N x J} matrix of the spatial random effects for each species. Only used to simulate data when \code{sp = TRUE}. If \code{factor.model = TRUE}, the first dimension is \code{n.factors}.}
  \item{psi}{a \eqn{N \times J}{N x J} matrix of the occurrence probabilities for each species at each site.}
  \item{z}{a \eqn{N \times J}{N x J} matrix of the latent occurrence states for each species at each site.}
  \item{p}{a \code{N x J x max(n.rep)} array of the detection probabilities for each species at each site and replicate combination. Sites with fewer than \code{max(n.rep)} replicates will contain \code{NA} values.}
  \item{y}{a \code{N x J x max(n.rep)} array of the raw detection-nondetection data for each species at each site and replicate combination. Sites with fewer than \code{max(n.rep)} replicates will contain \code{NA} values.}
  \item{X.p.re}{a three-dimensional numeric array containing the levels of any detection random effect included in the model. Only relevant when detection random effects are specified in \code{p.RE}.}
  \item{X.lambda.re}{a numeric matrix containing the levels of any occurrence random effect included in the model. Only relevant when occurrence random effects are specified in \code{psi.RE}.}
  \item{alpha.star}{a numeric matrix where each row contains the simulated detection random effects for each given level of the random effects included in the detection model. Only relevant when detection random effects are included in the model.}
  \item{beta.star}{a numeric matrix where each row contains the simulated occurrence random effects for each given level of the random effects included in the occurrence model. Only relevant when occurrence random effects are included in the model.}


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
psi.RE <- list(levels = c(10), 
               sigma.sq.psi = c(1.5))
p.RE <- list(levels = c(15), 
             sigma.sq.p = 0.8)
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

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, 
                alpha = alpha, psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, 
                cov.model = 'exponential', phi = phi, sigma.sq = sigma.sq)
}
