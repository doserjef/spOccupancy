\name{simTMsOcc}
\alias{simTMsOcc}
\title{Simulate Multi-Species Multi-Season Detection-Nondetection Data}

\description{
  The function \code{simTMsOcc} simulates multi-species multi-season detection-nondetection data for simulation studies, power assessments, or function testing. Data can be optionally simulated with a spatial Gaussian Process in the occurrence portion of the model, as well as an option to allow for species correlations using a factor modeling approach. Non-spatial random intercepts can also be included in the detection or occurrence portions of the occupancy model.
}

\usage{
simTMsOcc(J.x, J.y, n.time, n.rep, N, beta, alpha, sp.only = 0, 
	  trend = TRUE, psi.RE = list(), p.RE = list(), 
          sp = FALSE, svc.cols = 1, cov.model, 
	  sigma.sq, phi, nu, ar1 = FALSE, rho, sigma.sq.t, 
	  factor.model = FALSE, n.factors, range.probs, grid, ...)
}

\arguments{

\item{J.x}{a single numeric value indicating the number of sites to simulate detection-nondetection data along the horizontal axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{J.y}{a single numeric value indicating the number of sites to simulate detection-nondetection data along the vertical axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{n.time}{a single numeric value indicating the number of primary time periods (denoted T) over which sampling occurs.}

\item{n.rep}{a numeric matrix indicating the number of replicates at each site during each primary time period. The matrix must have \eqn{J = J.x \times J.y}{J = J.x x J.y} rows and T columns, where T is the number of primary time periods (e.g., years or seasons) over which sampling occurs.} 

\item{N}{a single numeric value indicating the number of species to simulate detection-nondetection data.}

\item{beta}{a numeric matrix with \eqn{N}{N} rows containing the intercept and regression coefficient parameters for the occurrence portion of the multi-species occupancy model. Each row corresponds to the regression coefficients for a given species.} 

\item{alpha}{a numeric matrix with \eqn{N}{N} rows containing the intercept and regression coefficient parameters for the detection portion of the multi-species occupancy model. Each row corresponds to the regression coefficients for a given species.}

\item{sp.only}{a numeric vector specifying which occurrence covariates should only vary over space and not over time. The numbers in the vector correspond to the elements in the vector of regression coefficients (\code{beta}). By default, all simulated occurrence covariates are assumed to vary over both space and time.}

\item{trend}{a logical value. If \code{TRUE}, a temporal trend will be used to simulate the detection-nondetection data and the second element of \code{beta} is assumed to be the trend parameter. If \code{FALSE} no trend is used to simulate the data and all elements of \code{beta} (except the first value which is the intercept) correspond to covariate effects.}

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

\item{ar1}{a logical value indicating whether to simulate a temporal random effect with an AR(1) process. By default, set to \code{FALSE}.}

\item{rho}{a vector of \code{N} values indicating the AR(1) temporal correlation parameter for each species. Ignored when \code{ar1 = FALSE}.}

\item{sigma.sq.t}{a vector of \code{N} values indicating the AR(1) temporal variance parameter for each species. Ignored when \code{ar1 = FALSE}.}

\item{factor.model}{a logical value indicating whether to simulate data following a factor modeling approach that explicitly incoporates species correlations. If \code{sp = TRUE}, the latent factors are simulated from independent spatial processes. If \code{sp = FALSE}, the latent factors are simulated from standard normal distributions.}

\item{n.factors}{a single numeric value specifying the number of latent factors to use to simulate the data if \code{factor.model = TRUE}.}

\item{range.probs}{a numeric vector of length \code{N} where each value should fall between 0 and 1, and indicates the probability that one of the \code{J} spatial locations simulated is within the simulated range of the given species. If set to 1, every species has the potential of being present at each location.}

\item{grid}{an atomic vector used to specify the grid across which to simulate the latent spatial processes. This argument is used to simulate the underlying spatial processes at a different resolution than the coordinates (e.g., if coordinates are distributed across a grid).}

\item{...}{currently no additional arguments}
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
}

\value{
  A list comprised of: 

  \item{X}{a \eqn{J \times T \times p.occ}{J x T x p.occ} numeric array containing the design matrix for the occurrence portion of the occupancy model.}
  \item{X.p}{a four-dimensional numeric array with dimensions corresponding to sites, primary time periods, repeat visits, and number of detection regression coefficients. This is the design matrix used for the detection portion of the occupancy model.}
  \item{coords}{a \eqn{J \times 2}{J x 2} numeric matrix of coordinates of each occupancy site. Required for spatial models.}
  \item{w}{a \eqn{N \times J}{N x J} matrix of the spatial random effects for each species. Only used to simulate data when \code{sp = TRUE}. If \code{factor.model = TRUE}, the first dimension is \code{n.factors}.}
  \item{psi}{a \eqn{N \times J \times T}{N x J x T} array of the occurrence probabilities for each species at each site during each primary time period.}
  \item{z}{a \eqn{N \times J \times T}{N x J x T} array of the latent occurrence status for each species at each site during each primary time period.}
  \item{p}{a \code{N x J x T x max(n.rep)} array of the detection probabilities for each species at each site, primary time period, and secondyary replicate combination. Sites with fewer than \code{max(n.rep)} replicates will contain \code{NA} values.}
  \item{y}{a \code{N x J x T x max(n.rep)} array of the raw detection-nondetection data for each species at each site, primary time period, and replicate combination. Sites with fewer than \code{max(n.rep)} replicates will contain \code{NA} values.}
  \item{X.p.re}{a four-dimensional numeric array containing the levels of any detection random effect included in the model. Only relevant when detection random effects are specified in \code{p.RE}.}
  \item{X.re}{a numeric matrix containing the levels of any occurrence random effect included in the model. Only relevant when occurrence random effects are specified in \code{psi.RE}.}
  \item{alpha.star}{a numeric matrix where each row contains the simulated detection random effects for each given level of the random effects included in the detection model. Only relevant when detection random effects are included in the model.}
  \item{beta.star}{a numeric matrix where each row contains the simulated occurrence random effects for each given level of the random effects included in the occurrence model. Only relevant when occurrence random effects are included in the model.}
  \item{eta}{a numeric matrix with each row corresponding to species and column corresponding to time period of the AR(1) temporal random effects.}

}

\examples{
# Simulate Data -----------------------------------------------------------
set.seed(500)
J.x <- 8
J.y <- 8
J <- J.x * J.y
# Years sampled
n.time <- sample(3:10, J, replace = TRUE)
# n.time <- rep(10, J)
n.time.max <- max(n.time)
# Replicates
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(2:4, n.time[j], replace = TRUE)
  # n.rep[j, 1:n.time[j]] <- rep(4, n.time[j])
}
N <- 7
# Community-level covariate effects
# Occurrence
beta.mean <- c(-3, -0.2, 0.5)
trend <- FALSE
sp.only <- 0
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 1.5, 1.4)
# Detection
alpha.mean <- c(0, 1.2, -1.5)
tau.sq.alpha <- c(1, 0.5, 2.3)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list()
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
sp <- TRUE
svc.cols <- c(1, 2)
p.svc <- length(svc.cols)
n.factors <- 3
phi <- runif(p.svc * n.factors, 3 / .9, 3 / .3)
factor.model <- TRUE
cov.model <- 'exponential'

dat <- simTMsOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, N = N,
		 beta = beta, alpha = alpha, sp.only = sp.only, trend = trend,
		 psi.RE = psi.RE, p.RE = p.RE, factor.model = factor.model,
                 svc.cols = svc.cols, n.factors = n.factors, phi = phi, sp = sp,
                 cov.model = cov.model)
str(dat)
}
