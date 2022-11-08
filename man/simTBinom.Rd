\name{simTBinom}
\alias{simTBinom}
\title{Simulate Multi-Season Single-Species Binomial Data}

\description{
  The function \code{simTBinom} simulates multi-season single-species binomial data for simulation studies, power assessments, or function testing. Data can be optionally simulated with a spatial Gaussian Process in the model. Non-spatial random intercepts can also be included in the model.
}

\usage{
simTBinom(J.x, J.y, n.time, weights, beta, sp.only = 0, 
          trend = TRUE, psi.RE = list(), sp = FALSE, 
          cov.model, sigma.sq, phi, nu, svc.cols = 1, 
          ar1 = FALSE, rho, sigma.sq.t, x.positive = FALSE, ...) 
}

\arguments{
\item{J.x}{a single numeric value indicating the number of sites to simulate data along the horizontal axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}.}

\item{J.y}{a single numeric value indicating the number of sites to simulate data along the vertical axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}.}

\item{n.time}{a single numeric value indicating the number of primary time periods (denoted T) over which sampling occurs.} 

\item{weights}{a numeric matrix with rows corresponding to sites and columns corresponding to primary time periods that indicates the number of Bernoulli trials at each of the site/time period combinations.}

\item{beta}{a numeric vector containing the intercept and regression coefficient parameters for the model.}

\item{sp.only}{a numeric vector specifying which occurrence covariates should only vary over space and not over time. The numbers in the vector correspond to the elements in the vector of regression coefficients (\code{beta}). By default, all simulated occurrence covariates are assumed to vary over both space and time.}

\item{trend}{a logical value. If \code{TRUE}, a temporal trend will be used to simulate the detection-nondetection data and the second element of \code{beta} is assumed to be the trend parameter. If \code{FALSE} no trend is used to simulate the data and all elements of \code{beta} (except the first value which is the intercept) correspond to covariate effects.}

\item{psi.RE}{a list used to specify the non-spatial random intercepts included in the model. The list must have two tags: \code{levels} and \code{sigma.sq.psi}. \code{levels} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the number of levels there are in each intercept. \code{sigma.sq.psi} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the variances for each random effect. If not specified, no random effects are included in the model.}

\item{sp}{a logical value indicating whether to simulate a spatially-explicit occupancy model with a Gaussian process. By default set to \code{FALSE}.}

\item{svc.cols}{a vector indicating the variables whose effects will be
  estimated as spatially-varying coefficients. \code{svc.cols} is an
  integer vector with values indicating the order of covariates specified
  in the model formula (with 1 being the intercept if specified).}

\item{cov.model}{a quoted keyword that specifies the covariance function used to model the spatial dependence structure among the latent occurrence values.  Supported covariance model key words are: \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and \code{"gaussian"}.}

\item{sigma.sq}{a numeric value indicating the spatial variance parameter. Ignored when \code{sp = FALSE}. If \code{svc.cols} has more than one value, there should be a distinct spatial variance parameter for each spatially-varying coefficient.}

\item{phi}{a numeric value indicating the spatial decay parameter. Ignored when \code{sp = FALSE}. If \code{svc.cols} has more than one value, there should be a distinct spatial decay parameter for each spatially-varying coefficient.}

\item{nu}{a numeric value indicating the spatial smoothness parameter. Only used when \code{sp = TRUE} and \code{cov.model = "matern"}. If \code{svc.cols} has more than one value, there should be a distinct spatial smoothness parameter for each spatially-varying coefficient.} 

\item{ar1}{a logical value indicating whether to simulate a temporal random effect with an AR(1) process. By default, set to \code{FALSE}.}

\item{rho}{a numeric value indicating the AR(1) temporal correlation parameter. Ignored when \code{ar1 = FALSE}.}

\item{sigma.sq.t}{a numeric value indicating the AR(1) temporal variance parameter. Ignored when \code{ar1 = FALSE}.}

\item{x.positive}{a logical value indicating whether the simulated covariates should be simulated as random standard normal covariates (\code{x.positive = FALSE}) or restricted to positive values (\code{x.positive = TRUE}). If \code{x.positive = TRUE}, covariates are simulated from a random normal and then the minimum value is added to each covariate value to ensure non-negative covariate values.}

\item{...}{currently no additional arguments}
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  A list comprised of: 

  \item{X}{a \eqn{J \times T \times p.occ}{J x T x p.occ} numeric array containing the design matrix for the model.}
  \item{coords}{a \eqn{J \times 2}{J x 2} numeric matrix of coordinates of each occupancy site. Required for spatial models.}
  \item{w}{a matrix of the spatial random effect values for each site. The number of columns is determined by the \code{svc.cols} argument (the number of spatially-varying coefficients).}
  \item{psi}{a \eqn{J \times T}{J x T} matrix of the occupancy probabilities for each site during each primary time period.}
  \item{z}{a  \eqn{J \times T}{J x T} matrix of the binomial data at each site during each primary time period.}
  \item{X.w}{a three dimensional array containing the covariate effects (including an intercept) whose effects are assumed to be spatially-varying. Dimensions correspond to sites, primary time periods, and covariate.}
  \item{X.re}{a numeric matrix containing the levels of any unstructured random effect included in the model. Only relevant when random effects are specified in \code{psi.RE}.}
  \item{beta.star}{a numeric vector that contains the simulated random effects for each given level of the random effects included in the model. Only relevant when random effects are included in the model.}
}

\examples{
set.seed(1000)
# Sites
J.x <- 15
J.y <- 15 
J <- J.x * J.y
# Years sampled
n.time <- sample(10, J, replace = TRUE)
# Binomial weights
weights <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  weights[j, 1:n.time[j]] <- sample(5, n.time[j], replace = TRUE)
}
# Occurrence --------------------------
beta <- c(-2, -0.5, -0.2, 0.75)
p.occ <- length(beta)
trend <- TRUE
sp.only <- 0
psi.RE <- list()
# Spatial parameters ------------------
sp <- TRUE
svc.cols <- c(1, 2, 3)
p.svc <- length(svc.cols)
cov.model <- "exponential"
sigma.sq <- runif(p.svc, 0.1, 1)
phi <- runif(p.svc, 3/1, 3/0.2)
# Temporal parameters -----------------
ar1 <- TRUE 
rho <- 0.8
sigma.sq.t <- 1

dat <- simTBinom(J.x = J.x, J.y = J.y, n.time = n.time, weights = weights, beta = beta, 
                 psi.RE = psi.RE, sp.only = sp.only, trend = trend, 
                 sp = sp, svc.cols = svc.cols, 
                 cov.model = cov.model, sigma.sq = sigma.sq, phi = phi,
                 rho = rho, sigma.sq.t = sigma.sq.t, ar1 = TRUE, x.positive = FALSE)
}
