\name{stMsPGOcc}
\alias{stMsPGOcc}
\title{Function for Fitting Multi-Species Multi-Season Spatial Occupancy Models}

\description{
  The function \code{stMsPGOcc} fits multi-species multi-season spatial occupancy models with species correlations (i.e., a spatially-explicit joint species distribution model with imperfect detection). We use Polya-Gamma latent variables and a spatial factor modeling approach. Models are implemented using a Nearest Neighbor Gaussian Process.
}

\usage{
stMsPGOcc(occ.formula, det.formula, data, inits, priors, tuning, 
          cov.model = 'exponential', NNGP = TRUE, 
          n.neighbors = 15, search.type = 'cb', 
          n.factors, n.batch, batch.length, 
          accept.rate = 0.43, n.omp.threads = 1, 
          verbose = TRUE, ar1 = FALSE, n.report = 100, 
          n.burn = round(.10 * n.batch * batch.length), n.thin = 1, 
          n.chains = 1, ...)
}

\arguments{

  \item{occ.formula}{a symbolic description of the model to be fit
    for the occurrence portion of the model using R's model syntax. Random
    intercepts are allowed using lme4 syntax (Bates et al. 2015). Only
    right-hand side of formula is specified. See example below.}
  
  \item{det.formula}{a symbolic description of the model to be fit
    for the detection portion of the model using R's model syntax. Only
    right-hand side of formula is specified. See example below. Random 
    intercepts are allowed using lme4 syntax (Bates et al. 2015).}

  \item{data}{a list containing data necessary for model fitting.
    Valid tags are \code{y}, \code{occ.covs}, \code{det.covs}, 
    \code{coords}, and \code{grid.index}. 
    \code{y} is a four-dimensional array with first dimension equal to the 
    number of species, second dimension equal to the number of sites, third
    dimension equal to the number of primary time periods, and 
    fourth dimension equal to the maximum number of secondary replicates at a given site. 
    \code{occ.covs} is a list of variables included in the occurrence portion of the model. Each
    list element is a different occurrence covariate, which can be site level
    or site/primary time period level. Site-level covariates are specified as a vector of 
    length \eqn{J}{J} while site/primary time period level covariates are specified as a matrix
    with rows corresponding to sites and columns correspond to primary time periods.
    Similarly, \code{det.covs} is a list of variables included in the detection
    portion of the model, with each list element corresponding to an individual variable.    
    In addition to site-level and/or site/primary time period-level, 
    detection covariates can also be observational-level. Observation-level covariates
    are specified as a three-dimensional array with first dimension corresponding to 
    sites, second dimension corresponding to primary time period, and third 
    dimension corresponding to replicate. \code{coords} is a matrix of the observation coordinates used
    to estimate the SVCs for each site. \code{coords} has two columns for the 
    easting and northing coordinate, respectively. Typically, each site in the data 
    set will have it's own coordinate, such that \code{coords} is a \eqn{J \times 2}{J x 2} 
    matrix and \code{grid.index} should not be specified. 
    If you desire to estimate SVCs at some larger spatial level, 
    e.g., if points fall within grid cells and you want to estimate an SVC for 
    each grid cell instead of each point, \code{coords} can be specified as the coordinate for 
    each grid cell. In such a case, \code{grid.index} is an indexing vector of length J, where each 
    value of \code{grid.index} indicates the corresponding row in \code{coords} that the given 
    site corresponds to. Note that \code{spOccupancy} assumes coordinates are specified 
    in a projected coordinate system.}

  \item{inits}{a list with each tag corresponding to a parameter name.
    Valid tags are \code{alpha.comm}, \code{beta.comm}, \code{beta}, 
    \code{alpha}, \code{tau.sq.beta}, \code{tau.sq.alpha}, 
    \code{sigma.sq.psi}, \code{sigma.sq.p}, \code{z}, 
    \code{phi}, \code{lambda}, \code{nu}, \code{sigma.sq.t}, 
    and \code{rho}. \code{nu} is only specified if 
    \code{cov.model = "matern"}, and \code{sigma.sq.psi} and \code{sigma.sq.p} are
    only specified if random effects are included in \code{occ.formula} or 
    \code{det.formula}, respectively. \code{sigma.sq.t} and \code{rho} are only 
    relevant when \code{ar1 = TRUE}. The value portion of each tag is 
    the parameter's initial value. See \code{priors} description for definition
    of each parameter name. Additionally, the tag \code{fix} can be set to \code{TRUE} 
    to fix the starting values across all chains. If \code{fix} is not specified
    (the default), starting values are varied randomly across chains.}

  \item{priors}{a list with each tag corresponding to a parameter name. 
  Valid tags are \code{beta.comm.normal}, \code{alpha.comm.normal}, 
  \code{tau.sq.beta.ig}, \code{tau.sq.alpha.ig}, \code{sigma.sq.psi}, 
  \code{sigma.sq.p}, \code{phi.unif}, \code{nu.unif}, \code{sigma.sq.t.ig}, and 
  \code{rho.unif}. Community-level occurrence (\code{beta.comm}) and detection 
  (\code{alpha.comm}) regression coefficients are assumed to follow a 
  normal distribution. The hyperparameters of the normal distribution
  are passed as a list of length two with the first and second elements 
  corresponding to the mean and variance of the normal distribution, 
  which are each specified as vectors of length equal to the number of 
  coefficients to be estimated or of length one if priors are the same for 
  all coefficients. If not specified, prior means are set 
  to 0 and prior variances set to 2.73. By default, community-level variance parameters 
  for occupancy (\code{tau.sq.beta}) and detection (\code{tau.sq.alpha}) are 
  assumed to follow an inverse Gamma distribution. The hyperparameters of 
  the inverse gamma distribution are passed as a list of length two with 
  the first and second elements corresponding to the shape and scale parameters,
  which are each specified as vectors of length equal to the number of 
  coefficients to be estimated or a single value if priors are the same for all 
  parameters. If not specified, prior shape and scale 
  parameters are set to 0.1. 
  The spatial factor model fits \code{n.factors} independent spatial processes. 
  The spatial decay \code{phi} and smoothness \code{nu} parameters 
  for each latent factor are assumed to follow Uniform distributions. 
  The hyperparameters of the Uniform are passed as a list with two elements, 
  with both elements being vectors of length \code{n.factors} 
  corresponding to the lower and 
  upper support, respectively, or as a single value if the same value is assigned
  for all factor combinations. The priors for the factor loadings 
  matrix \code{lambda} are fixed
  following the standard spatial factor model to ensure parameter 
  identifiability (Christensen and Amemlya 2002). The
  upper triangular elements of the \code{N x n.factors} matrix are fixed at 0 and the 
  diagonal elements are fixed at 1. The lower triangular elements are assigned a 
  standard normal prior (i.e., mean 0 and variance 1).
  \code{sigma.sq.psi} and \code{sigma.sq.p} are the random 
  effect variances for any occurrence or
  detection random effects, respectively, and are assumed to follow an inverse
  Gamma distribution. The hyperparameters of the inverse-Gamma distribution
  are passed as a list of length two with first and second elements corresponding
  to the shape and scale parameters, respectively, which are each specified as
  vectors of length equal to the number of random intercepts or of length one 
  if priors are the same for all random effect variances. 
  parameters are set to 0.1. \code{sigma.sq.t} and \code{rho} 
  are the AR(1) variance and correlation parameters for the AR(1) zero-mean
  temporal random effects, respectively. \code{sigma.sq.t} is assumed to follow an inverse-Gamma
  distribution, where the hyperparameters are specified as a list of length two with the 
  first and second elements corresponding to the shape and scale parameters, respectively, 
  which can each be specified as vector equal to the number of species in the model or a single value
  if the same prior is used for all species. \code{rho} is assumed to follow a 
  uniform distribution, where the hyperparameters are specified similarly as a list of length two
  with the first and second elements corresponding to the lower and upper bounds of the 
  uniform prior, which can each be specified as vector equal to the number of species in the 
  model or a single value if the same prior is used for all species. 
}

  \item{tuning}{a list with each tag corresponding to a parameter
    name. Valid tags are \code{phi}, \code{nu}, \code{rho}. The value portion of each
    tag defines the initial variance of the adaptive sampler. We assume the 
    initial variance of the adaptive sampler is the same for each species, 
    although the adaptive sampler will adjust the tuning variances separately
    for each species. See Roberts and Rosenthal (2009) for details.}

  \item{cov.model}{a quoted keyword that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.}
  
  \item{NNGP}{if \code{TRUE}, model is fit with an NNGP. If \code{FALSE}, 
    a full Gaussian process is used. See Datta et al. (2016) and 
    Finley et al. (2019) for more information. Only 
    \code{NNGP = TRUE} is currently supported.}
  
  \item{n.neighbors}{number of neighbors used in the NNGP. Only used if 
  \code{NNGP = TRUE}. Datta et al. (2016) showed that 15 neighbors is usually 
  sufficient, but that as few as 5 neighbors can be adequate for certain data
  sets, which can lead to even greater decreases in run time. We recommend
  starting with 15 neighbors (the default) and if additional gains in computation
  time are desired, subsequently compare the results with a smaller number
  of neighbors using WAIC.}
  
  \item{search.type}{a quoted keyword that specifies the type of nearest
    neighbor search algorithm. Supported method key words are: \code{"cb"} and
    \code{"brute"}. The \code{"cb"} should generally be much
    faster. If locations do not have identical coordinate values on
    the axis used for the nearest neighbor ordering then \code{"cb"} 
    and \code{"brute"} should produce identical neighbor sets. 
    However, if there are identical coordinate values on the axis used 
    for nearest neighbor ordering, then \code{"cb"} and \code{"brute"} 
    might produce different, but equally valid, neighbor sets, 
    e.g., if data are on a grid. }

  \item{n.factors}{the number of factors to use in the spatial factor model approach. 
    Typically, the number of factors is set to be small (e.g., 4-5) relative to the 
    total number of species in the community, which will lead to substantial 
    decreases in computation time. However, the value can be anywhere 
    between 1 and N (the number of species in the community).}

  \item{n.batch}{the number of MCMC batches  in each chain to run for the Adaptive MCMC 
    sampler. See Roberts and Rosenthal (2009) for details.}
  
  \item{batch.length}{the length of each MCMC batch to run for the Adaptive 
    MCMC sampler. See Roberts and Rosenthal (2009) for details.}
  
  \item{accept.rate}{target acceptance rate for Adaptive MCMC. Defaul is 
    0.43. See Roberts and Rosenthal (2009) for details.}
  
  \item{n.omp.threads}{a positive integer indicating
   the number of threads to use for SMP parallel processing within chains. The package must
   be compiled for OpenMP support. For most Intel-based machines, we
   recommend setting \code{n.omp.threads} up to the number of
   hyperthreaded cores. Note, \code{n.omp.threads} > 1 might not
   work on some systems.}
 
  \item{verbose}{if \code{TRUE}, messages about data preparation, 
    model specification, and progress of the sampler are printed to the screen. 
    Otherwise, no messages are printed.}

  \item{ar1}{logical value indicating whether to include an AR(1) zero-mean
    temporal random effect in the model. If \code{FALSE}, the model is
    fit without an AR(1) temporal autocovariance structure. If \code{TRUE},
    an AR(1) random effect is included in the model to account for temporal
    autocorrelation across the primary time periods.}
 
  \item{n.report}{the interval to report Metropolis sampler acceptance and
    MCMC progress. Note this is specified in terms of batches and not overall
    samples for spatial models.}

  \item{n.burn}{the number of samples out of the total \code{n.samples} to 
    discard as burn-in for each chain. By default, the first 10\% of samples is discarded.}
  
  \item{n.thin}{the thinning interval for collection of MCMC samples. The
    thinning occurs after the \code{n.burn} samples are discarded. Default 
    value is set to 1.}

  \item{n.chains}{the number of chains to run in sequence.}

  \item{...}{currently no additional arguments}
}

\note{
 Some of the underlying code used for generating random numbers from the
 Polya-Gamma distribution is taken from the \pkg{pgdraw} package written
 by Daniel F. Schmidt and Enes Makalic. Their code implements Algorithm
 6 in PhD thesis of Jesse Bennett Windle (2013) \url{https://repositories.lib.utexas.edu/handle/2152/21842}.
}

\references{
  Datta, A., S. Banerjee, A.O. Finley, and A.E. Gelfand. (2016)
  Hierarchical Nearest-Neighbor Gaussian process models for large
  geostatistical datasets. \emph{Journal of the American Statistical
    Association}, \doi{10.1080/01621459.2015.1044091}.

  Finley, A.O., A. Datta, B.D. Cook, D.C. Morton, H.E. Andersen, and
  S. Banerjee. (2019) Efficient algorithms for Bayesian Nearest Neighbor
  Gaussian Processes. \emph{Journal of Computational and Graphical
    Statistics}, \doi{10.1080/10618600.2018.1537924}.

  Finley, A. O., Datta, A., and Banerjee, S. (2020). spNNGP R package 
  for nearest neighbor Gaussian process models. \emph{arXiv} preprint arXiv:2001.09111.

  Polson, N.G., J.G. Scott, and J. Windle. (2013) Bayesian Inference for
  Logistic Models Using Polya-Gamma Latent Variables.
  \emph{Journal of the American Statistical Association}, 108:1339-1349.

  Roberts, G.O. and Rosenthal J.S. (2009) Examples  of adaptive MCMC.
  \emph{Journal of Computational and Graphical Statistics}, 18(2):349-367.

  Bates, Douglas, Martin Maechler, Ben Bolker, Steve Walker (2015).
  Fitting Linear Mixed-Effects Models Using lme4. Journal of
  Statistical Software, 67(1), 1-48. \doi{10.18637/jss.v067.i01}.

  Hooten, M. B., and Hobbs, N. T. (2015). A guide to Bayesian model 
  selection for ecologists. \emph{Ecological Monographs}, 85(1), 3-28.

  Christensen, W. F., and Amemiya, Y. (2002). Latent variable analysis 
  of multivariate spatial data. \emph{Journal of the American Statistical Association}, 
  97(457), 302-317.

}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{stMsPGOcc} that is a list comprised of: 

  \item{beta.comm.samples}{a \code{coda} object of posterior samples
    for the community level occurrence regression coefficients.}
  
  \item{alpha.comm.samples}{a \code{coda} object of posterior samples
    for the community level detection regression coefficients.}

  \item{tau.sq.beta.samples}{a \code{coda} object of posterior samples
    for the occurrence community variance parameters.}

  \item{tau.sq.alpha.samples}{a \code{coda} object of posterior samples
    for the detection community variance parameters.}

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the species level occurrence regression coefficients.}

  \item{alpha.samples}{a \code{coda} object of posterior samples
    for the species level detection regression coefficients.}

  \item{theta.samples}{a \code{coda} object of posterior samples
    for the species level correlation parameters and the species-level
    temporal autocorrelation parameters.}

  \item{lambda.samples}{a \code{coda} object of posterior samples
    for the latent spatial factor loadings.}

  \item{z.samples}{a four-dimensional array of posterior samples for
    the latent occurrence values for each species. Dimensions corresopnd to 
    MCMC sample, species, site, and primary time period.}

  \item{psi.samples}{a four-dimensional array of posterior samples for
    the latent occupancy probability values for each species. Dimensions correspond
    to MCMC sample, species, site, and primary time period.}

  \item{w.samples}{a three-dimensional array of posterior samples for 
    the latent spatial random effects for each spatial factor. 
    Dimensions correspond to MCMC sample, factor, and site.}

  \item{sigma.sq.psi.samples}{a \code{coda} object of posterior samples
    for variances of random intercepts included in the occurrence portion
    of the model. Only included if random intercepts are specified in 
    \code{occ.formula}.}

  \item{sigma.sq.p.samples}{a \code{coda} object of posterior samples
    for variances of random intercpets included in the detection portion 
    of the model. Only included if random intercepts are specified in 
    \code{det.formula}.}

  \item{beta.star.samples}{a \code{coda} object of posterior samples
    for the occurrence random effects. Only included if random intercepts
    are specified in \code{occ.formula}.}

  \item{alpha.star.samples}{a \code{coda} object of posterior samples
    for the detection random effects. Only included if random intercepts 
    are specified in \code{det.formula}.}

  \item{like.samples}{a four-dimensional array of posterior samples
    for the likelihood value used for calculating WAIC. Dimensions correspond 
    to MCMC sample, species, site, and time period.}

  \item{rhat}{a list of Gelman-Rubin diagnostic values for some of the model
    parameters.}

  \item{ESS}{a list of effective sample sizes for some of the model parameters.}

  \item{run.time}{MCMC sampler execution time reported using \code{proc.time()}.}

  The return object will include additional objects used for 
  subsequent prediction and/or model fit evaluation. Note that detection 
  probability estimated values are not included in the model object, but can
  be extracted using \code{fitted()}.
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
svc.cols <- c(1)
p.svc <- length(svc.cols)
n.factors <- 3
phi <- runif(p.svc * n.factors, 3 / .9, 3 / .3)
factor.model <- TRUE
cov.model <- 'exponential'
ar1 <- TRUE
sigma.sq.t <- runif(N, 0.05, 1)
rho <- runif(N, 0.1, 1)

dat <- simTMsOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, N = N,
                 beta = beta, alpha = alpha, sp.only = sp.only, trend = trend,
                 psi.RE = psi.RE, p.RE = p.RE, factor.model = factor.model,
                 svc.cols = svc.cols, n.factors = n.factors, phi = phi, sp = sp,
                 cov.model = cov.model, ar1 = ar1, sigma.sq.t = sigma.sq.t, rho = rho)

y <- dat$y
X <- dat$X
X.p <- dat$X.p
coords <- dat$coords
X.re <- dat$X.re
X.p.re <- dat$X.p.re

occ.covs <- list(occ.cov.1 = X[, , 2],
                 occ.cov.2 = X[, , 3])
det.covs <- list(det.cov.1 = X.p[, , , 2],
                 det.cov.2 = X.p[, , , 3])

data.list <- list(y = y, occ.covs = occ.covs,
                  det.covs = det.covs,
                  coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                   alpha.comm.normal = list(mean = 0, var = 2.72),
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                   tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                   rho.unif = list(a = -1, b = 1),
                   sigma.sq.t.ig = list(a = 0.1, b = 0.1),
                   phi.unif = list(a = 3 / .9, b = 3 / .1))
z.init <- apply(y, c(1, 2, 3), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha.comm = 0, beta.comm = 0, beta = 0,
                   alpha = 0, tau.sq.beta = 1, tau.sq.alpha = 1,
                   rho = 0.5, sigma.sq.t = 0.5,
                   phi = 3 / .5, z = z.init)
# Tuning
tuning.list <- list(phi = 1, rho = 0.5)

# Number of batches
n.batch <- 5
# Batch length
batch.length <- 25
n.burn <- 25
n.thin <- 1
n.samples <- n.batch * batch.length

# Note that this is just a test case and more iterations/chains may need to 
# be run to ensure convergence.
out <- stMsPGOcc(occ.formula = ~ occ.cov.1 + occ.cov.2,
                 det.formula = ~ det.cov.1 + det.cov.2,
                 data = data.list,
                 inits = inits.list,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 accept.rate = 0.43,
                 ar1 = TRUE,
                 NNGP = TRUE,
                 n.neighbors = 5,
                 n.factors = n.factors,
                 cov.model = 'exponential',
                 priors = prior.list,
                 tuning = tuning.list,
                 n.omp.threads = 1,
                 verbose = TRUE,
                 n.report = 1,
                 n.burn = n.burn,
                 n.thin = n.thin,
                 n.chains = 1)

summary(out)
}
