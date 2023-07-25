\name{tMsPGOcc}
\alias{tMsPGOcc}
\title{Function for Fitting Multi-Species Multi-Season Occupancy Models}

\usage{
tMsPGOcc(occ.formula, det.formula, data, inits, priors, tuning, 
         n.batch, batch.length, 
         accept.rate = 0.43, n.omp.threads = 1, 
         verbose = TRUE, ar1 = FALSE, n.report = 100, 
         n.burn = round(.10 * n.batch * batch.length), n.thin = 1, 
         n.chains = 1, ...)
}

\description{
  The function \code{tMsPGOcc} fits multi-species multi-season occupancy models using Polya-Gamma data augmentation.
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
    Valid tags are \code{y}, \code{occ.covs}, and \code{det.covs}. 
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
    dimension corresponding to replicate.}

  \item{inits}{a list with each tag corresponding to a parameter name.
    Valid tags are \code{alpha.comm}, \code{beta.comm}, \code{beta}, 
    \code{alpha}, \code{tau.sq.beta}, \code{tau.sq.alpha}, 
    \code{sigma.sq.psi}, \code{sigma.sq.p}, \code{z}, 
    \code{sigma.sq.t}, and \code{rho}. \code{sigma.sq.t} and \code{rho} are only
    relevant when \code{ar1 = TRUE}, and \code{sigma.sq.psi} and \code{sigma.sq.p} are
    only specified if random effects are included in \code{occ.formula} or 
    \code{det.formula}, respectively. The value portion of each tag is 
    the parameter's initial value. See \code{priors} description for definition
    of each parameter name. Additionally, the tag \code{fix} can be set to \code{TRUE} 
    to fix the starting values across all chains. If \code{fix} is not specified
    (the default), starting values are varied randomly across chains.}

  \item{priors}{a list with each tag corresponding to a parameter name. 
  Valid tags are \code{beta.comm.normal}, \code{alpha.comm.normal}, 
  \code{tau.sq.beta.ig}, \code{tau.sq.alpha.ig}, \code{sigma.sq.psi}, 
  \code{sigma.sq.p}, \code{sigma.sq.t.ig}, and \code{rho.unif}. 
  Community-level occurrence (\code{beta.comm}) and detection 
  (\code{alpha.comm}) regression coefficients are assumed to follow a 
  normal distribution. The hyperparameters of the normal distribution
  are passed as a list of length two with the first and second elements 
  corresponding to the mean and variance of the normal distribution, 
  which are each specified as vectors of length equal to the number of 
  coefficients to be estimated or of length one if priors are the same for 
  all coefficients. If not specified, prior means are set 
  to 0 and prior variances set to 2.72. By default, community-level variance parameters 
  for occupancy (\code{tau.sq.beta}) and detection (\code{tau.sq.alpha}) are 
  assumed to follow an inverse Gamma distribution. The hyperparameters of 
  the inverse gamma distribution are passed as a list of length two with 
  the first and second elements corresponding to the shape and scale parameters,
  which are each specified as vectors of length equal to the number of 
  coefficients to be estimated or a single value if priors are the same for all 
  parameters. If not specified, prior shape and scale 
  parameters are set to 0.1. \code{sigma.sq.t} and 
  \code{rho} are the AR(1) variance and correlation parameters for the AR(1) zero-mean
  temporal random effects, respectively. \code{sigma.sq.t} is assumed to follow an inverse-Gamma
  distribution, where the hyperparameters are specified as a list of length two with the 
  first and second elements corresponding to the shape and scale parameters, respectively, 
  which can each be specified as vector equal to the number of species in the model or a single value
  if the same prior is used for all species. \code{rho} is assumed to follow a 
  uniform distribution, where the hyperparameters are specified similarly as a list of length two
  with the first and second elements corresponding to the lower and upper bounds of the 
  uniform prior, which can each be specified as vector equal to the number of species in the 
  model or a single value if the same prior is used for all species. 
  \code{sigma.sq.psi} and \code{sigma.sq.p} are the random effect variances for any occurrence or
  detection random effects, respectively, and are assumed to follow an inverse
  Gamma distribution. The hyperparameters of the inverse-Gamma distribution
  are passed as a list of length two with first and second elements corresponding
  to the shape and scale parameters, respectively, which are each specified as
  vectors of length equal to the number of random intercepts or of length one 
  if priors are the same for all random effect variances.}

\item{tuning}{a list with each tag corresponding to a parameter
  name. Valid tags are \code{rho}. The value portion of each
  tag defines the initial tuning variance of the Adaptive sampler. See
  Roberts and Rosenthal (2009) for details.}

  \item{n.batch}{the number of MCMC batches  in each chain to run for the Adaptive MCMC 
    sampler. See Roberts and Rosenthal (2009) for details.}
  
  \item{batch.length}{the length of each MCMC batch to run for the Adaptive 
    MCMC sampler. See Roberts and Rosenthal (2009) for details.}
  
  \item{accept.rate}{target acceptance rate for Adaptive MCMC. Defaul is 
    0.43. See Roberts and Rosenthal (2009) for details.}
  
  \item{n.omp.threads}{a positive integer indicating
   the number of threads to use for SMP parallel processing. The package must
   be compiled for OpenMP support. For most Intel-based machines, we
   recommend setting \code{n.omp.threads} up to the number of
   hyperthreaded cores. Note, \code{n.omp.threads} > 1 might not
   work on some systems.}
 
  \item{verbose}{if \code{TRUE}, messages about data preparation, 
    model specification, and progress of the sampler are printed to the screen. 
    Otherwise, no messages are printed.}

  \item{ar1}{logical value indicating whether to include an AR(1) zero-mean 
    temporal random effect in the model for each species. If \code{FALSE}, the model is
    fit without an AR(1) temporal autocovariance structure. If \code{TRUE},
    a species-specific AR(1) random effect is included in the model to account for temporal 
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

  Polson, N.G., J.G. Scott, and J. Windle. (2013) Bayesian Inference for
  Logistic Models Using Polya-Gamma Latent Variables.
  \emph{Journal of the American Statistical Association}, 108:1339-1349.

  Roberts, G.O. and Rosenthal J.S. (2009) Examples  of adaptive MCMC.
  \emph{Journal of Computational and Graphical Statistics}, 18(2):349-367.

  Bates, Douglas, Martin Maechler, Ben Bolker, Steve Walker (2015).
  Fitting Linear Mixed-Effects Models Using lme4. Journal of
  Statistical Software, 67(1), 1-48. \doi{10.18637/jss.v067.i01}.

  Kery, M., & Royle, J. A. (2021). Applied hierarchical modeling in ecology: 
  Analysis of distribution, abundance and species richness in R and BUGS: 
  Volume 2: Dynamic and advanced models. Academic Press. Section 4.6.

}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{tMsPGOcc} that is a list comprised of: 

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
    for the species level AR(1) variance (\code{sigma.sq.t}) and 
    correlation (\code{rho}) parameters. Only included if \code{ar1 = TRUE}.}

  \item{eta.samples}{a three-dimensional array of posterior samples for 
    the species-specific AR(1) random effects for each primary time period.
    Dimensions correspond to MCMC sample, species, and primary time period.}

  \item{z.samples}{a four-dimensional array of posterior samples for
    the latent occurrence values for each species. Dimensions corresopnd to 
    MCMC sample, species, site, and primary time period.}

  \item{psi.samples}{a four-dimensional array of posterior samples for
    the latent occupancy probability values for each species. Dimensions correspond
    to MCMC sample, species, site, and primary time period.}

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
sp <- FALSE

dat <- simTMsOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, N = N,
                 beta = beta, alpha = alpha, sp.only = sp.only, trend = trend,
                 psi.RE = psi.RE, p.RE = p.RE, sp = sp)

y <- dat$y
X <- dat$X
X.p <- dat$X.p
X.re <- dat$X.re
X.p.re <- dat$X.p.re

occ.covs <- list(occ.cov.1 = X[, , 2],
                 occ.cov.2 = X[, , 3])
det.covs <- list(det.cov.1 = X.p[, , , 2],
                 det.cov.2 = X.p[, , , 3])

data.list <- list(y = y, occ.covs = occ.covs,
                  det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                   alpha.comm.normal = list(mean = 0, var = 2.72),
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
z.init <- apply(y, c(1, 2, 3), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha.comm = 0, beta.comm = 0, beta = 0,
                   alpha = 0, tau.sq.beta = 1, tau.sq.alpha = 1,
                   z = z.init)
# Tuning
tuning.list <- list(phi = 1)

# Number of batches
n.batch <- 5
# Batch length
batch.length <- 25
n.burn <- 25
n.thin <- 1
n.samples <- n.batch * batch.length

out <- tMsPGOcc(occ.formula = ~ occ.cov.1 + occ.cov.2,
                det.formula = ~ det.cov.1 + det.cov.2,
                data = data.list,
                inits = inits.list,
                n.batch = n.batch,
                batch.length = batch.length,
                accept.rate = 0.43,
                priors = prior.list,
                n.omp.threads = 1,
                verbose = TRUE,
                n.report = 1,
                n.burn = n.burn,
                n.thin = n.thin,
                n.chains = 1)

summary(out)
}

