\name{tIntPGOcc}
\alias{tIntPGOcc}
\title{Function for Fitting Multi-Season Single-Species Integrated Occupancy Models Using Polya-Gamma Latent Variables}

\usage{
tIntPGOcc(occ.formula, det.formula, data, inits, priors, tuning, 
          n.batch, batch.length, accept.rate = 0.43, n.omp.threads = 1, 
          verbose = TRUE, ar1 = FALSE, n.report = 100, 
          n.burn = round(.10 * n.batch * batch.length), n.thin = 1, n.chains = 1, 
          ...)
}

\description{
Function for fitting single-species multi-season integrated occupancy models using Polya-Gamma latent variables. Data integration is done using a joint likelihood framework, assuming distinct detection models for each data source that are each conditional on a single latent occurrence process. 
}

\arguments{
\item{occ.formula}{a symbolic description of the model to be fit
  for the occurrence portion of the model using R's model syntax. Only
  right-hand side of formula is specified. See example below. Random intercepts
  are allowed using lme4 syntax (Bates et al. 2015).}

 \item{det.formula}{a list of symbolic descriptions of the models to be fit
  for the detection portion of the model using R's model syntax for 
  each data set. Each element in the list is a formula for the detection 
  model of a given data set. Only right-hand side of formula is specified. 
  See example below. Random intercepts are allowed using lme4 syntax (Bates et al. 2015).}

  \item{data}{a list containing data necessary for model fitting.
  Valid tags are \code{y}, \code{occ.covs}, \code{det.covs}, \code{sites}, and \code{seasons}. 
  \code{y} is a list of three-dimensional arrays with first dimensional equal to the 
  number of sites surveyed in that data set, second dimension equal to the
  number of primary time periods (i.e., years or seasons), and third dimension 
  equal to the maximum number of replicate surveys at a site within a given season. 
  \code{occ.covs} is a list of variables included in the occurrence portion of the model. Each
  list element is a different occurrence covariate, which can be site level
  or site/primary time period level. Site-level covariates are specified as a vector of 
  length \eqn{J}{J} while site/primary time period level covariates are specified as a matrix
  with rows corresponding to sites and columns corresponding to primary time periods.
  \code{det.covs} is a list of variables included in the detection 
  portion of the model for each data source. \code{det.covs} should have the same 
  number of elements as \code{y}, where each element is itself a list. 
  Each element of the list for a given data source is a different detection covariate, 
  which can be site-level , site-season-level, or observation-level. Site-level covariates
  and site/primary time period level covariates are specified in the same manner as 
  \code{occ.covs}. Observation-level
  covariates are specified as a three-dimensional array with first dimension corresponding 
  to sites, second dimension corresponding to primary time period, and third dimension
  corresponding to replicate. \code{sites} is a list of site indices with number of 
  elements equal to the number of data sources being modeled. Each element contains a vector of length equal to the number of sites that specific data source contains. Each value in the vector indicates the corresponding site in \code{occ.covs} covariates that corresponds with the specific row of the detection-nondetection data for the data source. This is used to properly link sites across data sets. Similarly, \code{seasons} is a list of season indices with number of elements equal to the number of data sources being modeled. Each element contains a vector of length equal to the number of seasons that a specific data source is available for. This is used to properly link seasons across data sets. Each value in the vector indicates the corresponding season in \code{occ.covs} covariates that correspond with the specific column of the detection-nondetection data for the given data source. This is used to properly link seasons across data sets, which can have a differing number of seasons surveyed.}

\item{inits}{a list with each tag corresponding to a parameter name.
  Valid tags are \code{z}, \code{beta}, \code{alpha}, \code{sigma.sq.psi}, 
  \code{sigma.sq.p}, \code{sigma.sq.t}, and \code{rho}. The value portion of each tag is the 
  parameter's initial value. \code{sigma.sq.psi} and \code{sigma.sq.p} are 
  only relevant when including random effects in the occurrence and 
  detection portion of the occupancy model, respectively. \code{sigma.sq.t} and \code{rho}
  are only relevant when \code{ar1 = TRUE}. The tag \code{alpha} is a list comprised 
  of the initial values for the detection parameters for each data source. 
  Each element of the list should be a vector of initial values for all detection 
  parameters in the given data source or a single value for each data source to 
  assign all parameters for a given data source the same initial value. See \code{priors} 
  description for definition of each parameter name.
  Additionally, the tag \code{fix} can be set to \code{TRUE} 
  to fix the starting values across all chains. If \code{fix} is not specified
  (the default), starting values are varied randomly across chains.}

\item{priors}{a list with each tag corresponding to a parameter name. 
  Valid tags are \code{beta.normal}, \code{alpha.normal}, \code{sigma.sq.psi.ig}, 
  \code{sigma.sq.p.ig}, \code{sigma.sq.t.ig}, and \code{rho.unif}. 
  Occupancy (\code{beta}) and detection (\code{alpha}) 
  regression coefficients are assumed to follow a normal distribution.  
  For \code{beta} hyperparameters of the normal distribution are passed as a list of length two with the first and second elements corresponding to the mean and variance of the normal distribution, which are each specified as vectors of length equal to the number of coefficients to be estimated or of length one if priors are the same for all coefficients. For the detection coefficients \code{alpha}, the mean and variance hyperparameters are themselves passed in as lists, with each element of the list corresponding to the specific hyperparameters for the detection parameters in a given data source. If not specified, prior means are set to 0 and prior variances set to 2.72.
  \code{sigma.sq.psi} and \code{sigma.sq.p} are the random effect variances for any unstructured occurrence or 
  detection random effects, respectively, and are assumed to follow an inverse
  Gamma distribution. The hyperparameters of the inverse-Gamma distribution
  are passed as a list of length two with first and second elements corresponding
  to the shape and scale parameters, respectively, which are each specified as
  vectors of length equal to the number of random intercepts or of length one 
  if priors are the same for all random effect variances. \code{sigma.sq.t} and 
  \code{rho} are the AR(1) variance and correlation parameters for the AR(1) zero-mean
  temporal random effects, respectively. \code{sigma.sq.t} is assumed to follow an inverse-Gamma
  distribution, where the hyperparameters are specified as a vector with elements
  corresponding to the shape and scale parameters, respectively. \code{rho} is 
  assumed to follow a uniform distribution, where the hyperparameters are specified in 
  a vector of length two with elements corresponding to the lower and upper bounds of
  the uniform prior.}

\item{tuning}{a list with each tag corresponding to a parameter
  name. Valid tags are \code{rho}. The value portion of each
  tag defines the initial tuning variance of the Adaptive sampler. See
  Roberts and Rosenthal (2009) for details.}

\item{n.batch}{the number of MCMC batches in each chain to run for the Adaptive MCMC 
  sampler. See Roberts and Rosenthal (2009) for details.}

\item{batch.length}{the length of each MCMC batch in each chain to run for the Adaptive 
  MCMC sampler. See Roberts and Rosenthal (2009) for details.}

\item{accept.rate}{target acceptance rate for Adaptive MCMC. Default is 
  0.43. See Roberts and Rosenthal (2009) for details.}

\item{n.omp.threads}{a positive integer indicating the number of threads
  to use for SMP parallel processing within chains. This will have no impact on 
  model run times for non-spatial models. The package must be compiled for 
  OpenMP support. For most Intel-based machines, we recommend setting 
  \code{n.omp.threads} up to the number of hyperthreaded cores. Note, 
  \code{n.omp.threads} > 1 might not work on some systems. Currently only 
  relevant for spatial models.}

\item{verbose}{if \code{TRUE}, messages about data preparation, 
  model specification, and progress of the sampler are printed to the screen. 
  Otherwise, no messages are printed.}

\item{ar1}{logical value indicating whether to include an AR(1) zero-mean 
  temporal random effect in the model. If \code{FALSE}, the model is
  fit without an AR(1) temporal autocovariance structure. If \code{TRUE},
  an AR(1) random effect is included in the model to account for temporal 
  autocorrelation across the primary time periods.}

\item{n.report}{the interval to report MCMC progress. Note this is specified
  in terms of batches, not MCMC samples.}

\item{n.burn}{the number of samples out of the total \code{n.samples} to 
  discard as burn-in for each chain. By default, the first 10\% of samples is discarded.}

\item{n.thin}{the thinning interval for collection of MCMC samples. The
  thinning occurs after the \code{n.burn} samples are discarded. Default 
  value is set to 1.}

\item{n.chains}{the number of chains to run.}

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

  Bates, Douglas, Martin Maechler, Ben Bolker, Steve Walker (2015).
  Fitting Linear Mixed-Effects Models Using lme4. Journal of
  Statistical Software, 67(1), 1-48. \doi{10.18637/jss.v067.i01}.
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}
}

\value{
  An object of class \code{tIntPGOcc} that is a list comprised of: 

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the occupancy regression coefficients.}

  \item{alpha.samples}{a \code{coda} object of posterior samples
    for the detection regression coefficients for all data sources.}

  \item{z.samples}{a three-dimensional array of posterior samples 
    for the latent occupancy values, with dimensions corresponding to 
    posterior sample, site, and primary time period. Note this object
    will contain predicted occupancy values for sites/primary time periods
    that were not sampled.}

  \item{psi.samples}{a three-dimensional array of posterior samples
    for the latent occupancy probability values, with dimensions 
    corresponding to posterior sample, site, and primary time period. Note
    this object will contain predicted occupancy probabilities for
    sites/primary time periods that were not sampled.}

  \item{sigma.sq.psi.samples}{a \code{coda} object of posterior samples
    for variances of random intercepts included in the occupancy portion
    of the model. Only included if random intercepts are specified in 
    \code{occ.formula}.}

  \item{sigma.sq.p.samples}{a \code{coda} object of posterior samples
    for variances of random intercpets included in the detection portion 
    of the model. Includes random effect variances for all data sources. 
    Only included if random intercepts are specified in 
    \code{det.formula}.}

  \item{beta.star.samples}{a \code{coda} object of posterior samples
    for the occurrence random effects. Only included if random intercepts
    are specified in \code{occ.formula}.}

  \item{alpha.star.samples}{a \code{coda} object of posterior samples
    for the detection random effects in any of the data sources. 
    Only included if random intercepts are specified in at least one of the 
    individual data set detection formulas in \code{det.formula}.}

  \item{theta.samples}{a \code{coda} object of posterior samples for the 
    AR(1) variance (\code{sigma.sq.t}) and correlation (\code{rho}) 
    parameters. Only included if \code{ar1 = TRUE}.}

  \item{eta.samples}{a \code{coda} object of posterior samples for the 
    AR(1) random effects for each primary time period. Only included if 
    \code{ar1 = TRUE}.}
  
  \item{p.samples}{a list of four-dimensional arrays consisting of the posterior
    samples of detection probability for each data source. For each data source, the 
    dimensions of the four-dimensional array correspond to MCMC sample, site, season, 
    and replicate within season.}

  \item{like.samples}{a two-dimensional array of posterior samples
    for the likelihood values associated with each site and primary
    time period, for each individual data source. Used for calculating WAIC.}

  \item{rhat}{a list of Gelman-Rubin diagnostic values for some of the model
    parameters.}

  \item{ESS}{a list of effective sample sizes for some of the model parameters.}

  \item{run.time}{execution time reported using \code{proc.time()}.}

  The return object will include additional objects used for 
  subsequent prediction and/or model fit evaluation.   

}

\examples{
set.seed(332)

# Simulate Data -----------------------------------------------------------
# Number of locations in each direction. This is the total region of interest
# where some sites may or may not have a data source.
J.x <- 15 
J.y <- 15
J.all <- J.x * J.y
# Number of data sources.
n.data <- 3
# Sites for each data source.
J.obs <- sample(ceiling(0.2 * J.all):ceiling(0.4 * J.all), n.data, replace = TRUE)
# Maximum number of years for each data set
n.time.max <- c(4, 8, 10)
# Number of years each site in each data set is sampled
n.time <- list()
for (i in 1:n.data) {
  n.time[[i]] <- sample(1:n.time.max[i], J.obs[i], replace = TRUE)
}
# Replicates for each data source.
n.rep <- list()
for (i in 1:n.data) {
  n.rep[[i]] <- matrix(NA, J.obs[i], n.time.max[i])
  for (j in 1:J.obs[i]) {
    n.rep[[i]][j, sample(1:n.time.max[i], n.time[[i]][j], replace = FALSE)] <- 
      sample(1:4, n.time[[i]][j], replace = TRUE)
  }
}
# Total number of years across all data sets
n.time.total <- 10
# List denoting the specific years each data set was sampled during. 
data.seasons <- list()
for (i in 1:n.data) {
  data.seasons[[i]] <- sort(sample(1:n.time.total, n.time.max[i], replace = FALSE))
}

# Occupancy covariates
beta <- c(0, 0.4, 0.3)
trend <- TRUE
# Random occupancy covariates
psi.RE <- list(levels = c(20),
               sigma.sq.psi = c(0.6))
p.occ <- length(beta)
# Detection covariates
alpha <- list()
alpha[[1]] <- c(0, 0.2, -0.5)
alpha[[2]] <- c(-1, 0.5, 0.3, -0.8)
alpha[[3]] <- c(-0.5, 1)

p.RE <- list()
p.det.long <- sapply(alpha, length)
p.det <- sum(p.det.long)

# Simulate occupancy data.
dat <- simTIntOcc(n.data = n.data, J.x = J.x, J.y = J.y, J.obs = J.obs,
                  n.time = n.time, data.seasons = data.seasons, n.rep = n.rep,
                  beta = beta, alpha = alpha, trend = trend, 
                  psi.RE = psi.RE, p.RE = p.RE)

y <- dat$y
X <- dat$X.obs
X.re <- dat$X.re.obs
X.p <- dat$X.p
sites <- dat$sites

# Package all data into a list
occ.covs <- list(trend = X[, , 2], 
                 occ.cov.1 = X[, , 3], 
                 occ.factor.1 = X.re[, , 1])
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(det.cov.1.1 = X.p[[1]][, , , 2],
                      det.cov.1.2 = X.p[[1]][, , , 3])
det.covs[[2]] <- list(det.cov.2.1 = X.p[[2]][, , , 2],
                      det.cov.2.2 = X.p[[2]][, , , 3],
                      det.cov.2.3 = X.p[[2]][, , , 4])
det.covs[[3]] <- list(det.cov.3.1 = X.p[[3]][, , , 2])
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs,
                  sites = sites, seasons = data.seasons)

# Testing
occ.formula <- ~ trend + occ.cov.1 + (1 | occ.factor.1)
# Note that the names are not necessary.
det.formula <- list(f.1 = ~ det.cov.1.1 + det.cov.1.2,
                    f.2 = ~ det.cov.2.1 + det.cov.2.2 + det.cov.2.3,
                    f.3 = ~ det.cov.3.1)

# NOTE: this is a short run of the model, in reality we would run the 
#       model for much longer.
out <- tIntPGOcc(occ.formula = occ.formula,
                 det.formula = det.formula,
                 data = data.list,
                 n.batch = 3,
                 batch.length = 25, 
                 n.report = 1,
                 n.burn = 25,
                 n.thin = 1,
                 n.chains = 1)
summary(out)
}
