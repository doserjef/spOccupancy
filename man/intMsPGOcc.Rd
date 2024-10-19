\name{intMsPGOcc}
\alias{intMsPGOcc}
\title{Function for Fitting Integrated Multi-Species Occupancy Models Using Polya-Gamma Latent Variables}


\description{
Function for fitting integrated multi-species occupancy models using Polya-Gamma latent variables. 
}

\usage{
intMsPGOcc(occ.formula, det.formula, data, inits, priors, n.samples,
           n.omp.threads = 1, verbose = TRUE, n.report = 100, 
           n.burn = round(.10 * n.samples), n.thin = 1, n.chains = 1,
           ...)
}

\arguments{

\item{occ.formula}{a symbolic description of the model to be fit
  for the occurrence portion of the model using R's model syntax. Only
  right-hand side of formula is specified. See example below. Random 
  intercepts are allowed using \pkg{lme4} syntax (Bates et al. 2015).}

 \item{det.formula}{a list of symbolic descriptions of the models to be fit
  for the detection portion of the model using R's model syntax for 
  each data set. Each element in the list is a formula for the detection 
  model of a given data set. Only right-hand side of formula is specified. 
  Random effects are not currently supported. See example below.}

\item{data}{a list containing data necessary for model fitting.
  Valid tags are \code{y}, \code{occ.covs}, \code{det.covs}, \code{sites},
  and \code{species}. \code{y} is a list of three-dimensional arrays. Each 
  element of the list has first dimension equal to the number of species
  observed in that data source, second dimension equal to the number of 
  sites observed in that data source, and thir dimension equal to the 
  maximum number of replicates at a given site. 
  \code{occ.covs} is a matrix or data frame containing the variables 
  used in the occurrence portion of the model, with the number of rows
  being the number of sites with at least one data source for each 
  column (variable). \code{det.covs} is a list of variables included
  in the detection portion of the model for each data source. \code{det.covs}
  should have the same number of elements as \code{y}, where each element
  is itself a list. Each element of the list for a given data source is a
  different detection covariate, which can be site-level or observational-level.
  Site-level covariates are specified as a vector with length equal to the
  number of observed sites of that data source, while observational-level
  covariates are specified as a matrix or data frame with the number of 
  rows equal to the number of observed sites of that data source and number
  of columns equal to the maximum number of replicates at a given site. \code{sites}
  is a list of site indices with number of elements equal to the number of data sources being 
  modeled. Each element contains a vector of length equal to the number of 
  sites that specific data source contains. Each value in the vector indicates
  the row in \code{occ.covs} that corresponds with the specific row of the 
  detection-nondetection data for the data source. This is used to properly link sites across
  data sets. \code{species} is a list with number of data sources being modeled. Each
  element of the list is a vector of codes (these can be numeric or character) that
  indicate the species modeled in the specific data set.}
  
\item{inits}{a list with each tag corresponding to a parameter name.
  Valid tags are \code{alpha.comm}, \code{beta.comm}, \code{beta}, 
  \code{alpha}, \code{tau.sq.beta}, \code{tau.sq.alpha}, \code{sigma.sq.psi}, 
  and \code{z}. The value portion of each tag is the parameter's initial value. See 
  \code{priors} description for definition of each parameter name.
  Additionally, the tag \code{fix} can be set to \code{TRUE} 
  to fix the starting values across all chains. If \code{fix} is not specified
  (the default), starting values are varied randomly across chains.}

\item{priors}{a list with each tag corresponding to a parameter name. 
  Valid tags are \code{beta.comm.normal}, \code{alpha.comm.normal}, 
  \code{tau.sq.beta.ig}, \code{tau.sq.alpha.ig}, \code{sigma.sq.psi.ig}, and 
  \code{sigma.sq.p.ig}. Community-level occurrence 
  (\code{beta.comm}) regression coefficients are assumed to 
  follow a normal distribution. The hyperparameters of the 
  normal distribution are passed as a list of length two with the first
  and second elements corresponding to the mean and variance of the normal
  distribution, which are each specified as vectors of 
  length equal to the number of coefficients to be estimated or of length one
  if priors are the same for all coefficients. If not
  specified, prior means are set to 0 and prior variances set to 2.72. For the 
  community-level detection means (\code{alpha.comm}), the mean and variance
  hyperparameters are themselves passed in as lists, with each element of the list
  corresponding to the specific hyperparameters for the detection parameters in a 
  given data source. If not specified, prior means are set to 0 and prior 
  variances set to 2.72. Community-level variance parameters for occurrence (\code{tau.sq.beta}) and
  detection (\code{tau.sq.alpha}) are assumed to follow an inverse Gamma 
  distribution. For the occurrence parameters, the hyperparameters of the inverse gamma distribution are 
  passed as a list of length two with the first and second elements 
  corresponding to the shape and scale parameters, which are each specified
  as vectors of length equal to the number of coefficients to be estimated 
  or a single value if all parameters are assigned the same prior. 
  If not specified, prior shape and scale parameters are set to 0.1. For the detection
  community-level variance parameters (\code{tau.sq.alpha}), the shape and scale 
  parameters are passed in as lists, with each element of the list corresponding 
  to the specific hyperparameters for the detection variances in a given 
  data source. \code{sigma.sq.psi} and are the random effect variances for any occurrence
  random effects, respectively, and are assumed to follow an inverse
  Gamma distribution. The hyperparameters of the inverse-Gamma distribution
  are passed as a list of length two with first and second elements corresponding
  to the shape and scale parameters, respectively, which are each specified as
  vectors of length equal to the number of random intercepts or of length one 
  if priors are the same for all random effect variances.}

\item{n.samples}{the number of posterior samples to collect in each chain.}

\item{n.omp.threads}{a positive integer indicating the number of threads
  to use for SMP parallel processing within chains. This will have no impact on 
  model run times for non-spatial models. The package must be compiled for 
  OpenMP support. For most Intel-based machines, we recommend setting 
  \code{n.omp.threads} up to the number of hypterthreaded cores. Note, 
  \code{n.omp.threads} > 1 might not work on some systems. Currently only 
  relevant for spatial models.}

\item{verbose}{if \code{TRUE}, messages about data preparation, 
  model specification, and progress of the sampler are printed to the screen. 
  Otherwise, no messages are printed.}

\item{n.report}{the interval to report MCMC progress.}

\item{n.burn}{the number of samples out of the total \code{n.samples} to 
  discard as burn-in for each chain. By default, the first 10\% of samples is discarded.}

\item{n.thin}{the thinning interval for collection of MCMC samples. The
  thinning occurs after the \code{n.burn} samples are discarded. Default 
  value is set to 1.}

\item{n.chains}{the number of chains to run in sequence.}

\item{...}{currently no additional arguments}

}

\references{
  Polson, N.G., J.G. Scott, and J. Windle. (2013) Bayesian Inference for
  Logistic Models Using Polya-Gamma Latent Variables.
  \emph{Journal of the American Statistical Association}, 108:1339-1349.

  Bates, Douglas, Martin Maechler, Ben Bolker, Steve Walker (2015).
  Fitting Linear Mixed-Effects Models Using lme4. Journal of
  Statistical Software, 67(1), 1-48. \doi{10.18637/jss.v067.i01}.

  Dorazio, R. M., and Royle, J. A. (2005). Estimating size and composition 
  of biological communities by modeling the occurrence of species. 
  Journal of the American Statistical Association, 100(470), 389-398.

  Doser, J. W., Leuenberger, W., Sillett, T. S., Hallworth, M. T. 
  & Zipkin, E. F. (2022). Integrated community occupancy models: A framework 
  to assess occurrence and biodiversity dynamics using multiple data 
  sources. Methods in Ecology and Evolution, 00, 1-14. \doi{10.1111/2041-210X.13811}
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
}

\value{
  An object of class \code{intMsPGOcc} that is a list comprised of: 

  \item{beta.comm.samples}{a \code{coda} object of posterior samples
    for the community level occurrence regression coefficients.}
  
  \item{alpha.comm.samples}{a \code{coda} object of posterior samples
    for the community level detection regression coefficients for all 
    data sources.}

  \item{tau.sq.beta.samples}{a \code{coda} object of posterior samples
    for the occurrence community variance parameters.}

  \item{tau.sq.alpha.samples}{a \code{coda} object of posterior samples
    for the detection community variance parameters for all data sources.}

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the species level occurrence regression coefficients.}

  \item{alpha.samples}{a \code{coda} object of posterior samples
    for the species level detection regression coefficients for all 
    data sources.}

  \item{z.samples}{a three-dimensional array of posterior samples for
    the latent occurrence values for each species.}

  \item{psi.samples}{a three-dimensional array of posterior samples for
    the latent occurrence probability values for each species.}

  \item{sigma.sq.psi.samples}{a \code{coda} object of posterior samples
    for variances of random intercepts included in the occurrence portion
    of the model. Only included if random intercepts are specified in 
    \code{occ.formula}.}

  \item{beta.star.samples}{a \code{coda} object of posterior samples
    for the occurrence random effects. Only included if random intercepts
    are specified in \code{occ.formula}.}

  \item{like.samples}{a three-dimensional array of posterior samples
    for the likelihood value associated with each site and species. 
    Used for calculating WAIC.}

  \item{rhat}{a list of Gelman-Rubin diagnostic values for some of the model
    parameters.}

  \item{ESS}{a list of effective sample sizes for some of the model parameters.}

  \item{run.time}{MCMC sampler execution time reported using \code{proc.time()}.}

  The return object will include additional objects used for 
  subsequent prediction and/or model fit evaluation. 
}

\examples{
set.seed(91)
J.x <- 10
J.y <- 10
# Total number of data sources across the study region
J.all <- J.x * J.y
# Number of data sources.
n.data <- 2
# Sites for each data source.
J.obs <- sample(ceiling(0.2 * J.all):ceiling(0.5 * J.all), n.data, replace = TRUE)
n.rep <- list()
n.rep[[1]] <- rep(3, J.obs[1])
n.rep[[2]] <- rep(4, J.obs[2])

# Number of species observed in each data source
N <- c(8, 3)

# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 0.5)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.4, 0.3)
# Detection
# Detection covariates
alpha.mean <- list()
tau.sq.alpha <- list()
# Number of detection parameters in each data source
p.det.long <- c(4, 3)
for (i in 1:n.data) {
  alpha.mean[[i]] <- runif(p.det.long[i], -1, 1)
  tau.sq.alpha[[i]] <- runif(p.det.long[i], 0.1, 1)
}
# Random effects
psi.RE <- list()
p.RE <- list()
beta <- matrix(NA, nrow = max(N), ncol = p.occ)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(max(N), beta.mean[i], sqrt(tau.sq.beta[i]))
}
alpha <- list()
for (i in 1:n.data) {
  alpha[[i]] <- matrix(NA, nrow = N[i], ncol = p.det.long[i])
  for (t in 1:p.det.long[i]) {
    alpha[[i]][, t] <- rnorm(N[i], alpha.mean[[i]][t], sqrt(tau.sq.alpha[[i]])[t])
  }
}
sp <- FALSE
factor.model <- FALSE
# Simulate occupancy data
dat <- simIntMsOcc(n.data = n.data, J.x = J.x, J.y = J.y,
                   J.obs = J.obs, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
                   psi.RE = psi.RE, p.RE = p.RE, sp = sp, factor.model = factor.model,
                   n.factors = n.factors)
J <- nrow(dat$coords.obs)
y <- dat$y
X <- dat$X.obs
X.p <- dat$X.p
X.re <- dat$X.re.obs
X.p.re <- dat$X.p.re
sites <- dat$sites
species <- dat$species

# Package all data into a list
occ.covs <- cbind(X)
colnames(occ.covs) <- c('int', 'occ.cov.1')
#colnames(occ.covs) <- c('occ.cov')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(det.cov.1.1 = X.p[[1]][, , 2], 
                      det.cov.1.2 = X.p[[1]][, , 3], 
                      det.cov.1.3 = X.p[[1]][, , 4])
det.covs[[2]] <- list(det.cov.2.1 = X.p[[2]][, , 2], 
                      det.cov.2.2 = X.p[[2]][, , 3]) 

data.list <- list(y = y, 
                  occ.covs = occ.covs, 
                  det.covs = det.covs, 
                  sites = sites, 
                  species = species)
# Take a look at the data.list structure for integrated multi-species
# occupancy models.
# Priors 
prior.list <- list(beta.comm.normal = list(mean = 0,var = 2.73),
                   alpha.comm.normal = list(mean = list(0, 0),
                                            var = list(2.72, 2.72)), 
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                   tau.sq.alpha.ig = list(a = list(0.1, 0.1), 
                                          b = list(0.1, 0.1)))
inits.list <- list(alpha.comm = list(0, 0), 
                   beta.comm = 0, 
                   tau.sq.beta = 1, 
                   tau.sq.alpha = list(1, 1), 
                   alpha = list(a = matrix(rnorm(p.det.long[1] * N[1]), N[1], p.det.long[1]), 
                                b = matrix(rnorm(p.det.long[2] * N[2]), N[2], p.det.long[2])),
                   beta = 0)

# Fit the model. 
# Note that this is just a test case and more iterations/chains may need to 
# be run to ensure convergence.
out <- intMsPGOcc(occ.formula = ~ occ.cov.1,
                  det.formula = list(f.1 = ~ det.cov.1.1 + det.cov.1.2 + det.cov.1.3,
                                     f.2 = ~ det.cov.2.1 + det.cov.2.2),
                  inits = inits.list,
                  priors = prior.list,
                  data = data.list, 
                  n.samples = 100, 
                  n.omp.threads = 1, 
                  verbose = TRUE, 
                  n.report = 10, 
                  n.burn = 50, 
                  n.thin = 1, 
                  n.chains = 1) 
summary(out, level = 'community')
}
