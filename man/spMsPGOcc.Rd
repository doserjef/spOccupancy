\name{spMsPGOcc}
\alias{spMsPGOcc}
\title{Function for Fitting Multi-Species Spatial Occupancy Models Using Polya-Gamma Latent Variables}

\description{
  The function \code{spMsPGOcc} fits multi-species spatial occupancy models using Polya-Gamma latent variables. Models can be fit using either a full Gaussian process or a Nearest Neighbor Gaussian Process for large data sets. 
}

\usage{
spMsPGOcc(occ.formula, det.formula, data, inits, priors, tuning, 
          cov.model = 'exponential', NNGP = TRUE, 
          n.neighbors = 15, search.type = 'cb', n.batch, 
          batch.length, accept.rate = 0.43, n.omp.threads = 1, 
          verbose = TRUE, n.report = 100, 
          n.burn = round(.10 * n.batch * batch.length), n.thin = 1, 
          n.chains = 1, k.fold, k.fold.threads = 1, k.fold.seed, 
          k.fold.only = FALSE, ...)
}

\arguments{

  \item{occ.formula}{a symbolic description of the model to be fit
    for the occurrence portion of the model using R's model syntax. Only
    right-hand side of formula is specified. See example below. Random 
    intercepts are allowed using lme4 syntax (Bates et al. 2015).}
  
  \item{det.formula}{a symbolic description of the model to be fit
    for the detection portion of the model using R's model syntax. Only
    right-hand side of formula is specified. See example below. Random 
    intercepts are allowed using lme4 syntax (Bates et al. 2015).}

  \item{data}{a list containing data necessary for model fitting.
    Valid tags are \code{y}, \code{occ.covs}, \code{det.covs}, \code{coords}. 
    \code{y} is a three-dimensional array with first dimension equal to the 
    number of species, second dimension equal to the number of sites, and 
    third dimension equal to the maximum number of replicates at a given site. 
    \code{occ.covs} is a matrix or data frame containing the variables 
    used in the occurrence portion of the model, with \eqn{J}{J} rows for 
    each column (variable). \code{det.covs} is a list of variables included 
    in the detection portion of the model. Each list element is a different 
    detection covariate, which can be site-level or observational-level. 
    Site-level covariates are specified as a vector of length \eqn{J}{J} while 
    observation-level covariates are specified as a matrix or data frame 
    with the number of rows equal to \eqn{J}{J} and number of columns equal 
    to the maximum number of replicates at a given site. \code{coords} is a 
    \eqn{J \times 2}{J x 2} matrix of the observation coordinates. 
    Note that \code{spOccupancy} assumes coordinates are specified 
    in a projected coordinate system.}

  \item{inits}{a list with each tag corresponding to a parameter name.
    Valid tags are \code{alpha.comm}, \code{beta.comm}, \code{beta}, 
    \code{alpha}, \code{tau.sq.beta}, \code{tau.sq.alpha}, \code{sigma.sq.psi}, 
    \code{sigma.sq.p}, \code{z}, \code{sigma.sq}, 
    \code{phi}, \code{w}, and \code{nu}. \code{nu} is only specified if 
    \code{cov.model = "matern"}, \code{sigma.sq.psi} is only specified if
    there are random intercepts in \code{occ.formula}, and \code{sigma.sq.p} is
    only specified if there are random intercpets in \code{det.formula}. The value portion of each tag is 
    the parameter's initial value. See \code{priors} description for definition
    of each parameter name. Additionally, the tag \code{fix} can be set to \code{TRUE} 
    to fix the starting values across all chains. If \code{fix} is not specified
    (the default), starting values are varied randomly across chains.}

  \item{priors}{a list with each tag corresponding to a parameter name. 
  Valid tags are \code{beta.comm.normal}, \code{alpha.comm.normal}, 
  \code{tau.sq.beta.ig}, \code{tau.sq.alpha.ig}, \code{phi.unif}, \code{sigma.sq.ig}, 
  \code{sigma.sq.unif}, \code{nu.unif}, \code{sigma.sq.psi}, \code{sigma.sq.p}. 
  Community-level occurrence (\code{beta.comm}) and detection 
  (\code{alpha.comm}) regression coefficients are assumed to follow a 
  normal distribution. The hyperparameters of the normal distribution
  are passed as a list of length two with the first and second elements 
  corresponding to the mean and variance of the normal distribution, 
  which are each specified as vectors of length equal to the number of 
  coefficients to be estimated or of length one if priors are the same for 
  all coefficients. If not specified, prior means are set 
  to 0 and prior variances set to 2.73. Community-level variance parameters 
  for occupancy (\code{tau.sq.beta}) and detection (\code{tau.sq.alpha}) are 
  assumed to follow an inverse Gamma distribution. The hyperparameters of 
  the inverse gamma distribution are passed as a list of length two with 
  the first and second elements corresponding to the shape and scale parameters,
  which are each specified as vectors of length equal to the number of 
  coefficients to be estimated or a single value if priors are the same for all 
  parameters. If not specified, prior shape and scale 
  parameters are set to 0.1. The species-specific spatial variance parameter, 
  \code{sigma.sq}, is assumed to follow an inverse-Gamma distribution or a uniform
  distribution (default is inverse-Gamma). \code{sigma.sq} of all species
  can also be fixed at its initial value by setting the prior value to \code{"fixed"}.
  The spatial decay \code{phi} and smoothness \code{nu} parameters are assumed
  to follow Uniform distributions. The hyperparameters of the inverse-Gamma are
  passed as a list of length two, with the list elements being vectors of length 
  N corresponding to the species-specific shape and scale parameters, respectively, 
  or a single value if the same value is assigned for all species. 
  The hyperparameters of the Uniform are also passed as a list with two elements, 
  with both elements being vectors of length N corresponding to the lower and 
  upper support, respectively, or as a single value if the same value is assigned
  for all species. \code{sigma.sq.psi} and 
  \code{sigma.sq.p} are the random effect variances for any occurrence or 
  detection random effects, respectively, and are assumed to follow an inverse
  Gamma distribution. The hyperparameters of the inverse-Gamma distribution
  are passed as a list of length two with first and second elements corresponding
  to the shape and scale parameters, respectively, which are each specified as
  vectors of length equal to the number of random intercepts or of length one 
  if priors are the same for all random effect variances.}

  \item{tuning}{a list with each tag corresponding to a parameter
    name. Valid tags are \code{phi} and \code{nu}. The value portion of each
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
    Finley et al. (2019) for more information.}
  
  \item{n.neighbors}{number of neighbors used in the NNGP. Only used if 
  \code{NNGP = TRUE}. Datta et al. (2016) showed that 15 neighbors is usually 
  sufficient, but that as few as 5 neighbors can be adequate for certain data
  sets, which can lead to even greater decreases in run time. We recommend
  starting with 15 neighbors (the default) and if additional gains in computation
  time are desired, subsequently compare the results with a smaller number
  of neighbors using WAIC or k-fold cross-validation.}
  
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

  \item{n.batch}{the number of MCMC batches in each chain to run for the Adaptive MCMC 
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
 
  \item{n.report}{the interval to report Metropolis sampler acceptance and
    MCMC progress. Note this is specified in terms of batches and not overall
    samples for spatial models.}

  \item{n.burn}{the number of samples out of the total \code{n.samples} to 
    discard as burn-in for each chain. By default, the first 10\% of samples is discarded.}
  
  \item{n.thin}{the thinning interval for collection of MCMC samples. The
    thinning occurs after the \code{n.burn} samples are discarded. Default 
    value is set to 1.}

  \item{n.chains}{the number of chains to run in sequence.}

  \item{k.fold}{specifies the number of \emph{k} folds for cross-validation. 
    If not specified as an argument, then cross-validation is not performed
    and \code{k.fold.threads} and \code{k.fold.seed} are ignored. In \emph{k}-fold
    cross-validation, the data specified in \code{data} is randomly
    partitioned into \emph{k} equal sized subsamples. Of the \emph{k} subsamples, 
    \emph{k} - 1 subsamples are used to fit the model and the remaining \emph{k}
    samples are used for prediction. The cross-validation process is repeated 
    \emph{k} times (the folds). As a scoring rule, we use the model deviance 
    as described in Hooten and Hobbs (2015). Cross-validation is performed
    after the full model is fit using all the data. Cross-validation results
    are reported in the \code{k.fold.deviance} object in the return list.}
  
  \item{k.fold.threads}{number of threads to use for cross-validation. If 
    \code{k.fold.threads > 1} parallel processing is accomplished using the 
    \pkg{foreach} and \pkg{doParallel} packages. Ignored if \code{k.fold}
    is not specified.} 
  
  \item{k.fold.seed}{seed used to split data set into \code{k.fold} parts
    for k-fold cross-validation. Ignored if \code{k.fold} is not specified.}

  \item{k.fold.only}{a logical value indicating whether to only perform 
    cross-validation (\code{TRUE}) or perform cross-validation after fitting 
    the full model (\code{FALSE}). Default value is \code{FALSE}.} 
  
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

}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{spMsPGOcc} that is a list comprised of: 

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
    for the species level covariance parameters.}

  \item{z.samples}{a three-dimensional array of posterior samples for
    the latent occurrence values for each species.}

  \item{psi.samples}{a three-dimensional array of posterior samples for
    the latent occupancy probability values for each species.}

  \item{w.samples}{a three-dimensional array of posterior samples for 
    the latent spatial random effects for each species.}

  \item{sigma.sq.psi.samples}{a \code{coda} object of posterior samples
    for variances of random intercepts included in the occurrence portion
    of the model. Only included if random intercepts are specified in 
    \code{occ.formula}.}

  \item{sigma.sq.p.samples}{a \code{coda} object of posterior samples
    for variances of random intercpets included in the detection portion 
    of the model. Only included if random intercepts are specified in 
    \code{det.formula}.}

  \item{alpha.star.samples}{a \code{coda} object of posterior samples
    for the detection random effects. Only included if random intercepts 
    are specified in \code{det.formula}.}

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

  \item{k.fold.deviance}{vector of scoring rules (deviance) from k-fold cross-validation. 
    A separate value is reported for each species. 
    Only included if \code{k.fold} is specified in function call.}

  The return object will include additional objects used for 
  subsequent prediction and/or model fit evaluation. Note that detection probability
  estimated values are not included in the model object, but can be extracted using \code{fitted()}.
}

\examples{
set.seed(400)

# Simulate Data -----------------------------------------------------------
J.x <- 7
J.y <- 7
J <- J.x * J.y
n.rep <- sample(2:4, size = J, replace = TRUE)
N <- 5
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, -0.15)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 0.3)
# Detection
alpha.mean <- c(0.5, 0.2, -.2)
tau.sq.alpha <- c(0.2, 0.3, 0.8)
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
phi <- runif(N, 3/1, 3/.4)
sigma.sq <- runif(N, 0.3, 3)
sp <- TRUE

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
                phi = phi, sigma.sq = sigma.sq, sp = TRUE, cov.model = 'exponential')

# Number of batches
n.batch <- 30
# Batch length
batch.length <- 25
n.samples <- n.batch * batch.length

y <- dat$y
X <- dat$X
X.p <- dat$X.p
coords <- as.matrix(dat$coords)

# Package all data into a list
occ.covs <- X[, 2, drop = FALSE]
colnames(occ.covs) <- c('occ.cov')
det.covs <- list(det.cov.1 = X.p[, , 2], 
                 det.cov.2 = X.p[, , 3])
data.list <- list(y = y, 
                  occ.covs = occ.covs,
                  det.covs = det.covs, 
                  coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72), 
                   alpha.comm.normal = list(mean = 0, var = 2.72), 
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                   tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                   phi.unif = list(a = 3/1, b = 3/.1), 
                   sigma.sq.ig = list(a = 2, b = 2)) 
# Initial values
inits.list <- list(alpha.comm = 0, 
                   beta.comm = 0, 
                   beta = 0, 
                   alpha = 0,
                   tau.sq.beta = 1, 
                   tau.sq.alpha = 1, 
                   phi = 3 / .5, 
                   sigma.sq = 2,
                   w = matrix(0, nrow = N, ncol = nrow(X)),
                   z = apply(y, c(1, 2), max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1) 

out <- spMsPGOcc(occ.formula = ~ occ.cov, 
                 det.formula = ~ det.cov.1 + det.cov.2, 
                 data = data.list,
                 inits = inits.list, 
                 n.batch = n.batch, 
                 batch.length = batch.length, 
                 accept.rate = 0.43, 
                 priors = prior.list, 
                 cov.model = "exponential", 
                 tuning = tuning.list, 
                 n.omp.threads = 1, 
                 verbose = TRUE, 
                 NNGP = TRUE, 
                 n.neighbors = 5, 
                 search.type = 'cb', 
                 n.report = 10, 
                 n.burn = 500, 
                 n.thin = 1, 
                 n.chains = 1)

summary(out, level = 'both')
}
