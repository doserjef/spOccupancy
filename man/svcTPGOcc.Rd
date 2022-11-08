\name{svcTPGOcc}
\alias{svcTPGOcc}
\title{Function for Fitting Multi-Season Single-Species Spatially-Varying Coefficient Occupancy Models Using Polya-Gamma Latent Variables}

\usage{
svcTPGOcc(occ.formula, det.formula, data, inits, priors, 
          tuning, svc.cols = 1, cov.model = 'exponential', NNGP = TRUE, 
          n.neighbors = 15, search.type = 'cb', n.batch, 
          batch.length, accept.rate = 0.43, n.omp.threads = 1, 
          verbose = TRUE, ar1 = FALSE, n.report = 100, 
          n.burn = round(.10 * n.batch * batch.length), 
          n.thin = 1, n.chains = 1, k.fold, k.fold.threads = 1, 
          k.fold.seed = 100, k.fold.only = FALSE, ...)
}

\description{
  Function for fitting multi-season single-species spatially-varying coefficient occupancy models using Polya-Gamma latent variables. Models are fit using Nearest Neighbor Gaussian Processes. 
}

\arguments{
\item{occ.formula}{a symbolic description of the model to be fit
  for the occurrence portion of the model using R's model syntax. Only
  right-hand side of formula is specified. See example below. Random intercepts
  are allowed using lme4 syntax (Bates et al. 2015).}

\item{det.formula}{a symbolic description of the model to be fit
  for the detection portion of the model using R's model syntax. Only
  right-hand side of formula is specified. See example below. Random intercepts
  are allowed using lme4 syntax (Bates et al. 2015).}

\item{data}{a list containing data necessary for model fitting.
  Valid tags are \code{y}, \code{occ.covs}, \code{det.covs}, and \code{coords}. \code{y}
  is a three-dimensional array with first dimension equal to the number of 
  sites (\eqn{J}{J}), second dimension equal to the maximum number of 
  primary time periods (i.e., years or seasons), and third dimension equal 
  to the maximum number of replicates at a given site. \code{occ.covs} is a
  list of variables included in the occurrence portion of the model. Each
  list element is a different occurrence covariate, which can be site level
  or site/primary time period level. Site-level covariates are specified as a vector of 
  length \eqn{J}{J} while site/primary time period level covariates are specified as a matrix
  with rows corresponding to sites and columns correspond to primary time periods.
  Similarly, \code{det.covs} is a list of variables included in the detection
  portion of the model, with each list element corresponding to an 
  individual variable. In addition to site-level and/or site/primary time period-level, 
  detection covariates can also be observational-level. Observation-level covariates
  are specified as a three-dimensional array with first dimension corresponding to 
  sites, second dimension corresponding to primary time period, and third 
  dimension corresponding to replicate. \code{coords} is a \eqn{J \times 2}{J x 2} 
  matrix of the observation coordinates. Note that \code{spOccupancy} assumes 
  coordinates are specified in a projected coordinate system.}

\item{inits}{a list with each tag corresponding to a parameter name.
  Valid tags are \code{z}, \code{beta}, \code{alpha}, \code{sigma.sq}, \code{phi}, 
  \code{w}, \code{nu}, \code{sigma.sq.psi}, \code{sigma.sq.p}, \code{sigma.sq.t}, 
  \code{rho}. The value portion of each tag is the 
  parameter's initial value. \code{sigma.sq.psi} and \code{sigma.sq.p} are 
  only relevant when including random effects in the occurrence and 
  detection portion of the occupancy model, respectively. \code{nu} is only
  specified if \code{cov.model = "matern"}. \code{sigma.sq.t} and \code{rho} 
  are only relevant when \code{ar1 = TRUE}. See \code{priors} 
  description for definition of each parameter name.
  Additionally, the tag \code{fix} can be set to \code{TRUE} 
  to fix the starting values across all chains. If \code{fix} is not specified
  (the default), starting values are varied randomly across chains.}

\item{priors}{a list with each tag corresponding to a parameter name. 
  Valid tags are \code{beta.normal}, \code{alpha.normal}, \code{sigma.sq.psi.ig}, 
  \code{sigma.sq.p.ig}, \code{phi.unif}, \code{sigma.sq.ig}, \code{nu.unif}, 
  \code{sigma.sq.t.ig}, and \code{rho.unif}. 
  Occupancy (\code{beta}) and detection (\code{alpha}) 
  regression coefficients are assumed to follow a normal distribution. 
  The hyperparameters of the normal distribution are passed as a list of 
  length two with the first and second elements corresponding to the mean 
  and variance of the normal distribution, which are each specified as vectors of 
  length equal to the number of coefficients to be estimated or of length one
  if priors are the same for all coefficients. If not specified, prior means 
  are set to 0 and prior variances set to 2.72. \code{sigma.sq.psi} and 
  \code{sigma.sq.p} are the random effect variances for any occurrence or 
  detection random effects, respectively, and are assumed to follow an inverse
  Gamma distribution. The hyperparameters of the inverse-Gamma distribution
  are passed as a list of length two with first and second elements corresponding
  to the shape and scale parameters, respectively, which are each specified as
  vectors of length equal to the number of random intercepts or of length one 
  if priors are the same for all random effect variances. The spatial variance
  parameter, \code{sigma.sq}, is assumed to follow an inverse-Gamma distribution. 
  The spatial decay \code{phi} and smoothness \code{nu} parameters are assumed 
  to follow Uniform distributions. The hyperparameters of the inverse-Gamma for 
  \code{sigma.sq.ig} are passed as a list of length two, with the first and second
  elements corresponding to the shape and scale parameters, respectively, with each 
  element comprised of a vector equal to the number of spatially-varying coefficients to 
  be estimated or of length one if priors are the same for all coefficients. The 
  hyperparameters of the uniform are also passed as a list of length two
  with the first and second elements corresponding to the lower and upper support, 
  respectively, which can be passed as a vector equal to the number of spatially-varying coefficients
  to be estimated or of length one if priors are the same for all coefficients. \code{sigma.sq.t} and 
  \code{rho} are the AR(1) variance and correlation parameters for the AR(1) zero-mean
  temporal random effects, respectively. \code{sigma.sq.t} is assumed to follow an inverse-Gamma
  distribution, where the hyperparameters are specified as a vector with elements
  corresponding to the shape and scale parameters, respectively. \code{rho} is 
  assumed to follow a uniform distribution, where the hyperparameters are specified in 
  a vector of length two with elements corresponding to the lower and upper bounds of
  the uniform prior.}

  \item{tuning}{a list with each tag corresponding to a parameter
    name. Valid tags are \code{phi}, \code{sigma.sq}, \code{nu}, and \code{rho}. 
    The value portion of each tag defines the initial variance of the Adaptive sampler. See
    Roberts and Rosenthal (2009) for details.}

  \item{svc.cols}{a vector indicating the variables whose effects will be 
    estimated as spatially-varying coefficients. \code{svc.cols} can be an 
    integer vector with values indicating the order of covariates specified 
    in the model formula (with 1 being the intercept if specified), or it can 
    be specified as a character vector with names corresponding to variable 
    names in \code{occ.covs} (for the intercept, use '(Intercept)'). \code{svc.cols}
    default argument of 1 results in a spatial occupancy model analogous to 
    \code{stPGOcc} (assuming an intercept is included in the model).} 

\item{cov.model}{a quoted keyword that specifies the covariance
  function used to model the spatial dependence structure among the
  observations.  Supported covariance model key words are:
  \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
  \code{"gaussian"}.}

\item{NNGP}{if \code{TRUE}, model is fit with an NNGP. If \code{FALSE}, 
  a full Gaussian process is used. See Datta et al. (2016) and 
  Finley et al. (2019) for more information. Currently only \code{NNGP = TRUE} is
  supported for multi-season single-species occupancy models.}

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

\item{batch.length}{the length of each MCMC batch in each chain to run for the Adaptive 
  MCMC sampler. See Roberts and Rosenthal (2009) for details.}

\item{accept.rate}{target acceptance rate for Adaptive MCMC. Default is 
  0.43. See Roberts and Rosenthal (2009) for details.}

\item{n.omp.threads}{a positive integer indicating the number of threads
  to use for SMP parallel processing. The package must be compiled for 
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

\item{n.report}{the interval to report MCMC progress.}

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
  as described in Hooten and Hobbs (2015). For cross-validation in multi-season
  models, the data are split along the site dimension, such that each hold-out
  data set consists of a \code{J / k.fold} sites sampled over all primary time
  periods during which data are available at each given site. Cross-validation is performed
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
  Polson, N.G., J.G. Scott, and J. Windle. (2013) Bayesian Inference for
  Logistic Models Using Polya-Gamma Latent Variables. 
  \emph{Journal of the American Statistical Association}, 108:1339-1349.

  Bates, Douglas, Martin Maechler, Ben Bolker, Steve Walker (2015).
  Fitting Linear Mixed-Effects Models Using lme4. Journal of
  Statistical Software, 67(1), 1-48. \doi{10.18637/jss.v067.i01}.

  Kéry, M., & Royle, J. A. (2021). Applied hierarchical modeling in ecology: 
  Analysis of distribution, abundance and species richness in R and BUGS: 
  Volume 2: Dynamic and advanced models. Academic Press. Section 4.6.

  Finley, A. O., and  Banerjee, S. (2020). Bayesian spatially varying coefficient 
  models in the spBayes R package. \emph{Environmental Modelling and Software}, 
  125, 104608.

  Hooten, M. B., and Hobbs, N. T. (2015). A guide to Bayesian model 
  selection for ecologists. Ecological monographs, 85(1), 3-28.

  MacKenzie, D. I., J. D. Nichols, G. B. Lachman, S. Droege, 
  J. Andrew Royle, and C. A. Langtimm. 2002. Estimating Site Occupancy 
  Rates When Detection Probabilities Are Less Than One. Ecology 83: 2248-2255.
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{svcTPGOcc} that is a list comprised of: 

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the occupancy regression coefficients.}

  \item{alpha.samples}{a \code{coda} object of posterior samples
    for the detection regression coefficients.}

  \item{z.samples}{a three-dimensional array of posterior samples 
    for the latent occupancy values, with dimensions corresponding to 
    posterior sample, site, and primary time period.}

  \item{psi.samples}{a three-dimensional array of posterior samples
    for the latent occupancy probability values, with dimensions 
    corresponding to posterior sample, site, and primary time period.}

  \item{theta.samples}{a \code{coda} object of posterior samples
    for spatial covariance parameters and temporal covariance parameters
    if \code{ar1 = TRUE}.}

  \item{w.samples}{a three-dimensional array of posterior samples for the 
        latent spatial random effects for all spatially-varying coefficients. 
	Dimensions correspond to MCMC sample, coefficient, and sites.}

  \item{sigma.sq.psi.samples}{a \code{coda} object of posterior samples
    for variances of random intercepts included in the occupancy portion
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

  \item{eta.samples}{a \code{coda} object of posterior samples for the 
    AR(1) random effects for each primary time period. Only included if 
    \code{ar1 = TRUE}.}

  \item{like.samples}{a three-dimensional array of posterior samples
    for the likelihood values associated with each site and primary
    time period. Used for calculating WAIC.}

  \item{rhat}{a list of Gelman-Rubin diagnostic values for some of the model
    parameters.}

  \item{ESS}{a list of effective sample sizes for some of the model parameters.}

  \item{run.time}{execution time reported using \code{proc.time()}.}

  \item{k.fold.deviance}{scoring rule (deviance) from k-fold cross-validation. 
    Only included if \code{k.fold} is specified in function call.}

  The return object will include additional objects used for 
  subsequent prediction and/or model fit evaluation. Note that detection
  probability estimated values are not included in the model object, but can be 
  extracted using \code{fitted()}. Note that if \code{k.fold.only = TRUE}, the
  return list object will only contain \code{run.time} and \code{k.fold.deviance}.
}

\examples{
set.seed(1000)
# Sites
J.x <- 15
J.y <- 15
J <- J.x * J.y
# Years sampled
n.time <- sample(10, J, replace = TRUE)
n.time.max <- max(n.time)
# Replicates
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(4, n.time[j], replace = TRUE)
}
# Occurrence --------------------------
beta <- c(-2, -0.5, -0.2, 0.75)
trend <- TRUE 
sp.only <- 0
psi.RE <- list()
# Detection ---------------------------
alpha <- c(1, 0.7, -0.5)
p.RE <- list()
# Spatial parameters ------------------
sp <- TRUE
svc.cols <- c(1, 2, 3)
p.svc <- length(svc.cols)
cov.model <- "exponential"
sigma.sq <- runif(p.svc, 0.1, 1)
phi <- runif(p.svc, 3 / 1, 3 / 0.2)
rho <- 0.8
sigma.sq.t <- 1
ar1 <- TRUE	 
x.positive <- FALSE 

# Get all the data
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, 
               beta = beta, alpha = alpha, sp.only = sp.only, trend = trend, 
               psi.RE = psi.RE, p.RE = p.RE, 
               sp = sp, cov.model = cov.model, sigma.sq = sigma.sq, phi = phi, 
               svc.cols = svc.cols, ar1 = ar1, rho = rho, sigma.sq.t = sigma.sq.t, 
               x.positive = x.positive)

# Data summary ------------------------------------------------------------
apply(dat$psi, 2, mean)

# Prep the data for svcTPGOcc ---------------------------------------------
# Full data set 
y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re
coords <- dat$coords

# Package all data into a list
occ.covs <- list(int = X[, , 1], 
                 trend = X[, , 2], 
                 occ.cov.1 = X[, , 3], 
                 occ.cov.2 = X[, , 4]) 
# Detection
det.covs <- list(det.cov.1 = X.p[, , , 2], 
                 det.cov.2 = X.p[, , , 3]) 
# Data list bundle
data.list <- list(y = y, 
                  occ.covs = occ.covs,
                  det.covs = det.covs, 
                  coords = coords)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                   alpha.normal = list(mean = 0, var = 2.72),
                   phi.unif = list(a = 3/1, b = 3/.1)) 

# Starting values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(beta = 0, alpha = 0,
                   sigma.sq = 1, phi = 3 / 0.5,
                   z = z.init, nu = 1)
# Tuning
tuning.list <- list(phi = 0.4, nu = 0.3, rho = 0.5, sigma.sq = 0.5) 

# MCMC settings
n.batch <- 2 
n.burn <- 0 
n.thin <- 1

# Run the model
out <- svcTPGOcc(occ.formula = ~ trend + occ.cov.1 + occ.cov.2, 
                 det.formula = ~ det.cov.1 + det.cov.2, 
                 data = data.list,
                 inits = inits.list,
                 tuning = tuning.list,
                 priors = prior.list, 
                 cov.model = "exponential", 
                 svc.cols = svc.cols,
                 NNGP = TRUE, 
                 ar1 = TRUE,
                 n.neighbors = 5, 
                 n.batch = n.batch,
                 batch.length = 25,
                 verbose = TRUE, 
                 n.report = 25,
                 n.burn = n.burn, 
                 n.thin = n.thin,
                 n.chains = 1) 
}
