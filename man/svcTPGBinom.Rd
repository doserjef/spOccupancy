\name{svcTPGBinom}
\alias{svcTPGBinom}
\title{Function for Fitting Multi-Season Single-Species Spatially-Varying Coefficient Binomial Models Using Polya-Gamma Latent Variables}

\description{
  The function \code{svcTPGBinom} fits multi-season single-species spatially-varying coefficient binomial models using Polya-Gamma latent variables. Models are fit using Nearest Neighbor Gaussian Processes. 
}

\usage{
svcTPGBinom(formula, data, inits, priors, 
            tuning, svc.cols = 1, cov.model = 'exponential', NNGP = TRUE, 
            n.neighbors = 15, search.type = 'cb', n.batch, 
            batch.length, accept.rate = 0.43, n.omp.threads = 1, 
            verbose = TRUE, ar1 = FALSE, n.report = 100, 
            n.burn = round(.10 * n.batch * batch.length), 
            n.thin = 1, n.chains = 1, 
            k.fold, k.fold.threads = 1, k.fold.seed = 100, 
            k.fold.only = FALSE, ...)
}

\arguments{

  \item{formula}{a symbolic description of the model to be fit using R's model syntax. Only
    right-hand side of formula is specified. See example below. Random intercepts
    are allowed using lme4 syntax (Bates et al. 2015).}


\item{data}{a list containing data necessary for model fitting.
  Valid tags are \code{y}, \code{covs}, \code{weights}, and \code{coords}. \code{y}
  is a two-dimensional array with the rows corresponding to the number of
  sites (\eqn{J}{J}) and columns corresponding to the maximum number of
  primary time periods (i.e., years or seasons). \code{covs} is a
  list of variables included in the occurrence portion of the model. Each
  list element is a different occurrence covariate, which can be site level
  or site/primary time period level. Site-level covariates are specified as a vector of
  length \eqn{J}{J} while site/primary time period level covariates are specified as a matrix
  with rows corresponding to sites and columns correspond to primary time periods.
  \code{weights} is a site by time period matrix containing the binomial weights 
  (i.e., the total number of Bernoulli trials) at each site/time period combination. Note
  that missing values are allowed and should be specified as NA. 
  \code{coords} is a \eqn{J \times 2}{J x 2}
  matrix of the observation coordinates. Note that \code{spOccupancy} assumes
  coordinates are specified in a projected coordinate system.}

  \item{inits}{a list with each tag corresponding to a parameter name.
    Valid tags are \code{beta}, \code{sigma.sq}, 
    \code{phi}, \code{w}, \code{nu}, \code{sigma.sq.psi}, \code{sigma.sq.t}, and \code{rho}.  
    \code{nu} is only specified if \code{cov.model = "matern"}, and \code{sigma.sq.psi}
    is only specified if there are random effects in \code{formula}. 
    \code{sigma.sq.t} and \code{rho} are only relevant when \code{ar1 = TRUE}.
    The value portion of each tag is the parameter's initial value. See \code{priors}
    description for definition of each parameter name.
    Additionally, the tag \code{fix} can be set to \code{TRUE} 
    to fix the starting values across all chains. If \code{fix} is not specified
    (the default), starting values are varied randomly across chains.}
  
  \item{priors}{a list with each tag corresponding to a parameter name. 
    Valid tags are \code{beta.normal}, \code{phi.unif}, 
    \code{sigma.sq.ig}, \code{sigma.sq.unif}, \code{nu.unif},
    \code{sigma.sq.psi.ig}, \code{sigma.sq.t.ig}, and \code{rho.unif}. 
    Regression coefficients (\code{beta})
    are assumed to follow a normal distribution. The hyperparameters of the 
    normal distribution are passed as a list of length two with the first
    and second elements corresponding to the mean and variance of the normal
    distribution, which are each specified as vectors of 
    length equal to the number of coefficients to be estimated or of length
    one if priors are the same for all coefficients. If not
    specified, prior means are set to 0 and prior variances set to 2.73. The 
    spatial variance parameter, \code{sigma.sq}, for each spatially-varying 
    coefficient is assumed to follow an 
    inverse-Gamma distribution or a uniform distribution (default is 
    inverse-Gamma). The spatial decay \code{phi} and 
    smoothness \code{nu} parameters are assumed to follow Uniform 
    distributions. The hyperparameters of the inverse-Gamma for \code{sigma.sq}
    are passed as a list with two elements corresponding to the 
    shape and scale parametters, respetively, with each element comprised of 
    a vector equal to the number of spatially-varying coefficients to be estimated
    or of length one if priors are the same for all coefficients.
    The hyperparameters of any uniform priors are also passed as a list of 
    length two with the first and second elements corresponding to 
    the lower and upper support, respectively, which can be passed as
    a vector equal to the total number of spatially-varying coefficients to 
    be estimated or of length one if priors are the same for all coefficients.
    \code{sigma.sq.psi} are the random effect variances for any random effects, 
    respectively, and are assumed to follow an 
    inverse-Gamma distribution. The hyperparameters of
    the inverse-Gamma distribution are passed as a list of length two with the 
    first and second elements corresponding to the shape and scale parameters, 
    respectively, which are each specified as vectors of length equal to the 
    number of random intercepts or of length one if priors are the same for all
    random effect variances. \code{sigma.sq.t} and 
    \code{rho} are the AR(1) variance and correlation parameters for the AR(1) zero-mean
    temporal random effects, respectively. 
    \code{sigma.sq.t} is assumed to follow an inverse-Gamma
    distribution, where the hyperparameters are specified as a vector with elements
    corresponding to the shape and scale parameters, respectively. \code{rho} is 
    assumed to follow a uniform distribution, where the hyperparameters are specified in 
    a vector of length two with elements corresponding to the lower and upper bounds of
    the uniform prior.}

  \item{svc.cols}{a vector indicating the variables whose effects will be 
    estimated as spatially-varying coefficients. \code{svc.cols} can be an 
    integer vector with values indicating the order of covariates specified 
    in the model formula (with 1 being the intercept if specified), or it can 
    be specified as a character vector with names corresponding to variable 
    names in \code{covs} (for the intercept, use '(Intercept)').} 
 
  \item{cov.model}{a quoted keyword that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.}

  \item{tuning}{a list with each tag corresponding to a parameter
    name. Valid tags are \code{phi}, \code{sigma.sq}, \code{nu}, and \code{rho}. 
    The value portion of each tag defines the initial variance of the Adaptive sampler. See
    Roberts and Rosenthal (2009) for details.}
  
  \item{NNGP}{if \code{TRUE}, model is fit with an NNGP. If \code{FALSE}, 
    a full Gaussian process is used. See Datta et al. (2016) and 
    Finley et al. (2019) for more information. Currently, only \code{NNGP = TRUE} is
    supported for multi-season occupancy models.}
  
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
    MCMC progress.}

  \item{n.burn}{the number of samples out of the total \code{n.batch * batch.length} 
    samples in each chain to discard as burn-in. By default, the first 
    10\% of samples is discarded.}
  
  \item{n.thin}{the thinning interval for collection of MCMC samples. The
    thinning occurs after the \code{n.burn} samples are discarded. Default 
    value is set to 1.}

  \item{n.chains}{the number of MCMC chains to run.}

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

\references{
  Doser, J. W., Finley, A. O., Saunders, S. P., Kery, M., Weed, A. S., &
  Zipkin, E. F. (2024A). Modeling complex species-environment
  relationships through spatially-varying coefficient occupancy models.
  Journal of Agricultural, Biological and Environmental Statistics.
  \doi{10.1007/s13253-023-00595-6}.

  Doser, J. W., Kery, M., Saunders, S. P., Finley, A. O., Bateman, B. L.,
  Grand, J., Reault, S., Weed, A. S., & Zipkin, E. F. (2024B). Guidelines
  for the use of spatially varying coefficients in species distribution
  models. Global Ecology and Biogeography, 33(4), e13814.
  \doi{10.1111/geb.13814}.

  Bates, Douglas, Martin Maechler, Ben Bolker, Steve Walker (2015).
  Fitting Linear Mixed-Effects Models Using lme4. Journal of
  Statistical Software, 67(1), 1-48. \doi{10.18637/jss.v067.i01}.

  Datta, A., S. Banerjee, A.O. Finley, and A.E. Gelfand. (2016)
  Hierarchical Nearest-Neighbor Gaussian process models for large
  geostatistical datasets. \emph{Journal of the American Statistical
    Association}, \doi{10.1080/01621459.2015.1044091}.
  
  Finley, A.O., A. Datta, B.D. Cook, D.C. Morton, H.E. Andersen, and
  S. Banerjee. (2019) Efficient algorithms for Bayesian Nearest Neighbor
  Gaussian Processes. \emph{Journal of Computational and Graphical
    Statistics}, \doi{10.1080/10618600.2018.1537924}.

  Polson, N.G., J.G. Scott, and J. Windle. (2013) Bayesian Inference for
  Logistic Models Using Polya-Gamma Latent Variables. 
  \emph{Journal of the American Statistical Association}, 108:1339-1349.

  Roberts, G.O. and Rosenthal J.S. (2009) Examples  of adaptive MCMC. 
  \emph{Journal of Computational and Graphical Statistics}, 18(2):349-367.
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{svcTPGBinom} that is a list comprised of: 

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the regression coefficients.}

  \item{y.rep.samples}{a three-dimensional array of posterior samples 
    for the fitted data values, with dimensions corresponding to posterior
    sample, site, and primary time period.}
  
  \item{psi.samples}{a three-dimensional array of posterior samples 
    for the occurrence probability values, with dimensions corresponding to posterior
    sample, site, and primary time period.}

  \item{theta.samples}{a \code{coda} object of posterior samples
    for spatial covariance parameters and temporal covariance parameters
    if \code{ar1 = TRUE}.}

  \item{w.samples}{a three-dimensional array of posterior samples for the 
        latent spatial random effects for all spatially-varying coefficients. 
	Dimensions correspond to MCMC sample, coefficient, and sites.}

  \item{sigma.sq.psi.samples}{a \code{coda} object of posterior samples
    for variances of unstructured random intercepts included in the 
    model. Only included if random intercepts are specified in 
    \code{formula}.}

  \item{beta.star.samples}{a \code{coda} object of posterior samples
    for the unstructured random effects. Only included if random intercepts
    are specified in \code{formula}.}

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

  \item{k.fold.deviance}{soring rule (deviance) from k-fold cross-validation. 
    Only included if \code{k.fold} is specified in function call.}

  The return object will include additional objects used for 
  subsequent prediction and/or model fit evaluation. 
  Note that if \code{k.fold.only = TRUE}, the
  return list object will only contain \code{run.time} and \code{k.fold.deviance}
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

# Get all the data
dat <- simTBinom(J.x = J.x, J.y = J.y, n.time = n.time, weights = weights, beta = beta, 
                 psi.RE = psi.RE, sp.only = sp.only, trend = trend, 
                 sp = sp, svc.cols = svc.cols, 
                 cov.model = cov.model, sigma.sq = sigma.sq, phi = phi,
                 rho = rho, sigma.sq.t = sigma.sq.t, ar1 = TRUE, x.positive = FALSE)

# Prep the data for spOccupancy -------------------------------------------
y <- dat$y
X <- dat$X
X.re <- dat$X.re
coords <- dat$coords

# Package all data into a list
covs <- list(int = X[, , 1],
             trend = X[, , 2],
             cov.1 = X[, , 3], 
             cov.2 = X[, , 4])
# Data list bundle
data.list <- list(y = y,
                  covs = covs,
                  weights = weights, 
                  coords = coords)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   sigma.sq.ig = list(a = 2, b = 1),
                   phi.unif = list(a = 3/1, b = 3/.1), 
                   sigma.sq.t.ig = c(2, 0.5), 
                   rho.unif = c(-1, 1))

# Starting values
inits.list <- list(beta = beta, alpha = 0,
                   sigma.sq = 1, phi = 3 / 0.5, 
                   sigma.sq.t = 0.5, rho = 0)
# Tuning
tuning.list <- list(phi = 0.4, nu = 0.3, rho = 0.2)

# MCMC settings
n.batch <- 2
n.burn <- 0
n.thin <- 1

# Note that this is just a test case and more iterations/chains may need to 
# be run to ensure convergence.
out <- svcTPGBinom(formula = ~ trend + cov.1 + cov.2, 
                   svc.cols = svc.cols,
                   data = data.list, 
                   n.batch = n.batch, 
                   batch.length = 25, 
                   inits = inits.list, 
                   priors = prior.list,
                   accept.rate = 0.43, 
                   cov.model = "exponential", 
                   ar1 = TRUE,
                   tuning = tuning.list, 
                   n.omp.threads = 1, 
                   verbose = TRUE, 
                   NNGP = TRUE, 
                   n.neighbors = 5,
                   n.report = 1, 
                   n.burn = n.burn, 
                   n.thin = n.thin, 
                   n.chains = 1) 
}
