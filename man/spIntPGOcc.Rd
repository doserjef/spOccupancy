\name{spIntPGOcc}
\alias{spIntPGOcc}
\title{Function for Fitting Single-Species Integrated Spatial Occupancy Models Using Polya-Gamma Latent Variables}

\description{
  The function \code{spIntPGOcc} fits single-species integrated spatial occupancy models using Polya-Gamma latent variables. Models can be fit using either a full Gaussian process or a Nearest Neighbor Gaussian Process for large data sets. Data integration is done using a joint likelihood framework, assuming distinct detection models for each data source that are each conditional on a single latent occupancy process. 
}

\usage{
spIntPGOcc(occ.formula, det.formula, data, inits, priors, 
           tuning, cov.model = "exponential", NNGP = TRUE, 
           n.neighbors = 15, search.type = 'cb', n.batch, 
           batch.length, accept.rate = 0.43, n.omp.threads = 1, 
           verbose = TRUE, n.report = 100, 
           n.burn = round(.10 * n.batch * batch.length), 
           n.thin = 1, n.chains = 1, k.fold, 
           k.fold.threads = 1, k.fold.seed, k.fold.data, 
           k.fold.only = FALSE, ...)
}

\arguments{

  \item{occ.formula}{a symbolic description of the model to be fit
  for the occurrence portion of the model using R's model syntax. Only
  right-hand side of formula is specified. See example below. Random 
  intercepts are allowed using lme4 syntax (Bates et al. 2015).}

 \item{det.formula}{a list of symbolic descriptions of the models to be fit
  for the detection portion of the model using R's model syntax for 
  each data set. Each element in the list is a formula for the detection 
  model of a given data set. Only right-hand side of formula is specified. 
  See example below. Random intercepts are allowed using lme4 syntax (Bates et al. 2015).}
 
  \item{data}{a list containing data necessary for model fitting.
  Valid tags are \code{y}, \code{occ.covs}, \code{det.covs}, \code{sites} and \code{coords}. 
  \code{y} is a list of matrices or data frames for each data set used in 
  the integrated model. Each element of the list has first dimension equal 
  to the number of sites with that data source and second dimension equal 
  to the maximum number of replicates at a given site. \code{occ.covs} is 
  a matrix or data frame containing the variables used in the occurrence 
  portion of the model, with the number of rows being the number of sites 
  with at least one data source for each column (variable). \code{det.covs} is a list of
  variables included in the detection portion of the model for each data source. 
  \code{det.covs} should have the same number of elements as \code{y}, where
  each element is itself a list. Each element of the list for a given data
  source is a different detection covariate, which can be site-level or 
  observational-level. Site-level covariates are specified as a vector 
  with length equal to the number of observed sites of that data source,
  while observation-level covariates are specified as a matrix or data frame 
  with the number of rows equal to the number of observed sites of that data source
  and number of columns equal to the maximum number of replicates at a given site. 
  \code{coords} is a matrix of the observation site coordinates. Note that \code{spOccupancy} assumes coordinates are specified 
  in a projected coordinate system. \code{sites}
  is a list of site indices with number of elements equal to the number of data sources being 
  modeled. Each element contains a vector of length equal to the number of 
  sites that specific data source contains. Each value in the vector indicates
  the row in \code{occ.covs} that corresponds with the specific row of the 
  detection-nondetection data for the data source. This is used to properly link sites across
  data sets.}

  \item{inits}{a list with each tag corresponding to a parameter name.
  Valid tags are \code{z}, \code{beta}, \code{alpha}, \code{sigma.sq}, \code{phi}, 
  \code{w}, \code{nu}, \code{sigma.sq.psi}, \code{sigma.sq.p}. 
  The value portion of all tags except \code{alpha} 
  is the parameter's initial value. \code{sigma.sq.psi} and \code{sigma.sq.p}
  are only relevant when including random effects in the occurrence and 
  detection portion of the occupancy model, respectively. The tag \code{alpha} is a list comprised 
  of the initial values for the detection parameters for each data source. 
  Each element of the list should be a vector of initial values for all detection 
  parameters in the given data source or a single value for each data source to 
  assign all parameters for a given data source the same initial value. See 
  \code{priors} description for definition of each parameter name.
  Additionally, the tag \code{fix} can be set to \code{TRUE} 
  to fix the starting values across all chains. If \code{fix} is not specified
  (the default), starting values are varied randomly across chains.}
  
  \item{priors}{a list with each tag corresponding to a parameter name. 
    Valid tags are \code{beta.normal}, \code{alpha.normal}, \code{phi.unif}, 
    \code{sigma.sq.ig}, \code{sigma.sq.unif}, \code{nu.unif}, \code{sigma.sq.psi.ig}, 
    and \code{sigma.sq.p.ig}. Occurrence (\code{beta}) and detection 
    (\code{alpha}) regression coefficients are assumed to follow a normal 
    distribution. For \code{beta} hyperparameters of the normal distribution 
    are passed as a list of length two with the first and second elements 
    corresponding to the mean and variance of the normal distribution, which are 
    each specified as vectors of length equal to the number of coefficients 
    to be estimated or of length one if priors are the same for all coefficients. 
    For the detection coefficients \code{alpha}, the mean 
    and variance hyperparameters are themselves passed in as lists, with each 
    element of the list corresponding to the specific hyperparameters for the 
    detection parameters in a given data source. If not specified, prior means 
    are set to 0 and prior variances set to 2.72 for normal priors. The 
    spatial variance parameter, \code{sigma.sq}, is assumed to follow an 
    inverse-Gamma distribution or a uniform distribution (default is 
    inverse-Gamma). \code{sigma.sq} can also be fixed at its 
    initial value by setting the prior value to 
    \code{"fixed"}. The spatial decay \code{phi} and 
    smoothness \code{nu} parameters are assumed to follow Uniform 
    distributions. The hyperparameters of the inverse-Gamma are passed as 
    a vector of length two, with the first and second elements corresponding
    to the \emph{shape} and \emph{scale}, respectively. The hyperparameters
    of the Uniform are also passed as a vector of length two with the first 
    and second elements corresponding to the lower and upper support, 
    respectively. \code{sigma.sq.psi} and 
    \code{sigma.sq.p} are the random effect variances for any occurrence or 
    detection random effects, respectively, and are assumed to follow an inverse
    Gamma distribution. The hyperparameters of the inverse-Gamma distribution
    are passed as a list of length two with first and second elements corresponding
    to the shape and scale parameters, respectively, which are each specified as
    vectors of length equal to the number of random intercepts or of length one 
    if priors are the same for all random effect variances.}

  \item{tuning}{a list with each tag corresponding to a parameter
    name. Valid tags are \code{phi} and \code{nu}. The value portion of each
    tag defines the initial variance of the Adaptive sampler. See
    Roberts and Rosenthal (2009) for details.}
  
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

  \item{n.batch}{the number of MCMC batches to run for each chain for the Adaptive MCMC 
    sampler. See Roberts and Rosenthal (2009) for details.}
  
  \item{batch.length}{the length of each MCMC batch to run for the Adaptive 
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
  
  \item{n.report}{the interval to report Metropolis sampler acceptance and
    MCMC progress. Note this is specified in terms of batches and not overall
    samples for spatial models.}

  \item{n.burn}{the number of samples out of the total \code{n.batch * batch.length} 
    samples to discard as burn-in. By default, the first 10\% of samples is discarded.}
  
  \item{n.thin}{the thinning interval for collection of MCMC samples. The
    thinning occurs after the \code{n.burn} samples are discarded. Default 
    value is set to 1.}

  \item{n.chains}{the number of chains to run.}

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

  \item{k.fold.data}{an integer specifying the specific data set to hold out
    values from. If not specified, data from all data set locations will 
    be incorporated into the k-fold cross-validation.}

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

  Hooten, M. B., and Hobbs, N. T. (2015). A guide to Bayesian model 
  selection for ecologists. \emph{Ecological Monographs}, 85(1), 3-28.

  Hooten, M. B., and Hefley, T. J. (2019). Bringing Bayesian models to life. 
  \emph{CRC Press}.

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
  An object of class \code{spIntPGOcc} that is a list comprised of: 

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the occurrence regression coefficients.}

  \item{alpha.samples}{a \code{coda} object of posterior samples
    for the detection regression coefficients for all data sources.}

  \item{z.samples}{a \code{coda} object of posterior samples 
    for the latent occurrence values}

  \item{psi.samples}{a \code{coda} object of posterior samples
    for the latent occurrence probability values}

  \item{theta.samples}{a \code{coda} object of posterior samples
    for covariance parameters.}

  \item{w.samples}{a \code{coda} object of posterior samples
    for latent spatial random effects.}

  \item{rhat}{a list of Gelman-Rubin diagnostic values for some of the model
    parameters.}

  \item{ESS}{a list of effective sample sizes for some of the model parameters.}

  \item{run.time}{execution time reported using \code{proc.time()}.}

  \item{k.fold.deviance}{scoring rule (deviance) from k-fold cross-validation. A 
    separate deviance value is returned for each data source. Only included if 
    \code{k.fold} is specified in function call. Only a single value is returned
    if \code{k.fold.data} is specified.}

  The return object will include additional objects used for 
  subsequent prediction and/or model fit evaluation. Note that detection
  probability estimated values are not included in the model object, but can be 
  extracted using \code{fitted()}.

}

\examples{
set.seed(400)

# Simulate Data -----------------------------------------------------------
# Number of locations in each direction. This is the total region of interest
# where some sites may or may not have a data source. 
J.x <- 8
J.y <- 8
J.all <- J.x * J.y
# Number of data sources.
n.data <- 4
# Sites for each data source. 
J.obs <- sample(ceiling(0.2 * J.all):ceiling(0.5 * J.all), n.data, replace = TRUE)
# Replicates for each data source.
n.rep <- list()
for (i in 1:n.data) {
  n.rep[[i]] <- sample(1:4, size = J.obs[i], replace = TRUE)
}
# Occupancy covariates
beta <- c(0.5, 0.5)
p.occ <- length(beta)
# Detection covariates
alpha <- list()
alpha[[1]] <- runif(2, 0, 1)
alpha[[2]] <- runif(3, 0, 1)
alpha[[3]] <- runif(2, -1, 1)
alpha[[4]] <- runif(4, -1, 1)
p.det.long <- sapply(alpha, length)
p.det <- sum(p.det.long)
sigma.sq <- 2
phi <- 3 / .5
sp <- TRUE

# Simulate occupancy data from multiple data sources. 
dat <- simIntOcc(n.data = n.data, J.x = J.x, J.y = J.y, J.obs = J.obs, 
                 n.rep = n.rep, beta = beta, alpha = alpha, sp = sp, 
                 sigma.sq = sigma.sq, phi = phi, cov.model = 'exponential')

y <- dat$y
X <- dat$X.obs
X.p <- dat$X.p
sites <- dat$sites
X.0 <- dat$X.pred
psi.0 <- dat$psi.pred
coords <- as.matrix(dat$coords.obs)
coords.0 <- as.matrix(dat$coords.pred)

# Package all data into a list
occ.covs <- X[, 2, drop = FALSE]
colnames(occ.covs) <- c('occ.cov')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(det.cov.1.1 = X.p[[1]][, , 2])
det.covs[[2]] <- list(det.cov.2.1 = X.p[[2]][, , 2], 
                      det.cov.2.2 = X.p[[2]][, , 3])
det.covs[[3]] <- list(det.cov.3.1 = X.p[[3]][, , 2])
det.covs[[4]] <- list(det.cov.4.1 = X.p[[4]][, , 2], 
                      det.cov.4.2 = X.p[[4]][, , 3], 
                      det.cov.4.3 = X.p[[4]][, , 4])
data.list <- list(y = y, 
                  occ.covs = occ.covs,
                  det.covs = det.covs, 
                  sites = sites, 
                  coords = coords)

J <- length(dat$z.obs)

# Initial values
inits.list <- list(alpha = list(0, 0, 0, 0), 
                   beta = 0, 
                   phi = 3 / .5, 
                   sigma.sq = 2, 
                   w = rep(0, J), 
                   z = rep(1, J))
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                   alpha.normal = list(mean = list(0, 0, 0, 0), 
                                       var = list(2.72, 2.72, 2.72, 2.72)),
                   phi.unif = c(3/1, 3/.1), 
                   sigma.sq.ig = c(2, 2))
# Tuning
tuning.list <- list(phi = 0.3) 

# Number of batches
n.batch <- 2
# Batch length
batch.length <- 25

# Note that this is just a test case and more iterations/chains may need to 
# be run to ensure convergence.
out <- spIntPGOcc(occ.formula = ~ occ.cov, 
                  det.formula = list(f.1 = ~ det.cov.1.1, 
                                     f.2 = ~ det.cov.2.1 + det.cov.2.2, 
                                     f.3 = ~ det.cov.3.1, 
                                     f.4 = ~ det.cov.4.1 + det.cov.4.2 + det.cov.4.3), 
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
                  NNGP = FALSE, 
                  n.report = 10, 
                  n.burn = 10, 
                  n.thin = 1)

summary(out)
}
