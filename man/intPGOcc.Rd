\name{intPGOcc}
\alias{intPGOcc}
\title{Function for Fitting Single-Species Integrated Occupancy Models Using Polya-Gamma Latent Variables}


\description{
Function for fitting single-species integrated occupancy models using Polya-Gamma latent variables. Data integration is done using a joint likelihood framework, assuming distinct detection models for each data source that are each conditional on a single latent occurrence process. 
}

\usage{
intPGOcc(occ.formula, det.formula, data, inits, priors, n.samples, 
         n.omp.threads = 1, verbose = TRUE, n.report = 1000, 
         n.burn = round(.10 * n.samples), n.thin = 1, n.chains = 1,
         k.fold, k.fold.threads = 1, k.fold.seed, 
         k.fold.data, ...)
}
\arguments{
  \item{occ.formula}{a symbolic description of the model to be fit
  for the occurrence portion of the model using R's model syntax. Only
  right-hand side of formula is specified. See example below.}

 \item{det.formula}{a list of symbolic descriptions of the models to be fit
  for the detection portion of the model using R's model syntax for 
  each data set. Each element in the list is a formula for the detection 
  model of a given data set. Only right-hand side of formula is specified. 
  See example below.}
 
  \item{data}{a list containing data necessary for model fitting.
  Valid tags are \code{y}, \code{occ.covs}, \code{det.covs}, and \code{sites}. 
  \code{y} is a list of matrices or data frames for each data set used in 
  the integrated model. Each element of the list has first dimension equal 
  to the number of sites with that data source and second dimension equal 
  to the maximum number of replicates at a given site. \code{occ.covs} is 
  a matrix or data frame containing the variables used in the occupancy 
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
  and number of columns equal to the maximum number of replicates at a given site.}

  \item{inits}{a list with each tag corresponding to a parameter name.
  Valid tags are \code{z}, \code{beta}, and \code{alpha}. The value portion
  of tags \code{z} and \code{beta} is the parameter's initial value. The tag
  \code{alpha} is a list comprised of the initial values for the detection 
  parameters for each data source. Each element of the list should be a vector of 
  initial values for all detection parameters in the given data source or a single value
  for each data source to assign all parameters for a given data source
  the same initial value. See \code{priors} description for definition of 
  each parameter name. Additionally, the tag \code{fix} can be set to \code{TRUE} 
  to fix the starting values across all chains. If \code{fix} is not specified
  (the default), starting values are varied randomly across chains.}

  \item{priors}{a list with each tag corresponding to a parameter name. 
    Valid tags are \code{beta.normal} and \code{alpha.normal}. Occurrence 
    (\code{beta}) and detection (\code{alpha}) regression coefficients 
    are assumed to follow a normal distribution. For \code{beta} hyperparameters of the 
    normal distribution are passed as a list of length two with the first
    and second elements corresponding to the mean and variance of the normal
    distribution, which are each specified as vectors of 
    length equal to the number of coefficients to be estimated or of 
    length one if priors are the same for all coefficients. For 
    the detection coefficients \code{alpha}, the mean and variance 
    hyperparameters are themselves passed in as lists, with each element 
    of the list corresponding to the specific hyperparameters for the 
    detection parameters in a given data source. If not specified, prior means 
    are set to 0 and prior variances set to 2.72.}
  
  \item{n.samples}{the number of posterior samples to collect in each chain.}

  \item{n.omp.threads}{a positive integer indicating the number of threads
    to use for SMP parallel processing. The package must be compiled for 
    OpenMP support. For most Intel-based machines, we recommend setting 
    \code{n.omp.threads} up to the number of hypterthreaded cores. Note, 
    \code{n.omp.threads} > 1 might not work on some systems.}
  
  \item{verbose}{if \code{TRUE}, messages about data preparation, 
    model specification, and progress of the sampler are printed to the screen. 
    Otherwise, no messages are printed.}
  
  \item{n.report}{the interval to report MCMC progress.}

  \item{n.burn}{the number of samples out of the total \code{n.samples} to 
    discard as burn-in. By default, the first 10\% of samples is discarded.}
  
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

  \item{k.fold.data}{an integer specifying the specific data set to hold out
    values from. If not specified, data from all data set locations will 
    be incorporated into the k-fold cross-validation.}

  
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

  Hooten, M. B., and Hobbs, N. T. (2015). A guide to Bayesian model 
  selection for ecologists. Ecological monographs, 85(1), 3-28.

  Finley, A. O., Datta, A., and Banerjee, S. (2020). spNNGP R 
  package for nearest neighbor Gaussian process models. 
  arXiv preprint arXiv:2001.09111.

}
\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{intPGOcc} that is a list comprised of: 

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the occupancy regression coefficients.}

  \item{alpha.samples}{a \code{coda} object of posterior samples
    for the detection regression coefficients for all data sources.}

  \item{z.samples}{a \code{coda} object of posterior samples 
    for the latent occupancy values}

  \item{psi.samples}{a \code{coda} object of posterior samples
    for the latent occupancy probability values}

  \item{rhat}{a list of Gelman-Rubin diagnostic values for some of the model
    parameters.}

  \item{ESS}{a list of effective sample sizes for some of the model parameters.}

  \item{y.rep.samples}{a list of three-dimensional arrays of fitted values for 
    each data source for use in Goodness of Fit assessments.}

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
set.seed(1008)

# Simulate Data -----------------------------------------------------------
J.x <- 15
J.y <- 15
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
beta <- c(0.5, 1)
p.occ <- length(beta)
# Detection covariates
alpha <- list()
for (i in 1:n.data) {
  alpha[[i]] <- runif(2, -1, 1)
}
p.det.long <- sapply(alpha, length)
p.det <- sum(p.det.long)

# Simulate occupancy data. 
dat <- simIntOcc(n.data = n.data, J.x = J.x, J.y = J.y, J.obs = J.obs, 
                 n.rep = n.rep, beta = beta, alpha = alpha, sp = FALSE)

y <- dat$y
X <- dat$X.obs
X.p <- dat$X.p
sites <- dat$sites

# Package all data into a list
occ.covs <- X[, 2, drop = FALSE]
colnames(occ.covs) <- c('occ.cov')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(det.cov.1.1 = X.p[[1]][, , 2]) 
det.covs[[2]] <- list(det.cov.2.1 = X.p[[2]][, , 2]) 
det.covs[[3]] <- list(det.cov.3.1 = X.p[[3]][, , 2]) 
det.covs[[4]] <- list(det.cov.4.1 = X.p[[4]][, , 2]) 
data.list <- list(y = y, 
                  occ.covs = occ.covs,
                  det.covs = det.covs, 
                  sites = sites)

J <- length(dat$z.obs)
# Initial values
inits.list <- list(alpha = list(0, 0, 0, 0), 
                   beta = 0, 
                   z = rep(1, J))
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                   alpha.normal = list(mean = list(0, 0, 0, 0), 
                                       var = list(2.72, 2.72, 2.72, 2.72)))
n.samples <- 5000
out <- intPGOcc(occ.formula = ~ occ.cov, 
                det.formula = list(f.1 = ~ det.cov.1.1, 
                                   f.2 = ~ det.cov.2.1, 
                                   f.3 = ~ det.cov.3.1, 
                                   f.4 = ~ det.cov.4.1), 
                data = data.list,
                inits = inits.list,
                n.samples = n.samples, 
                priors = prior.list, 
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = 1000, 
                n.thin = 1, 
                n.chains = 1)

summary(out)
}
