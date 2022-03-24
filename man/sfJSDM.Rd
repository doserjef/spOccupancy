\name{sfJSDM}
\alias{sfJSDM}
\title{Function for Fitting a Spatial Factor Joint Species Distribution Model}

\description{
  The function \code{sfJSDM} fits a spatially-explicit joint species distribution model. This model does not explicitly account for imperfect detection (see \code{sfMsPGOcc()}). We use Polya-Gamma latent variables and a spatial factor modeling approach. Currently, models are implemented using a Nearest Neighbor Gaussian Process. Future development will allow for running the models using a full Gaussian Process.}
  
\usage{
sfJSDM(formula, data, inits, priors, tuning, 
       cov.model = 'exponential', NNGP = TRUE, 
       n.neighbors = 15, search.type = 'cb', n.factors, n.batch, 
       batch.length, accept.rate = 0.43, n.omp.threads = 1, 
       verbose = TRUE, n.report = 100, 
       n.burn = round(.10 * n.batch * batch.length), n.thin = 1, 
       n.chains = 1, k.fold, k.fold.threads = 1, k.fold.seed, ...)
}

\arguments{

  \item{formula}{a symbolic description of the model to be fit
    for the model using R's model syntax. Only right-hand side of formula 
    is specified. See example below. Random intercepts are allowed 
    using \pkg{lme4} syntax (Bates et al. 2015).}
  
  \item{data}{a list containing data necessary for model fitting.
    Valid tags are \code{y}, \code{covs}, and \code{coords}. \code{y}
    is a two-dimensional array with first dimension equal to the number 
    of species and second dimension equal to the number of sites. Note how this differs
    from other \code{spOccupancy} functions in that \code{y} does not have any 
    replicate surveys. This is because \code{sfJSDM} does not account for imperfect
    detection. \code{covs} is a matrix or data frame containing the variables 
    used in the model, with \eqn{J}{J} rows for 
    each column (variable). \code{coords} is a matrix with \eqn{J}{J} 
    rows and 2 columns consisting of the spatial coordinates of each site in 
    the data. Note that \code{spOccupancy} assumes coordinates are specified 
  in a projected coordinate system.}

  \item{inits}{a list with each tag corresponding to a parameter name.
    Valid tags are \code{beta.comm}, \code{beta},  \code{tau.sq.beta}, 
    \code{phi}, \code{lambda}, \code{sigma.sq.psi}, and \code{nu}. 
    \code{nu} is only specified if \code{cov.model = "matern"}. 
    \code{sigma.sq.psi} is only specified if
    random intercepts are included in \code{formula}. The value portion of each tag is 
    the parameter's initial value. See \code{priors} description for definition
    of each parameter name. Additionally, the tag \code{fix} can be set to \code{TRUE} 
    to fix the starting values across all chains. If \code{fix} is not specified
    (the default), starting values are varied randomly across chains.}

  \item{priors}{a list with each tag corresponding to a parameter name. 
  Valid tags are \code{beta.comm.normal}, \code{tau.sq.beta.ig}, \code{phi.unif}, 
  \code{nu.unif}, and \code{sigma.sq.psi.ig}. Community-level occurrence 
  (\code{beta.comm}) regression coefficients are assumed to follow a 
  normal distribution. The hyperparameters of the normal distribution
  are passed as a list of length two with the first and second elements 
  corresponding to the mean and variance of the normal distribution, 
  which are each specified as vectors of length equal to the number of 
  coefficients to be estimated or of length one if priors are the same for 
  all coefficients. If not specified, prior means are set 
  to 0 and prior variances set to 2.73. Community-level variance parameters 
  (\code{tau.sq.beta}) are assumed to follow an 
  inverse Gamma distribution. The hyperparameters of 
  the inverse gamma distribution are passed as a list of length two with 
  the first and second elements corresponding to the shape and scale parameters,
  which are each specified as vectors of length equal to the number of 
  coefficients to be estimated or a single value if priors are the same for all 
  parameters. If not specified, prior shape and scale 
  parameters are set to 0.1. The spatial factor model fits \code{n.factors} independent
  spatial processes. The spatial decay \code{phi} and smoothness \code{nu} parameters 
  for each latent factor are assumed to follow Uniform distributions. 
  The hyperparameters of the Uniform are passed as a list with two elements, 
  with both elements being vectors of length \code{n.factors} corresponding to the lower and 
  upper support, respectively, or as a single value if the same value is assigned
  for all factors. The priors for the factor loadings matrix \code{lambda} are fixed
  following the standard spatial factor model to ensure parameter 
  identifiability (Christensen and Amemlya 2002). The
  upper triangular elements of the \code{N x n.factors} matrix are fixed at 0 and the 
  diagonal elements are fixed at 1. The lower triangular elements are assigned a 
  standard normal prior (i.e., mean 0 and variance 1).
  \code{sigma.sq.psi} is the random effect variance for any random effects, and is assumed
  to follow an inverse Gamma distribution. The hyperparameters of the inverse-Gamma distribution
  are passed as a list of length two with first and second elements corresponding to the shape
  and scale parameters, respectively, which are each specified as vectors of length equal to 
  the number of random intercepts or of length one if priors are the same for all random 
  effect variances.}

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
    Finley et al. (2019) for more information. For spatial factor models, only 
    \code{NNGP = TRUE} is currently supported.}
  
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
  An object of class \code{sfJSDM} that is a list comprised of: 

  \item{beta.comm.samples}{a \code{coda} object of posterior samples
    for the community level occurrence regression coefficients.}
  
  \item{tau.sq.beta.samples}{a \code{coda} object of posterior samples
    for the occurrence community variance parameters.}

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the species level occurrence regression coefficients.}

  \item{theta.samples}{a \code{coda} object of posterior samples
    for the species level correlation parameters.}

  \item{lambda.samples}{a \code{coda} object of posterior samples
    for the latent spatial factor loadings.}

  \item{psi.samples}{a three-dimensional array of posterior samples for
    the latent occurrence probability values for each species.}

  \item{w.samples}{a three-dimensional array of posterior samples for 
    the latent spatial random effects for each latent factor.}

  \item{sigma.sq.psi.samples}{a \code{coda} object of posterior samples
    for variances of random intercepts included in the occurrence portion
    of the model. Only included if random intercepts are specified in 
    \code{formula}.}

  \item{beta.star.samples}{a \code{coda} object of posterior samples
    for the occurrence random effects. Only included if random intercepts
    are specified in \code{formula}.}

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
  estimated values are not included in the model object, but can be extracted
  using \code{fitted()}.
}

\examples{
set.seed(408)
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- rep(1, J) 
N <- 10
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 1.3, 0.5)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, .5, 1.6)
# Detection
# Fix this to be constant and really close to 1. 
alpha.mean <- c(9)
tau.sq.alpha <- c(0.05)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list()
# Include a single random effect
psi.RE <- list(levels = c(50), 
               sigma.sq.psi = c(1.5))
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
alpha.true <- alpha
# Factor model
factor.model <- TRUE
n.factors <- 4
sigma.sq <- rep(1, n.factors)
phi <- rep(3/.4, n.factors)
sp <- TRUE
cov.model <- "exponential"

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
                psi.RE = psi.RE, p.RE = p.RE, sp = sp, sigma.sq = sigma.sq, 
                phi = phi, cov.model = cov.model, factor.model = TRUE, n.factors = n.factors)

X <- dat$X
y <- dat$y
coords <- dat$coords
X.re <- dat$X.re

occ.covs <- data.frame(X, X.re)
occ.covs$X.re <- occ.covs$X.re
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2', 'occ.re.1')
data.list <- list(y = y[, , 1], 
                  covs = occ.covs, 
                  coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                   phi.unif = list(a = 3 / 1, b = 3 / .1))
inits.list <- list(beta.comm = 0, 
                   beta = 0, 
                   tau.sq.beta = 1, 
                   phi = 3 / .5) 
tuning.list <- list(phi = 1)

out <- sfJSDM(formula = ~ occ.cov.1 + occ.cov.2 + (1 | occ.re.1), 
              data = data.list, 
              inits = inits.list, 
              priors = prior.list, 
              tuning = tuning.list, 
              cov.model = 'spherical', 
              NNGP = TRUE, 
              n.neighbors = 5,
              n.factors = 4, 
              n.batch = 10, 
              batch.length = 25, 
              n.report = 250, 
              n.burn = 50,
              n.thin = 2,
              n.chains = 1) 
summary(out)
}
