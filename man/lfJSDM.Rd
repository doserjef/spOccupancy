\name{lfJSDM}
\alias{lfJSDM}
\title{Function for Fitting a Latent Factor Joint Species Distribution Model}


\description{
Function for fitting a joint species distribution model with species correlations. This model does not explicitly account for imperfect detection (see \code{lfMsPGOcc()}). We use Polya-gamma latent variables and a factor modeling approach.}

\usage{
lfJSDM(formula, data, inits, priors, n.factors, 
       n.samples, n.omp.threads = 1, verbose = TRUE, n.report = 100, 
       n.burn = round(.10 * n.samples), n.thin = 1, n.chains = 1,
       k.fold, k.fold.threads = 1, k.fold.seed, k.fold.only = FALSE, ...)
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
  replicate surveys. This is because \code{lfJSDM} does not account for imperfect
  detection. \code{covs} is a matrix or data frame containing the variables 
  used in the model, with \eqn{J}{J} rows for 
  each column (variable). \code{coords} is a matrix with \eqn{J}{J} 
  rows and 2 columns consisting of the spatial coordinates of each site in 
  the data. Note that \code{spOccupancy} assumes coordinates are specified 
  in a projected coordinate system.}

\item{inits}{a list with each tag corresponding to a parameter name.
  Valid tags are \code{beta.comm}, \code{beta}, \code{tau.sq.beta}, \code{sigma.sq.psi},
  \code{lambda}. The value portion of each tag is the parameter's initial value. See 
  \code{priors} description for definition of each parameter name.
  Additionally, the tag \code{fix} can be set to \code{TRUE} 
  to fix the starting values across all chains. If \code{fix} is not specified
  (the default), starting values are varied randomly across chains.}

\item{priors}{a list with each tag corresponding to a parameter name. 
  Valid tags are \code{beta.comm.normal}, \code{tau.sq.beta.ig}, and 
  \code{sigma.sq.psi.ig}. Community-level (\code{beta.comm}) regression coefficients 
  are assumed to follow a normal distribution. The hyperparameters of the 
  normal distribution are passed as a list of length two with the first
  and second elements corresponding to the mean and variance of the normal
  distribution, which are each specified as vectors of 
  length equal to the number of coefficients to be estimated or of length one
  if priors are the same for all coefficients. If not
  specified, prior means are set to 0 and prior variances set to 2.72. 
  Community-level variance parameters (\code{tau.sq.beta})
  are assumed to follow an inverse Gamma distribution. 
  The hyperparameters of the inverse gamma distribution are 
  passed as a list of length two with the first and second elements 
  corresponding to the shape and scale parameters, which are each specified
  as vectors of length equal to the number of coefficients to be estimated 
  or a single value if all parameters are assigned the same prior. 
  If not specified, prior shape and scale parameters are set to 0.1. The 
  factor model fits \code{n.factors} independent latent factors. The priors for the 
  factor loadings matrix \code{lambda} are fixed following standard approaches to 
  ensure parameter identifiability. The upper triangular elements of the \code{N x n.factors} 
  matrix are fixed at 0 and the diagonal elements are fixed at 1. The lower 
  triangular elements are assigned a standard normal prior (i.e., mean 0 and variance 1). 
  \code{sigma.sq.psi} is the random effect variance for any random effects, and is assumed
  to follow an inverse Gamma distribution. The hyperparameters of the inverse-Gamma distribution
  are passed as a list of length two with first and second elements corresponding to the shape
  and scale parameters, respectively, which are each specified as vectors of length equal to 
  the number of random intercepts or of length one if priors are the same for all random 
  effect variances.}

\item{n.factors}{the number of factors to use in the latent factor model approach.
    Typically, the number of factors is set to be small (e.g., 4-5) relative to the
    total number of species in the community, which will lead to substantial
    decreases in computation time. However, the value can be anywhere
    between 0 and N (the number of species in the community). When set to 0, the model
    assumes there are no residual species correlations, which is equivalent to the 
    \code{msPGOcc()} function but without imperfect detection.}

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
  Polson, N.G., J.G. Scott, and J. Windle. (2013) Bayesian Inference for
  Logistic Models Using Polya-Gamma Latent Variables.
  \emph{Journal of the American Statistical Association}, 108:1339-1349.

  Bates, Douglas, Martin Maechler, Ben Bolker, Steve Walker (2015).
  Fitting Linear Mixed-Effects Models Using lme4. Journal of
  Statistical Software, 67(1), 1-48. \doi{10.18637/jss.v067.i01}.

  Hooten, M. B., and Hobbs, N. T. (2015). A guide to Bayesian model 
  selection for ecologists. Ecological monographs, 85(1), 3-28.
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{lfJSDM} that is a list comprised of: 

  \item{beta.comm.samples}{a \code{coda} object of posterior samples
    for the community level occurrence regression coefficients.}
  
  \item{tau.sq.beta.samples}{a \code{coda} object of posterior samples
    for the occurrence community variance parameters.}

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the species level occurrence regression coefficients.}

  \item{lambda.samples}{a \code{coda} object of posterior samples
    for the latent factor loadings.}

  \item{psi.samples}{a three-dimensional array of posterior samples for
    the latent probability of occurrence/detection values for each species.}

  \item{sigma.sq.psi.samples}{a \code{coda} object of posterior samples
    for variances of random intercepts included in the occurrence portion
    of the model. Only included if random intercepts are specified in 
    \code{occ.formula}.}

  \item{w.samples}{a three-dimensional array of posterior samples for 
    the latent effects for each latent factor.}

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
  estimated values are not included in the model object, but can be extracted
  using \code{fitted()}. 
}

\examples{
set.seed(400)
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- rep(1, J)
N <- 10
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 0.6, 1.5)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 1.2, 1.7)
# Detection
# Fix this to be constant and really close to 1. 
alpha.mean <- c(9)
tau.sq.alpha <- c(0.05)
p.det <- length(alpha.mean)
# Random effects
# Include a single random effect
psi.RE <- list(levels = c(20), 
               sigma.sq.psi = c(2))
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

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
                psi.RE = psi.RE, p.RE = p.RE, sp = FALSE,
                factor.model = TRUE, n.factors = 4)

X <- dat$X
y <- dat$y
X.re <- dat$X.re
coords <- dat$coords
occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2', 'occ.re.1')
data.list <- list(y = y[, , 1], 
                  covs = occ.covs, 
                  coords = coords) 
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1)) 
inits.list <- list(beta.comm = 0, beta = 0, tau.sq.beta = 1) 
out <- lfJSDM(formula = ~ occ.cov.1 + occ.cov.2 + (1 | occ.re.1), 
              data = data.list, 
              inits = inits.list, 
              priors = prior.list, 
              n.factors = 4, 
              n.samples = 1000,
              n.report = 500, 
              n.burn = 500,
              n.thin = 2,
              n.chains = 1) 
summary(out)
}
