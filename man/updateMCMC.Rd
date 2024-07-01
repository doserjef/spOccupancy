\name{updateMCMC}
\alias{updateMCMC}
\title{Update a spOccupancy or spAbundance model run with more MCMC iterations}

\usage{
updateMCMC(object, n.batch, n.samples, n.burn = 0, n.thin, 
           keep.orig = TRUE, verbose = TRUE, n.report = 100, 
           save.fitted = TRUE, ...)
}

\description{
  Function for updating a previously run spOccupancy or spAbundance model with 
  additional MCMC iterations. This function is useful for situations where a model
  is run for a long time but convergence/adequate mixing of the MCMC chains is not 
  reached. Instead of re-running the entire model again, this function allows you 
  to pick up where you left off. This function is currently in development, and only currently works 
  with the following spOccupancy and spAbundance model objects: msAbund, sfJSDM, lfJSDM. Note that
  cross-validation is not possible when updating the model.
}

\arguments{
  \item{object}{a \code{spOccupancy} or \code{spAbundance} model object. Currently 
    supports objects of class \code{msAbund} and \code{sfJSDM}.}

  \item{n.batch}{the number of additional MCMC batches in each chain 
    to run for the adaptive MCMC sampler. Only valid for model types fit with 
    an adaptive MCMC sampler}

  \item{n.samples}{the number of posterior samples to collect in each chain. Only 
    valid for model types that are run with a fully Gibbs sampler and have 
    \code{n.samples} as an argument in the original model fitting function.}

  \item{n.burn}{the number of samples out of the total \code{n.batch * batchlength} to 
  discard as burn-in for each chain from the updated samples. Note this argument
  does not discard samples from the previous model run, and rather only applies
  to the samples in the updated run of the model. Defaults to 0}

  \item{n.thin}{the thinning interval for collection of MCMC samples in 
    the updated model run. The thinning occurs after the \code{n.burn} 
    samples are discarded. Default value is set to 1.}

  \item{keep.orig}{A logical value indicating whether or not the samples from 
    the original run of the model should be kept or discarded.}

  \item{verbose}{if \code{TRUE}, messages about data preparation, 
    model specification, and progress of the sampler are printed to the screen. 
    Otherwise, no messages are printed.}

  \item{n.report}{the interval to report Metropolis sampler acceptance and
    MCMC progress.}
    
  \item{save.fitted}{logical value indicating whether or not fitted values 
    and likelihood values should be saved in the resulting model object. This is only
    relevant for models of class \code{msAbund}. If \code{save.fitted = FALSE}, 
    the components \code{y.rep.samples}, \code{mu.samples}, and \code{like.samples} 
    will not be included in the model object, and subsequent functions for 
    calculating WAIC, fitted values, and posterior predictive checks will not work, 
    although they all can be calculated manually if desired. Setting 
    \code{save.fitted = FALSE} can be useful when working with very large 
    data sets to minimize the amount of RAM needed when fitting and storing 
    the model object in memory.}

  \item{...}{currently no additional arguments}
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
}

\value{
An object of the same class as the original model fit provided in the \code{argument} object. See the manual page for the original model type for complete details.
}

\examples{
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6)
# Detection
alpha.mean <- c(0)
tau.sq.alpha <- c(1)
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
alpha.true <- alpha
n.factors <- 3
phi <- rep(3 / .7, n.factors)
sigma.sq <- rep(2, n.factors)
nu <- rep(2, n.factors)

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
                psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq,
                phi = phi, nu = nu, cov.model = 'matern', factor.model = TRUE,
                n.factors = n.factors)

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , drop = FALSE]

y <- apply(y, c(1, 2), max, na.rm = TRUE)
data.list <- list(y = y, coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                   nu.unif = list(0.5, 2.5))
# Starting values
inits.list <- list(beta.comm = 0,
                   beta = 0,
                   fix = TRUE,
                   tau.sq.beta = 1)
# Tuning
tuning.list <- list(phi = 1, nu = 0.25)

batch.length <- 25
n.batch <- 2
n.report <- 100
formula <- ~ 1

out <- sfJSDM(formula = formula,
              data = data.list,
              inits = inits.list,
              n.batch = n.batch,
              batch.length = batch.length,
              accept.rate = 0.43,
              priors = prior.list,
              cov.model = "matern",
              tuning = tuning.list,
              n.factors = 3,
              n.omp.threads = 1,
              verbose = TRUE,
              NNGP = TRUE,
              n.neighbors = 5,
              search.type = 'cb',
              n.report = 10,
              n.burn = 0,
              n.thin = 1,
              n.chains = 2)
summary(out)

# Update the initial model fit
out.new <- updateMCMC(out, n.batch = 1, keep.orig = TRUE, 
		     verbose = TRUE, n.report = 1) 
summary(out.new)
}
