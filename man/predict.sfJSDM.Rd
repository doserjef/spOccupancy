\name{predict.sfJSDM}
\alias{predict.sfJSDM}
\title{Function for prediction at new locations for spatial factor joint species distribution model}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new locations given an object of class `sfJSDM`.
}

\usage{
\method{predict}{sfJSDM}(object, X.0, coords.0, n.omp.threads = 1, verbose = TRUE, 
        n.report = 100, ignore.RE = FALSE, ...)
}

\arguments{

  \item{object}{an object of class sfJSDM}

  \item{X.0}{the design matrix of covariates at the prediction locations. This should include a column of 1s for the intercept if an intercept is included in the model. If random effects are included in the model, the levels of the random effects at the new locations should be included as a column in the design matrix. The ordering of the levels should match the ordering used to fit the data in \code{sfJSDM}. Columns should correspond to the order of how covariates were specified in the formula argument of \code{sfJSDM}. Column names of the random effects must match the name of the random effects, if specified in the formula argument of \code{sfJSDM}.}

  \item{coords.0}{the spatial coordinates corresponding to \code{X.0}. Note that \code{spOccupancy} assumes coordinates are specified 
  in a projected coordinate system.}

  \item{n.omp.threads}{a positive integer indicating
    the number of threads to use for SMP parallel processing. The package must
    be compiled for OpenMP support. For most Intel-based machines, we recommend setting
    \code{n.omp.threads} up to the number of hyperthreaded cores. 
    Note, \code{n.omp.threads} > 1 might not work on some systems.}
 
  \item{verbose}{if \code{TRUE}, model specification and progress of the 
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.}

  \item{n.report}{the interval to report sampling progress.}

  \item{ignore.RE}{a logical value indicating whether to include unstructured random 
    effects for prediction. If TRUE, unstructured random effects will be ignored and 
    prediction will only use the fixed effects and the spatial random effects. If FALSE, 
    random effects will be included in the prediction for both observed and unobserved 
    levels of the unstructured random effects.}

  \item{...}{currently no additional arguments}
}

\note{
  When \code{ignore.RE = FALSE}, both sampled levels and non-sampled levels of random effects are supported for prediction. For sampled levels, the posterior distribution for the random intercept corresponding to that level of the random effect will be used in the prediction. For non-sampled levels, random values are drawn from a normal distribution using the posterior samples of the random effect variance, which results in fully propagated uncertainty in predictions with models that incorporate random effects.
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An list object of class \code{predict.sfJSDM} that consists of:  

  \item{psi.0.samples}{a three-dimensional array of posterior predictive samples for the
    latent occurrence probability values.}

  \item{z.0.samples}{a three-dimensional array of posterior predictive samples for the
    latent occurrence values.}

  \item{w.0.samples}{a three-dimensional array of posterior predictive samples
    for the latent spatial factors.}

  \item{run.time}{execution time reported using \code{proc.time()}.}

  The return object will include additional objects used for standard 
  extractor functions. 
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
n.factors <- 3
phi <- runif(n.factors, 3/1, 3/.4)
sp <- TRUE

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
                phi = phi, sigma.sq = sigma.sq, sp = TRUE, cov.model = 'exponential', 
                factor.model = TRUE, n.factors = n.factors)

# Number of batches
n.batch <- 10
# Batch length
batch.length <- 25
n.samples <- n.batch * batch.length

# Split into fitting and prediction data set
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
# Summarize the multiple replicates into a single value for use in a JSDM
y <- apply(dat$y[, -pred.indx, ], c(1, 2), max, na.rm = TRUE)
# Occupancy covariates
X <- dat$X[-pred.indx, ]
# Coordinates
coords <- as.matrix(dat$coords[-pred.indx, ])
# Prediction values
X.0 <- dat$X[pred.indx, ]
coords.0 <- as.matrix(dat$coords[pred.indx, ])
psi.0 <- dat$psi[, pred.indx]

# Package all data into a list
covs <- X[, 2, drop = FALSE]
colnames(covs) <- c('occ.cov')
data.list <- list(y = y, 
                  covs = covs,
                  coords = coords)

# Priors 
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72), 
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                   phi.unif = list(a = 3/1, b = 3/.1)) 
# Starting values
lambda.inits <- matrix(0, N, n.factors)
diag(lambda.inits) <- 1
lambda.inits[lower.tri(lambda.inits)] <- rnorm(sum(lower.tri(lambda.inits)))
inits.list <- list(beta.comm = 0, 
                   beta = 0, 
                   tau.sq.beta = 1, 
                   phi = 3 / .5, 
                   sigma.sq = 2,
                   lambda = lambda.inits)
# Tuning
tuning.list <- list(phi = 1) 

out <- sfJSDM(formula = ~ occ.cov, 
              data = data.list,
              inits = inits.list, 
              n.batch = n.batch, 
              batch.length = batch.length, 
              accept.rate = 0.43, 
              n.factors = 3,
              priors = prior.list, 
              cov.model = "exponential", 
              tuning = tuning.list, 
              n.omp.threads = 1, 
              verbose = TRUE, 
              NNGP = TRUE, 
              n.neighbors = 5, 
              search.type = 'cb', 
              n.report = 10, 
              n.burn = 100, 
              n.thin = 1)

summary(out, level = 'both')

# Predict at new locations ------------------------------------------------
out.pred <- predict(out, X.0, coords.0, verbose = FALSE)
}
