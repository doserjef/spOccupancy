\name{predict.svcPGBinom}
\alias{predict.svcPGBinom}
\title{Function for prediction at new locations for single-species spatially-varying coefficient Binomial models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new 
  locations given an object of class `svcPGBinom`.
}

\usage{
\method{predict}{svcPGBinom}(object, X.0, coords.0, weights.0, n.omp.threads = 1, verbose = TRUE, 
        n.report = 100, ignore.RE = FALSE, ...)
}

\arguments{

  \item{object}{an object of class \code{svcPGBinom}}

  \item{X.0}{the design matrix of covariates at the prediction locations. Note that for spatially-varying coefficients models the order of covariates in \code{X.0} must be the same as the order of covariates specified in the model formula. This should include a column of 1s for the intercept if an intercept is included in the model. If unstructured random effects are included in the model, the levels of the random effects at the new locations should be included as a column in the design matrix. The ordering of the levels should match the ordering used to fit the data in \code{svcPGBinom}. Columns should correspond to the order of how covariates were specified in the corresponding formula argument of \code{svcPGBinom}. Column names of the random effects must match the name of the random effects, if specified in the corresponding formula argument of \code{svcPGBinom}.}

  \item{coords.0}{the spatial coordinates corresponding to \code{X.0}. Note that \code{spOccupancy} assumes coordinates are specified in a projected coordinate system.}

  \item{weights.0}{a numeric vector containing the binomial weights (i.e., the total number of 
    Bernoulli trials) at each site. If \code{weights.0} is not specified, 
    we assume 1 trial at each site (i.e., presence/absence).}

   \item{n.omp.threads}{a positive integer indicating
    the number of threads to use for SMP parallel processing. The package must
    be compiled for OpenMP support. For most Intel-based machines, we recommend setting
    \code{n.omp.threads} up to the number of hyperthreaded cores. 
    Note, \code{n.omp.threads} > 1 might not work on some systems.}
 
  \item{verbose}{if \code{TRUE}, model specification and progress of the 
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.}

  \item{ignore.RE}{a logical value indicating whether to include unstructured random 
    effects for prediction. If TRUE, unstructured random effects will be ignored and 
    prediction will only use the fixed effects and the spatial random effects. If FALSE, 
    random effects will be included in the prediction for both observed and unobserved 
    levels of the unstructured random effects.}

  \item{n.report}{the interval to report sampling progress.}

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
  A list object of class \code{predict.svcPGBinom} consisting of:  

  \item{psi.0.samples}{a \code{coda} object of posterior predictive samples 
    for the binomial probability values.}
  
  \item{y.0.samples}{a \code{coda} object of posterior predictive samples 
    for the binomial data.}
  
  \item{w.0.samples}{a three-dimensional array of posterior predictive samples
    for the spatial random effects, with dimensions corresponding to MCMC iteration, 
    coefficient, and site.}

  \item{run.time}{execution time reported using \code{proc.time()}.}
}

\examples{
set.seed(1000)
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
# Binomial weights
weights <- sample(10, J, replace = TRUE)
beta <- c(0, 0.5, -0.2, 0.75)
p <- length(beta)
# No unstructured random effects
psi.RE <- list()
# Spatial parameters
sp <- TRUE
# Two spatially-varying covariates. 
svc.cols <- c(1, 2)
p.svc <- length(svc.cols)
cov.model <- "exponential"
sigma.sq <- runif(p.svc, 0.4, 1.5)
phi <- runif(p.svc, 3/1, 3/0.2)

# Simulate the data  
dat <- simBinom(J.x = J.x, J.y = J.y, weights = weights, beta = beta, 
                psi.RE = psi.RE, sp = sp, svc.cols = svc.cols, 
                cov.model = cov.model, sigma.sq = sigma.sq, phi = phi)

# Binomial data
y <- dat$y
# Covariates
X <- dat$X
# Spatial coordinates
coords <- dat$coords

# Subset data for prediction if desired
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y.0 <- y[pred.indx, drop = FALSE]
X.0 <- X[pred.indx, , drop = FALSE]
coords.0 <- coords[pred.indx, ]
y <- y[-pred.indx, drop = FALSE]
X <- X[-pred.indx, , drop = FALSE]
coords <- coords[-pred.indx, ]
weights.0 <- weights[pred.indx]
weights <- weights[-pred.indx]

# Package all data into a list
# Covariates
covs <- cbind(X)
colnames(covs) <- c('int', 'cov.1', 'cov.2', 'cov.3')

# Data list bundle
data.list <- list(y = y, 
                  covs = covs,
                  coords = coords, 
                  weights = weights)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                   sigma.sq.ig = list(a = 2, b = 1), 
                   phi.unif = list(a = 3 / 1, b = 3 / 0.1)) 

# Starting values
inits.list <- list(beta = 0, alpha = 0,
                   sigma.sq = 1, phi = phi)
# Tuning
tuning.list <- list(phi = 1) 

n.batch <- 10
batch.length <- 25
n.burn <- 100
n.thin <- 1

out <- svcPGBinom(formula = ~ cov.1 + cov.2 + cov.3, 
                  svc.cols = c(1, 2),
                  data = data.list, 
                  n.batch = n.batch, 
                  batch.length = batch.length, 
                  inits = inits.list, 
                  priors = prior.list,
                  accept.rate = 0.43, 
                  cov.model = "exponential", 
                  tuning = tuning.list, 
                  n.omp.threads = 1, 
                  verbose = TRUE, 
                  NNGP = TRUE, 
                  n.neighbors = 5,
                  n.report = 2, 
                  n.burn = n.burn, 
                  n.thin = n.thin, 
                  n.chains = 1) 

summary(out)

# Predict at new locations ------------------------------------------------
out.pred <- predict(out, X.0, coords.0, weights.0, verbose = FALSE)
str(out.pred)
}
