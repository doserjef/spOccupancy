\name{predict.svcTPGBinom}
\alias{predict.svcTPGBinom}
\title{Function for prediction at new locations for multi-season single-species spatially-varying coefficient binomial models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new locations given an object of class `svcTPGBinom`. Prediction is possible for both the latent occupancy state as well as detection. Predictions are currently only possible for sampled primary time periods.  
}

\usage{
\method{predict}{svcTPGBinom}(object, X.0, coords.0, t.cols, weights.0,  n.omp.threads = 1, 
        verbose = TRUE, n.report = 100, ignore.RE = FALSE, ...)
}

\arguments{

  \item{object}{an object of class svcTPGBinom}

  \item{X.0}{the design matrix of covariates at the prediction locations. This should be a three-dimensional array, with dimensions corresponding to site, primary time period, and covariate, respectively. Note that the first covariate should consist of all 1s for the intercept if an intercept is included in the model. If random effects are included in the occupancy (or detection if \code{type = 'detection'}) portion of the model, the levels of the random effects at the new locations/time periods should be included as an element of the three-dimensional array. The ordering of the levels should match the ordering used to fit the data in \code{svcTPGBinom}. The covariates should be organized in the same order as they were specified in the corresponding formula argument of \code{svcTPGBinom}. Names of the third dimension (covariates) of any random effects in X.0 must match the name of the random effects used to fit the model, if specified in the corresponding formula argument of \code{svcTPGBinom}. See example below.}

  \item{coords.0}{the spatial coordinates corresponding to \code{X.0}. Note that \code{spOccupancy} assumes coordinates are specified 
  in a projected coordinate system.}

  \item{weights.0}{a numeric site by primary time period matrix containing the binomial weights (i.e., the total number of 
    Bernoulli trials) at each site and primary time period. If \code{weights.0} is not specified, 
    we assume 1 trial at each site/primary time period (i.e., presence/absence).}

  \item{t.cols}{an indexing vector used to denote which primary time periods are contained in the design matrix of covariates at the prediction locations (\code{X.0}). The values should denote the specific primary time periods used to fit the model. The values should indicate the columns in \code{data$y} used to fit the model for which prediction is desired. See example below.} 
   \item{n.omp.threads}{a positive integer indicating
    the number of threads to use for SMP parallel processing. The package must
    be compiled for OpenMP support. For most Intel-based machines, we recommend setting
    \code{n.omp.threads} up to the number of hyperthreaded cores. 
    Note, \code{n.omp.threads} > 1 might not work on some systems.}
 
  \item{verbose}{if \code{TRUE}, model specification and progress of the 
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.}

  \item{ignore.RE}{logical value that specifies whether or not to remove random unstructured occurrence (or detection if \code{type = 'detection'}) effects from the subsequent predictions. If \code{TRUE}, random effects will be included. If \code{FALSE}, unstructured random effects will be set to 0 and predictions will only be generated from the fixed effects, the spatial random effects, and AR(1) random effects if the model was fit with \code{ar1 = TRUE}.}

  \item{n.report}{the interval to report sampling progress.}

  \item{...}{currently no additional arguments}
}

\note{
  When \code{ignore.RE = FALSE}, both sampled levels and non-sampled levels of unstructured random effects are supported for prediction. For sampled levels, the posterior distribution for the random intercept corresponding to that level of the random effect will be used in the prediction. For non-sampled levels, random values are drawn from a normal distribution using the posterior samples of the random effect variance, which results in fully propagated uncertainty in predictions with models that incorporate random effects. 

  Occurrence predictions at sites that are only sampled for a subset of the total number of primary time periods are obtained directly when fitting the model. See the \code{psi.samples} and \code{y.rep.samples} portions of the output list from the model object of class \code{svcTPGBinom}. 
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  A list object of class \code{predict.svcTPGBinom} that consists of: 

  \item{psi.0.samples}{a three-dimensional object of posterior predictive samples for the
    occurrence probability values with dimensions corresponding to posterior predictive
    sample, site, and primary time period.}

  \item{y.0.samples}{a three-dimensional object of posterior predictive samples for the
    predicted binomial data with dimensions corresponding to posterior predictive sample, site, 
    and primary time period.}

  \item{w.0.samples}{a three-dimensional array of posterior predictive samples
    for the spatial random effects, with dimensions corresponding to MCMC iteration, 
    coefficient, and site.}

  \item{run.time}{execution time reported using \code{proc.time()}.}

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
# Subset data for prediction
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
y.0 <- dat$y[pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , , drop = FALSE]
# Prediction covariates
X.0 <- dat$X[pred.indx, , , drop = FALSE]
# Spatial coordinates
coords <- as.matrix(dat$coords[-pred.indx, ])
coords.0 <- as.matrix(dat$coords[pred.indx, ])
psi.0 <- dat$psi[pred.indx, ]
w.0 <- dat$w[pred.indx, ]
weights.0 <- weights[pred.indx, ]
weights <- weights[-pred.indx, ]

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
                   phi.unif = list(a = 3/1, b = 3/.1))

# Starting values
inits.list <- list(beta = beta, alpha = 0,
                   sigma.sq = 1, phi = 3 / 0.5, nu = 1)
# Tuning
tuning.list <- list(phi = 0.4, nu = 0.3, rho = 0.2)

# MCMC information
n.batch <- 2
n.burn <- 0
n.thin <- 1


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
                   n.report = 25,
                   n.burn = n.burn,
                   n.thin = n.thin,
                   n.chains = 1)
# Predict at new locations ------------------------------------------------
out.pred <- predict(out, X.0, coords.0, t.cols = 1:max(n.time), 
                    weights = weights.0, n.report = 10)
str(out.pred)
}
