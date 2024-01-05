\name{predict.svcTMsPGOcc}
\alias{predict.svcTMsPGOcc}
\title{Function for prediction at new locations for multi-season multi-species spatially-varying coefficient occupancy models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new locations given an object of class `svcTMsPGOcc`. Prediction is possible for both the latent occupancy state as well as detection. Predictions are currently only possible for sampled primary time periods.  
}

\usage{
\method{predict}{svcTMsPGOcc}(object, X.0, coords.0, t.cols, n.omp.threads = 1, 
                          verbose = TRUE, n.report = 100, 
                          ignore.RE = FALSE, type = 'occupancy', grid.index.0, ...)
}

\arguments{

  \item{object}{an object of class svcTMsPGOcc}

  \item{X.0}{the design matrix of covariates at the prediction locations. This should be a three-dimensional array, with dimensions corresponding to site, primary time period, and covariate, respectively. Note that the first covariate should consist of all 1s for the intercept if an intercept is included in the model. If random effects are included in the occupancy (or detection if \code{type = 'detection'}) portion of the model, the levels of the random effects at the new locations/time periods should be included as an element of the three-dimensional array. The ordering of the levels should match the ordering used to fit the data in \code{svcTMsPGOcc}. The covariates should be organized in the same order as they were specified in the corresponding formula argument of \code{svcTMsPGOcc}. Names of the third dimension (covariates) of any random effects in X.0 must match the name of the random effects used to fit the model, if specified in the corresponding formula argument of \code{svcTMsPGOcc}. See example below.}

  \item{coords.0}{the spatial coordinates corresponding to \code{X.0}. Note that \code{spOccupancy} assumes coordinates are specified 
  in a projected coordinate system.}

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

  \item{type}{a quoted keyword indicating what type of prediction to produce. Valid keywords are 'occupancy' to predict latent occupancy probability and latent occupancy values (this is the default), or 'detection' to predict detection probability given new values of detection covariates.}

  \item{grid.index.0}{an indexing vector used to specify how each row in \code{X.0} corresponds to the coordinates specified in \code{coords.0}. Only relevant if the spatial random effect was estimated at a higher spatial resolution (e.g., grid cells) than point locations.}

  \item{...}{currently no additional arguments}
}

\note{
  When \code{ignore.RE = FALSE}, both sampled levels and non-sampled levels of unstructured random effects are supported for prediction. For sampled levels, the posterior distribution for the random intercept corresponding to that level of the random effect will be used in the prediction. For non-sampled levels, random values are drawn from a normal distribution using the posterior samples of the random effect variance, which results in fully propagated uncertainty in predictions with models that incorporate random effects. 

  Occurrence predictions at sites that are only sampled for a subset of the total number of primary time periods are obtained directly when fitting the model. See the \code{psi.samples} and \code{z.samples} portions of the output list from the model object of class \code{svcTMsPGOcc}. 
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  A list object of class \code{predict.svcTMsPGOcc}. When \code{type = 'occupancy'}, the list consists of: 

  \item{psi.0.samples}{a four-dimensional object of posterior predictive samples for the
    latent occupancy probability values with dimensions corresponding to posterior predictive
    sample, species, site, and primary time period.}

  \item{z.0.samples}{a three-dimensional object of posterior predictive samples for the
    latent occupancy values with dimensions corresponding to posterior predictive sample, species, 
    site, and primary time period.}

  \item{w.0.samples}{a four-dimensional array of posterior predictive samples
    for the latent spatial factors with dimensions correpsonding to MCMC sample, latent factor, 
    site, and spatially-varying coefficient.}

  When \code{type = 'detection'}, the list consists of: 

  \item{p.0.samples}{a four-dimensional object of posterior predictive samples for the 
    detection probability values with dimensions corresponding to posterior predictive 
    sample, site, and primary time period.}

  The return object will include additional objects used for standard 
  extractor functions. 
}

\examples{
# Simulate Data -----------------------------------------------------------
set.seed(500)
J.x <- 8
J.y <- 8
J <- J.x * J.y
# Years sampled
n.time <- sample(3:10, J, replace = TRUE)
# n.time <- rep(10, J)
n.time.max <- max(n.time)
# Replicates
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(2:4, n.time[j], replace = TRUE)
  # n.rep[j, 1:n.time[j]] <- rep(4, n.time[j])
}
N <- 7
# Community-level covariate effects
# Occurrence
beta.mean <- c(-3, -0.2, 0.5)
trend <- FALSE
sp.only <- 0
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 1.5, 1.4)
# Detection
alpha.mean <- c(0, 1.2, -1.5)
tau.sq.alpha <- c(1, 0.5, 2.3)
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
sp <- TRUE
svc.cols <- c(1, 2)
p.svc <- length(svc.cols)
n.factors <- 2
phi <- runif(p.svc * n.factors, 3 / .9, 3 / .3)
factor.model <- TRUE
cov.model <- 'exponential'
ar1 <- TRUE
sigma.sq.t <- runif(N, 0.05, 1)
rho <- runif(N, 0.1, 1)

dat <- simTMsOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, N = N,
		 beta = beta, alpha = alpha, sp.only = sp.only, trend = trend,
		 psi.RE = psi.RE, p.RE = p.RE, factor.model = factor.model,
                 svc.cols = svc.cols, n.factors = n.factors, phi = phi, sp = sp,
                 cov.model = cov.model, ar1 = ar1, sigma.sq.t = sigma.sq.t, rho = rho)

# Subset data for prediction
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, , , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , , drop = FALSE]
# Prediction covariates
X.0 <- dat$X[pred.indx, , , drop = FALSE]
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , , drop = FALSE]
# Coordinates
coords <- dat$coords[-pred.indx, ]
coords.0 <- dat$coords[pred.indx, ]

occ.covs <- list(occ.cov.1 = X[, , 2],
		 occ.cov.2 = X[, , 3])
det.covs <- list(det.cov.1 = X.p[, , , 2],
		 det.cov.2 = X.p[, , , 3])

data.list <- list(y = y, occ.covs = occ.covs,
                  det.covs = det.covs,
                  coords = coords)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
		   rho.unif = list(a = -1, b = 1),
		   sigma.sq.t.ig = list(a = 0.1, b = 0.1),
                   phi.unif = list(a = 3 / .9, b = 3 / .1))
z.init <- apply(y, c(1, 2, 3), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(alpha.comm = 0, beta.comm = 0, beta = 0,
		   alpha = 0, tau.sq.beta = 1, tau.sq.alpha = 1,
		   rho = 0.5, sigma.sq.t = 0.5,
		   phi = 3 / .5, z = z.init)
# Tuning
tuning.list <- list(phi = 1, rho = 0.5)

# Number of batches
n.batch <- 5
# Batch length
batch.length <- 25
n.burn <- 25
n.thin <- 1
n.samples <- n.batch * batch.length

out <- svcTMsPGOcc(occ.formula = ~ occ.cov.1 + occ.cov.2,
                det.formula = ~ det.cov.1 + det.cov.2,
                data = data.list,
                inits = inits.list,
                n.batch = n.batch,
                batch.length = batch.length,
                accept.rate = 0.43,
		ar1 = TRUE,
		svc.cols = svc.cols,
		NNGP = TRUE,
		n.neighbors = 5,
		n.factors = n.factors,
		cov.model = 'exponential',
                priors = prior.list,
                tuning = tuning.list,
                n.omp.threads = 1,
                verbose = TRUE,
                n.report = 1,
                n.burn = n.burn,
		n.thin = n.thin,
		n.chains = 1)

summary(out)

# Predict at new sites across all n.max.years
# Take a look at array of covariates for prediction
str(X.0)
# Subset to only grab time periods 1, 2, and 5
t.cols <- c(1, 2, 5)
X.pred <- X.0[, t.cols, ]
out.pred <- predict(out, X.pred, coords.0, t.cols = t.cols, type = 'occupancy')
str(out.pred)

# Extract SVC samples for each species at prediction locations
svc.samples <- getSVCSamples(out, out.pred)
str(svc.samples)
}
