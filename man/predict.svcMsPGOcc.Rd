\name{predict.svcMsPGOcc}
\alias{predict.svcMsPGOcc}
\title{Function for prediction at new locations for spatially varying coefficient multi-species occupancy models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new locations given an object of class `svcMsPGOcc`. Prediction is possible for both the latent occupancy state as well as detection.
}

\usage{
\method{predict}{svcMsPGOcc}(object, X.0, coords.0, n.omp.threads = 1, verbose = TRUE, 
        n.report = 100, ignore.RE = FALSE, type = 'occupancy', ...)
}

\arguments{

  \item{object}{an object of class svcMsPGOcc}

  \item{X.0}{the design matrix of covariates at the prediction locations. This should include a column of 1s for the intercept if an intercept is included in the model. If random effects are included in the occupancy (or detection if \code{type = 'detection'}) portion of the model, the levels of the random effects at the new locations should be included as a column in the design matrix. The ordering of the levels should match the ordering used to fit the data in \code{svcMsPGOcc}. Columns should correspond to the order of how covariates were specified in the corresponding formula argument of \code{svcMsPGOcc}. Column names of the random effects must match the name of the random effects, if specified in the corresponding formula argument of \code{svcMsPGOcc}.}

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

  \item{type}{a quoted keyword indicating what type of prediction to produce. Valid keywords are 'occupancy' to predict latent occupancy probability and latent occupancy values (this is the default), or 'detection' to predict detection probability given new values of detection covariates.}

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
  An list object of class \code{predict.svcMsPGOcc}. When \code{type = 'occupancy'}, the list consists of:  

  \item{psi.0.samples}{a three-dimensional array of posterior predictive samples for the
    latent occurrence probability values.}

  \item{z.0.samples}{a three-dimensional array of posterior predictive samples for the
    latent occurrence values.}

  \item{w.0.samples}{a four-dimensional array of posterior predictive samples
    for the spatially-varying coefficients, with dimensions corresponding to 
    MCMC sample, spatial factor, site, and spatially varying coefficient.}

  \item{run.time}{execution time reported using \code{proc.time()}.}

  When \code{type = 'detection'}, the list consists of: 

  \item{p.0.samples}{a three-dimensional array of posterior predictive samples for the 
    detection probability values.}
  
  \item{run.time}{execution time reported using \code{proc.time()}.}

  The return object will include additional objects used for standard 
  extractor functions. 
}

\examples{
set.seed(400)

# Simulate Data -----------------------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(5, size = J, replace = TRUE)
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, -0.2, 0.3, -0.1, 0.4)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 1.5, 0.4, 0.5, 0.3)
# Detection
alpha.mean <- c(0, 1.2, -0.5)
tau.sq.alpha <- c(1, 0.5, 1.3)
p.det <- length(alpha.mean)
# No random effects
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
# Number of spatial factors for each SVC
n.factors <- 2
# The intercept and first two covariates have spatially-varying effects
svc.cols <- c(1, 2, 3)
p.svc <- length(svc.cols)
q.p.svc <- n.factors * p.svc
# Spatial decay parameters
phi <- runif(q.p.svc, 3 / 0.9, 3 / 0.1)
# A length N vector indicating the proportion of simulated locations
# that are within the range for a given species.
range.probs <- runif(N, 1, 1)
factor.model <- TRUE
cov.model <- 'spherical'
sp <- TRUE

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
		psi.RE = psi.RE, p.RE = p.RE, phi = phi, sp = sp, svc.cols = svc.cols,
		cov.model = cov.model, n.factors = n.factors,
		factor.model = factor.model, range.probs = range.probs)

# Split into fitting and prediction data set
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, ]
# Occupancy covariates
X <- dat$X[-pred.indx, ]
# Coordinates
coords <- as.matrix(dat$coords[-pred.indx, ])
# Detection covariates
X.p <- dat$X.p[-pred.indx, , ]
# Prediction values
X.0 <- dat$X[pred.indx, ]
coords.0 <- as.matrix(dat$coords[pred.indx, ])

# Prep data for spOccupancy -----------------------------------------------
# Occurrence covariates
occ.covs <- cbind(X)
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2', 'occ.cov.3',
			'occ.cov.4')
# Detection covariates
det.covs <- list(det.cov.1 = X.p[, , 2],
		 det.cov.2 = X.p[, , 3])
# Data list
data.list <- list(y = y, coords = coords, occ.covs = occ.covs,
                  det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                   phi.unif = list(a = 3 / 1, b = 3 / .1))
inits.list <- list(alpha.comm = 0,
		   beta.comm = 0,
		   beta = 0,
		   alpha = 0,
		   tau.sq.beta = 1,
		   tau.sq.alpha = 1,
		   z = apply(y, c(1, 2), max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1)

# Number of batches
n.batch <- 2
# Batch length
batch.length <- 25
n.burn <- 0
n.thin <- 1
n.samples <- n.batch * batch.length

out <- svcMsPGOcc(occ.formula = ~ occ.cov.1 + occ.cov.2 + occ.cov.3 +
                                  occ.cov.4,
                  det.formula = ~ det.cov.1 + det.cov.2,
                  data = data.list,
                  inits = inits.list,
                  n.batch = n.batch,
		  n.factors = n.factors,
                  batch.length = batch.length,
		  std.by.sp = TRUE,
                  accept.rate = 0.43,
                  priors = prior.list,
		  svc.cols = svc.cols,
                  cov.model = "spherical",
                  tuning = tuning.list,
                  n.omp.threads = 1,
                  verbose = TRUE,
                  NNGP = TRUE,
                  n.neighbors = 5,
                  search.type = 'cb',
                  n.report = 10,
                  n.burn = n.burn,
		  n.thin = n.thin,
		  n.chains = 1)

summary(out)
# Predict at new locations ------------------------------------------------
out.pred <- predict(out, X.0, coords.0, verbose = FALSE)

# Get SVC samples for each species at prediction locations
svc.samples <- getSVCSamples(out, out.pred)
}
