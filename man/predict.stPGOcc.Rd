\name{predict.stPGOcc}
\alias{predict.stPGOcc}
\title{Function for prediction at new locations for multi-season single-species spatial occupancy models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new locations given an object of class `stPGOcc`. Prediction is possible for both the latent occupancy state as well as detection. Predictions are currently only possible for sampled primary time periods.  
}

\usage{
\method{predict}{stPGOcc}(object, X.0, coords.0, t.cols, n.omp.threads = 1, 
                          verbose = TRUE, n.report = 100, 
                          ignore.RE = FALSE, type = 'occupancy', 
                          forecast = FALSE, grid.index.0, ...)
}

\arguments{

  \item{object}{an object of class stPGOcc}

  \item{X.0}{the design matrix of covariates at the prediction locations. This should be a three-dimensional array, with dimensions corresponding to site, primary time period, and covariate, respectively. Note that the first covariate should consist of all 1s for the intercept if an intercept is included in the model. If random effects are included in the occupancy (or detection if \code{type = 'detection'}) portion of the model, the levels of the random effects at the new locations/time periods should be included as an element of the three-dimensional array. The ordering of the levels should match the ordering used to fit the data in \code{stPGOcc}. The covariates should be organized in the same order as they were specified in the corresponding formula argument of \code{stPGOcc}. Names of the third dimension (covariates) of any random effects in X.0 must match the name of the random effects used to fit the model, if specified in the corresponding formula argument of \code{stPGOcc}. See example below.}

  \item{coords.0}{the spatial coordinates corresponding to \code{X.0}. Note that \code{spOccupancy} assumes coordinates are specified 
  in a projected coordinate system.}

  \item{t.cols}{an indexing vector used to denote which primary time periods are contained in the design matrix of covariates at the prediction locations (\code{X.0}). The values should denote the specific primary time periods used to fit the model. The values should indicate the columns in \code{data$y} used to fit the model for which prediction is desired. See example below. Not required when \code{forecast = TRUE}.} 

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

  \item{forecast}{a logical value indicating whether prediction is occurring at non-sampled primary time periods (e.g., forecasting).}

  \item{grid.index.0}{an indexing vector used to specify how each row in \code{X.0} corresponds to the coordinates specified in \code{coords.0}. Only relevant if the spatial random effect was estimated at a higher spatial resolution (e.g., grid cells) than point locations.}

  \item{...}{currently no additional arguments}
}

\note{
  When \code{ignore.RE = FALSE}, both sampled levels and non-sampled levels of unstructured random effects are supported for prediction. For sampled levels, the posterior distribution for the random intercept corresponding to that level of the random effect will be used in the prediction. For non-sampled levels, random values are drawn from a normal distribution using the posterior samples of the random effect variance, which results in fully propagated uncertainty in predictions with models that incorporate random effects. 

  Occurrence predictions at sites that are only sampled for a subset of the total number of primary time periods are obtained directly when fitting the model. See the \code{psi.samples} and \code{z.samples} portions of the output list from the model object of class \code{stPGOcc}. 
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  A list object of class \code{predict.stPGOcc}. When \code{type = 'occupancy'}, the list consists of: 

  \item{psi.0.samples}{a three-dimensional object of posterior predictive samples for the
    latent occupancy probability values with dimensions corresponding to posterior predictive
    sample, site, and primary time period.}

  \item{z.0.samples}{a three-dimensional object of posterior predictive samples for the
    latent occupancy values with dimensions corresponding to posterior predictive sample, site, 
    and primary time period.}

  \item{w.0.samples}{a \code{coda} object of posterior predictive samples 
    for the latent spatial random effects.}

  When \code{type = 'detection'}, the list consists of: 

  \item{p.0.samples}{a three-dimensional object of posterior predictive samples for the 
    detection probability values with dimensions corresponding to posterior predictive 
    sample, site, and primary time period.}

  The return object will include additional objects used for standard 
  extractor functions. 
}

\examples{
set.seed(500)
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
# Primary time periods
n.time <- sample(10, J, replace = TRUE)
n.time.max <- max(n.time)
# Replicates
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
# Occurrence --------------------------
beta <- c(0.4, 0.5, -0.9)
trend <- TRUE 
sp.only <- 0
psi.RE <- list()
# Detection ---------------------------
alpha <- c(-1, 0.7, -0.5)
p.RE <- list()
# Spatial -----------------------------
sp <- TRUE
cov.model <- "exponential"
sigma.sq <- 2
phi <- 3 / .4

# Get all the data
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, 
               beta = beta, alpha = alpha, sp.only = sp.only, trend = trend, 
               psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, sigma.sq = sigma.sq, 
               phi = phi, cov.model = cov.model, ar1 = FALSE)

# Subset data for prediction
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , , drop = FALSE]
# Prediction covariates
X.0 <- dat$X[pred.indx, , , drop = FALSE]
# Detection covariates
X.p <- dat$X.p[-pred.indx, , , , drop = FALSE]
psi.0 <- dat$psi[pred.indx, ]
# Coordinates
coords <- dat$coords[-pred.indx, ]
coords.0 <- dat$coords[pred.indx, ]

# Package all data into a list
# Occurrence
occ.covs <- list(int = X[, , 1], 
                 trend = X[, , 2], 
                 occ.cov.1 = X[, , 3]) 
# Detection
det.covs <- list(det.cov.1 = X.p[, , , 2], 
                 det.cov.2 = X.p[, , , 3]) 
# Data list bundle
data.list <- list(y = y, 
                  occ.covs = occ.covs,
                  det.covs = det.covs, 
                  coords = coords) 
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                   alpha.normal = list(mean = 0, var = 2.72), 
                   sigma.sq.ig = c(2, 2), 
                   phi.unif = c(3 / 1, 3 / 0.1))

# Initial values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(beta = 0, alpha = 0, z = z.init, phi = 3 / .5, sigma.sq = 2, 
                   w = rep(0, J))
# Tuning
tuning.list <- list(phi = 1)
# Number of batches
n.batch <- 10
# Batch length
batch.length <- 25
n.iter <- n.batch * batch.length

# Run the model
out <- stPGOcc(occ.formula = ~ trend + occ.cov.1, 
               det.formula = ~ det.cov.1 + det.cov.2, 
               data = data.list, 
               inits = inits.list, 
               n.batch = n.batch, 
               batch.length = batch.length, 
               priors = prior.list,
               cov.model = "exponential", 
               tuning = tuning.list, 
               NNGP = TRUE, 
               ar1 = FALSE,
               n.neighbors = 5, 
               search.type = 'cb', 
               n.report = 10, 
               n.burn = 50, 
               n.chains = 1)

summary(out)

# Predict at new sites across all n.max.years
# Take a look at array of covariates for prediction
str(X.0)
# Subset to only grab time periods 1, 2, and 5
t.cols <- c(1, 2, 5)
X.pred <- X.0[, t.cols, ]
out.pred <- predict(out, X.0, coords.0, t.cols = t.cols, type = 'occupancy')
str(out.pred)
}
