\name{predict.tIntPGOcc}
\alias{predict.tIntPGOcc}
\title{Function for prediction at new locations for multi-season single-species integrated occupancy models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new locations given an object of class `tIntPGOcc`. Prediction is currently only possible for the latent occupancy state. Predictions are currently only possible for sampled primary time periods. 
}

\usage{
\method{predict}{tIntPGOcc}(object, X.0, t.cols, ignore.RE = FALSE, type = 'occupancy', ...)
}

\arguments{

  \item{object}{an object of class tIntPGOcc}

  \item{X.0}{the design matrix of covariates at the prediction locations. This should be a three-dimensional array, with dimensions corresponding to site, primary time period, and covariate, respectively. Note that the first covariate should consist of all 1s for the intercept if an intercept is included in the model. If random effects are included in the occupancy (or detection if \code{type = 'detection'}) portion of the model, the levels of the random effects at the new locations/time periods should be included as an element of the three-dimensional array. The ordering of the levels should match the ordering used to fit the data in \code{tIntPGOcc}. The covariates should be organized in the same order as they were specified in the corresponding formula argument of \code{tIntPGOcc}. Names of the third dimension (covariates) of any random effects in X.0 must match the name of the random effects used to fit the model, if specified in the corresponding formula argument of \code{tIntPGOcc}. See example below.}

  \item{t.cols}{an indexing vector used to denote which primary time periods are contained in the design matrix of covariates at the prediction locations (\code{X.0}). The values should denote the specific primary time periods used to fit the model. The values should indicate the columns in \code{data$y} used to fit the model for which prediction is desired. See example below.} 

  \item{ignore.RE}{logical value that specifies whether or not to remove random unstructured occurrence (or detection if \code{type = 'detection'}) effects from the subsequent predictions. If \code{TRUE}, unstructured random effects will be included. If \code{FALSE}, unstructured random effects will be set to 0 and predictions will only be generated from the fixed effects and AR(1) random effects if the model was fit with \code{ar1 = TRUE}.}

  \item{type}{a quoted keyword indicating what type of prediction to produce. Valid keywords are 'occupancy' to predict latent occupancy probability and latent occupancy values (this is the default), or 'detection' to predict detection probability given new values of detection covariates. Detection prediction is not currently supported for integrated models.}

  \item{...}{currently no additional arguments}
}

\note{
  When \code{ignore.RE = FALSE}, both sampled levels and non-sampled levels of unstructured random effects are supported for prediction. For sampled levels, the posterior distribution for the random intercept corresponding to that level of the random effect will be used in the prediction. For non-sampled levels, random values are drawn from a normal distribution using the posterior samples of the random effect variance, which results in fully propagated uncertainty in predictions with models that incorporate random effects. 

  Occurrence predictions at sites that are only sampled for a subset of the total number of primary time periods are obtained directly when fitting the model. See the \code{psi.samples} and \code{z.samples} portions of the output list from the model object of class \code{tIntPGOcc}. 
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}
}

\value{
  A list object of class \code{predict.tIntPGOcc}. When \code{type = 'occupancy'}, the list consists of: 

  \item{psi.0.samples}{a three-dimensional object of posterior predictive samples for the
    latent occupancy probability values with dimensions corresponding to posterior predictive
    sample, site, and primary time period.}

  \item{z.0.samples}{a three-dimensional object of posterior predictive samples for the
    latent occupancy values with dimensions corresponding to posterior predictive sample, site, 
    and primary time period.}

  The return object will include additional objects used for standard 
  extractor functions. 
}

\examples{
set.seed(332)

# Simulate Data -----------------------------------------------------------
# Number of locations in each direction. This is the total region of interest
# where some sites may or may not have a data source.
J.x <- 15 
J.y <- 15
J.all <- J.x * J.y
# Number of data sources.
n.data <- 3
# Sites for each data source.
J.obs <- sample(ceiling(0.2 * J.all):ceiling(0.4 * J.all), n.data, replace = TRUE)
# Maximum number of years for each data set
n.time.max <- c(4, 8, 10)
# Number of years each site in each data set is sampled
n.time <- list()
for (i in 1:n.data) {
  n.time[[i]] <- sample(1:n.time.max[i], J.obs[i], replace = TRUE)
}
# Replicates for each data source.
n.rep <- list()
for (i in 1:n.data) {
  n.rep[[i]] <- matrix(NA, J.obs[i], n.time.max[i])
  for (j in 1:J.obs[i]) {
    n.rep[[i]][j, sample(1:n.time.max[i], n.time[[i]][j], replace = FALSE)] <- 
      sample(1:4, n.time[[i]][j], replace = TRUE)
  }
}
# Total number of years across all data sets
n.time.total <- 10
# List denoting the specific years each data set was sampled during. 
data.seasons <- list()
for (i in 1:n.data) {
  data.seasons[[i]] <- sort(sample(1:n.time.total, n.time.max[i], replace = FALSE))
}

# Occupancy covariates
beta <- c(0, 0.4, 0.3)
trend <- TRUE
# Random occupancy covariates
psi.RE <- list()
p.occ <- length(beta)
# Detection covariates
alpha <- list()
alpha[[1]] <- c(0, 0.2, -0.5)
alpha[[2]] <- c(-1, 0.5, 0.3, -0.8)
alpha[[3]] <- c(-0.5, 1)

p.RE <- list()
p.det.long <- sapply(alpha, length)
p.det <- sum(p.det.long)

# Simulate occupancy data.
dat <- simTIntOcc(n.data = n.data, J.x = J.x, J.y = J.y, J.obs = J.obs,
                  n.time = n.time, data.seasons = data.seasons, n.rep = n.rep,
                  beta = beta, alpha = alpha, trend = trend, 
                  psi.RE = psi.RE, p.RE = p.RE)

y <- dat$y
X <- dat$X.obs
X.p <- dat$X.p
sites <- dat$sites

# Package all data into a list
occ.covs <- list(trend = X[, , 2], 
                 occ.cov.1 = X[, , 3])
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(det.cov.1.1 = X.p[[1]][, , , 2],
                      det.cov.1.2 = X.p[[1]][, , , 3])
det.covs[[2]] <- list(det.cov.2.1 = X.p[[2]][, , , 2],
                      det.cov.2.2 = X.p[[2]][, , , 3],
                      det.cov.2.3 = X.p[[2]][, , , 4])
det.covs[[3]] <- list(det.cov.3.1 = X.p[[3]][, , , 2])
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs,
                  sites = sites, seasons = data.seasons)

# Testing
occ.formula <- ~ trend + occ.cov.1
# Note that the names are not necessary.
det.formula <- list(f.1 = ~ det.cov.1.1 + det.cov.1.2,
                    f.2 = ~ det.cov.2.1 + det.cov.2.2 + det.cov.2.3,
                    f.3 = ~ det.cov.3.1)

# NOTE: this is a short run of the model, in reality we would run the 
#       model for much longer.
out <- tIntPGOcc(occ.formula = occ.formula,
                 det.formula = det.formula,
                 data = data.list,
                 n.batch = 3,
                 batch.length = 25, 
                 n.report = 1,
                 n.burn = 25,
                 n.thin = 1,
                 n.chains = 1)

t.cols <- 1:n.time.total
out.pred <- predict(out, X.0 = dat$X.pred, t.cols = t.cols, 
                    type = 'occupancy')
str(out.pred)
}
