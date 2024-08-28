\name{predict.tPGOcc}
\alias{predict.tPGOcc}
\title{Function for prediction at new locations for multi-season single-species occupancy models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new locations given an object of class `tPGOcc`. Prediction is possible for both the latent occupancy state as well as detection. Predictions are currently only possible for sampled primary time periods. 
}

\usage{
\method{predict}{tPGOcc}(object, X.0, t.cols, ignore.RE = FALSE, type = 'occupancy', ...)
}

\arguments{

  \item{object}{an object of class tPGOcc}

  \item{X.0}{the design matrix of covariates at the prediction locations. This should be a three-dimensional array, with dimensions corresponding to site, primary time period, and covariate, respectively. Note that the first covariate should consist of all 1s for the intercept if an intercept is included in the model. If random effects are included in the occupancy (or detection if \code{type = 'detection'}) portion of the model, the levels of the random effects at the new locations/time periods should be included as an element of the three-dimensional array. The ordering of the levels should match the ordering used to fit the data in \code{tPGOcc}. The covariates should be organized in the same order as they were specified in the corresponding formula argument of \code{tPGOcc}. Names of the third dimension (covariates) of any random effects in X.0 must match the name of the random effects used to fit the model, if specified in the corresponding formula argument of \code{tPGOcc}. See example below.}

  \item{t.cols}{an indexing vector used to denote which primary time periods are contained in the design matrix of covariates at the prediction locations (\code{X.0}). The values should denote the specific primary time periods used to fit the model. The values should indicate the columns in \code{data$y} used to fit the model for which prediction is desired. See example below.} 

  \item{ignore.RE}{logical value that specifies whether or not to remove random unstructured occurrence (or detection if \code{type = 'detection'}) effects from the subsequent predictions. If \code{TRUE}, unstructured random effects will be included. If \code{FALSE}, unstructured random effects will be set to 0 and predictions will only be generated from the fixed effects and AR(1) random effects if the model was fit with \code{ar1 = TRUE}.}

  \item{type}{a quoted keyword indicating what type of prediction to produce. Valid keywords are 'occupancy' to predict latent occupancy probability and latent occupancy values (this is the default), or 'detection' to predict detection probability given new values of detection covariates.}

  \item{...}{currently no additional arguments}
}

\note{
  When \code{ignore.RE = FALSE}, both sampled levels and non-sampled levels of unstructured random effects are supported for prediction. For sampled levels, the posterior distribution for the random intercept corresponding to that level of the random effect will be used in the prediction. For non-sampled levels, random values are drawn from a normal distribution using the posterior samples of the random effect variance, which results in fully propagated uncertainty in predictions with models that incorporate random effects. 

  Occurrence predictions at sites that are only sampled for a subset of the total number of primary time periods are obtained directly when fitting the model. See the \code{psi.samples} and \code{z.samples} portions of the output list from the model object of class \code{tPGOcc}. 
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  A list object of class \code{predict.tPGOcc}. When \code{type = 'occupancy'}, the list consists of: 

  \item{psi.0.samples}{a three-dimensional object of posterior predictive samples for the
    latent occupancy probability values with dimensions corresponding to posterior predictive
    sample, site, and primary time period.}

  \item{z.0.samples}{a three-dimensional object of posterior predictive samples for the
    latent occupancy values with dimensions corresponding to posterior predictive sample, site, 
    and primary time period.}

  When \code{type = 'detection'}, the list consists of: 

  \item{p.0.samples}{a three-dimensional object of posterior predictive samples for the 
    detection probability values with dimensions corresponding to posterior predictive 
    sample, site, and primary time period.}

  The return object will include additional objects used for standard 
  extractor functions. 
}

\examples{
set.seed(990)
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

# Get all the data
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, 
               beta = beta, alpha = alpha, sp.only = sp.only, trend = trend, 
               psi.RE = psi.RE, p.RE = p.RE, sp = FALSE, ar1 = FALSE)

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
                  det.covs = det.covs) 
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                   alpha.normal = list(mean = 0, var = 2.72))

# Starting values
z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(beta = 0, alpha = 0, z = z.init)

n.batch <- 100
batch.length <- 25
n.burn <- 2000
n.thin <- 1

# Run the model
# Note that this is just a test case and more iterations/chains may need to 
# be run to ensure convergence.
out <- tPGOcc(occ.formula = ~ trend + occ.cov.1, 
              det.formula = ~ det.cov.1 + det.cov.2, 
              data = data.list,
              inits = inits.list,
              priors = prior.list, 
              n.batch = n.batch,
              batch.length = batch.length,
              ar1 = FALSE,
              verbose = TRUE, 
              n.report = 500,
              n.burn = n.burn, 
              n.thin = n.thin,
              n.chains = 1) 

# Predict at new sites across during time periods 1, 2, and 5
# Take a look at array of covariates for prediction
str(X.0)
# Subset to only grab time periods 1, 2, and 5
t.cols <- c(1, 2, 5)
X.pred <- X.0[, t.cols, ]
out.pred <- predict(out, X.pred, t.cols = t.cols, type = 'occupancy')
str(out.pred)
}
