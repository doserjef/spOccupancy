\name{predict.msPGOcc}
\alias{predict.msPGOcc}
\title{Function for prediction at new locations for multi-species occupancy models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new locations given an object of class `msPGOcc`. Prediction is possible for both the latent occupancy state as well as detection.
}

\usage{
\method{predict}{msPGOcc}(object, X.0, ignore.RE = FALSE, type = 'occupancy', ...)
}

\arguments{

  \item{object}{an object of class msPGOcc}

  \item{X.0}{the design matrix of covariates at the prediction locations. This should include a column of 1s for the intercept if an intercept is included in the model. If random effects are included in the occupancy (or detection if \code{type = 'detection'}) portion of the model, the levels of the random effects at the new locations should be included as a column in the design matrix. The ordering of the levels should match the ordering used to fit the data in \code{msPGOcc}. Columns should correspond to the order of how covariates were specified in the corresponding formula argument of \code{msPGOcc}. Column names of the random effects must match the name of the random effects, if specified in the corresponding formula argument of \code{msPGOcc}.}

  \item{ignore.RE}{a logical value indicating whether to include unstructured random effects for prediction. If TRUE, random effects will be ignored and prediction will only use the fixed effects. If FALSE, random effects will be included in the prediction for both observed and unobserved levels of the random effect.}

  \item{...}{currently no additional arguments}

  \item{type}{a quoted keyword indicating what type of prediction to produce. Valid keywords are 'occupancy' to predict latent occupancy probability and latent occupancy values (this is the default), or 'detection' to predict detection probability given new values of detection covariates.}

}

\note{
  When \code{ignore.RE = FALSE}, both sampled levels and non-sampled levels of random effects are supported for prediction. For sampled levels, the posterior distribution for the random intercept corresponding to that level of the random effect will be used in the prediction. For non-sampled levels, random values are drawn from a normal distribution using the posterior samples of the random effect variance, which results in fully propagated uncertainty in predictions with models that incorporate random effects.
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  A list object of class \code{predict.msPGOcc}. When \code{type = 'occupancy'}, the list consists of:  

  \item{psi.0.samples}{a three-dimensional array of posterior predictive samples for the
    latent occurrence probability values.}

  \item{z.0.samples}{a three-dimensional array of posterior predictive samples for the
    latent occurrence values.}

  When \code{type = 'detection'}, the list consists of: 

  \item{p.0.samples}{a three-dimensional array of posterior predictive samples for the 
    detection probability values.}

  The return object will include additional objects used for standard 
  extractor functions. 
}

\examples{
set.seed(400)
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 0.5)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 0.3)
# Detection
alpha.mean <- c(0.5, 0.2, -0.1)
tau.sq.alpha <- c(0.2, 0.3, 1)
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

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
                sp = FALSE)
n.samples <- 5000
# Split into fitting and prediction data set
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, ]
# Occupancy covariates
X <- dat$X[-pred.indx, ]
# Detection covariates
X.p <- dat$X.p[-pred.indx, , ]
# Prediction values
X.0 <- dat$X[pred.indx, ]
psi.0 <- dat$psi[, pred.indx]
# Package all data into a list
occ.covs <- X[, 2, drop = FALSE]
colnames(occ.covs) <- c('occ.cov')
det.covs <- list(det.cov.1 = X.p[, , 2], 
                 det.cov.2 = X.p[, , 3])
data.list <- list(y = y, 
                  occ.covs = occ.covs,
                  det.covs = det.covs)

# Occupancy initial values
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72), 
                   alpha.comm.normal = list(mean = 0, var = 2.72), 
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Initial values
inits.list <- list(alpha.comm = 0, 
                   beta.comm = 0, 
                   beta = 0, 
                   alpha = 0,
                   tau.sq.beta = 1, 
                   tau.sq.alpha = 1, 
                   z = apply(y, c(1, 2), max, na.rm = TRUE))

out <- msPGOcc(occ.formula = ~ occ.cov, 
               det.formula = ~ det.cov.1 + det.cov.2, 
               data = data.list, 
               inits = inits.list, 
               n.samples = n.samples, 
               priors = prior.list, 
               n.omp.threads = 1, 
               verbose = TRUE, 
               n.report = 1000, 
               n.burn = 4000)

summary(out, level = 'community')

# Predict at new locations ------------------------------------------------
out.pred <- predict(out, X.0)
}
