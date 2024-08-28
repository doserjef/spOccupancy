\name{predict.intMsPGOcc}
\alias{predict.intMsPGOcc}
\title{Function for prediction at new locations for integrated multi-species occupancy models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new locations given an object of class `intMsPGOcc`. Prediction is currently possible only for the latent occupancy state.
}

\usage{
\method{predict}{intMsPGOcc}(object, X.0, ignore.RE = FALSE, ...)
}

\arguments{

  \item{object}{an object of class intMsPGOcc}

  \item{X.0}{the design matrix of covariates at the prediction locations. This should include a column of 1s for the intercept if an intercept is included in the model. If random effects are included in the occupancy (or detection if \code{type = 'detection'}) portion of the model, the levels of the random effects at the new locations should be included as a column in the design matrix. The ordering of the levels should match the ordering used to fit the data in \code{intMsPGOcc}. Columns should correspond to the order of how covariates were specified in the corresponding formula argument of \code{intMsPGOcc}. Column names of the random effects must match the name of the random effects, if specified in the corresponding formula argument of \code{intMsPGOcc}.}

  \item{ignore.RE}{a logical value indicating whether to include unstructured random effects for prediction. If TRUE, random effects will be ignored and prediction will only use the fixed effects. If FALSE, random effects will be included in the prediction for both observed and unobserved levels of the random effect.}

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
  A list object of class \code{predict.intMsPGOcc} consisting of:   

  \item{psi.0.samples}{a three-dimensional array of posterior predictive samples for the
    latent occurrence probability values.}

  \item{z.0.samples}{a three-dimensional array of posterior predictive samples for the
    latent occurrence values.}

  The return object will include additional objects used for standard 
  extractor functions. 
}

\examples{
set.seed(91)
J.x <- 10
J.y <- 10
# Total number of data sources across the study region
J.all <- J.x * J.y
# Number of data sources.
n.data <- 2
# Sites for each data source.
J.obs <- sample(ceiling(0.2 * J.all):ceiling(0.5 * J.all), n.data, replace = TRUE)
n.rep <- list()
n.rep[[1]] <- rep(3, J.obs[1])
n.rep[[2]] <- rep(4, J.obs[2])

# Number of species observed in each data source
N <- c(8, 3)

# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 0.5)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.4, 0.3)
# Detection
# Detection covariates
alpha.mean <- list()
tau.sq.alpha <- list()
# Number of detection parameters in each data source
p.det.long <- c(4, 3)
for (i in 1:n.data) {
  alpha.mean[[i]] <- runif(p.det.long[i], -1, 1)
  tau.sq.alpha[[i]] <- runif(p.det.long[i], 0.1, 1)
}
# Random effects
psi.RE <- list()
p.RE <- list()
beta <- matrix(NA, nrow = max(N), ncol = p.occ)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(max(N), beta.mean[i], sqrt(tau.sq.beta[i]))
}
alpha <- list()
for (i in 1:n.data) {
  alpha[[i]] <- matrix(NA, nrow = N[i], ncol = p.det.long[i])
  for (t in 1:p.det.long[i]) {
    alpha[[i]][, t] <- rnorm(N[i], alpha.mean[[i]][t], sqrt(tau.sq.alpha[[i]])[t])
  }
}
sp <- FALSE
factor.model <- FALSE
# Simulate occupancy data
dat <- simIntMsOcc(n.data = n.data, J.x = J.x, J.y = J.y,
                   J.obs = J.obs, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
                   psi.RE = psi.RE, p.RE = p.RE, sp = sp, factor.model = factor.model,
                   n.factors = n.factors)
J <- nrow(dat$coords.obs)
y <- dat$y
X <- dat$X.obs
X.p <- dat$X.p
X.re <- dat$X.re.obs
X.p.re <- dat$X.p.re
sites <- dat$sites
species <- dat$species

# Package all data into a list
occ.covs <- cbind(X)
colnames(occ.covs) <- c('int', 'occ.cov.1')
#colnames(occ.covs) <- c('occ.cov')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(det.cov.1.1 = X.p[[1]][, , 2],
                      det.cov.1.2 = X.p[[1]][, , 3],
                      det.cov.1.3 = X.p[[1]][, , 4])
det.covs[[2]] <- list(det.cov.2.1 = X.p[[2]][, , 2],
                      det.cov.2.2 = X.p[[2]][, , 3])

data.list <- list(y = y,
                  occ.covs = occ.covs,
                  det.covs = det.covs,
                  sites = sites,
                  species = species)
# Take a look at the data.list structure for integrated multi-species
# occupancy models.
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.73),
                   alpha.comm.normal = list(mean = list(0, 0),
                                            var = list(2.72, 2.72)),
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                   tau.sq.alpha.ig = list(a = list(0.1, 0.1),
                                          b = list(0.1, 0.1)))
inits.list <- list(alpha.comm = list(0, 0),
                   beta.comm = 0,
                   tau.sq.beta = 1,
                   tau.sq.alpha = list(1, 1),
                   alpha = list(a = matrix(rnorm(p.det.long[1] * N[1]), N[1], p.det.long[1]),
                                b = matrix(rnorm(p.det.long[2] * N[2]), N[2], p.det.long[2])),
                   beta = 0)

# Fit the model. 
# Note that this is just a test case and more iterations/chains may need to 
# be run to ensure convergence.
out <- intMsPGOcc(occ.formula = ~ occ.cov.1,
                  det.formula = list(f.1 = ~ det.cov.1.1 + det.cov.1.2 + det.cov.1.3,
                                     f.2 = ~ det.cov.2.1 + det.cov.2.2),
                  inits = inits.list,
                  priors = prior.list,
                  data = data.list,
                  n.samples = 100,
                  n.omp.threads = 1,
                  verbose = TRUE,
                  n.report = 10,
                  n.burn = 50,
                  n.thin = 1,
                  n.chains = 1)
#Predict at new locations. 
X.0 <- dat$X.pred
psi.0 <- dat$psi.pred
out.pred <- predict(out, X.0, ignore.RE = TRUE)

# Create prediction for one species. 
curr.sp <- 2
psi.hat.quants <- apply(out.pred$psi.0.samples[,curr.sp, ], 
                        2, quantile, c(0.025, 0.5, 0.975))
plot(psi.0[curr.sp, ], psi.hat.quants[2, ], pch = 19, xlab = 'True',
     ylab = 'Predicted', ylim = c(min(psi.hat.quants), max(psi.hat.quants)), 
     main = paste("Species ", curr.sp, sep = ''))
segments(psi.0[curr.sp, ], psi.hat.quants[1, ], psi.0[curr.sp, ], psi.hat.quants[3, ])
lines(psi.0[curr.sp, ], psi.0[curr.sp, ])
}
