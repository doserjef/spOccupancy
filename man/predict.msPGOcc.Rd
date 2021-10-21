\name{predict.msPGOcc}
\alias{predict.msPGOcc}
\title{Function for prediction at new locations for multi-species occupancy models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new locations given an object of class `msPGOcc`. 
}

\usage{
\method{predict}{msPGOcc}(object, X.0, ...)
}

\arguments{

  \item{object}{an object of class msPGOcc}

  \item{X.0}{the design matrix for prediction locations. This should include a column of 1s for the intercept.}

  \item{...}{currently no additional arguments}
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{predict.msPGOcc} that is a list comprised of: 

  \item{psi.0.samples}{a three-dimensional array of posterior predictive samples for the
    latent occupancy probability values.}

  \item{z.0.samples}{a three-dimensional array of posterior predictive samples for the
    latent occupancy values.}

  The return object will include additional objects used for standard 
  extractor functions. 
}

\examples{
\dontrun{
set.seed(400)
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 10
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
		 det.cov.2 = X.p[, , 3]
		 )
data.list <- list(y = y, 
		  occ.covs = occ.covs,
		  det.covs = det.covs)

# Occupancy initial values
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72), 
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Starting values
starting.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE))

out <- msPGOcc(occ.formula = ~ occ.cov, 
	       det.formula = ~ det.cov.1 + det.cov.2, 
	       data = data.list, 
	       starting = starting.list, 
	       n.samples = n.samples, 
	       priors = prior.list, 
               n.omp.threads = 1, 
	       verbose = TRUE, 
	       n.report = 1000)

summary(out, level = 'community')

# Predict at new locations ------------------------------------------------

out.pred <- predict(out, X.0)
# Create a plot one species at a time: 
for (i in 1:N) {
  psi.hat.quants <- apply(out.pred$psi.0.samples[,i, ], 
  			2, quantile, c(0.025, 0.5, 0.975))
  plot(psi.0[i, ], psi.hat.quants[2, ], pch = 19, xlab = 'True',
       ylab = 'Predicted', ylim = c(min(psi.hat.quants), max(psi.hat.quants)), 
       main = paste("Species ", i, sep = ''))
  segments(psi.0[i, ], psi.hat.quants[1, ], psi.0[i, ], psi.hat.quants[3, ])
  lines(psi.0[i, ], psi.0[i, ])
  Sys.sleep(1)
}
}
}
