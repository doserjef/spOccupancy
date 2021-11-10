\name{predict.PGOcc}
\alias{predict.PGOcc}
\title{Function for prediction at new locations for single species occupancy models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new locations given an object of class `PGOcc`. When non-spatial random effects are included in the occurrence portion of the model, \code{predict} only allows prediction at random effect levels that are included in the fitted data.
}

\usage{
\method{predict}{PGOcc}(object, X.0, ...)
}

\arguments{

  \item{object}{an object of class PGOcc}

  \item{X.0}{the design matrix for prediction locations. This should include a column of 1s for the intercept. If random effects are included in the occurrence portion of the model, the levels of the random effects at the new locations should be included as a column in the design matrix. The ordering of the levels should match the ordering used to fit the data in \code{PGOcc}. Covariates should have the same column names as those used when fitting the model with \code{PGOcc}.}

  \item{...}{currently no additional arguments}
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{predict.PGOcc} that is a list comprised of: 

  \item{psi.0.samples}{a \code{coda} object of posterior predictive samples for the
    latent occupancy probability values.}

  \item{z.0.samples}{a \code{coda} object of posterior predictive samples for the
    latent occupancy values.}

  The return object will include additional objects used for standard 
  extractor functions. 
}

\examples{
set.seed(400)
# Simulate Data -----------------------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.5, 2)
p.occ <- length(beta)
alpha <- c(0, 1)
p.det <- length(alpha)
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      sp = FALSE)
# Split into fitting and prediction data set
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, ]
# Occupancy covariates
X <- dat$X[-pred.indx, ]
# Prediction covariates
X.0 <- dat$X[pred.indx, ]
# Detection covariates
X.p <- dat$X.p[-pred.indx, , ]

# Package all data into a list
occ.covs <- X[, 2, drop = FALSE]
colnames(occ.covs) <- c('occ.cov')
det.covs <- list(det.cov = X.p[, , 2])
data.list <- list(y = y, 
		  occ.covs = occ.covs,
		  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = rep(0, p.occ),
				      var = rep(2.72, p.occ)),
		   alpha.normal = list(mean = rep(0, p.det),
			               var = rep(2.72, p.det)))
# Initial values
inits.list <- list(alpha = rep(0, p.det),
		      beta = rep(0, p.occ),
		      z = apply(y, 1, max, na.rm = TRUE))

n.samples <- 5000
n.report <- 1000

out <- PGOcc(occ.formula = ~ occ.cov, 
	     det.formula = ~ det.cov,
	     data = data.list, 
	     inits = inits.list,
	     n.samples = n.samples,
	     priors = prior.list,
	     n.omp.threads = 1,
	     verbose = TRUE,
	     n.report = n.report, 
	     n.burn = 4000, 
	     n.thin = 1)

summary(out)

# Predict at new locations ------------------------------------------------

out.pred <- predict(out, X.0)
psi.0.quants <- apply(out.pred$psi.0.samples, 2, quantile, c(0.025, 0.5, 0.975))
plot(dat$psi[pred.indx], psi.0.quants[2, ], pch = 19, xlab = 'True', 
     ylab = 'Fitted', ylim = c(min(psi.0.quants), max(psi.0.quants)))
segments(dat$psi[pred.indx], psi.0.quants[1, ], dat$psi[pred.indx], psi.0.quants[3, ])
lines(dat$psi[pred.indx], dat$psi[pred.indx])
}
