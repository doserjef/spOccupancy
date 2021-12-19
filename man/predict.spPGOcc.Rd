\name{predict.spPGOcc}
\alias{predict.spPGOcc}
\title{Function for prediction at new locations for single-species spatial occupancy models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new 
  locations given an object of class `spPGOcc`. 
}

\usage{
\method{predict}{spPGOcc}(object, X.0, coords.0, n.omp.threads = 1, verbose = TRUE, 
                          n.report = 100, ...)
}

\arguments{

  \item{object}{an object of class \code{spPGOcc}}

  \item{X.0}{the design matrix for prediction locations. This should 
    include a column of 1s for the intercept if an intercept is included in the original model fit.
    Covariates should have the same column names as those used when fitting the model with 
    \code{spPGOcc}.}

  \item{coords.0}{the spatial coordinates corresponding to \code{X.0}.}

   \item{n.omp.threads}{a positive integer indicating
    the number of threads to use for SMP parallel processing. The package must
    be compiled for OpenMP support. For most Intel-based machines, we recommend setting
    \code{n.omp.threads} up to the number of hyperthreaded cores. 
    Note, \code{n.omp.threads} > 1 might not work on some systems.}
 
  \item{verbose}{if \code{TRUE}, model specification and progress of the 
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.}

  \item{n.report}{the interval to report sampling progress.}

  \item{...}{currently no additional arguments}
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\references{
  Hooten, M. B., and Hefley, T. J. (2019). Bringing Bayesian models to life. CRC Press.
}

\value{
  An object of class \code{predict.spPGOcc} that is a list comprised of: 

  \item{psi.0.samples}{a \code{coda} object of posterior predictive samples 
    for the latent occurrence probability values.}
  
  \item{z.0.samples}{a \code{coda} object of posterior predictive samples 
    for the latent occurrence values.}
  
  \item{w.0.samples}{a \code{coda} object of posterior predictive samples 
    for the latent spatial random effects.}

  \item{run.time}{execution time reported using \code{proc.time()}.}

  The return object will include additional objects used for standard 
  extractor functions. 
}

\examples{
set.seed(400)
# Simulate Data -----------------------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.5, 2)
p.occ <- length(beta)
alpha <- c(0, 1)
p.det <- length(alpha)
phi <- 3 / .6
sigma.sq <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha, 
              sigma.sq = sigma.sq, phi = phi, sp = TRUE, cov.model = 'exponential')
# Split into fitting and prediction data set
pred.indx <- sample(1:J, round(J * .5), replace = FALSE)
y <- dat$y[-pred.indx, ]
# Occupancy covariates
X <- dat$X[-pred.indx, ]
# Prediction covariates
X.0 <- dat$X[pred.indx, ]
# Detection covariates
X.p <- dat$X.p[-pred.indx, , ]
coords <- as.matrix(dat$coords[-pred.indx, ])
coords.0 <- as.matrix(dat$coords[pred.indx, ])
psi.0 <- dat$psi[pred.indx]
w.0 <- dat$w[pred.indx]

# Package all data into a list
occ.covs <- X[, -1, drop = FALSE]
colnames(occ.covs) <- c('occ.cov')
det.covs <- list(det.cov.1 = X.p[, , 2])
data.list <- list(y = y, 
		  occ.covs = occ.covs, 
		  det.covs = det.covs, 
		  coords = coords)

# Number of batches
n.batch <- 50
# Batch length
batch.length <- 25
n.iter <- n.batch * batch.length
# Priors 
prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
		   alpha.normal = list(mean = 0, var = 2.72),
		   sigma.sq.ig = c(2, 2), 
		   phi.unif = c(3/1, 3/.1)) 
# Initial values
inits.list <- list(alpha = 0, beta = 0,
		      phi = 3 / .5, 
		      sigma.sq = 2,
		      w = rep(0, nrow(X)),
		      z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1) 

out <- spPGOcc(occ.formula = ~ occ.cov, 
	       det.formula = ~ det.cov.1, 
	       data = data.list, 
	       inits = inits.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length, 
	       accept.rate = 0.43, 
	       priors = prior.list,
	       cov.model = 'exponential', 
	       tuning = tuning.list, 
	       n.omp.threads = 1, 
	       verbose = TRUE, 
	       NNGP = TRUE, 
	       n.neighbors = 15, 
	       search.type = 'cb', 
	       n.report = 10, 
	       n.burn = 500, 
	       n.thin = 1)

summary(out) 

# Predict at new locations ------------------------------------------------
out.pred <- predict(out, X.0, coords.0, verbose = FALSE)
}
