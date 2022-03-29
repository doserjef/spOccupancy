\name{predict.spIntPGOcc}
\alias{predict.spIntPGOcc}
\title{Function for prediction at new locations for single-species integrated spatial occupancy models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new locations given an object of class `spIntPGOcc`. 
}

\usage{
\method{predict}{spIntPGOcc}(object, X.0, coords.0, n.omp.threads = 1, verbose = TRUE, 
        n.report = 100, ...)
}

\arguments{

  \item{object}{an object of class \code{spIntPGOcc}.}

  \item{X.0}{the design matrix for prediction locations. This should include a column of 1s for the intercept. Covariates should have the same column names as those used when fitting the model with \code{spIntPGOcc}.}

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

  \item{...}{currently no additional arguments}
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{predict.spIntPGOcc} that is a list comprised of: 

  \item{psi.0.samples}{a \code{coda} object of posterior predictive samples for the
    latent occurrence probability values.}

  \item{z.0.samples}{a \code{coda} object of posterior predictive samples for the
    latent occurrence values.}

  The return object will include additional objects used for standard 
  extractor functions. 
}

\references{
  Hooten, M. B., and Hefley, T. J. (2019). Bringing Bayesian models to life. CRC Press.
}

\examples{
set.seed(400)

# Simulate Data -----------------------------------------------------------
# Number of locations in each direction. This is the total region of interest
# where some sites may or may not have a data source. 
J.x <- 8
J.y <- 8
J.all <- J.x * J.y
# Number of data sources.
n.data <- 4
# Sites for each data source. 
J.obs <- sample(ceiling(0.2 * J.all):ceiling(0.5 * J.all), n.data, replace = TRUE)
# Replicates for each data source.
n.rep <- list()
for (i in 1:n.data) {
  n.rep[[i]] <- sample(1:4, size = J.obs[i], replace = TRUE)
}
# Occupancy covariates
beta <- c(0.5, 0.5)
p.occ <- length(beta)
# Detection covariates
alpha <- list()
alpha[[1]] <- runif(2, 0, 1)
alpha[[2]] <- runif(3, 0, 1)
alpha[[3]] <- runif(2, -1, 1)
alpha[[4]] <- runif(4, -1, 1)
p.det.long <- sapply(alpha, length)
p.det <- sum(p.det.long)
sigma.sq <- 2
phi <- 3 / .5
sp <- TRUE

# Simulate occupancy data. 
dat <- simIntOcc(n.data = n.data, J.x = J.x, J.y = J.y, J.obs = J.obs, 
                 n.rep = n.rep, beta = beta, alpha = alpha, sp = sp, 
                 phi = phi, sigma.sq = sigma.sq, cov.model = 'spherical')

y <- dat$y
X <- dat$X.obs
X.p <- dat$X.p
sites <- dat$sites
X.0 <- dat$X.pred
psi.0 <- dat$psi.pred
coords <- as.matrix(dat$coords.obs)
coords.0 <- as.matrix(dat$coords.pred)

# Package all data into a list
occ.covs <- X[, 2, drop = FALSE]
colnames(occ.covs) <- c('occ.cov')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(det.cov.1.1 = X.p[[1]][, , 2])
det.covs[[2]] <- list(det.cov.2.1 = X.p[[2]][, , 2], 
                      det.cov.2.2 = X.p[[2]][, , 3])
det.covs[[3]] <- list(det.cov.3.1 = X.p[[3]][, , 2])
det.covs[[4]] <- list(det.cov.4.1 = X.p[[4]][, , 2], 
                      det.cov.4.2 = X.p[[4]][, , 3], 
                      det.cov.4.3 = X.p[[4]][, , 4])
data.list <- list(y = y, 
                  occ.covs = occ.covs,
                  det.covs = det.covs, 
                  sites = sites, 
                  coords = coords)

J <- length(dat$z.obs)

# Initial values
inits.list <- list(alpha = list(0, 0, 0, 0), 
                   beta = 0, 
                   phi = 3 / .5, 
                   sigma.sq = 2, 
                   w = rep(0, J), 
                   z = rep(1, J))
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                   alpha.normal = list(mean = list(0, 0, 0, 0), 
                                       var = list(2.72, 2.72, 2.72, 2.72)),
                   phi.unif = c(3/1, 3/.1), 
                   sigma.sq.ig = c(2, 2))
# Tuning
tuning.list <- list(phi = 1) 

# Number of batches
n.batch <- 40
# Batch length
batch.length <- 25

out <- spIntPGOcc(occ.formula = ~ occ.cov, 
                  det.formula = list(f.1 = ~ det.cov.1.1, 
                                     f.2 = ~ det.cov.2.1 + det.cov.2.2, 
                                     f.3 = ~ det.cov.3.1, 
                                     f.4 = ~ det.cov.4.1 + det.cov.4.2 + det.cov.4.3), 
                  data = data.list,  
                  inits = inits.list, 
                  n.batch = n.batch, 
                  batch.length = batch.length, 
                  accept.rate = 0.43, 
                  priors = prior.list, 
                  cov.model = "spherical", 
                  tuning = tuning.list, 
                  n.omp.threads = 1, 
                  verbose = TRUE, 
                  NNGP = TRUE, 
                  n.neighbors = 5, 
                  search.type = 'cb', 
                  n.report = 10, 
                  n.burn = 500, 
                  n.thin = 1)
summary(out)

# Predict at new locations ------------------------------------------------
out.pred <- predict(out, X.0, coords.0, verbose = FALSE)
}
