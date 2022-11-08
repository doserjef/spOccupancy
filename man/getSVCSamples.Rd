\name{getSVCSamples}
\alias{getSVCSamples}
\title{Extract spatially-varying coefficient MCMC samples}

\usage{
getSVCSamples(object, pred.object, ...)
}

\description{
  Function for extracting the full spatially-varying coefficient MCMC samples
  from an spOccupancy model object. 
}

\arguments{
  \item{object}{an object of class \code{svcPGOcc}, \code{svcPGBinom}, 
    \code{svcTPGOcc}, \code{svcTPGBinom}.}

  \item{pred.object}{a prediction object from a spatially-varying coefficient 
    model fit using spOccupancy. Should be of class \code{predict.svcPGOcc}, 
    \code{predict.svcPGBinom}, \code{predict.svcTPGOcc}, or \code{predict.svcTPGBinom}. 
    If specified, SVC samples are extracted at the prediction locations.}

  \item{...}{currently no additional arguments}

}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
}

\value{
  A list of \code{coda::mcmc} objects of the spatially-varying coefficient MCMC samples
  for all spatially-varying coefficients estimated in the model (including the
  intercept if specified). Note these values correspond to the sum of the estimated
  spatial and non-spatial effect to give the overall effect of the covariate at
  each location. Each element of the list is a two-dimensional matrix 
  where dimensions correspond to MCMC sample and site. If \code{pred.object} is specified, 
  values are returned for the prediction locations instead of the sampled locations. 
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
phi <- c(3 / .6, 3 / .8)
sigma.sq <- c(1.2, 0.7)
svc.cols <- c(1, 2)
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha, 
              sigma.sq = sigma.sq, phi = phi, sp = TRUE, cov.model = 'exponential', 
              svc.cols = svc.cols)
# Detection-nondetection data
y <- dat$y
# Occupancy covariates
X <- dat$X
# Detection covarites
X.p <- dat$X.p
# Spatial coordinates
coords <- dat$coords

# Package all data into a list
occ.covs <- X[, -1, drop = FALSE]
colnames(occ.covs) <- c('occ.cov')
det.covs <- list(det.cov.1 = X.p[, , 2])
data.list <- list(y = y, 
                  occ.covs = occ.covs, 
                  det.covs = det.covs, 
                  coords = coords)

# Number of batches
n.batch <- 10
# Batch length
batch.length <- 25
n.iter <- n.batch * batch.length
# Priors 
prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                   alpha.normal = list(mean = 0, var = 2.72),
                   sigma.sq.ig = list(a = 2, b = 1), 
                   phi.unif = list(a = 3/1, b = 3/.1)) 
# Initial values
inits.list <- list(alpha = 0, beta = 0,
                   phi = 3 / .5, 
                   sigma.sq = 2,
                   w = matrix(0, nrow = length(svc.cols), ncol = nrow(X)),
                   z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1) 

out <- svcPGOcc(occ.formula = ~ occ.cov, 
                det.formula = ~ det.cov.1, 
                data = data.list, 
                inits = inits.list, 
                n.batch = n.batch, 
                batch.length = batch.length, 
                accept.rate = 0.43, 
                priors = prior.list,
                cov.model = 'exponential', 
                svc.cols = c(1, 2),
                tuning = tuning.list, 
                n.omp.threads = 1, 
                verbose = TRUE, 
                NNGP = TRUE, 
                n.neighbors = 5, 
                search.type = 'cb', 
                n.report = 10, 
                n.burn = 50, 
                n.thin = 1)

svc.samples <- getSVCSamples(out)
str(svc.samples)
}

