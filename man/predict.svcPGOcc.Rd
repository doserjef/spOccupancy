\name{predict.svcPGOcc}
\alias{predict.svcPGOcc}
\title{Function for prediction at new locations for single-species spatially-varying coefficient occupancy models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new 
  locations given an object of class `svcPGOcc`. Prediction is possible for both the 
  latent occupancy state as well as detection.
}

\usage{

\method{predict}{svcPGOcc}(object, X.0, coords.0, weights.0, n.omp.threads = 1, verbose = TRUE, 
        n.report = 100, ignore.RE = FALSE, type = 'occupancy', grid.index.0, ...)
}

\arguments{

  \item{object}{an object of class \code{svcPGOcc}}

  \item{X.0}{the design matrix of covariates at the prediction locations. This should include a column of 1s for the intercept if an intercept is included in the model. If random effects are included in the occupancy (or detection if \code{type = 'detection'}) portion of the model, the levels of the random effects at the new locations should be included as a column in the design matrix. The ordering of the levels should match the ordering used to fit the data in \code{svcPGOcc}. Columns should correspond to the order of how covariates were specified in the corresponding formula argument of \code{svcPGOcc}. Column names of the random effects must match the name of the random effects, if specified in the corresponding formula argument of \code{svcPGOcc}.}

  \item{coords.0}{the spatial coordinates corresponding to \code{X.0}. Note that \code{spOccupancy} assumes coordinates are specified 
  in a projected coordinate system.}

  \item{weights.0}{not used for objects of class \code{svcTPGOcc}. Used when 
    calling other functions.}

   \item{n.omp.threads}{a positive integer indicating
    the number of threads to use for SMP parallel processing. The package must
    be compiled for OpenMP support. For most Intel-based machines, we recommend setting
    \code{n.omp.threads} up to the number of hyperthreaded cores. 
    Note, \code{n.omp.threads} > 1 might not work on some systems.}
 
  \item{verbose}{if \code{TRUE}, model specification and progress of the 
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.}

  \item{ignore.RE}{a logical value indicating whether to include unstructured random 
    effects for prediction. If TRUE, unstructured random effects will be ignored and 
    prediction will only use the fixed effects and the spatial random effects. If FALSE, 
    random effects will be included in the prediction for both observed and unobserved 
    levels of the unstructured random effects.}

  \item{n.report}{the interval to report sampling progress.}

  \item{type}{a quoted keyword indicating what type of prediction to produce. Valid keywords are 'occupancy' to predict latent occupancy probability and latent occupancy values (this is the default), or 'detection' to predict detection probability given new values of detection covariates.}

  \item{grid.index.0}{an indexing vector used to specify how each row in \code{X.0} corresponds to the coordinates specified in \code{coords.0}. Only relevant if the SVCs were estimated at a higher spatial resolution (e.g., grid cells) than point locations.}

  \item{...}{currently no additional arguments}
}

\note{
  When \code{ignore.RE = FALSE}, both sampled levels and non-sampled levels of random effects are supported for prediction. For sampled levels, the posterior distribution for the random intercept corresponding to that level of the random effect will be used in the prediction. For non-sampled levels, random values are drawn from a normal distribution using the posterior samples of the random effect variance, which results in fully propagated uncertainty in predictions with models that incorporate random effects.
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\references{
  Hooten, M. B., and Hefley, T. J. (2019). Bringing Bayesian models to life. CRC Press.
}

\value{
  A list object of class \code{predict.svcPGOcc}. When \code{type = 'occupancy'}, the list consists of:  

  \item{psi.0.samples}{a \code{coda} object of posterior predictive samples 
    for the latent occurrence probability values.}
  
  \item{z.0.samples}{a \code{coda} object of posterior predictive samples 
    for the latent occurrence values.}
  
  \item{w.0.samples}{a three-dimensional array of posterior predictive samples
    for the spatial random effects, with dimensions corresponding to MCMC iteration, 
    coefficient, and site.}

  \item{run.time}{execution time reported using \code{proc.time()}.}

  When \code{type = 'detection'}, the list consists of: 

  \item{p.0.samples}{a \code{coda} object of posterior predictive samples for the 
    detection probability values.}

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
phi <- c(3 / .6, 3 / .8)
sigma.sq <- c(0.5, 0.9)
svc.cols <- c(1, 2)
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha, 
              sigma.sq = sigma.sq, phi = phi, sp = TRUE, cov.model = 'exponential', 
              svc.cols = svc.cols)
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
w.0 <- dat$w[pred.indx, , drop = FALSE]

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
                   sigma.sq.ig = list(a = 2, b = 0.5), 
                   phi.unif = list(a = 3/1, b = 3/.1)) 
# Initial values
inits.list <- list(alpha = 0, beta = 0,
                   phi = 3 / .5, 
                   sigma.sq = 0.5,
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
                tuning = tuning.list, 
                n.omp.threads = 1, 
                verbose = TRUE, 
                NNGP = TRUE, 
                svc.cols = c(1, 2),
                n.neighbors = 15, 
                search.type = 'cb', 
                n.report = 10, 
                n.burn = 50, 
                n.thin = 1)

summary(out) 

# Predict at new locations ------------------------------------------------
out.pred <- predict(out, X.0, coords.0, verbose = FALSE)
}
