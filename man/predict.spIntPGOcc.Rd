\name{predict.spIntPGOcc}
\alias{predict.spIntPGOcc}
\title{Function for prediction at new locations for single species integrated spatial occupancy models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new locations given an object of class `spIntPGOcc`. 
}

\usage{
\method{predict}{spIntPGOcc}(object, X.0, coords.0, sub.sample, 
	                     n.omp.threads = 1, verbose = TRUE, 
			     n.report = 100, ...)
}

\arguments{

  \item{object}{an object of class spIntPGOcc}

  \item{X.0}{the design matrix for prediction locations. This should include a column of 1s for the intercept.}

  \item{coords.0}{the spatial coordinates corresponding to \code{X.0}.}

  \item{sub.sample}{an optional list that specifies the samples to include in the composition sampling algorithm. Valid tags are \code{start}, \code{end}, and \code{thin}. Given the values associated with tags, the sample subset is selected using \code{seq(as.integer(start), as.integer(end), by = as.integer(thin))}}

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
    latent occupancy probability values.}

  \item{z.0.samples}{a \code{coda} object of posterior predictive samples for the
    latent occupancy values.}

  The return object will include additional objects used for standard 
  extractor functions. 
}

\examples{
\dontrun{
set.seed(400)

# Simulate Data -----------------------------------------------------------
# Number of locations in each direction. This is the total region of interest
# where some sites may or may not have a data source. 
J.x <- 15
J.y <- 15
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
		 n.rep = n.rep, beta = beta, alpha = alpha, sp = sp)

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

# Starting values
alpha.start <- list()
for (q in 1:n.data) {
  alpha.start[[q]] <- rep(0, p.det.long[q])
}
starting.list <- list(alpha = alpha.start, 
		      beta = rep(0, p.occ), 
		      phi = 3 / .5, 
		      sigma.sq = 2, 
		      w = rep(0, J), 
		      z = rep(1, J))
# Priors
alpha.mu <- list()
alpha.var <- list()
for (q in 1:n.data) {
  alpha.mu[[q]] <- rep(0, p.det.long[q])
  alpha.var[[q]] <- rep(2.72, p.det.long[q])
}
prior.list <- list(beta.normal = list(mean = rep(0, p.occ), 
				      var = rep(2.72, p.occ)), 
		   alpha.normal = list(mean = alpha.mu, 
			               var = alpha.var),
		   phi.unif = c(3/1, 3/.1), 
		   sigma.sq.ig = c(2, 2))
# Tuning
tuning.list <- list(phi = 1, 
		    nu = 1)

# Number of batches
n.batch <- 200
# Batch length
batch.length <- 25

out <- spIntPGOcc(occ.formula = ~ occ.cov, 
		  det.formula = list(f.1 = ~ det.cov.1.1, 
				     f.2 = ~ det.cov.2.1 + det.cov.2.2, 
				     f.3 = ~ det.cov.3.1, 
				     f.4 = ~ det.cov.4.1 + det.cov.4.2 + det.cov.4.3), 
		  data = data.list,  
		  starting = starting.list, 
		  n.batch = n.batch, 
		  batch.length = batch.length, 
		  accept.rate = 0.43, 
		  priors = prior.list, 
		  cov.model = "exponential", 
		  tuning = tuning.list, 
		  n.omp.threads = 1, 
		  verbose = TRUE, 
		  NNGP = TRUE, 
		  n.neighbors = 15, 
		  search.type = 'cb', 
		  n.report = 50)

summary(out)

# Predict at new locations ------------------------------------------------
sub.sample <- list(start = 4001, end = 5000, thin = 1)

out.pred <- predict(out, X.0, coords.0, sub.sample)
psi.hat.quants <- apply(out.pred$psi.0.samples, 2, quantile, c(0.025, 0.5, 0.975))

sq.cex=2.1
par(cex.lab=1.5,cex.main=1.5,xaxt="n",yaxt="n",bty="n",mar=c(5,5,4,2))
layout(matrix(1:2,1,2))
plot(coords.0,type="n",xlab="",ylab="",asp=TRUE,main="Simulated")
points(coords.0,pch=15,cex=sq.cex,col=rgb(0,0,0,(psi.0-min(psi.0))/diff(range(psi.0))))
plot(coords.0,type="n",xlab="",ylab="",asp=TRUE,main="Spatial")
points(coords.0,pch=15,cex=sq.cex,col=rgb(0,0,0,(psi.hat.quants[2,]-min(psi.hat.quants[2,]))/diff(range(psi.hat.quants[2,]))))
}
}
