\name{predict.intPGOcc}
\alias{predict.intPGOcc}
\title{Function for prediction at new locations for single species integrated occupancy models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new locations given an object of class `intPGOcc`. 
}

\usage{
\method{predict}{intPGOcc}(object, X.0, sub.sample, ...)
}

\arguments{

  \item{object}{an object of class intPGOcc}

  \item{X.0}{the design matrix for prediction locations. This should include a column of 1s for the intercept.}

  \item{sub.sample}{an optional list that specifies the samples to include in the composition sampling algorithm. Valid tags are \code{start}, \code{end}, and \code{thin}. Given the values associated with tags, the sample subset is selected using \code{seq(as.integer(start), as.integer(end), by = as.integer(thin))}}

  \item{...}{currently no additional arguments}
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{predict.intPGOcc} that is a list comprised of: 

  \item{psi.0.samples}{a \code{coda} object of posterior predictive samples for the
    latent occupancy probability values.}

  \item{z.0.samples}{a \code{coda} object of posterior predictive samples for the
    latent occupancy values.}

  The return object will include additional objects used for standard 
  extractor functions. 
}

\examples{
\dontrun{
set.seed(1008)

# Simulate Data -----------------------------------------------------------
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
beta <- c(0.5, 1)
p.occ <- length(beta)
# Detection covariates
alpha <- list()
for (i in 1:n.data) {
  alpha[[i]] <- runif(sample(1:4, 1), -1, 1)
}
p.det.long <- sapply(alpha, length)
p.det <- sum(p.det.long)

# Simulate occupancy data. 
dat <- simIntOcc(n.data = n.data, J.x = J.x, J.y = J.y, J.obs = J.obs, 
		 n.rep = n.rep, beta = beta, alpha = alpha, sp = FALSE)

y <- dat$y
X <- dat$X.obs
X.p <- dat$X.p
sites <- dat$sites

# Package all data into a list
occ.covs <- X[, 2, drop = FALSE]
colnames(occ.covs) <- c('occ.cov')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(det.cov.1.1 = X.p[[1]][, , 2], 
		      det.cov.1.2 = X.p[[1]][, , 3], 
		      det.cov.1.3 = X.p[[1]][, , 4])
det.covs[[2]] <- list(det.cov.2.1 = X.p[[2]][, , 2])
det.covs[[3]] <- list(det.cov.3.1 = X.p[[3]][, , 2], 
		      det.cov.3.2 = X.p[[3]][, , 3], 
		      det.cov.3.3 = X.p[[3]][, , 4])
det.covs[[4]] <- list(det.cov.4.1 = X.p[[4]][, , 2], 
		      det.cov.4.2 = X.p[[4]][, , 3], 
		      det.cov.4.3 = X.p[[4]][, , 4])
data.list <- list(y = y, 
		  occ.covs = occ.covs,
		  det.covs = det.covs, 
		  sites = sites)

J <- length(dat$z.obs)
# Starting values
alpha.start <- list()
for (q in 1:n.data) {
  alpha.start[[q]] <- rep(0, p.det.long[q])
}
starting.list <- list(alpha = alpha.start, 
		      beta = rep(0, p.occ), 
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
			               var = alpha.var))

n.samples <- 5000
out <- intPGOcc(occ.formula = ~ occ.cov, 
                det.formula = list(f.1 = ~ det.cov.1.1 + det.cov.1.2 + det.cov.1.3, 
		              f.2 = ~ det.cov.2.1, 
		              f.3 = ~ det.cov.3.1 + det.cov.3.2 + det.cov.3.3, 
		              f.4 = ~ det.cov.4.1 + det.cov.4.2 + det.cov.4.3), 
		data = data.list,
		starting = starting.list,
		n.samples = n.samples, 
		priors = prior.list, 
		n.omp.threads = 1, 
		verbose = TRUE, 
		n.report = 1000)

summary(out)

# Prediction
X.0 <- dat$X.pred
psi.0 <- dat$psi.pred
sub.sample <- list(start = 4001, end = 5000, thin = 1)

out.pred <- predict(out, X.0, sub.sample)
psi.hat.quants <- apply(out.pred$psi.0.samples, 2, quantile, c(0.025, 0.5, 0.975))
plot(psi.0, psi.hat.quants[2, ], pch = 19, xlab = 'True', 
     ylab = 'Fitted', ylim = c(min(psi.hat.quants), max(psi.hat.quants)))
segments(psi.0, psi.hat.quants[1, ], psi.0, psi.hat.quants[3, ])
lines(psi.0, psi.0)
}
}
