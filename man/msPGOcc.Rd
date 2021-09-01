\name{msPGOcc}
\alias{msPGOcc}
\title{Function for Fitting Multi-Species Occupancy Models Using Polya-Gamma Latent Variables}


\description{
Function for fitting multi-species occupancy models using Polya-Gamma latent variables. 
}

\usage{
msPGOcc(occ.formula, det.formula, data, starting, n.samples, priors, 
	n.omp.threads = 1, verbose = TRUE, n.report = 100, ...)
}

\arguments{

\item{occ.formula}{a symbolic description of the model to be fit
  for the occupancy portion of the model using R's model syntax. Only
  right-hand side of formula is specified. See example below.}

\item{det.formula}{a symbolic description of the model to be fit
  for the detection portion of the model using R's model syntax. Only
  right-hand side of formula is specified. See example below.}

\item{data}{a list containing data necessary for model fitting.
  Valid tags are \code{y}, \code{occ.covs}, and \code{det.covs}. \code{y}
  is a three-dimensional array with first dimension equal to the number 
  of species, second dimension equal to the number of sites, and third 
  dimension equal to the maximum number of replicates at a given site. 
  \code{occ.covs} is a matrix or data frame containing the variables 
  used in the occupancy portion of the model, with \eqn{J}{J} rows for 
  each column (variable). \code{det.covs} is a list of variables included 
  in the detection portion of the model. Each list element is a different 
  detection covariate, which can be site-level or observational-level. 
  Site-level covariates are specified as a vector of length \eqn{J}{J} while 
  observation-level covariates are specified as a matrix or data frame 
  with the number of rows equal to \eqn{J}{J} and number of columns equal 
  to the maximum number of replicates at a given site.}

\item{starting}{a list with each tag corresponding to a parameter name.
  Valid tags are \code{alpha.comm}, \code{beta.comm}, \code{beta}, 
  \code{alpha}, \code{tau.beta}, \code{tau.alpha}, \code{z}. 
  The value portion of each tag is the parameter's starting value.}

\item{n.samples}{the number of posterior samples to collect.}

\item{priors}{a list with each tag corresponding to a parameter name. 
  Valid tags are \code{beta.comm.normal}, \code{alpha.comm.normal}, 
  \code{tau.beta.ig}, and \code{tau.alpha.ig}. Community-level occupancy 
  (\code{beta.comm}) and detection (\code{alpha.comm}) regression coefficients 
  are assumed to follow a normal distribution. The hyperparameters of the 
  normal distribution are passed as a list of length two with the first
  and second elements corresponding to the mean and variance of the normal
  distribution, which are each specified as vectors of 
  length equal to the number of coefficients to be estimated. If not
  specified, prior means are set to 0 and prior variances set to 2.73. 
  Community-level variance parameters for occupancy (\code{tau.beta}) and
  detection (\code{tau.alpha}) are assumed to follow an inverse Gamma 
  distribution. The hyperparameters of the inverse gamma distribution are 
  passed as a list of length two with the first and second elements 
  corresponding to the shape and scale parameters, which are each specified
  as vectors of length equal to the number of coefficients to be estimated. 
  If not specified, prior shape and scale parameters are set to 0.1.}

\item{n.omp.threads}{a positive integer indicating the number of threads
  to use for SMP parallel processing. The package must be compiled for 
  OpenMP support. For most Intel-based machines, we recommend setting 
  \code{n.omp.threads} up to the number of hypterthreaded cores. Note, 
  \code{n.omp.threads} > 1 might not work on some systems.}

\item{verbose}{if \code{TRUE}, model specification and progress of the 
  sampler is printed to the screen. Otherwise, only messages in data 
  preparation phase are printed.}

\item{n.report}{the interval to report MCMC progress.}

\item{...}{currently no additional arguments}

}

\note{
 Some of the underlying code used for generating random numbers from the
 Polya-Gamma distribution is taken from the \pkg{pgdraw} package written
 by Daniel F. Schmidt and Enes Makalic. Their code implements Algorithm
 6 in PhD thesis of Jesse Bennett Windle (2013) \url{https://repositories.lib.utexas.edu/handle/2152/21842}.
}

\references{
  Polson, N.G., J.G. Scott, and J. Windle. (2013) Bayesian Inference for
  Logistic Models Using Polya-Gamma Latent Variables.
  \emph{Journal of the American Statistical Association}, 108:1339-1349.
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{msPGOcc} that is a list comprised of: 

  \item{beta.comm.samples}{a \code{coda} object of posterior samples
    for the community level occupancy regression coefficients.}
  
  \item{alpha.comm.samples}{a \code{coda} object of posterior samples
    for the community level detection regression coefficients.}

  \item{tau.beta.samples}{a \code{coda} object of posterior samples
    for the occupancy community variance parameters.}

  \item{tau.alpha.samples}{a \code{coda} object of posterior samples
    for the detection community variance parameters.}

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the species level occupancy regression coefficients.}

  \item{alpha.samples}{a \code{coda} object of posterior samples
    for the species level detection regression coefficients.}

  \item{z.samples}{a three-dimensional array of posterior samples for
    the latent occupancy values for each species.}

  \item{psi.samples}{a three-dimensional array of posterior samples for
    the latent occupancy probability values for each species.}

  \item{y.rep.samples}{a four-dimensional array of fitted values for use
    in Goodness of Fit assessments.}

  \item{run.time}{MCMC sampler execution time reported using \code{proc.time()}.}

  The return object will include additional objects used for 
  subsequent prediction and/or model fit evaluation. 
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
# These could potentially be treated as random effects, but maybe not necessary.
phi <- runif(N, 3/1, 3/.1)
sigma.sq <- runif(N, 0.3, 3)
sp = FALSE

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
		sp = FALSE)
n.samples <- 5000
y <- dat$y
X <- dat$X
X.p <- dat$X.p
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
prior.list <- list(beta.comm.normal = list(mean = rep(0, p.occ), 
					   var = rep(2.72, p.occ)), 
		   alpha.comm.normal = list(mean = rep(0, p.det), 
					    var = rep(2.72, p.det)), 
		   tau.beta.ig = list(a = rep(0.1, p.occ), 
					  b = rep(0.1, p.occ)), 
		   tau.alpha.ig = list(a = rep(0.1, p.det), 
					   b = rep(0.1, p.det)))
# Starting values
starting.list <- list(alpha.comm = rep(0, p.det), 
		      beta.comm = rep(0, p.occ), 
		      beta = matrix(0, N, p.occ), 
		      alpha = matrix(0, N, p.det),
		      tau.beta = rep(1, p.occ), 
		      tau.alpha = rep(1, p.det), 
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
}
}
