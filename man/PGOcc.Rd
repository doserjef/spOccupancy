\name{PGOcc}
\alias{PGOcc}
\title{Function for Fitting Single Species Occupancy Models Using Polya-Gamma Latent Variables}

\usage{
PGOcc(y, X, X.p, starting, n.rep, n.samples, priors,
      n.omp.threads = 1, verbose = TRUE, 
      n.report = 100, ...
)
}

\description{
  Function for fitting single species occupancy models using Polya-Gamma latent variables. 
}

\arguments{
\item{y}{a matrix or data frame containing the detection-nondetection 
  data with number of rows equal to the number of sites (\eqn{J}{J}) 
  and number of columns equal to the maximum number of replicates 
  at a given site (\code{max(n.rep)})}

\item{X}{a matrix or data frame containing the design matrix for the 
  occupancy portion of the model with number of rows equal to the number 
  of sites (\eqn{J}{J}) and number of columns equal to the number of 
  occupancy regression coefficients.}

\item{X.p}{a three-dimensional array containing the design matrix for the 
  detection portion of the model with first dimension equal to the number 
  of sites (\eqn{J}{J}), second dimension equal to the maximum number of 
  replicates at a given site (\code{max(n.rep)}), and third dimension equal 
  to the number of detection regression coefficients.} 

\item{starting}{a list with each tag corresponding to a parameter name.
  Valid tags are \code{z}, \code{beta}, and \code{alpha}. The value portion
  of each tag is the parameter's starting value.}

\item{n.rep}{a vector of length \eqn{J}{J} consisting of the number of 
  replicate detection-nondetection surveys at each data point.}

\item{n.samples}{the number of posterior samples to collect.}

\item{priors}{a list with each tag corresponding to a parameter name. 
  Valid tags are \code{beta.normal} and \code{alpha.normal}. Occupancy 
  (\code{beta}) and detection (\code{alpha}) regression coefficients 
  are assumed to follow a normal distribution. The hyperparameters of the 
  normal distribution are passed as a list of length two with the first
  and second elements corresponding to the mean and variance of the normal
  distribution, which are each specified as vectors of 
  length equal to the number of coefficients to be estimated. If not
  specified, prior means are set to 0 and prior variances set to 2.73.}

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
  An object of class \code{PGOcc} that is a list comprised of: 

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the occupancy regression coefficients.}

  \item{alpha.samples}{a \code{coda} object of posterior samples
    for the detection regression coefficients.}

  \item{z.samples}{a \code{coda} object of posterior samples 
    for the latent occupancy values}

  \item{psi.samples}{a \code{coda} object of posterior samples
    for the latent occupancy probability values}

  \item{y.rep.samples}{a three-dimensional array of fitted values
    for use in Goodness of Fit assessments.}

  \item{run.time}{MCMC sampler execution time reported using \code{proc.time()}.}

  The return object will include additional objects used for 
  subsequent prediction and/or model fit evaluation. 

}

\examples{
set.seed(400)
# Simulate Data -----------------------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.5, -0.15)
p.occ <- length(beta)
alpha <- c(0.7, 0.4)
p.det <- length(alpha)
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      sp = FALSE)
# Detection-nondetection data
y <- dat$y
# Occupancy covariates
X <- dat$X.psi
# Detection covariates
X.p <- dat$X.p

# Priors
prior.list <- list(beta.normal = list(mean = rep(0, p.occ),
				      var = rep(2.72, p.occ)),
		   alpha.normal = list(mean = rep(0, p.det),
			               var = rep(2.72, p.det)))
# Starting values
starting.list <- list(alpha = rep(0, p.det),
		      beta = rep(0, p.occ),
		      z = apply(y, 1, max, na.rm = TRUE))

n.samples <- 5000
n.report <- 1000

out <- PGOcc(y = y,
	     X = as.data.frame(X),
	     X.p = X.p,
	     starting = starting.list,
	     n.rep = n.rep,
	     n.samples = n.samples,
	     priors = prior.list,
	     n.omp.threads = 1,
	     verbose = TRUE,
	     n.report = n.report)
summary(out)
}
