\name{spPGOcc}
\alias{spPGOcc}
\title{Function for Fitting Single Species Spatial Occupancy Models Using Polya-Gamma Latent Variables}

\description{
  The function \code{spPGOcc} fits single species spatial occupancy models using Polya-Gamma latent variables. Models can be fit using either a full Gaussian process or a Nearest Neighbor Gaussian Process for large data sets. 
}

\usage{
spPGOcc(occ.formula, det.formula, data, starting, n.batch, 
	batch.length, accept.rate = 0.43, priors, 
	cov.model = "exponential", tuning, n.omp.threads = 1, 
	verbose = TRUE, NNGP = FALSE, n.neighbors = 15, 
	search.type = "cb", n.report = 100, ...)
}

\arguments{

  \item{occ.formula}{a symbolic description of the model to be fit
    for the occupancy portion of the model using R's model syntax. Only
    right-hand side of formula is specified. See example below.}
  
  \item{det.formula}{a symbolic description of the model to be fit
  for the detection portion of the model using R's model syntax. Only
  right-hand side of formula is specified. See example below.}

  \item{data}{a list containing data necessary for model fitting.
    Valid tags are \code{y}, \code{occ.covs}, \code{det.covs}, and \code{coords}. 
    \code{y} is the detection-nondetection data matrix or data frame with 
    first dimension equal to the number of sites (\eqn{J}{J}) and second 
    dimension equal to the maximum number of replicates at a given site. 
    \code{occ.covs} is a matrix or data frame containing the variables used 
    in the occupancy portion of the model, with \eqn{J}{J} rows for each column 
    (variable). \code{det.covs} is a list of variables included in the 
    detection portion of the model. Each list element is a different detection 
    covariate, which can be site-level or observational-level. Site-level 
    covariates are specified as a vector of length \eqn{J}{J} while 
    observation-level covariates are specified 
    as a matrix or data frame with the number of rows equal to \eqn{J}{J} and 
    number of columns equal to the maximum number of replicates at a given site. 
    \code{coords} is a \eqn{J \times 2}{J x 2} matrix of the observation coordinates.}
  
  \item{starting}{a list with each tag corresponding to a parameter name.
    Valid tags are \code{z}, \code{beta}, \code{alpha}, \code{sigma.sq}, 
    \code{phi}, \code{w}, and \code{nu}. \code{nu} is only specified if 
    \code{cov.model = "matern"}. The value portion of each tag is the 
    parameter's starting value.}
 
  \item{n.batch}{the number of MCMC batches to run for the Adaptive MCMC 
    sampler. See Roberts and Rosenthal (2009) for details.}
  
  \item{batch.length}{the length of each MCMC batch to run for the Adaptive 
    MCMC sampler. See Roberts and Rosenthal (2009) for details.}
  
  \item{accept.rate}{target acceptance rate for Adaptive MCMC. Defaul is 
    0.43. See Roberts and Rosenthal (2009) for details.}
  
  \item{priors}{a list with each tag corresponding to a parameter name. 
    Valid tags are \code{beta.normal}, \code{alpha.normal}, \code{phi.unif}, 
    \code{sigma.sq.ig}, and \code{nu.unif}. Occupancy 
    (\code{beta}) and detection (\code{alpha}) regression coefficients 
    are assumed to follow a normal distribution. The hyperparameters of the 
    normal distribution are passed as a list of length two with the first
    and second elements corresponding to the mean and variance of the normal
    distribution, which are each specified as vectors of 
    length equal to the number of coefficients to be estimated. If not
    specified, prior means are set to 0 and prior variances set to 2.73. The 
    spatial variance parameter, \code{sigma.sq}, is assumed to follow an 
    inverse-Gamma distribution, whereas the spatial decay \code{phi} and 
    smoothness \code{nu} parameters are assumed to follow Uniform 
    distributions. The hyperparameters of the inverse-Gamma are passed as 
    a vector of length two, with the first and second elements corresponding
    to the \emph{shape} and \emph{scale}, respectively. The hyperparameters
    of the Uniform are also passed as a vector of length two with the first 
    and second elements corresponding to the lower and upper support, 
    respectively.}
 
  \item{cov.model}{a quoted keyword that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.}

  \item{tuning}{a list with each tag corresponding to a parameter
    name. Valid tags are \code{phi} and \code{nu}. The value portion of each
    tag defines the initial variance of the Adaptive sampler. See
    Roberts and Rosenthal (2009) for details.}
  
  \item{n.omp.threads}{a positive integer indicating
   the number of threads to use for SMP parallel processing. The package must
   be compiled for OpenMP support. For most Intel-based machines, we
   recommend setting \code{n.omp.threads} up to the number of
   hyperthreaded cores. Note, \code{n.omp.threads} > 1 might not
   work on some systems.}
 
  \item{verbose}{if \code{TRUE}, model specification and progress of the sampler 
    is printed on the screen. Otherwise, only messages in data preparation
    phase are printed.}
  
  \item{NNGP}{if \code{TRUE}, model is fit with an NNGP. If \code{FALSE}, 
    a full Gaussian process is used. See Datta et al. (2016) and 
    Finley et al. (2019) for more information.}
  
  \item{n.neighbors}{number of neighbors used in the NNGP. Only used if 
  \code{NNGP = TRUE}.}
  
  \item{search.type}{a quoted keyword that specifies the type of nearest
    neighbor search algorithm. Supported method key words are: \code{"cb"} and
    \code{"brute"}. The \code{"cb"} should generally be much
    faster. If locations do not have identical coordinate values on
    the axis used for the nearest neighbor ordering then \code{"cb"} 
    and \code{"brute"} should produce identical neighbor sets. 
    However, if there are identical coordinate values on the axis used 
    for nearest neighbor ordering, then \code{"cb"} and \code{"brute"} 
    might produce different, but equally valid, neighbor sets, 
    e.g., if data are on a grid. }
 
  \item{n.report}{the interval to report Metropolis sampler acceptance and
    MCMC progress.}
  
  \item{...}{currently no additional arguments}
}

\references{
  Datta, A., S. Banerjee, A.O. Finley, and A.E. Gelfand. (2016)
  Hierarchical Nearest-Neighbor Gaussian process models for large
  geostatistical datasets. \emph{Journal of the American Statistical
    Association}, \doi{10.1080/01621459.2015.1044091}.
  
  Finley, A.O., A. Datta, B.D. Cook, D.C. Morton, H.E. Andersen, and
  S. Banerjee. (2019) Efficient algorithms for Bayesian Nearest Neighbor
  Gaussian Processes. \emph{Journal of Computational and Graphical
    Statistics}, \doi{10.1080/10618600.2018.1537924}.

  Polson, N.G., J.G. Scott, and J. Windle. (2013) Bayesian Inference for
  Logistic Models Using Polya-Gamma Latent Variables. 
  \emph{Journal of the American Statistical Association}, 108:1339-1349.

  Roberts, G.O. and Rosenthal J.S. (2009) Examples  of adaptive MCMC. 
  \emph{Journal of Computational and Graphical Statistics}, 18(2):349-367.
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{spPGOcc} that is a list comprised of: 

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

  \item{theta.samples}{a \code{coda} object of posterior samples
    for covariance parameters.}

  \item{w.samples}{a \code{coda} object of posterior samples
    for latent spatial random effects.}

  \item{run.time}{MCMC sampler execution time reported using \code{proc.time()}. 
    This time does not include nearest neighbor search time for building
    the neighbor set when \code{NNGP = TRUE}.}

  The return object will include additional objects used for 
  subsequent prediction and/or model fit evaluation. 
}

\examples{
\dontrun{
set.seed(400)
# Simulate Data -----------------------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.5, -0.15)
p.occ <- length(beta)
alpha <- c(0.7, 0.4, -0.2)
p.det <- length(alpha)
phi <- 3 / .6
sigma.sq <- 2
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha, 
              sigma.sq = sigma.sq, phi = phi, sp = TRUE)
y <- dat$y
X <- dat$X
X.p <- dat$X.p
coords <- as.matrix(dat$coords)

# Package all data into a list
occ.covs <- X[, -1, drop = FALSE]
colnames(occ.covs) <- c('occ.cov')
det.covs <- list(det.cov.1 = X.p[, , 2], 
		 det.cov.2 = X.p[, , 3])
data.list <- list(y = y, 
		  occ.covs = occ.covs, 
		  det.covs = det.covs, 
		  coords = coords)

# Number of batches
n.batch <- 200
# Batch length
batch.length <- 25
n.iter <- n.batch * batch.length
# Priors 
prior.list <- list(beta.normal = list(mean = rep(0, p.occ), 
				      var = rep(2.72, p.occ)), 
		   alpha.normal = list(mean = rep(0, p.det), 
			               var = rep(2.72, p.det)), 
		   sigma.sq.ig = c(2, 2), 
		   phi.unif = c(3/1, 3/.1), 
		   nu.unif = c(0.5, 10))
# Starting values
starting.list <- list(alpha = rep(0, p.det), 
		      beta = rep(0, p.occ), 
		      phi = 3 / .5, 
		      sigma.sq = 2,
		      nu = 1, 
		      w = rep(0, nrow(X)),
		      z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1, 
		    nu = 1)


out <- spPGOcc(occ.formula = ~ occ.cov, 
	       det.formula = ~ det.cov.1 + det.cov.2, 
	       data = data.list, 
	       starting = starting.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length, 
	       priors = prior.list,
	       cov.model = "exponential", 
	       tuning = tuning.list, 
	       NNGP = TRUE, 
	       n.neighbors = 15, 
	       search.type = 'cb', 
	       n.report = 50)

summary(out)

# Plot of spatial GP. 
w.summary <- summary(out$w.samples)$quantiles[, c(3, 1, 5)]
w.true <- dat$w

sq.cex=2.1
par(cex.lab=1.5, cex.main=1.5, xaxt="n", yaxt="n", bty="n", mar=c(5,5,4,2))
layout(matrix(1:2,1,2))
plot(coords, type="n", xlab="", ylab="", asp=TRUE, main="TRUE")
points(coords,pch=15,cex=sq.cex,col=rgb(0,0,0,(w.true-min(w.true))/
					diff(range(w.true))))
plot(coords,type="n",xlab="",ylab="",asp=TRUE,main="FIT")
points(coords,pch=15,cex=sq.cex,
       col=rgb(0,0,0,(w.summary[, 1]-min(w.summary[, 1]))/
	       diff(range(w.summary[, 1]))))
}
}
