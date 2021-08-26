\name{ppcOcc}
\alias{ppcOcc}
\title{Function for performing posterior predictive checks}

\usage{
ppcOcc(object, fit.stat, sub.sample, group, ...)
}

\description{
  Function for performing posterior predictive checks on \code{spOccupancy} model objects.
}

\arguments{
\item{object}{an object of class \code{PGOcc}, \code{spPGOcc}, \code{msPGOcc}, 
\code{spMsPGOcc}, \code{intPGOcc}, or \code{spIntPGOcc}.}

\item{fit.stat}{a quoted keyword that specifies the fit statistic
to use in the posterior predictive check. Supported fit statistics are
\code{"freeman-tukey"} and \code{"chi-square"}.}

\item{sub.sample}{an optional list that specifies the samples to include in
the posterior predictive check. Valid tags are \code{start},
\code{end}, and \code{thin}. Given the value associated with the tags,
the sample subset is selected using \code{seq(as.integer(start),
as.integer(end), by=as.integer(thin))}. If not specified, the entire posterior
sample is used.}

\item{group}{a positive integer indicating the way to group the
detection-nondetection data for the posterior predictive check. Value 
1 will group values by row (site) and value 2 will group values by 
column (replicate).} 

\item{...}{currently no additional arguments}
}

\details{
  Standard GoF assessments are not valid for binary data, and posterior
  predictive checks must be performed on some sort of binned data. 
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{ppcOcc} that is a list comprised of: 

  \item{fit.y}{a \code{coda} object of posterior samples for the 
  fit statistic calculated on the observed data.}

  \item{fit.y.rep}{a \code{coda} object of posterior samples for the 
  fit statistic calculated on a replicate data set generated from the 
  model.}

  \item{fit.y.group.quants}{a matrix consisting of posterior quantiles
  for the fit statistic using the observed data for each unique element
  the fit statistic is calculated for (i.e., sites when group = 1, 
  replicates when group = 2).}

  \item{fit.y.rep.group.quants}{a matrix consisting of posterior quantiles
  for the fit statistic using the model replicated data for each unique element
  the fit statistic is calculated for (i.e., sites when group = 1, 
  replicates when group = 2).}

  The return object will include additional objects used for 
  standard extractor functions. 

}

\examples{
set.seed(400)
# Simulate Data -----------------------------------------------------------
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.5, -0.15)
p.occ <- length(beta)
alpha <- c(0.7, 0.4)
p.det <- length(alpha)
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      sp = FALSE)
occ.covs <- dat$X[, 2, drop = FALSE]
colnames(occ.covs) <- c('occ.cov')
det.covs <- list(det.cov = dat$X.p[, , 2])
# Data bundle
data.list <- list(y = dat$y, 
		  occ.covs = occ.covs, 
		  det.covs = det.covs)

# Priors
prior.list <- list(beta.normal = list(mean = rep(0, p.occ),
				      var = rep(2.72, p.occ)),
		   alpha.normal = list(mean = rep(0, p.det),
			               var = rep(2.72, p.det)))
# Starting values
starting.list <- list(alpha = rep(0, p.det),
		      beta = rep(0, p.occ),
		      z = apply(data.list$y, 1, max, na.rm = TRUE))

n.samples <- 5000
n.report <- 1000

out <- PGOcc(occ.formula = ~ occ.cov, 
	     det.formula = ~ det.cov, 
	     data = data.list, 
	     starting = starting.list,
	     n.samples = n.samples,
	     priors = prior.list,
	     n.omp.threads = 1,
	     verbose = TRUE,
	     n.report = n.report)

# Posterior predictive check
ppc.out <- ppcOcc(out, 'chi-square', list(start = 4001, end = 5000, thin = 1), 1)
summary(ppc.out)
}
