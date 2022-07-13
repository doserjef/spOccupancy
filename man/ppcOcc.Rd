\name{ppcOcc}
\alias{ppcOcc}
\title{Function for performing posterior predictive checks}

\usage{
ppcOcc(object, fit.stat, group, ...)
}

\description{
  Function for performing posterior predictive checks on \code{spOccupancy} model objects.
}

\arguments{
\item{object}{an object of class \code{PGOcc}, \code{spPGOcc}, \code{msPGOcc}, 
\code{spMsPGOcc}, \code{intPGOcc}, \code{spIntPGOcc}, \code{lfMsPGOcc}, \code{sfMsPGOcc}, 
\code{tPGOcc}, or \code{stPGOcc}.}

\item{fit.stat}{a quoted keyword that specifies the fit statistic
to use in the posterior predictive check. Supported fit statistics are
\code{"freeman-tukey"} and \code{"chi-squared"}.}

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

  \item{fit.y}{a numeric vector of posterior samples for the 
  fit statistic calculated on the observed data when \code{object} is 
  of class \code{PGOcc} or \code{spPGOcc}. When \code{object} is of 
  class \code{msPGOcc}, \code{spMsPGOcc}, \code{lfMsPGOcc}, or \code{sfMsPGOcc}, 
  this is a numeric matrix with rows corresponding to posterior samples
  and columns corresponding to species. When \code{object} is of class
  \code{intPGOcc} or \code{spIntPGOcc}, this is a list, with each element 
  of the list being a vector of posterior samples for each data set. 
  When \code{object} is of class \code{tPGOcc} or \code{stPGOcc}, this is a
  numeric matrix with rows corresponding to posterior samples and columns 
  corresponding to primary sampling periods.}

  \item{fit.y.rep}{a numeric vector of posterior samples for the 
  fit statistic calculated on a replicate data set generated from the 
  model when \code{object} is of class \code{PGOcc} or \code{spPGOcc}. When \code{object} is of 
  class \code{msPGOcc}, \code{spMsPGOcc}, \code{lfMsPGOcc}, or \code{sfMsPGOcc}, 
  this is a numeric matrix with rows corresponding to posterior samples
  and columns corresponding to species. When \code{object} is of class
  \code{intPGOcc} or \code{spIntPGOcc}, this is a list, with each element 
  of the list being a vector of posterior samples for each data set. 
  When \code{object} is of class \code{tPGOcc} or \code{stPGOcc}, this is a
  numeric matrix with rows corresponding to posterior samples and columns 
  corresponding to primary sampling periods.}

  \item{fit.y.group.quants}{a matrix consisting of posterior quantiles
  for the fit statistic using the observed data for each unique element
  the fit statistic is calculated for (i.e., sites when group = 1, 
  replicates when group = 2) when \code{object} is of class \code{PGOcc} or \code{spPGOcc}. 
  When \code{object} is of class \code{msPGOcc}, 
  \code{spMsPGOcc}, \code{lfMsPGOcc}, or \code{sfMsPGOcc}, this is a three-dimensional
  array with the additional dimension corresponding to species. When \code{object} is 
  of class \code{intPGOcc} or \code{spIntPGOcc}, this is a list, with each element 
  consisting of the posterior quantile matrix for each data set. When \code{object} 
  is of class \code{tPGOcc} or \code{stPGOcc}, this is a three-dimensional array 
  with the additional dimension corresponding to primary sampling periods.}

  \item{fit.y.rep.group.quants}{a matrix consisting of posterior quantiles
  for the fit statistic using the model replicated data for each unique element
  the fit statistic is calculated for (i.e., sites when group = 1, 
  replicates when group = 2) when \code{object} is of class \code{PGOcc} or \code{spPGOcc}. 
  When \code{object} is of class \code{msPGOcc}, 
  \code{spMsPGOcc}, \code{lfMsPGOcc}, or \code{sfMsPGOcc}, this is a three-dimensional
  array with the additional dimension corresponding to species. When \code{object} is 
  of class \code{intPGOcc} or \code{spIntPGOcc}, this is a list, with each element 
  consisting of the posterior quantile matrix for each data set. When \code{object} 
  is of class \code{tPGOcc} or \code{stPGOcc}, this is a three-dimensional array 
  with the additional dimension corresponding to primary sampling periods.}

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
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   alpha.normal = list(mean = 0, var = 2.72))
# Initial values
inits.list <- list(alpha = 0, beta = 0,
                   z = apply(data.list$y, 1, max, na.rm = TRUE))

n.samples <- 5000
n.report <- 1000

out <- PGOcc(occ.formula = ~ occ.cov, 
             det.formula = ~ det.cov, 
             data = data.list, 
             inits = inits.list,
             n.samples = n.samples,
             priors = prior.list,
             n.omp.threads = 1,
             verbose = TRUE,
             n.report = n.report, 
             n.burn = 4000, 
             n.thin = 1)

# Posterior predictive check
ppc.out <- ppcOcc(out, fit.stat = 'chi-squared', group = 1)
summary(ppc.out)
}
