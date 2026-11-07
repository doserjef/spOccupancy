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
\code{tPGOcc}, \code{stPGOcc}, \code{svcPGOcc}, \code{svcMsPGOcc}, \code{tMsPGOcc}, 
\code{stMsPGOcc}, \code{svcTMsPGOcc}, \code{tIntPGOcc}, \code{stIntPGOcc}, 
\code{svcTIntPGOcc}, \code{fpPGOcc}, \code{spFpPGOcc}, \code{svcFpPGocc}.}

\item{fit.stat}{a quoted keyword that specifies the fit statistic
to use in the posterior predictive check. Supported fit statistics are
\code{"freeman-tukey"} and \code{"chi-squared"}.}

\item{group}{a positive integer indicating the way to group the
detection-nondetection data for the posterior predictive check. Value 
1 will group values by row (site) and value 2 will group values by 
column (replicate). For false-positive model objects only, value 0 performs
the check on the individual site/replicate counts without grouping.} 

\item{...}{currently no additional arguments}
}

\details{
  Standard GoF assessments are not valid for binary data, and posterior
  predictive checks must be performed on some sort of binned data. This is
  what is done for all occupancy models that input binary data.

  For false positive model objects, two separate checks are performed: one for 
  the detection counts \eqn{y}, with expected value \eqn{z\mu + \omega}, and 
  one for the validation data \eqn{d}, with expected value 
  \eqn{n \cdot z\mu / (z\mu + \omega)}, where \eqn{n} is the observed number of 
  validated signals. Because the detection data are counts, they do not 
  need to be binned, and \code{group = 0} computes the fit statistic on each site/replicate count.
}

\author{
  Jeffrey W. Doser \email{jwdoser@ncsu.edu}, \cr
}

\value{

  An object of class \code{ppcOcc} that is a list comprised of: 

  \item{fit.y}{a numeric vector of posterior samples for the 
  fit statistic calculated on the observed data when \code{object} is 
  of class \code{PGOcc}, \code{spPGOcc}, or \code{svcPGOcc}. When \code{object} is of 
  class \code{msPGOcc}, \code{spMsPGOcc}, \code{lfMsPGOcc}, 
  \code{sfMsPGOcc}, or \code{svcMsPGOcc} 
  this is a numeric matrix with rows corresponding to posterior samples
  and columns corresponding to species. When \code{object} is of class
  \code{intPGOcc} or \code{spIntPGOcc}, this is a list, with each element 
  of the list being a vector of posterior samples for each data set. 
  When \code{object} is of class \code{tPGOcc} or \code{stPGOcc}, this is a
  numeric matrix with rows corresponding to posterior samples and columns 
  corresponding to primary sampling periods. When \code{object} is of class \code{tMsPGOcc}, 
  \code{stMsPGOcc}, or \code{svcTMsPGOcc}, this is a three-dimensional array with dimensions
  corresponding to MCMC sample, species, and primary time period. }

  \item{fit.y.rep}{a numeric vector of posterior samples for the 
  fit statistic calculated on a replicate data set generated from the 
  model when \code{object} is of class \code{PGOcc}, \code{spPGOcc}, or \code{svcPGOcc}. 
  When \code{object} is of class \code{msPGOcc}, \code{spMsPGOcc}, 
  \code{lfMsPGOcc}, \code{sfMsPGOcc}, or \code{svcMsPGOcc}
  this is a numeric matrix with rows corresponding to posterior samples
  and columns corresponding to species. When \code{object} is of class
  \code{intPGOcc} or \code{spIntPGOcc}, this is a list, with each element 
  of the list being a vector of posterior samples for each data set. 
  When \code{object} is of class \code{tPGOcc} or \code{stPGOcc}, this is a
  numeric matrix with rows corresponding to posterior samples and columns 
  corresponding to primary sampling periods. When \code{object} is of class \code{tMsPGOcc}, 
  \code{stMsPGOcc}, or \code{svcTMsPGOcc}, this is a three-dimensional array with dimensions
  corresponding to MCMC sample, species, and primary time period. }

  \item{fit.y.group.quants}{a matrix consisting of posterior quantiles
  for the fit statistic using the observed data for each unique element
  the fit statistic is calculated for (i.e., sites when group = 1, 
  replicates when group = 2) when \code{object} is of class \code{PGOcc}, \code{spPGOcc}, or \code{svcPGOcc}. 
  When \code{object} is of class \code{msPGOcc}, \code{spMsPGOcc}, 
  \code{lfMsPGOcc}, \code{sfMsPGOcc}, \code{svcMsPGOcc}, this is a three-dimensional
  array with the additional dimension corresponding to species. When \code{object} is 
  of class \code{intPGOcc} or \code{spIntPGOcc}, this is a list, with each element 
  consisting of the posterior quantile matrix for each data set. When \code{object} 
  is of class \code{tPGOcc} or \code{stPGOcc}, this is a three-dimensional array 
  with the additional dimension corresponding to primary sampling periods. When \code{object}
  is of class \code{tMsPGOcc}, \code{stMsPGOcc}, \code{svcTMsPGOcc}, this is a four-dimensional array
  with dimensions corresponding to quantile, species, grouping element, and primary time period.}

  \item{fit.y.rep.group.quants}{a matrix consisting of posterior quantiles
  for the fit statistic using the model replicated data for each unique element
  the fit statistic is calculated for (i.e., sites when group = 1, 
  replicates when group = 2) when \code{object} is of class \code{PGOcc}, \code{spPGOcc}, \code{svcPGOcc}. 
  When \code{object} is of class \code{msPGOcc}, 
  \code{spMsPGOcc}, \code{lfMsPGOcc}, \code{sfMsPGOcc}, or \code{svcMsPGOcc}, 
  this is a three-dimensional
  array with the additional dimension corresponding to species. When \code{object} is 
  of class \code{intPGOcc} or \code{spIntPGOcc}, this is a list, with each element 
  consisting of the posterior quantile matrix for each data set. When \code{object} 
  is of class \code{tPGOcc} or \code{stPGOcc}, this is a three-dimensional array 
  with the additional dimension corresponding to primary sampling periods. When \code{object}
  is of class \code{tMsPGOcc}, \code{stMsPGOcc}, \code{svcTMsPGOcc}, this is a four-dimensional array
  with dimensions corresponding to quantile, species, grouping element, and primary time period.}

  When \code{object} is a model fit from a false positive model, \code{fit.y}, \code{fit.y.rep}, 
  \code{fit.y.group.quants}, and \code{fit.y.rep.group.quants} correspond to 
  the detection count data \eqn{y}. When \code{group = 0}, the group quantiles
  are three-dimensional arrays with dimensions corresponding to quantile, site, 
  and replicate, with \code{NA} values for missing site/replicate combinations.
  The following additional elements are returned for the validation data \eqn{d}: 

  \item{fit.d}{a numeric vector of posterior samples for the fit statistic 
  calculated on the observed validation data.}

  \item{fit.d.rep}{a numeric vector of posterior samples for the fit statistic 
  calculated on replicate validation data generated from the model.}

  \item{fit.d.group.quants}{posterior quantiles of the fit statistic using
  the observed validation data for each unique element the fit statistic is 
  calculated for, in the same format as \code{fit.y.group.quants}.}

  \item{fit.d.rep.group.quants}{posterior quantiles of the fit statistic using
  the replicate validation data for each unique element the fit statistic is 
  calculated for, in the same format as \code{fit.y.group.quants}.}

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
