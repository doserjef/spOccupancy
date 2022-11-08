\name{waicOcc}
\alias{waicOcc}
\title{Compute Widely Applicable Information Criterion for spOccupancy Model Objects}

\usage{
waicOcc(object, ...)
}

\description{
  Function for computing the Widely Applicable Information Criterion 
  (WAIC; Watanabe 2010) for \code{spOccupancy} model objects. 
}

\arguments{
  \item{object}{an object of class \code{PGOcc}, \code{spPGOcc}, \code{msPGOcc}, 
  \code{spMsPGOcc}, \code{intPGOcc}, \code{spIntPGOcc}, \code{lfJSDM}, 
  \code{sfJSDM}, \code{lfMsPGOcc}, \code{sfMsPGOcc}, \code{tPGOcc}, \code{stPGOcc}, 
  \code{svcPGBinom}, \code{svcPGOcc}, \code{svcTPGBinom}, or \code{svcTPGOcc}.}

\item{...}{currently no additional arguments}

}

\references{
  Watanabe, S. (2010). Asymptotic equivalence of Bayes cross
  validation and widely applicable information criterion in 
  singular learning theory. \emph{Journal of Machine Learning Research}, 
  11:3571-3594.

  Gelman, A., J. B. Carlin, H. S. Stern, D. B. Dunson, A. Vehtari, 
  and D. B. Rubin. (2013). Bayesian Data Analysis. 3rd edition. 
  CRC Press, Taylor and Francis Group

  Gelman, A., J. Hwang, and A. Vehtari (2014). Understanding predictive
  information criteria for Bayesian models. \emph{Statistics and Computing}, 
  24:997-1016.
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  When \code{object} is of class \code{PGOcc}, \code{spPGOcc}, \code{msPGOcc}, 
  \code{spMsPGOcc}, \code{lfJSDM}, \code{sfJSDM}, \code{lfMsPGOcc}, \code{sfMsPGOcc},
  \code{tPGOcc}, \code{stPGOcc}, \code{svcPGBinom}, or \code{svcPGOcc}, 
  returns a vector with three elements corresponding to 
  estimates of the expected log pointwise predictive density (elpd), the 
  effective number of parameters (pD), and the WAIC. When \code{object} is 
  of class \code{intPGOcc} or \code{spIntPGOcc}, returns a data frame with 
  columns elpd, pD, and WAIC, with each row corresponding to the estimated 
  values for each data source in the integrated model.
}

\details{
  The effective number of parameters is calculated following the recommendations
  of Gelman et al. (2014).
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
# Initial values
inits.list <- list(alpha = rep(0, p.det),
                   beta = rep(0, p.occ),
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

# Calculate WAIC
waicOcc(out)
}

