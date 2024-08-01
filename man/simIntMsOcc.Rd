\name{simIntMsOcc}
\alias{simIntMsOcc}
\title{Simulate Multi-Species Detection-Nondetection Data from Multiple Data Sources}

\description{
  The function \code{simIntMsOcc} simulates multi-species detection-nondetection data from multiple data sources for simulation studies, power assessments, or function testing of integrated occupancy models. Data can optionally be simulated with a spatial Gaussian Process on the occurrence process.
}

\usage{
simIntMsOcc(n.data, J.x, J.y, J.obs, n.rep, n.rep.max, N, beta, alpha, psi.RE = list(),
            p.RE = list(), sp = FALSE, svc.cols = 1, cov.model, sigma.sq, phi, nu,
            factor.model = FALSE, n.factors, range.probs, ...)
}
\arguments{

\item{n.data}{an integer indicating the number of detection-nondetection data sources to simulate.}

\item{J.x}{a single numeric value indicating the number of sites across the region of interest along the horizontal axis. Total number of sites across the simulated region of interest is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{J.y}{a single numeric value indicating the number of sites across the region of interest along the vertical axis. Total number of sites across the simulated region of interest is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{J.obs}{a numeric vector of length \code{n.data} containing the number of sites to simulate each data source at. Data sources can be obtained at completely different sites, the same sites, or anywhere inbetween. Maximum number of sites a given data source is available at is equal to \eqn{J = J.x \times J.y}{J = J.x x J.y}.}

\item{n.rep}{a list of length \code{n.data}. Each element is a numeric vector with length corresponding to the number of sites that given data source is observed at (in \code{J.obs}). Each vector indicates the number of repeat visits at each of the sites for a given data source.}

\item{n.rep.max}{a vector of numeric values indicating the maximum number of replicate surveys for each data set. This is an optional argument, with its default value set to \code{max(n.rep)} for each data set. This can be used to generate data sets with different types of missingness (e.g., simulate data across 20 days (replicate surveys) but sites are only sampled a maximum of ten times each).}

\item{N}{a numeric vector of length \code{N} containing the number of species each data source samples. These can be the same if both data sets sample the same species, or can be different.}

\item{beta}{a numeric matrix with \code{max(N)} rows containing the intercept and regression coefficient parameters for the occurrence portion of the multi-species occupancy model. Each row corresponds to the regression coefficients for a given species.} 

\item{alpha}{a list of length \code{n.data}. Each element is a numeric matrix with the rows corresponding to the number of species that data source contains and columns corresponding to the regression coefficients for each data source.}

\item{psi.RE}{a list used to specify the non-spatial random intercepts included in the occurrence portion of the model. The list must have two tags: \code{levels} and \code{sigma.sq.psi}. \code{levels} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the number of levels there are in each intercept. \code{sigma.sq.psi} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the variances for each random effect. If not specified, no random effects are included in the occurrence portion of the model.}

\item{p.RE}{this argument is not currently supported. In a later version, this argument will allow for simulating data with detection random effects in the different data sources.}

\item{sp}{a logical value indicating whether to simulate a spatially-explicit occupancy model with a Gaussian process. By default set to \code{FALSE}.}

\item{svc.cols}{a vector indicating the variables whose effects will be
  estimated as spatially-varying coefficients. \code{svc.cols} is an
  integer vector with values indicating the order of covariates specified
  in the model formula (with 1 being the intercept if specified).}

\item{cov.model}{a quoted keyword that specifies the covariance function used to model the spatial dependence structure among the latent occurrence values.  Supported covariance model key words are: \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and \code{"gaussian"}.}

\item{sigma.sq}{a numeric vector of length \code{max(N)} containing the spatial variance parameter for each species. Ignored when \code{sp = FALSE} or when \code{factor.model = TRUE}.}

\item{phi}{a numeric vector of length \code{max(N)} containing the spatial decay parameter for each species. Ignored when \code{sp = FALSE}. If \code{factor.model = TRUE}, this should be of length \code{n.factors}.}

\item{nu}{a numeric vector of length \code{max(N)} containing the spatial smoothness parameter for each species. Only used when \code{sp = TRUE} and \code{cov.model = 'matern'}. If \code{factor.model = TRUE}, this should be of length \code{n.factors}.} 

\item{factor.model}{a logical value indicating whether to simulate data following a factor modeling approach that explicitly incoporates species correlations. If \code{sp = TRUE}, the latent factors are simulated from independent spatial processes. If \code{sp = FALSE}, the latent factors are simulated from standard normal distributions.}

\item{n.factors}{a single numeric value specifying the number of latent factors to use to simulate the data if \code{factor.model = TRUE}.}

\item{range.probs}{a numeric vector of length \code{N} where each value should fall between 0 and 1, and indicates the probability that one of the \code{J} spatial locations simulated is within the simulated range of the given species. If set to 1, every species has the potential of being present at each location.}

\item{...}{currently no additional arguments}
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
}

\references{
 Doser, J. W., Leuenberger, W., Sillett, T. S., Hallworth, M. T. & Zipkin, E. F. (2022). 
 Integrated community occupancy models: A framework to assess occurrence and biodiversity 
 dynamics using multiple data sources. Methods in Ecology and Evolution, 00, 1-14. \doi{10.1111/2041-210X.13811}
}


\value{
  A list comprised of: 

  \item{X.obs}{a numeric design matrix for the occurrence portion of the model. This matrix contains the intercept and regression coefficients for only the observed sites.}
  \item{X.pred}{a numeric design matrix for the occurrence portion of the model at sites where there are no observed data sources.} 
  \item{X.p}{a list of design matrices for the detection portions of the integrated multi-species occupancy model. Each element in the list is a design matrix of detection covariates for each data source.}
  \item{coords.obs}{a numeric matrix of coordinates of each observed site. Required for spatial models.}
  \item{coords.pred}{a numeric matrix of coordinates of each site in the study region without any data sources. Only used for spatial models.}
  \item{w}{a species (or factor) x site matrix of the spatial random effects for each species. Only used to simulate data when \code{sp = TRUE}. If \code{factor.model = TRUE}, the first dimension is \code{n.factors}.}
  \item{w.pred}{a matrix of the spatial random random effects for each species (or factor) at locations without any observation.}
  \item{psi.obs}{a species x site matrix of the occurrence probabilities for each species at the observed sites. Note that values are provided for all species, even if some species are only monitored at a subset of these points.}
  \item{psi.pred}{a species x site matrix of the occurrence probabilities for sites without any observations.}
  \item{z.obs}{a species x site matrix of the latent occurrence states at each observed site. Note that values are provided for all species, even if some species are only monitored at a subset of these points.}
  \item{z.pred}{a species x site matrix of the latent occurrence states at each site without any observations.}
  \item{p}{a list of detection probability arrays for each of the \code{n.data} data sources. Each array has dimensions corresponding to species, site, and replicate, respectively.}
  \item{y}{a list of arrays of the raw detection-nondetection data for each site and replicate combination for each species in the data set. Each array has dimensions corresponding to species, site, and replicate, respectively.}

}

\examples{
set.seed(91)
J.x <- 10
J.y <- 10
# Total number of data sources across the study region
J.all <- J.x * J.y
# Number of data sources.
n.data <- 2
# Sites for each data source.
J.obs <- sample(ceiling(0.2 * J.all):ceiling(0.5 * J.all), n.data, replace = TRUE)
n.rep <- list()
n.rep[[1]] <- rep(3, J.obs[1])
n.rep[[2]] <- rep(4, J.obs[2])

# Number of species observed in each data source
N <- c(8, 3)

# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 0.5)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.4, 0.3)
# Detection
# Detection covariates
alpha.mean <- list()
tau.sq.alpha <- list()
# Number of detection parameters in each data source
p.det.long <- c(4, 3)
for (i in 1:n.data) {
  alpha.mean[[i]] <- runif(p.det.long[i], -1, 1)
  tau.sq.alpha[[i]] <- runif(p.det.long[i], 0.1, 1)
}
# Random effects
psi.RE <- list()
p.RE <- list()
beta <- matrix(NA, nrow = max(N), ncol = p.occ)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(max(N), beta.mean[i], sqrt(tau.sq.beta[i]))
}
alpha <- list()
for (i in 1:n.data) {
  alpha[[i]] <- matrix(NA, nrow = N[i], ncol = p.det.long[i])
  for (t in 1:p.det.long[i]) {
    alpha[[i]][, t] <- rnorm(N[i], alpha.mean[[i]][t], sqrt(tau.sq.alpha[[i]])[t])
  }
}
sp <- FALSE
factor.model <- FALSE
# Simulate occupancy data
dat <- simIntMsOcc(n.data = n.data, J.x = J.x, J.y = J.y,
		   J.obs = J.obs, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
	           psi.RE = psi.RE, p.RE = p.RE, sp = sp, factor.model = factor.model,
                   n.factors = n.factors)
str(dat)
}
