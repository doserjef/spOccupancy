\name{simTIntOcc}
\alias{simTIntOcc}
\title{Simulate Single-Species Multi-Season Detection-Nondetection Data from Multiple Data Sources}

\description{
  The function \code{simTIntOcc} simulates single-species detection-nondetection data from multiple data sources over multiple seasons for simulation studies, power assessments, or function testing of integrated multi-season occupancy models. Data can optionally be simulated with a spatial Gaussian Process on the occurrence process. Non-spatial random intercepts can be included in the detection or occurrence portions of the model.
}

\usage{
simTIntOcc(n.data, J.x, J.y, J.obs, n.time, data.seasons, n.rep, n.rep.max, 
           beta, alpha, sp.only = 0, trend = TRUE, psi.RE = list(), 
           p.RE = list(), sp = FALSE, svc.cols = 1, cov.model, 
           sigma.sq, phi, nu, ar1 = FALSE, rho, sigma.sq.t, 
           x.positive = FALSE, ...)
}
\arguments{

\item{n.data}{an integer indicating the number of detection-nondetection data sources to simulate.}

\item{J.x}{a single numeric value indicating the number of sites across the region of interest along the horizontal axis. Total number of sites across the simulated region of interest is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{J.y}{a single numeric value indicating the number of sites across the region of interest along the vertical axis. Total number of sites across the simulated region of interest is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{J.obs}{a numeric vector of length \code{n.data} containing the number of sites to simulate each data source at. Data sources can be obtained at completely different sites, the same sites, or anywhere inbetween. Maximum number of sites a given data source is available at is equal to \eqn{J = J.x \times J.y}{J = J.x x J.y}.}

\item{n.time}{a numeric vector of lencth \code{n.data} indicating the number of primary time periods (denoted T) over which sampling occurs for each site within each data source. Data sources can be simulated over differing numbers of primary time periods, and within a given data source sites can be sampled for a differing number of years.} 

\item{data.seasons}{a list of length \code{n.data} where each list element denotes the specific overall years that the given data source is simulated for. The length of vector should be equal to the maximum number of seasons any one given site in a given data source is sampled as specified in \code{n.time}.}

\item{n.rep}{a list of length \code{n.data}. Each element is a numeric matrix with rows equal to the number of sites for the given data set and columns equal number of primary time periods over which sampling occurs for the given data set. The value in cell indicates the number of repeat visits (secondary sampling events) for each site within a given primary time period.}

\item{n.rep.max}{a vector of numeric values indicating the maximum number of replicate surveys for each data set. This is an optional argument, with its default value set to \code{max(n.rep)} for each data set. This can be used to generate data sets with different types of missingness (e.g., simulate data across 20 days (replicate surveys) but sites are only sampled a maximum of ten times each).}

\item{beta}{a numeric vector containing the intercept and regression coefficient parameters for the occupancy portion of the model. Note that if \code{trend = TRUE}, the second value in the vector corresponds to the estimated occurrence trend.}

\item{alpha}{a list of length \code{n.data}. Each element is a numeric vector containing the intercept and regression coefficient parameters for the detection portion of the single-species occupancy model for each data source.}		

\item{sp.only}{a numeric vector specifying which occurrence covariates should only vary over space and not over time. The numbers in the vector correspond to the elements in the vector of regression coefficients (\code{beta}). By default, all simulated occurrence covariates are assumed to vary over both space and time.}

\item{trend}{a logical value. If \code{TRUE}, a temporal trend will be used to simulate the detection-nondetection data and the second element of \code{beta} is assumed to be the trend parameter. If \code{FALSE} no trend is used to simulate the data and all elements of \code{beta} (except the first value which is the intercept) correspond to covariate effects.}

\item{psi.RE}{a list used to specify the non-spatial random intercepts included in the occupancy portion of the model. The list must have two tags: \code{levels} and \code{sigma.sq.psi}. \code{levels} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the number of levels there are in each intercept. \code{sigma.sq.psi} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the variances for each random effect. If not specified, no random effects are included in the occupancy portion of the model.}

\item{p.RE}{a list used to specify the non-spatial random intercepts included in the detection portion of the model. The list must be a list of lists, where the individual lists contain the detection coefficients for each data set in the integrated  model. Each of the lists must have two tags: \code{levels} and \code{sigma.sq.p}. \code{levels} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the number of levels there are in each intercept. \code{sigma.sq.p} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the variances for each random effect. If not specified, no random effects are included in the detection portion of the model.}

\item{sp}{a logical value indicating whether to simulate a spatially-explicit occupancy model with a Gaussian process. By default set to \code{FALSE}.}

\item{svc.cols}{a vector indicating the variables whose effects will be
  estimated as spatially-varying coefficients. \code{svc.cols} is an
  integer vector with values indicating the order of covariates specified
  in the model formula (with 1 being the intercept if specified).}

\item{cov.model}{a quoted keyword that specifies the covariance function used to model the spatial dependence structure among the latent occurrence values.  Supported covariance model key words are: \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and \code{"gaussian"}.}

\item{sigma.sq}{a numeric value indicating the spatial variance parameter. Ignored when \code{sp = FALSE}. When \code{svc.cols} is specified with more than one SVC, \code{sigma.sq} must be of length \code{length(svc.cols)}.}

\item{phi}{a numeric value indicating the spatial range parameter. Ignored when \code{sp = FALSE}. When \code{svc.cols} is specified with more than one SVC, \code{phi} must be of length \code{length(svc.cols)}.}

\item{nu}{a numeric value indicating the spatial smoothness parameter. Only used when \code{sp = TRUE} and \code{cov.model = "matern"}. When \code{svc.cols} is specified with more than one SVC, \code{nu} must be of length \code{length(svc.cols)}.}

\item{ar1}{a logical value indicating whether to simulate a temporal random effect with an AR(1) process. By default, set to \code{FALSE}.}

\item{rho}{a numeric value indicating the AR(1) temporal correlation parameter. Ignored when \code{ar1 = FALSE}.}

\item{sigma.sq.t}{a numeric value indicating the AR(1) temporal variance parameter. Ignored when \code{ar1 = FALSE}.}

\item{x.positive}{a logical value indicating whether the simulated covariates should be simulated as random standard normal covariates (\code{x.positive = FALSE}) or restricted to positive values (\code{x.positive = TRUE}). If \code{x.positive = TRUE}, covariates are simulated from a random normal and then the minimum value is added to each covariate value to ensure non-negative covariate values.}

\item{...}{currently no additional arguments}

}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}
}

\value{
  A list comprised of: 

  \item{X.obs}{a three-dimensional numeric array with dimensions corresponding to sites, primary time periods, and occurrence covariate containing the design matrix for the occurrence portion of the occupancy model. This matrix contains the intercept and regression coefficients for only the observed sites.}
  \item{X.pred}{a three-dimensional numeric array with dimensions corresponding to sites, primary time periods, and occurrence covariate containing the design matrix for the occurrence portion of the occupancy model. This matrix contains the intercept and regression coefficients for the sites in the study region where there are no observed data sources.}
  \item{X.pred}{a numeric design matrix for the occurrence portion of the model at sites where there are no observed data sources.} 
  \item{X.p}{a list of design matrices for the detection portions of the integrated occupancy model. Each element in the list is a design matrix of detection covariates for each data source. Each design matrix is formatted as a four-dimensional array with dimensions corresponding to sites, primary time period, secondary time period, and covariate.}
  \item{coords.obs}{a numeric matrix of coordinates of each observed site. Required for spatial models.}
  \item{coords.pred}{a numeric matrix of coordinates of each site in the study region without any data sources. Only used for spatial models.}
  \item{w.obs}{a matrix of the spatial random effects at observed locations. Only used to simulate data when \code{sp = TRUE}}. 
  \item{w.pred}{a matrix of the spatial random random effects at locations without any observation.}
  \item{psi.obs}{a matrix of the occurrence probabilities for each observed site and primary time period.}
  \item{psi.pred}{a matrix of the occurrence probabilities for sites without any observations.}
  \item{z.obs}{a matrix of the latent occurrence states at each observed site and primary time period.}
  \item{z.pred}{a matrix of the latent occurrence states at each site without any observations.}
  \item{p}{a list of detection probability arrays for each of the \code{n.data} data sources. Each array has three dimensions corresponding to site, primary time period, and secondary time period.}
  \item{y}{a list of arrays of the raw detection-nondetection data for each site, primary time period, and replicate combination.}
  \item{X.p.re}{a list of four-dimensional numeric arrays containing the levels of any detection random effect included in the model for each data source. Only relevant when detection random effects are specified in \code{p.RE}. Dimensions of each array correspond to site, primary time period, secondary time period, and random effect.}
  \item{X.re.obs}{a numeric array containing the levels of any occurrence random effect included in the model at the sites where there is at least one data source. Dimensions correspond to site, primary time period, and parameter. Only relevant when occurrence random effects are specified in \code{psi.RE}.}
  \item{X.re.pred}{a numeric array containing the levels of any occurrence random effect included in the model at the sites where there are no data sources sampled. Dimensions correspond to site, primary time period, and parameter. Only relevant when occurrence random effects are specified in \code{psi.RE}.}
  \item{alpha.star}{a list of numeric vectors that contains the simulated detection random effects for each given level of the random effects included in the detection model for each data set. Only relevant when detection random effects are included in the model.}
  \item{beta.star}{a numeric vector that contains the simulated occurrence random effects for each given level of the random effects included in the occurrence model. Only relevant when occurrence random effects are included in the model.}
  \item{eta}{a \eqn{T \times 1}{T x 1} matrix of the latent AR(1) random effects. Only included when \code{ar1 = TRUE}.}

}

\examples{
# Number of locations in each direction. This is the total region of interest
# where some sites may or may not have a data source.
J.x <- 15 
J.y <- 15
J.all <- J.x * J.y
# Number of data sources.
n.data <- 3
# Sites for each data source.
J.obs <- sample(ceiling(0.2 * J.all):ceiling(0.4 * J.all), n.data, replace = TRUE)
# Maximum number of years for each data set
n.time.max <- c(4, 8, 10)
# Number of years each site in each data set is sampled
n.time <- list()
for (i in 1:n.data) {
  n.time[[i]] <- sample(1:n.time.max[i], J.obs[i], replace = TRUE)
}
# Replicates for each data source.
n.rep <- list()
for (i in 1:n.data) {
  n.rep[[i]] <- matrix(NA, J.obs[i], n.time.max[i])
  for (j in 1:J.obs[i]) {
    n.rep[[i]][j, sample(1:n.time.max[i], n.time[[i]][j], replace = FALSE)] <- 
      sample(1:4, n.time[[i]][j], replace = TRUE)
  }
}
# Total number of years across all data sets
n.time.total <- 10
# List denoting the specific years each data set was sampled during. 
data.seasons <- list()
for (i in 1:n.data) {
  data.seasons[[i]] <- sort(sample(1:n.time.total, n.time.max[i], replace = FALSE))
}

# Occupancy covariates
beta <- c(0, 0.4, 0.3)
# Random occupancy effects
psi.RE <- list(levels = c(20), sigma.sq.psi = c(0.6))
# Detection covariates
alpha <- list()
for (i in 1:n.data) {
  alpha[[i]] <- runif(3, 0, 1)
}
# Detection random effects
p.RE <- list()
p.RE[[1]] <- list(levels = c(35), sigma.sq.p = c(0.5))
p.RE[[2]] <- list(levels = c(20, 10), sigma.sq.p = c(0.7, 0.3))
p.RE[[3]] <- list(levels = c(20),  sigma.sq.p = c(0.6))
p.det.long <- sapply(alpha, length)
p.det <- sum(p.det.long)
# Spatial components
sigma.sq <- 2
phi <- 3 / .5
nu <- 1
sp <- TRUE
# Temporal parameters
ar1 <- TRUE 
rho <- 0.9
sigma.sq.t <- 1.5
svc.cols <- c(1)
n.rep.max <- sapply(n.rep, max, na.rm = TRUE)

# Simulate occupancy data.
dat <- simTIntOcc(n.data = n.data, J.x = J.x, J.y = J.y, J.obs = J.obs,
                  n.time = n.time, data.seasons = data.seasons, 
                  n.rep = n.rep, n.rep.max = n.rep.max, 
                  beta = beta, alpha = alpha, trend = TRUE, 
                  psi.RE = psi.RE, p.RE = p.RE, sp = sp, svc.cols = svc.cols, 
                  cov.model = 'exponential', sigma.sq = sigma.sq, phi = phi, 
                  nu = nu, ar1 = ar1, rho = rho, sigma.sq.t = sigma.sq.t)
}