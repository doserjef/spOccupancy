\name{simTOcc}
\alias{simTOcc}
\title{Simulate Multi-Season Single-Species Detection-Nondetection Data}

\description{
  The function \code{simTOcc} simulates multi-season single-species occurrence data for simulation studies, power assessments, or function testing. Data can be optionally simulated with a spatial Gaussian Process in the occurrence portion of the model. Non-spatial random intercepts can also be included in the detection or occurrence portions of the occupancy model.
}

\usage{
simTOcc(J.x, J.y, n.time, n.rep, beta, alpha, sp.only = 0, trend = TRUE, 
        psi.RE = list(), p.RE = list(), sp = FALSE, cov.model, 
        sigma.sq, phi, nu, ar1 = FALSE, rho, sigma.sq.t, ...)
}

\arguments{
\item{J.x}{a single numeric value indicating the number of sites to simulate detection-nondetection data along the horizontal axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{J.y}{a single numeric value indicating the number of sites to simulate detection-nondetection data along the vertical axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{n.time}{a single numeric value indicating the number of primary time periods (denoted T) over which sampling occurs.} 

\item{n.rep}{a numeric matrix indicating the number of replicates at each site during each primary time period. The matrix must have \eqn{J = J.x \times J.y}{J = J.x x J.y} rows and T columns, where T is the number of primary time periods (e.g., years or seasons) over which sampling occurs.} 

\item{beta}{a numeric vector containing the intercept and regression coefficient parameters for the occupancy portion of the single-species occupancy model. Note that if \code{trend = TRUE}, the second value in the vector corresponds to the estimated occurrence trend.}

\item{alpha}{a numeric vector containing the intercept and regression coefficient parameters for the detection portion of the single-species occupancy model.}

\item{sp.only}{a numeric vector specifying which occurrence covariates should only vary over space and not over time. The numbers in the vector correspond to the elements in the vector of regression coefficients (\code{beta}). By default, all simulated occurrence covariates are assumed to vary over both space and time.}

\item{trend}{a logical value. If \code{TRUE}, a temporal trend will be used to simulate the detection-nondetection data and the second element of \code{beta} is assumed to be the trend parameter. If \code{FALSE} no trend is used to simulate the data and all elements of \code{beta} (except the first value which is the intercept) correspond to covariate effects.}

\item{psi.RE}{a list used to specify the unstructured random intercepts included in the occupancy portion of the model. The list must have two tags: \code{levels} and \code{sigma.sq.psi}. \code{levels} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the number of levels there are in each intercept. \code{sigma.sq.psi} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the variances for each random effect. An additional tag \code{site.RE} can be set to \code{TRUE} to simulate data with a site-specific non-spatial random effect on occurrence. If not specified, no random effects are included in the occupancy portion of the model.}

\item{p.RE}{a list used to specify the unstructured random intercepts included in the detection portion of the model. The list must have two tags: \code{levels} and \code{sigma.sq.p}. \code{levels} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the number of levels there are in each intercept. \code{sigma.sq.p} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the variances for each random effect. If not specified, no random effects are included in the detection portion of the model.}

\item{sp}{a logical value indicating whether to simulate a spatially-explicit occupancy model with a Gaussian process. By default set to \code{FALSE}.}

\item{cov.model}{a quoted keyword that specifies the covariance function used to model the spatial dependence structure among the latent occurrence values.  Supported covariance model key words are: \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and \code{"gaussian"}.}

\item{sigma.sq}{a numeric value indicating the spatial variance parameter. Ignored when \code{sp = FALSE}.}

\item{phi}{a numeric value indicating the spatial decay parameter. Ignored when \code{sp = FALSE}. }

\item{nu}{a numeric value indicating the spatial smoothness parameter. Only used when \code{sp = TRUE} and \code{cov.model = "matern"}.} 

\item{ar1}{a logical value indicating whether to simulate a temporal random effect with an AR(1) process. By default, set to \code{FALSE}.}

\item{rho}{a numeric value indicating the AR(1) temporal correlation parameter. Ignored when \code{ar1 = FALSE}.}

\item{sigma.sq.t}{a numeric value indicating the AR(1) temporal variance parameter. Ignored when \code{ar1 = FALSE}.}

\item{...}{currently no additional arguments}
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
}

\value{
  A list comprised of: 

  \item{X}{a \eqn{J \times T \times p.occ}{J x T x p.occ} numeric array containing the design matrix for the occurrence portion of the occupancy model.}
  \item{X.p}{a four-dimensional numeric array with dimensions corresponding to sites, primary time periods, repeat visits, and number of detection regression coefficients. This is the design matrix used for the detection portion of the occupancy model.}
  \item{coords}{a \eqn{J \times 2}{J x 2} numeric matrix of coordinates of each occupancy site. Required for spatial models.}
  \item{w}{a \eqn{J \times 1}{J x 1} matrix of the spatial random effects. Only used to simulate data when \code{sp = TRUE}.}
  \item{psi}{a \eqn{J \times T}{J x T} matrix of the occupancy probabilities for each site during each primary time period.}
  \item{z}{a  \eqn{J \times T}{J x T} matrix of the latent occupancy states at each site during each primary time period.}
  \item{p}{a \code{J x T x max(n.rep)} array of the detection probabilities for each site, primary time period, and replicate combination. Site/time periods with fewer than \code{max(n.rep)} replicates will contain \code{NA} values.}
  \item{y}{a \code{J x T x max(n.rep)} array of the raw detection-nondetection data for each sit, primary time period, and replicate combination.}
  \item{X.p.re}{a four-dimensional numeric array containing the levels of any detection random effect included in the model. Only relevant when detection random effects are specified in \code{p.RE}.}
  \item{X.re}{a numeric matrix containing the levels of any occurrence random effect included in the model. Only relevant when occurrence random effects are specified in \code{psi.RE}.}
  \item{alpha.star}{a numeric vector that contains the simulated detection random effects for each given level of the random effects included in the detection model. Only relevant when detection random effects are included in the model.}
  \item{beta.star}{a numeric vector that contains the simulated occurrence random effects for each given level of the random effects included in the occurrence model. Only relevant when occurrence random effects are included in the model.}
  \item{eta}{a \eqn{T \times 1}{T x 1} matrix of the latent AR(1) random effects. Only included when \code{ar1 = TRUE}.}
}

\examples{
J.x <- 10
J.y <- 10
J <- J.x * J.y
# Number of time periods sampled
n.time <- sample(10, J, replace = TRUE)
n.time.max <- max(n.time)
# Replicates
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
# Occurrence --------------------------
# Fixed
beta <- c(0.4, 0.5, -0.9)
trend <- TRUE 
sp.only <- 0
psi.RE <- list(levels = c(10), 
               sigma.sq.psi = c(1))
# Detection ---------------------------
alpha <- c(-1, 0.7, -0.5)
p.RE <- list(levels = c(10), 
             sigma.sq.p = c(0.5))
# Spatial parameters ------------------
sp <- TRUE
cov.model <- "exponential"
sigma.sq <- 2
phi <- 3 / .4
nu <- 1
# Temporal parameters -----------------
ar1 <- TRUE
rho <- 0.5
sigma.sq.t <- 0.8
# Get all the data
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep, 
               beta = beta, alpha = alpha, sp.only = sp.only, trend = trend, 
               psi.RE = psi.RE, p.RE = p.RE, 
               sp = sp, cov.model = cov.model, sigma.sq = sigma.sq, phi = phi, 
               ar1 = ar1, rho = rho, sigma.sq.t = sigma.sq.t)
str(dat)
}
