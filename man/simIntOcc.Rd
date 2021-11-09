\name{simIntOcc}
\alias{simIntOcc}
\title{Simulate Single-Species Detection-Nondetection Data from Multiple Data Sources}

\description{
  The function \code{simIntOcc} simulates single-species detection-nondetection data from multiple data sources for simulation studies, power assessments, or function testing of integrated occupancy models. Data can optionally be simulated with a spatial Gaussian Process on the occurrence process.
}

\usage{
simIntOcc(n.data, J.x, J.y, J.obs, n.rep, beta, alpha,
          sp = FALSE, cov.model, sigma.sq, phi, nu, ...)
}
\arguments{

\item{n.data}{an integer indicating the number of detection-nondetection data sources to simulate.}

\item{J.x}{a single numeric value indicating the number of sites across the region of interest along the horizontal axis. Total number of sites across the simulated region of interest is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{J.y}{a single numeric value indicating the number of sites across the region of interest along the vertical axis. Total number of sites across the simulated region of interest is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{J.obs}{a numeric vector of length \code{n.data} containing the number of sites to simulate each data source at. Data sources can be obtained at completely different sites, the same sites, or anywhere inbetween. Maximum number of sites a given data source is available at is equal to \eqn{J = J.x \times J.y}{J = J.x x J.y}.}

\item{n.rep}{a list of length \code{n.data}. Each element is a numeric vector with length corresponding to the number of sites that given data source is observed at (in \code{J.obs}). Each vector indicates the number of repeat visits at each of the sites for a given data source.}

\item{beta}{a numeric vector containing the intercept and regression coefficient parameters for the occurrence portion of the single-species occupancy model.}

\item{alpha}{a list of length \code{n.data}. Each element is a numeric vector containing the intercept and regression coefficient parameters for the detection portion of the single-species occupancy model for each data source.}		

\item{sp}{a logical value indicating whether to simulate a spatially-explicit occupancy model with a Gaussian process and exponential correlation function. By default set to \code{FALSE}.}

\item{cov.model}{a quoted keyword that specifies the covariance function used to model the spatial dependence structure among the latent occurrence values.  Supported covariance model key words are: \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and \code{"gaussian"}.}

\item{sigma.sq}{a numeric value indicating the spatial variance parameter. Ignored when \code{sp = FALSE}.}

\item{phi}{a numeric value indicating the spatial range parameter. Ignored when \code{sp = FALSE}. }

\item{nu}{a numeric value indicating the spatial smoothness parameter. Only used when \code{sp = TRUE} and \code{cov.model = "matern"}.}

\item{...}{currently no additional arguments}
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  A list comprised of: 

  \item{X.obs}{a numeric design matrix for the occurrence portion of the model. This matrix contains the intercept and regression coefficients for only the observed sites.}
  \item{X.pred}{a numeric design matrix for the occurrence portion of the model at sites where there are no observed data sources.} 
  \item{X.p}{a list of design matrices for the detection portions of the integrated occupancy model. Each element in the list is a design matrix of detection covariates for each data source.}
  \item{coords.obs}{a numeric matrix of coordinates of each observed site. Required for spatial models.}
  \item{coords.pred}{a numeric matrix of coordinates of each site in the study region without any data sources. Only used for spatial models.}
  \item{D.obs}{a distance matrix of observed sites. Only used for spatial models.}
  \item{D.pred}{a distance matrix of sites in the study region without any observed data. Only used for spatial models.}
  \item{w.obs}{a matrix of the spatial random effects at observed locations. Only used to simulate data when \code{sp = TRUE}}. 
  \item{w.pred}{a matrix of the spatial random random effects at locations without any observation.}
  \item{psi.obs}{a matrix of the occurrence probabilities for each observed site.}
  \item{psi.pred}{a matrix of the occurrence probabilities for sites without any observations.}
  \item{z.obs}{a vector of the latent occurrence states at each observed site.}
  \item{z.pred}{a vector of the latent occurrence states at each site without any observations.}
  \item{p}{a list of detection probability matrices for each of the \code{n.data} data sources.}
  \item{y}{a list of matrices of the raw detection-nondetection data for each site and replicate combination.}

}

\examples{
set.seed(400)

# Simulate Data -----------------------------------------------------------
J.x <- 15
J.y <- 15
J.all <- J.x * J.y
# Number of data sources.
n.data <- 4
# Sites for each data source. 
J.obs <- sample(ceiling(0.2 * J.all):ceiling(0.5 * J.all), n.data, replace = TRUE)
# Replicates for each data source.
n.rep <- list()
for (i in 1:n.data) {
  n.rep[[i]] <- sample(1:4, size = J.obs[i], replace = TRUE)
}
# Occupancy covariates
beta <- c(0.5, 1, -3)
p.occ <- length(beta)
# Detection covariates
alpha <- list()
for (i in 1:n.data) {
  alpha[[i]] <- runif(sample(1:4, 1), -1, 1)
}
p.det.long <- sapply(alpha, length)
p.det <- sum(p.det.long)
sigma.sq <- 2
phi <- 3 / .5
sp <- TRUE

# Simulate occupancy data. 
dat <- simIntOcc(n.data = n.data, J.x = J.x, J.y = J.y, J.obs = J.obs, 
		 n.rep = n.rep, beta = beta, alpha = alpha, sp = TRUE, 
		 cov.model = 'gaussian', sigma.sq = sigma.sq, phi = phi)
}
