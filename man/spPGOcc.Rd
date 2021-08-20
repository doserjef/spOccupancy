\name{spPGOcc}
\alias{spPGOcc}
\title{Function for Fitting Single Species Spatial Occupancy Models Using Polya-Gamma Latent Variables}

\description{
  The function \code{spPGOcc} fits single species spatial occupancy models using Polya-Gamma latent variables. Models can be fit using either a full Gaussian process or a Nearest Neighbor Gaussian Process for large data sets. 
}

\usage{
spPGOcc(y, X, X.p, starting, coords, n.rep, n.batch, 
	batch.length, accept.rate = 0.43, priors, 
	cov.model = "exponential", tuning, n.omp.threads = 1, 
	verbose = TRUE, NNGP = FALSE, n.neighbors = 15, 
	search.type = "cb", n.report = 100, ...)
}

\arguments{
  \item{y}{Data}
  
  \item{X}{Occupancy covariates}
  
  \item{X.p}{Detection covariates}
  
  \item{starting}{Starting values}
  
  \item{coords}{Coordinates}
  
  \item{n.rep}{Temporal replicates}
  
  \item{n.batch}{Number of MCMC batches}
  
  \item{batch.length}{Length of MCMC batch}
  
  \item{accept.rate}{Target acceptance rates for adaptive Metropolis}
  
  \item{priors}{Priors}
  
  \item{cov.model}{Covariance model}
  
  \item{tuning}{Tuning parameters}
  
  \item{n.omp.threads}{Number of threads}
  
  \item{verbose}{Be wordy}
  
  \item{NNGP}{Should we fit an NNGP?}
  
  \item{n.neighbors}{Number of nearest neighbors for NNGP.}
  
  \item{search.type}{Search type algorithm for nearest neighbors}
  
  \item{n.report}{Value to report}
  
  \item{...}{currently no additional arguments}
}
\value{
stuff
}

