\name{spIntPGOcc}
\alias{spIntPGOcc}
\title{Function for Fitting Single Species Integrated Spatial Occupancy Models Using Polya-Gamma Latent Variables}

\description{
  The function \code{spIntPGOcc} fits single species integrated spatial occupancy models using Polya-Gamma latent variables. Models can be fit using either a full Gaussian process or a Nearest Neighbor Gaussian Process for large data sets. Data integration is done using a joint likelihood framework, assuming distinct detection models for each data source that are each conditional on a single latent occupancy process. 
}

\usage{
spIntPGOcc(y, X, X.p, sites, coords, n.data, starting, n.rep, 
	   cov.model = "exponential", n.batch, batch.length, 
	   accept.rate = 0.43, priors, tuning, n.omp.threads = 1, 
	   verbose = TRUE, NNGP = FALSE, n.neighbors = 15, 
	   search.type = "cb", n.report = 100, ...)
}

\arguments{
  \item{y}{Data}
  
  \item{X}{Occupancy covariates}
  
  \item{X.p}{Detection covariates}

  \item{sites}{Indices for latent occupancy sites}
  
  \item{coords}{Coordinates}

  \item{n.data}{Number of data sources}
  
  \item{starting}{Starting values}
  
  \item{n.rep}{Temporal replicates}
  
  \item{cov.model}{Covariance model}
  
  \item{n.batch}{Number of MCMC batches}
  
  \item{batch.length}{Length of MCMC batch}
  
  \item{accept.rate}{Target acceptance rates for adaptive Metropolis}
  
  \item{priors}{Priors}
  
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
