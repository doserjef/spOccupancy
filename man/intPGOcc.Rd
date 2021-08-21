\name{intPGOcc}
\alias{intPGOcc}
\title{Function for Fitting Single Species Integrated Occupancy Models Using Polya-Gamma Latent Variables}


\description{
Function for fitting single species integrated occupancy models using Polya-Gamma latent variables. Data integration is done using a joint likelihood framework, assuming distinct detection models for each data source that are each connected to a single latent occupancy process. 
}

\usage{
intPGOcc(y, X, X.p, sites, n.data, starting, n.rep, n.samples, priors, 
	 n.omp.threads = 1, verbose = TRUE, n.report = 1000, ...
)
}
\arguments{
  \item{y}{Data in list format}
  
  \item{X}{Occupancy covariates}
  
  \item{X.p}{Detection covariates}
  
  \item{sites}{Indices for latent occupancy sites}
  
  \item{n.data}{Number of data sources}
  
  \item{starting}{Starting values}
  
  \item{n.rep}{Temporal replicates}
  
  \item{n.samples}{MCMC samples}
  
  \item{priors}{Priors}
  
  \item{n.omp.threads}{Number of threads}
  
  \item{verbose}{Be wordy}
  
  \item{n.report}{Value to report}
  
  \item{...}{currently no additional arguments}
}
\value{
stuff
}
