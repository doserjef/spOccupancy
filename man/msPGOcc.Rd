\name{msPGOcc}
\alias{msPGOcc}
\title{Function for Fitting Multi-Species Occupancy Models Using Polya-Gamma Latent Variables}


\description{
Function for fitting multi-species occupancy models using Polya-Gamma latent variables. 
}

\usage{
msPGOcc(
  y,
  X,
  X.p,
  starting,
  n.rep,
  n.samples,
  priors,
  n.omp.threads = 1,
  verbose = TRUE,
  n.report = 100,
  ...
)
}
\arguments{
\item{y}{Data}

\item{X}{Occupancy covariates}

\item{X.p}{Detection covariates}

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
