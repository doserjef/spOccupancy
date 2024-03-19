\name{updateMCMC}
\alias{updateMCMC}
\title{Update a spOccupancy or spAbundance model run with more MCMC iterations}

\usage{
updateMCMC(object, n.batch, n.samples, n.burn = 0, n.thin, 
           keep.orig = TRUE, verbose = TRUE, n.report = 100, ...)
}

\description{
  Function for updating a previously run spOccupancy or spAbundance model with 
  additional MCMC iterations. This function is useful for situations where a model
  is run for a long time but convergence/adequate mixing of the MCMC chains is not 
  reached. Instead of re-running the entire model again, this function allows you 
  to pick up where you left off. This function is currently in development, and only functions
  with the following spOccupancy and spAbundance model objects: msAbund and sfJSDM. 
}

\arguments{
  \item{object}{a \code{spOccupancy} or \code{spAbundance} model object. Currently 
    supports objects of class \code{msAbund} and \code{sfJSDM}.}

  \item{n.batch}{the number of additional MCMC batches in each chain 
    to run for the adaptive MCMC sampler. Only valid for model types fit with 
    an adaptive MCMC sampler}

  \item{n.samples}{the number of posterior samples to collect in each chain. Only 
    valid for model types that are run with a fully Gibbs sampler and have 
    \code{n.samples} as an argument in the original model fitting function.}

  \item{n.burn}{the number of samples out of the total \code{n.batch * batchlength} to 
  discard as burn-in for each chain from the updated samples. Note this argument
  does not discard samples from the previous model run, and rather only applies
  to the samples in the updated run of the model. Defaults to 0}

  \item{n.thin}{the thinning interval for collection of MCMC samples in 
    the updated model run. The thinning occurs after the \code{n.burn} 
    samples are discarded. Default value is set to 1.}

  \item{keep.orig}{A logical value indicating whether or not the samples from 
    the original run of the model should be kept or discarded.}

  \item{verbose}{if \code{TRUE}, messages about data preparation, 
    model specification, and progress of the sampler are printed to the screen. 
    Otherwise, no messages are printed.}

  \item{n.report}{the interval to report Metropolis sampler acceptance and
    MCMC progress.}

  \item{...}{currently no additional arguments}
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
}

\value{
An object of the same class as the original model fit provided in the \code{argument} object. See the manual page for the original model type for complete details.
}

\examples{
}
