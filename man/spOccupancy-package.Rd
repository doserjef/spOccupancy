\name{spOccupancy-package}
\alias{spOccupancy-package}
\alias{spOccupancy}
\keyword{package}
\title{Single-Species, Multi-Species, and Integrated Spatial Occupancy Models}

\description{
Fits single-species, multi-species, and integrated non-spatial and spatial 
occupancy models using Markov Chain Monte Carlo (MCMC). Models are fit using 
Polya-Gamma data augmentation detailed in Polson, Scott, and Windle (2013). 
Spatial models are fit using either Gaussian processes or Nearest Neighbor 
Gaussian Processes (NNGP) for large spatial datasets. Details on NNGPs are 
given in Datta, Banerjee, Finley, and Gelfand (2016). Provides functionality 
for data integration of multiple single-species occupancy data sets using a 
joint likelihood framework. Details on data integration are given in 
Miller, Pacifici, Sanderlin, and Reich (2019). Details on single-species and 
multi-species models are found in MacKenzie et al. (2002) and Dorazio and Royle (2005), 
respectively. Details on the package functionality is given in Doser et al. (2022) and 
Doser, Finley, Banerjee (2022). See \code{citation('spOccupancy')} for how to 
cite spOccupancy in publications.

\strong{Model Fitting Functions}

\code{\link{PGOcc}} fits single-species occupancy models.

\code{\link{spPGOcc}} fits single-species spatial occupancy models.

\code{\link{msPGOcc}} fits multi-species occupancy models.

\code{\link{spMsPGOcc}} fits multi-species spatial occupancy models.

\code{\link{intPGOcc}} fits single-species integrated occupancy models (i.e., an occupancy model with multiple data sources).

\code{\link{spIntPGOcc}} fits single-species integrated spatial occupancy models.

\code{\link{lfJSDM}} fits a joint species distribution model without imperfect detection. 

\code{\link{sfJSDM}} fits a spatial joint species distribution model without imperfect detection.

\code{\link{lfMsPGOcc}} fits a joint species distribution model with imperfect detection (i.e., a multi-species occupancy model with residual species correlations).

\code{\link{sfMsPGOcc}} fits a spatial joint species distribution model with imperfect detection.

\code{\link{tPGOcc}} fits a multi-season single-species occupancy model.

\code{\link{stPGOcc}} fits a multi-season single-species spatial occupancy model.

\code{\link{svcPGBinom}} fits a single-species spatially-varying coefficient GLM.

\code{\link{svcPGOcc}} fits a single-species spatially-varying coefficient occupancy model.

\code{\link{svcTPGBinom}} fits a single-species spatially-varying coefficient multi-season GLM.

\code{\link{svcTPGOcc}} fits a single-species spatially-varying coefficient multi-season occupancy model.


\strong{Goodness of Fit and Model Assessment Functions}

\code{\link{ppcOcc}} performs posterior predictive checks.

\code{\link{waicOcc}} computes the Widely Applicable Information Criterion for spOccupancy model objects.


\strong{Data Simulation Functions}

\code{\link{simOcc}} simulates single-species occupancy data.

\code{\link{simTOcc}} simulates single-species multi-season occupancy data.

\code{\link{simBinom}} simulates detection-nondetection data with perfect detection.


\code{\link{simTBinom}} simulates multi-season detection-nondetection data with perfect detection.

\code{\link{simMsOcc}} simulates multi-species occupancy data.

\code{\link{simIntOcc}} simulates single-species occupancy data from multiple data sources. 


All objects from model-fitting functions have support with the \code{summary} function for 
displaying a concise summary of model results, the \code{fitted} function for extracting 
model fitted values, and the \code{predict} function for predicting occupancy and/or detection
across an area of interest.
}

\references{

  Doser, J. W., Finley, A. O., Kéry, M., & Zipkin, E. F. (2022). 
  spOccupancy: An R package for single‐species, multi‐species, and 
  integrated spatial occupancy models. Methods in Ecology and Evolution.

  Doser, J. W., Finley, A. O., & Banerjee, S. (2022). Joint species 
  distribution models with imperfect detection for high-dimensional 
  spatial data. arXiv preprint arXiv:2204.02707.
}

\author{
  Jeffrey W. Doser, Andrew O. Finley, Marc Kéry
}

\docType{package}
