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
for data integration of multiple occupancy data sets using a 
joint likelihood framework. Details on data integration are given in 
Miller, Pacifici, Sanderlin, and Reich (2019). Details on single-species and 
multi-species models are found in MacKenzie et al. (2002) and Dorazio and Royle (2005), 
respectively. Details on the package functionality is given in Doser et al. (2022), 
Doser, Finley, Banerjee (2023), Doser et al. (2024a,b). 
See \code{citation('spOccupancy')} for how to cite spOccupancy in publications.

\strong{Single-species models}

\code{\link{PGOcc}} fits single-species occupancy models.

\code{\link{spPGOcc}} fits single-species spatial occupancy models.

\code{\link{intPGOcc}} fits single-species integrated occupancy models (i.e., an occupancy model with multiple data sources).

\code{\link{spIntPGOcc}} fits single-species integrated spatial occupancy models.

\code{\link{tPGOcc}} fits a multi-season single-species occupancy model.

\code{\link{stPGOcc}} fits a multi-season single-species spatial occupancy model.

\code{\link{svcPGBinom}} fits a single-species spatially-varying coefficient GLM.

\code{\link{svcPGOcc}} fits a single-species spatially-varying coefficient occupancy model.

\code{\link{svcTPGBinom}} fits a single-species spatially-varying coefficient multi-season GLM.

\code{\link{svcTPGOcc}} fits a single-species spatially-varying coefficient multi-season occupancy model.


\strong{Multi-species models}

\code{\link{msPGOcc}} fits multi-species occupancy models.

\code{\link{spMsPGOcc}} fits multi-species spatial occupancy models.

\code{\link{lfJSDM}} fits a joint species distribution model without imperfect detection. 

\code{\link{sfJSDM}} fits a spatial joint species distribution model without imperfect detection.

\code{\link{lfMsPGOcc}} fits a joint species distribution model with imperfect detection (i.e., a multi-species occupancy model with residual species correlations).

\code{\link{sfMsPGOcc}} fits a spatial joint species distribution model with imperfect detection.

\code{\link{svcMsPGOcc}} fits a multi-species spatially-varying coefficient occupancy model.

\code{\link{tMsPGOcc}} fits a multi-season multi-species occupancy model.

\code{\link{stMsPGOcc}} fits a multi-season multi-species spatial occupancy model.

\code{\link{svcTMsPGOcc}} fits a multi-season multi-species spatially-varying coefficient occupancy model.


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

\code{\link{simTMsOcc}} simulates multi-species multi-season occupancy data from multiple data sources. 

\strong{Miscellaneous}

\code{\link{postHocLM}} fits post-hoc linear (mixed) models.

\code{\link{getSVCSamples}} extracts spatially varying coefficient MCMC samples.

\code{\link{updateMCMC}} updates a spOccupancy or spAbundance model object with more MCMC iterations.

All objects from model-fitting functions have support with the \code{summary} function for 
displaying a concise summary of model results, the \code{fitted} function for extracting 
model fitted values, and the \code{predict} function for predicting occupancy and/or detection
across an area of interest.
}

\references{

  Doser, J. W., Finley, A. O., Kery, M., & Zipkin, E. F. (2022). 
  spOccupancy: An R package for single-species, multi-species, and 
  integrated spatial occupancy models. 
  Methods in Ecology and Evolution, 13, 1670-1678. \doi{10.1111/2041-210X.13897}.

  Doser, J. W., Finley, A. O., & Banerjee, S. (2023). Joint species 
  distribution models with imperfect detection for high-dimensional 
  spatial data. Ecology, 104(9), e4137. \doi{10.1002/ecy.4137}.

  Doser, J. W., Finley, A. O., Saunders, S. P., Kery, M., Weed, A. S., &
  Zipkin, E. F. (2024A). Modeling complex species-environment
  relationships through spatially-varying coefficient occupancy models.
  Journal of Agricultural, Biological and Environmental Statistics.
  \doi{10.1007/s13253-023-00595-6}.

  Doser, J. W., Kery, M., Saunders, S. P., Finley, A. O., Bateman, B. L.,
  Grand, J., Reault, S., Weed, A. S., & Zipkin, E. F. (2024B). Guidelines
  for the use of spatially varying coefficients in species distribution
  models. Global Ecology and Biogeography, 33, e13814.
  \doi{10.1111/geb.13814}.
 
}

\author{
  Jeffrey W. Doser, Andrew O. Finley, Marc Kery
}

\docType{package}
