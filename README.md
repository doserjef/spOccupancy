# spOccupancy: An R Package for single species, multispecies, and integrated spatial occupancy models

Fits single-species and multi-species non-spatial and spatial occupancy models using Markov Chain Monte Carlo (MCMC). Models are fit using Polya-Gamma data augmentation. Spatial models are fit using either Gaussian processes or Nearest Neighbor Gaussian Processes (NNGP) for large spatial datasets. Provides functionality for data integration of multiple single species occupancy data sets using a joint likelihood framework. 

## Core Functions 

+ `PGOcc`: function for fitting single species occupancy models using Polya-Gamma latent variables. 
+ `spPGOcc`: function for fitting single species spatial occupancy models using Polya-Gamma latent variables. Models can be fit using either a full Gaussian process or a NNGP for large data sets.
+ `msPGOcc`: function for fitting multispecies occupancy models using Polya-Gamma latent variables.
+ `spMsPGOcc`: function for fitting multi-species spatial occupancy models using Polya-Gamma latent variables. Models can be fit using either a full Gaussian process or a NNGP for large data sets.
+ `intPGOcc`: function for fitting single species integrated occupancy models using Polya-Gamma latent variables. Data integration is done using a joint likelihood framework, assuming distinct detection models for each data source that are all connected to a single latent occupancy process.
+ `spIntPGOcc`: function for fitting single species integrated spatial occupancy models using Polya-Gamma latent variables. Models can be fit using either a full Gaussian process or a NNGP for large data sets. Data integration is done using a joint likelihood framework, assuming distinct detection models for each data source that are all connected to a single latent occupancy process

## Additional Functionality

+ `ppcOcc`: function for performing posterior predictive checks on `spOccupancy` model objects.
+ `waicOcc`: function for computing the Widely Applicable Information Criterion on `spOccupancy` model objects for model comparison and selection.
+ `simOcc`: simulates single-species occupancy data. Data can be optionally simulated with a spatial Gaussian process following an exponential correlation function. 
+ `simMsOcc`: simulates multi-species occupancy data. Data can be optionally simulated with a spatial Gaussian process following an exponential correlation function.
+ `simIntOcc`: simulates single-species occupancy data from multiple data sources. Data can be optionally simulated with a spatial Gaussian process following an exponential correlation function.
+ All model objects from the core functions have a set of standard extractor functions including the following: `summary`, `print`, `fitted`. 
+ Out of sample predictions are possible for all model objects using `predict`. 
+ All `spOccupancy` model fitting functions have options for k-fold cross-validation with options for parallelization.
