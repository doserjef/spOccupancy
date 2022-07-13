# spOccupancy 0.4.0

+ Major new functionality for fitting multi-season (i.e., spatio-temporal) single-species occupancy models using the functions `tPGOcc()` and `stPGOcc()`. 
+ Fixed a bug in calculation of the detection probability values in `fitted()` functions for all spOccupancy model objects. See [this Github issue](https://github.com/doserjef/spOccupancy/issues/10) for more details. 
+ Fixed an error that occurred when predicting for multi-species models and setting `ignore.RE = TRUE`.  
+ Fixed other small bugs that caused model fitting functions to break under specific circumstances.

# spOccupancy 0.3.2

+ Fixed a bug in `waicOcc()` for integrated models (`intPGOcc()` and `spIntPGOcc()`) that sometimes resulted in incorrect estimates of WAIC for data sets other than the first data set. We strongly encourage users who have used `waicOcc()` with an integrated model to rerun their analyses using v0.3.2. 
+ Fixed a bug introduced in v0.3.0 that sometimes resulted in incorrect predictions from a spatially-explicit model with non-spatial random effects in the occurrence portion of the model. We strongly encourage users who have used `predict()` on a spatially-explicit model with non-spatial random effects in the occurrence portion of the model to rerun their analyses using v0.3.2.
+ Users can now specify a uniform prior on the spatial variance parameter instead of an inverse-Gamma prior. We also allow users to fix the value of the spatial variance parameter at the initial value. See the reference pages of spatially-explicit functions for more details. 
+ Slight changes in the information printed when fitting spatially-explicit models. 
+ Removed dependency on spBayes to pass CRAN checks. 

# spOccupancy 0.3.1

+ Fixed two small problems with `intPGOcc()` and `spIntPGOcc()` that were accidentally introduced in v0.3.0. See [this Github issue](https://github.com/doserjef/spOccupancy/issues/8) for more details.
+ Adapted C/C++ code to properly handle characters strings when calling Fortran BLAS/LAPACK routines following the new requirements for R 4.2.0.  

# spOccupancy 0.3.0

spOccupancy Version 0.3.0 contains numerous substantial updates that provide new functionality, improved computational performance for model fitting and subsequent model checking/comparison, and minor bug fixes. The changes include: 

+ Additional functionality for fitting spatial and non-spatial multi-species occupancy models with residual species correlations (i.e., joint species distribution models with imperfect detection). See documentation for `lfMsPGOcc()` and `sfMsPGOcc()`. We also included the functions `lfJSDM()` and `sfJSDM()` which are more typical joint species distribution models that fail to explicitly account for imperfect detection.
+ All single-species and multi-species models allow for unstructured random intercepts in both the occurrence and detection portions of the occupancy model. Prior to this version, random intercepts were not supported in the occurrence portion of spatially-explicit models. 
+ `predict()` functions for single-species and multi-species models now include the argument `type`, which allows for prediction of detection probability (`type = 'detection'`) at a set of covariate values as well as predictions of occurrence (`type = 'occupancy'`). 
+ All models are substantially faster than version 0.2.1. We improved performance by implementing a change in how we sample the latent Polya-Gamma variables in the detection component of the model. This results in substantial increases in speed for models where the number of replicates varies across sites. We additionally updated how non-spatial random effects were sampled, which also contributes to improved computational performance.
+ All model fitting functions now include the object `like.samples` in the resulting model object, which contains model likelihood values needed for calculation of WAIC. This leads to much shorter run times for `waicOcc()` compared to previous versions.
+ All `fitted.*()` functions now return both the fitted values and the estimated detection probability samples from a fitted `spOccupancy` model. 
+ Improved error handling for models with missing values and random effects.
+ Added the argument `ignore.RE` to all `predict()` functions. If non-spatial random intercepts are included when fitting the model, setting `ignore.RE = TRUE` will yield predictions that ignore the values of the random effects. If `ignore.RE = FALSE`, the model will predict new values using the random intercepts for both sampled and non-sampled levels of the effects.  
+ Fixed a bug in the cross-validation component of all `spOccupancy` model fitting functions that occurred when random effects were included in the occurrence and/or detection component of the model.
+ Fixed minor bug in `simOcc()` and `simMsOcc()` that prevented simulating data with multiple random intercepts on detection. 
+ Fixed minor bug in spatially-explicit models that resulted in an error when setting `NNGP = FALSE` and not specifying initial values for the spatial range parameter `phi`. 
+ Fixed a bug in the `predict()` functions for `spMsPGOcc` and `spPGOcc` objects that resulted in potentially inaccurate predictions when `n.omp.threads` > 1. 

# spOccupancy 0.2.1

+ Minor changes related to arguments in C++ code in header files to pass CRAN additional issues.

# spOccupancy 0.2.0

+ Added an `n.chains` argument to all model-fitting functions for running multiple chains in sequence.
+ Added posterior means, standard deviations, Gelman-Rubin diagnostic (Rhat) and Effective Sample Size (ESS) to `summary` displays for each model-fitting function.
+ Fixed spatially-explicit `predict` functions to return occurrence probabilities at sampled sites instead of NAs.

# spOccupancy 0.1.3

+ Minor bug fixes related to memory allocation in C++ code.

# spOccupancy 0.1.2

* This is the first release of `spOccupancy`.
