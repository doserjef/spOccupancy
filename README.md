
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spOccupancy <a href='https://www.doserlab.com/files/spoccupancy-web/'><img src="man/figures/logo.png" align="right" height="139" width="120"/></a>

[![](http://cranlogs.r-pkg.org/badges/grand-total/spOccupancy?color=blue)](https://CRAN.R-project.org/package=spOccupancy)
[![CRAN](https://www.r-pkg.org/badges/version/spOccupancy)](https://CRAN.R-project.org/package=spOccupancy)
[![Codecov test
coverage](https://codecov.io/gh/doserjef/spOccupancy/branch/main/graph/badge.svg)](https://codecov.io/gh/doserjef/spOccupancy?branch=main)

spOccupancy fits single-species, multi-species, and integrated spatial
occupancy models using Markov chain Monte Carlo (MCMC). Models are fit
using Pólya-Gamma data augmentation. Spatial models are fit using either
Gaussian processes or Nearest Neighbor Gaussian Processes (NNGP) for
large spatial datasets. The package provides functionality for data
integration of multiple single-species occupancy data sets using a joint
likelihood framework. For multi-species models, spOccupancy provides
functions to account for residual species correlations in a joint
species distribution model framework while accounting for imperfect
detection. `spOccupancy` also provides functions for multi-season (i.e.,
spatio-temporal) single-species occupancy models. Below we give a very
brief introduction to some of the package’s functionality, and
illustrate just one of the model fitting functions. For more
information, see the resources referenced at the bottom of this page.

## Installation

You can install the released version of `spOccupancy` from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("spOccupancy")
```

## Functionality

| `spOccupancy` Function | Description                                                               |
|------------------------|---------------------------------------------------------------------------|
| `PGOcc()`              | Single-species occupancy model                                            |
| `spPGOcc()`            | Single-species spatial occupancy model                                    |
| `intPGOcc()`           | Single-species occupancy model with multiple data sources                 |
| `spIntPGOcc()`         | Single-species spatial occupancy model with multiple data sources         |
| `msPGOcc()`            | Multi-species occupancy model                                             |
| `spMsPGOcc()`          | Multi-species spatial occupancy model                                     |
| `lfJSDM()`             | Joint species distribution model without imperfect detection              |
| `sfJSDM()`             | Spatial joint species distribution model without imperfect detection      |
| `lfMsPGOcc()`          | Multi-species occupancy model with species correlations                   |
| `sfMsPGOcc()`          | Multi-species spatial occupancy model with species correlations           |
| `intMsPGOcc()`         | Multi-species occupancy model with multiple data sources                  |
| `tPGOcc()`             | Single-species multi-season occupancy model                               |
| `stPGOcc()`            | Single-species multi-season spatio-temporal occupancy model               |
| `svcPGBinom()`         | Single-species spatially-varying coefficient GLM                          |
| `svcPGOcc()`           | Single-species spatially-varying coefficient occupancy model              |
| `svcTPGBinom()`        | Single-species spatially-varying coefficient multi-season GLM             |
| `svcTPGOcc()`          | Single-species spatially-varying coefficient multi-season occupancy model |
| `svcMsPGOcc()`         | Multi-species spatially-varying coefficient occupancy model               |
| `tMsPGOcc()`           | Multi-species, multi-season occupancy model                               |
| `stMsPGOcc()`          | Multi-species, multi-season spatial occupancy model                       |
| `svcTMsPGOcc()`        | Multi-species, multi-season spatially-varying coefficient occupancy model |
| `tIntPGOcc()`          | Multi-season occupancy model with multiple data sources                   |
| `stIntPGOcc()`         | Spatial multi-season occupancy model with multiple data sources           |
| `svcTIntPGOcc()`       | SVC multi-season occupancy model with multiple data sources               |
| `postHocLM()`          | Fit a linear (mixed) model using estimates from a previous model fit      |
| `ppcOcc()`             | Posterior predictive check using Bayesian p-values                        |
| `waicOcc()`            | Compute Widely Applicable Information Criterion (WAIC)                    |
| `updateMCMC()`         | Update an existing model object with more MCMC samples (in development)   |
| `simOcc()`             | Simulate single-species occupancy data                                    |
| `simTOcc()`            | Simulate single-species multi-season occupancy data                       |
| `simBinom()`           | Simulate detection-nondetection data with perfect detection               |
| `simTBinom()`          | Simulate multi-season detection-nondetection data with perfect detection  |
| `simMsOcc()`           | Simulate multi-species occupancy data                                     |
| `simTMsOcc()`          | Simulate multi-species, multi-season occupancy data                       |
| `simIntOcc()`          | Simulate single-species occupancy data from multiple data sources         |
| `simIntMsOcc()`        | Simulate multi-species occupancy data from multiple data sources          |
| `simTIntOcc()`         | Simulate multi-season occupancy data from multiple data sources           |

## Example usage

### Load package and data

To get started with `spOccupancy` we load the package and an example
data set. We use data on twelve foliage-gleaning birds from the [Hubbard
Brook Experimental Forest](https://hubbardbrook.org/), which is
available in the `spOccupancy` package as the `hbef2015` object. Here we
will only work with one bird species, the black-throated blue warbler
(BTBW), and so we subset the `hbef2015` object to only include this
species.

``` r
library(spOccupancy)
data(hbef2015)
sp.names <- dimnames(hbef2015$y)[[1]]
btbwHBEF <- hbef2015
btbwHBEF$y <- btbwHBEF$y[sp.names == "BTBW", , ]
```

### Fit a spatial occupancy model using `spPGOcc()`

Below we fit a single-species spatial occupancy model to the BTBW data
using a Nearest Neighbor Gaussian Process. We use the default priors and
initial values for the occurrence (`beta`) and detection (`alpha`)
coefficients, the spatial variance (`sigma.sq`), the spatial decay
parameter (`phi`), the spatial random effects (`w`), and the latent
occurrence values (`z`). We assume occurrence is a function of linear
and quadratic elevation along with a spatial random intercept. We model
detection as a function of linear and quadratic day of survey and linear
time of day the survey occurred.

``` r
# Specify model formulas
btbw.occ.formula <- ~ scale(Elevation) + I(scale(Elevation)^2)
btbw.det.formula <- ~ scale(day) + scale(tod) + I(scale(day)^2)
```

We run the model using an adaptive MCMC sampler with a target acceptance
rate of 0.43. We run 3 chains of the model in parallel each for 20,000
iterations split into 800 batches each of length 25. For each chain, we
discard the first 8000 iterations as burn-in and use a thinning rate of
4 for a resulting 9000 samples from the joint posterior. We fit the
model using 5 nearest neighbors and an exponential correlation function.
We also specify the `k.fold` argument to perform 2-fold cross-validation
after fitting the full model. Run `?spPGOcc` for more detailed
information on all function arguments.

``` r
# Run the model
out <- spPGOcc(occ.formula = btbw.occ.formula,
               det.formula = btbw.det.formula,
               data = btbwHBEF, n.batch = 800, batch.length = 25,
               accept.rate = 0.43, cov.model = "exponential", 
               NNGP = TRUE, n.neighbors = 5, n.burn = 8000, 
               n.thin = 4, n.chains = 3, verbose = FALSE, 
               parallel.chains = TRUE, k.fold = 2, k.fold.threads = 2)
```

This will produce a large output object, and you can use `str(out)` to
get an overview of what’s in there. Here we use the `summary()` function
to print a concise but informative summary of the model fit.

``` r
summary(out)
#> 
#> Call:
#> spPGOcc(occ.formula = btbw.occ.formula, det.formula = btbw.det.formula, 
#>     data = btbwHBEF, cov.model = "exponential", NNGP = TRUE, 
#>     n.neighbors = 5, n.batch = 800, batch.length = 25, accept.rate = 0.43, 
#>     verbose = FALSE, n.burn = 8000, n.thin = 4, n.chains = 3, 
#>     parallel.chains = TRUE, k.fold = 2, k.fold.threads = 2)
#> 
#> Samples per Chain: 20000
#> Burn-in: 8000
#> Thinning Rate: 4
#> Number of Chains: 3
#> Total Posterior Samples: 9000
#> Run Time (min): 0.648
#> 
#> Occurrence (logit scale): 
#>                          Mean     SD    2.5%     50%   97.5%   Rhat  ESS
#> (Intercept)            3.9819 0.5888  3.0184  3.9144  5.2920 1.0199  298
#> scale(Elevation)      -0.5278 0.2167 -0.9750 -0.5218 -0.1174 1.0042 2148
#> I(scale(Elevation)^2) -1.1603 0.2145 -1.6386 -1.1399 -0.7961 1.0133  447
#> 
#> Detection (logit scale): 
#>                    Mean     SD    2.5%     50%  97.5%   Rhat  ESS
#> (Intercept)      0.6626 0.1140  0.4390  0.6626 0.8868 0.9999 8178
#> scale(day)       0.2902 0.0708  0.1528  0.2902 0.4278 1.0002 9054
#> scale(tod)      -0.0326 0.0692 -0.1669 -0.0331 0.1061 1.0003 9897
#> I(scale(day)^2) -0.0754 0.0862 -0.2451 -0.0758 0.0965 1.0004 9000
#> 
#> Spatial Covariance: 
#>            Mean     SD   2.5%    50%  97.5%   Rhat ESS
#> sigma.sq 1.1206 1.0169 0.1908 0.8124 3.8735 1.0230 137
#> phi      0.0090 0.0083 0.0007 0.0055 0.0281 1.0456  66
```

### Posterior predictive check

The function `ppcOcc` performs a posterior predictive check on the
resulting list from the call to `spPGOcc`. For binary data, we need to
perform Goodness of Fit assessments on some binned form of the data
rather than the raw binary data. Below we perform a posterior predictive
check on the data grouped by site with a Freeman-Tukey fit statistic,
and then use the `summary` function to summarize the check with a
Bayesian p-value.

``` r
ppc.out <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out)
#> 
#> Call:
#> ppcOcc(object = out, fit.stat = "freeman-tukey", group = 1)
#> 
#> Samples per Chain: 20000
#> Burn-in: 8000
#> Thinning Rate: 4
#> Number of Chains: 3
#> Total Posterior Samples: 9000
#> 
#> Bayesian p-value:  0.4939 
#> Fit statistic:  freeman-tukey
```

### Model selection using WAIC and k-fold cross-validation

The `waicOcc` function computes the Widely Applicable Information
Criterion (WAIC) for use in model selection and assessment (note that
due to Monte Carlo error your results will differ slightly).

``` r
waicOcc(out)
#>       elpd         pD       WAIC 
#> -681.77267   21.27419 1406.09372
```

Alternatively, we can perform k-fold cross-validation (CV) directly in
our call to `spPGOcc` using the `k.fold` argument and compare models
using a deviance scoring rule. We fit the model with `k.fold = 2` and so
below we access the deviance scoring rule from the 2-fold
cross-validation. If we have additional candidate models to compare this
model with, then we might select for inference the one with the lowest
value of this CV score.

``` r
out$k.fold.deviance
#> [1] 1415.012
```

### Prediction

Prediction is possible using the `predict` function, a set of occurrence
covariates at the new locations, and the spatial coordinates of the new
locations. The object `hbefElev` contains elevation data across the
entire Hubbard Brook Experimental Forest. Below we predict BTBW
occurrence across the forest, which are stored in the `out.pred` object.

``` r
# First standardize elevation using mean and sd from fitted model
elev.pred <- (hbefElev$val - mean(btbwHBEF$occ.covs[, 1])) / sd(btbwHBEF$occ.covs[, 1])
coords.0 <- as.matrix(hbefElev[, c('Easting', 'Northing')])
X.0 <- cbind(1, elev.pred, elev.pred^2)
out.pred <- predict(out, X.0, coords.0, verbose = FALSE)
```

## Learn more

The `vignette("modelFitting")` provides a more detailed description and
tutorial of the core functions in `spOccupancy`. For full statistical
details on the MCMC samplers for core functions in `spOccupancy`, see
`vignette("mcmcSamplers")`. In addition, see [the introductory
spOccupancy paper](https://doi.org/10.1111/2041-210X.13897) that
describes the package in more detail (Doser et al. 2022). For a detailed
description and tutorial of joint species distribution models in
`spOccupancy` that account for residual species correlations, see
`vignette("factorModels")`, `vignette("mcmcFactorModels")`, and our
[open-access paper](https://doi.org/10.1002/ecy.4137) (Doser et
al. 2023). For a description and tutorial of multi-season
(spatio-temporal) occupancy models in `spOccupancy`, see
`vignette("spaceTimeModels")`. For a tutorial on spatially-varying
coefficient models in `spOccupancy`, see `vignette("svcModels")` and
`vignette(mcmcSVCModels)` and our associated papers that describe the
[methods](https://www.doserlab.com/files/pubs/doser2024JABES.pdf) (Doser
et al. 2024A) and [applications to
ecology](https://onlinelibrary.wiley.com/doi/epdf/10.1111/geb.13814)
(Doser et al. 2024B) in much more detail.

## References

Doser, J. W., Finley, A. O., Kery, M., and Zipkin, E. F. (2022).
spOccupancy: An R package for single-species, multi-species, and
integrated spatial occupancy models. Methods in Ecology and Evolution.
13(8) 1670-1678. <https://doi.org/10.1111/2041-210X.13897>.

Doser, J. W., Finley, A. O., and Banerjee, S. (2023). Joint species
distribution models with imperfect detection for high-dimensional
spatial data. Ecology, 104(9), e4137.
<https://doi.org/10.1002/ecy.4137>.

Doser, J. W., Finley, A. O., Saunders, S. P., Kéry, M., Weed, A. S., &
Zipkin, E. F. (2024A). Modeling complex species-environment
relationships through spatially-varying coefficient occupancy models.
Journal of Agricultural, Biological and Environmental Statistics.
<https://doi.org/10.1007/s13253-023-00595-6>.

Doser, J. W., Kéry, M., Saunders, S. P., Finley, A. O., Bateman, B. L.,
Grand, J., Reault, S., Weed, A. S., & Zipkin, E. F. (2024B). Guidelines
for the use of spatially varying coefficients in species distribution
models. Global Ecology and Biogeography, 33, e13814.
<https://doi.org/10.1111/geb.13814>
