
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spOccupancy <a href='https://www.jeffdoser.com/files/spoccupancy-web/'><img src="man/figures/logo.png" align="right" height="139" width="120"/></a>

[![](http://cranlogs.r-pkg.org/badges/grand-total/spOccupancy?color=blue)](https://CRAN.R-project.org/package=spOccupancy)
[![CRAN](https://www.r-pkg.org/badges/version/spOccupancy)](https://CRAN.R-project.org/package=spOccupancy)
[![Codecov test
coverage](https://codecov.io/gh/doserjef/spOccupancy/branch/main/graph/badge.svg)](https://app.codecov.io/gh/doserjef/spOccupancy?branch=main)

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
rate of 0.43. We run 3 chains of the model each for 10,000 iterations
split into 400 batches each of length 25. For each chain, we discard the
first 2000 iterations as burn-in and use a thinning rate of 4 for a
resulting 6000 samples from the joint posterior. We fit the model using
5 nearest neighbors and an exponential correlation function. We also
specify the `k.fold` argument to perform 2-fold cross-validation after
fitting the full model. Run `?spPGOcc` for more detailed information on
all function arguments.

``` r
# Run the model
out <- spPGOcc(occ.formula = btbw.occ.formula,
               det.formula = btbw.det.formula,
               data = btbwHBEF, n.batch = 400, batch.length = 25,
               accept.rate = 0.43, cov.model = "exponential", 
               NNGP = TRUE, n.neighbors = 5, n.burn = 2000, 
               n.thin = 4, n.chains = 3, verbose = FALSE, k.fold = 2)
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
#>     n.neighbors = 5, n.batch = 400, batch.length = 25, accept.rate = 0.43, 
#>     verbose = FALSE, n.burn = 2000, n.thin = 4, n.chains = 3, 
#>     k.fold = 2)
#> 
#> Samples per Chain: 10000
#> Burn-in: 2000
#> Thinning Rate: 4
#> Number of Chains: 3
#> Total Posterior Samples: 6000
#> Run Time (min): 0.7627
#> 
#> Occurrence (logit scale): 
#>                          Mean     SD    2.5%     50%   97.5%   Rhat ESS
#> (Intercept)            4.0762 0.6183  3.0681  4.0071  5.4855 1.0019 261
#> scale(Elevation)      -0.5203 0.2249 -0.9772 -0.5132 -0.0825 1.0039 984
#> I(scale(Elevation)^2) -1.1806 0.2257 -1.6958 -1.1593 -0.7969 1.0033 244
#> 
#> Detection (logit scale): 
#>                    Mean     SD    2.5%     50%  97.5%   Rhat  ESS
#> (Intercept)      0.6647 0.1155  0.4398  0.6623 0.8955 1.0012 5605
#> scale(day)       0.2899 0.0709  0.1507  0.2900 0.4273 1.0020 6000
#> scale(tod)      -0.0319 0.0697 -0.1680 -0.0327 0.1055 0.9999 6000
#> I(scale(day)^2) -0.0756 0.0867 -0.2476 -0.0752 0.0926 1.0000 6000
#> 
#> Spatial Covariance: 
#>            Mean     SD   2.5%    50%  97.5%   Rhat ESS
#> sigma.sq 1.2610 1.0218 0.2063 0.9694 4.0863 1.0130 101
#> phi      0.0093 0.0085 0.0009 0.0056 0.0294 1.0763  45
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
#> Samples per Chain: 10000
#> Burn-in: 2000
#> Thinning Rate: 4
#> Number of Chains: 3
#> Total Posterior Samples: 6000
#> 
#> Bayesian p-value:  0.4828 
#> Fit statistic:  freeman-tukey
```

### Model selection using WAIC and k-fold cross-validation

The `waicOcc` function computes the Widely Applicable Information
Criterion (WAIC) for use in model selection and assessment (note that
due to Monte Carlo error your results will differ slightly).

``` r
waicOcc(out)
#>       elpd         pD       WAIC 
#> -679.88774   22.96844 1405.71235
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
#> [1] 1414.417
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
[methods](https://www.jeffdoser.com/files/pubs/doser2024JABES.pdf)
(Doser et al. 2024A) and [applications to
ecology](https://onlinelibrary.wiley.com/doi/epdf/10.1111/geb.13814)
(Doser et al. 2024B) in much more detail.

## References

Doser, J. W., Finley, A. O., Kery, M., and Zipkin, E. F. (2022a).
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
