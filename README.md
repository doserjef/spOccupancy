
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spOccupancy <a href='https://www.jeffdoser.com/files/spoccupancy-web/'><img src="man/figures/logo.png" align="right" height="139"/></a>

[![](https://www.r-pkg.org/badges/version/spOccupancy?color=green)](https://cran.r-project.org/package=spOccupancy)
[![](http://cranlogs.r-pkg.org/badges/grand-total/spOccupancy?color=blue)](https://cran.r-project.org/package=spOccupancy)
[![](https://travis-ci.org/doserjef/spOccupancy.svg?branch=main)](https://travis-ci.org/doserjef/spOccupancy)

[![Codecov test
coverage](https://codecov.io/gh/doserjef/spOccupancy/branch/main/graph/badge.svg)](https://app.codecov.io/gh/doserjef/spOccupancy?branch=main)

spOccupancy fits single species, multispecies, and integrated spatial
occupancy models using Markov Chain Monte Carlo (MCMC). Models are fit
using Polya-Gamma data augmentation. Spatial models are fit using either
Gaussian processes or Nearest Neighbor Gaussian Processes (NNGP) for
large spatial datasets. Provides functionality for data integration of
multiple single species occupancy data sets using a joint likelihood
framework.

## Installation

You can install the released version of `spOccupancy` from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("spOccupancy")
```

## Functionality

| `spOccupancy` Function | Description                                                       |
| ---------------------- | ----------------------------------------------------------------- |
| `PGOcc`                | Single species occupancy model                                    |
| `spPGOcc`              | Single species spatial occupancy model                            |
| `intPGOcc`             | Single species occupancy model with multiple data sources         |
| `spIntPGOcc`           | Single species spatial occupancy model with multiple data sources |
| `msPGOcc`              | Multispecies occupancy model                                      |
| `spMsPGOcc`            | Multispecies spatial occupancy model                              |
| `ppcOcc`               | Posterior predictive check using Bayesian p-values                |
| `waicOcc`              | Compute Widely Applicable Information Criterion                   |
| `simOcc`               | Simulate single species occupancy data                            |
| `simMsOcc`             | Simulate multispecies occupancy data                              |
| `simIntOcc`            | Simulate single species occupancy data from multiple data sources |

## Example usage

### Load package and data

To get started with `spOccupancy` we load the package and an example
data set. We use data on twelve foliage-gleaning birds from the [Hubbard
Brook Experimental Forest](https://hubbardbrook.org/), which is
available in the `spOccupancy` package as the `hbef2015` object. Here we
will only work with one bird species, the Black-throated Blue Warbler
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

Below we fit a single species spatial occupancy model to the BTBW data
using a Nearest Neighbor Gaussian Process. We use the default priors and
initial values for the occurrence (`beta`) and regression (`alpha`)
coefficients, the spatial variance (`sigma.sq`), the spatial range
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

We run the model using an Adaptive MCMC sampler with a target acceptance
rate of 0.43. We run the model for 5000 iterations split into 200
batches each of length 25. We discard the first 2000 iterations as
burn-in and use a thinning rate of 3 for a resulting 1000 samples from
the joint posterior. We fit the model using 5 nearest neighbors and an
exponential correlation function. We also specify the `k.fold` argument
to perform 2-fold cross-validation after fitting the full model.

``` r
# Run the model
out <- spPGOcc(occ.formula = btbw.occ.formula,
               det.formula = btbw.det.formula,
               data = btbwHBEF, n.batch = 200, batch.length = 25,
               accept.rate = 0.43, cov.model = "exponential", 
               NNGP = TRUE, n.neighbors = 5, n.burn = 2000, 
               n.thin = 3, verbose = FALSE, k.fold = 2)
summary(out)
#> 
#> Call:
#> spPGOcc(occ.formula = btbw.occ.formula, det.formula = btbw.det.formula, 
#>     data = btbwHBEF, cov.model = "exponential", NNGP = TRUE, 
#>     n.neighbors = 5, n.batch = 200, batch.length = 25, accept.rate = 0.43, 
#>     verbose = FALSE, n.burn = 2000, n.thin = 3, k.fold = 2)
#> 
#> Chain Information:
#> Total samples: 5000
#> Burn-in: 2000
#> Thin: 3
#> Total Posterior Samples: 1000
#> 
#> Occurrence: 
#>                          2.5%     25%     50%     75%   97.5%
#> (Intercept)            3.2174  3.7268  4.0711  4.4382  5.1944
#> scale(Elevation)      -0.9919 -0.6673 -0.5048 -0.3682 -0.1265
#> I(scale(Elevation)^2) -1.5978 -1.3035 -1.1701 -1.0291 -0.8359
#> 
#> Detection: 
#>                    2.5%     25%     50%     75%  97.5%
#> (Intercept)      0.4494  0.5846  0.6680  0.7466 0.8972
#> scale(day)       0.1495  0.2402  0.2863  0.3335 0.4290
#> scale(tod)      -0.1685 -0.0777 -0.0311  0.0141 0.1028
#> I(scale(day)^2) -0.2375 -0.1345 -0.0784 -0.0163 0.0923
#> 
#> Covariance: 
#>            2.5%    25%    50%    75%  97.5%
#> sigma.sq 0.3997 0.8595 1.3070 2.0027 3.5470
#> phi      0.0016 0.0037 0.0076 0.0138 0.0262
```

### Posterior predictive check

The function `ppcOcc` performs a posterior predictive check on the
resulting list from the call to `spPGOcc`. For binary data, we first
need to perform Goodness of Fit assessments on some binned form of the
data rather than the raw binary data. Below we perform a posterior
predictive check on the data grouped by site with a Freeman-Tukey fit
statistic, and use the `summary` function to summarize the check with a
Bayesian p-value.

``` r
ppc.out <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out)
#> 
#> Call:
#> ppcOcc(object = out, fit.stat = "freeman-tukey", group = 1)
#> 
#> Chain Information:
#> Total samples: 5000
#> Burn-in: 2000
#> Thin: 3
#> Total Posterior Samples: 1000
#> 
#> Bayesian p-value:  0.408 
#> Fit statistic:  freeman-tukey
```

### Model selection using WAIC and k-fold cross-validation

The `waicOcc` function computes the Widely Applicable Information
Criterion (WAIC) for use in model selection and assessment.

``` r
waicOcc(out)
#>       elpd         pD       WAIC 
#> -675.88847   26.65018 1405.07730
```

Alternatively, we can perform k-fold cross-validation directly in our
call to `spPGOcc` using the `k.fold` argument and compare models using a
deviance scoring rule. We fit the model with `k.fold = 2` and so below
we access the devaince scoring rule from the 2-fold cross-validation.

``` r
out$k.fold.deviance
#> [1] 1506.063
```

### Prediction

Out-of-sample prediction is possible using the `predict` function, a set
of occurrence covariates at the new locations, and the spatial
coordinates of the new locations. The object `hbefElev` contains
elevation data across the entire Hubbard Brook Experimental Forest.
Below we predict BTBW occurrence across the forest, which are stored in
the `out.pred` object.

``` r
# First standardize elevation using mean and sd from fitted model
elev.pred <- (hbefElev$val - mean(btbwHBEF$occ.covs[, 1])) / sd(btbwHBEF$occ.covs[, 1])
coords.0 <- as.matrix(hbefElev[, c('Easting', 'Northing')])
X.0 <- cbind(1, elev.pred, elev.pred^2)
out.pred <- predict(out, X.0, coords.0, verbose = FALSE)
```

## Learn more

The `vignette("modelFitting")` provides a more detailed description and
tutorial of all functions in `spOccupancy`. For full statistical details
on the MCMC samplers used in `spOccupancy`, see
`vignette("mcmcSamplers")`.
