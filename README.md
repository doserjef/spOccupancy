
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spOccupancy: An R Package for single species, multispecies, and integrated spatial occupancy models

Fits single-species and multi-species non-spatial and spatial occupancy
models using Markov Chain Monte Carlo (MCMC). Models are fit using
Polya-Gamma data augmentation. Spatial models are fit using either
Gaussian processes or Nearest Neighbor Gaussian Processes (NNGP) for
large spatial datasets. Provides functionality for data integration of
multiple single species occupancy data sets using a joint likelihood
framework.

## Installation

Coming soon\!

<!--You can install the released version of `spOccupancy` from [CRAN](https://CRAN.R-project.org) with:-->

<!--Alternatively, you can install the development version from GitHub, although note this will only work if you have the proper compilation tools, as the `C/C++` code must be compiled:-->

## Functionality

| `spOccupancy` Function | Description                                                        |
| ---------------------- | ------------------------------------------------------------------ |
| `PGOcc`                | Single species occupancy model                                     |
| `spPGOcc`              | Single species spatial occupancy model                             |
| `intPGOcc`             | Single species occupancy model with multiple data sources          |
| `spIntPGOcc`           | Single species spatial occupancy model with multiplde data sources |
| `msPGOcc`              | Multispecies occupancy model                                       |
| `spMsPGOcc`            | Multispecies spatial occupancy model                               |
| `ppcOcc`               | Posterior predictive check using Bayesian p-values                 |
| `waicOcc`              | Compute Widely Applicable Information Criterion                    |
| `simOcc`               | Simulate single species occupancy data                             |
| `simMsOcc`             | Simulate multispecies occupancy data                               |
| `simIntOcc`            | Simulate single species occupancy data from multiple data sources  |

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

### Fit a spatial occupancy model using `spPGOcc`

Below we fit a single species spatial occupancy model to the BTBW data
using a Nearest Neighbor Gaussian Process. We use the default priors and
starting values for the occurrence (`beta`) and regression (`alpha`)
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
rate of 0.43. We run the model for 10,000 iterations split into 400
batches each of length 25. We discard the first 2000 iterations as
burn-in and use a thinning rate of 8 for a resulting 1000 samples from
the joint posterior. We fit the model using 5 nearest neighbors and an
exponential correlation function. We also specify the `k.fold` argument
to perform 4-fold cross-validation after fitting the full model.

``` r
# Run the model
out <- spPGOcc(occ.formula = btbw.occ.formula,
               det.formula = btbw.det.formula,
               data = btbwHBEF, n.batch = 400, batch.length = 25,
               accept.rate = 0.43, cov.model = "exponential", 
               NNGP = TRUE, n.neighbors = 5, n.burn = 2000, 
               n.thin = 8, verbose = FALSE, k.fold = 4)
summary(out)
#> 
#> Call:
#> spPGOcc(occ.formula = btbw.occ.formula, det.formula = btbw.det.formula, 
#>     data = btbwHBEF, cov.model = "exponential", NNGP = TRUE, 
#>     n.neighbors = 5, n.batch = 400, batch.length = 25, accept.rate = 0.43, 
#>     verbose = FALSE, n.burn = 2000, n.thin = 8, k.fold = 4)
#> 
#> Chain Information:
#> Total samples: 10000
#> Burn-in: 2000
#> Thin: 8
#> Total Posterior Samples: 1000
#> 
#> Occurrence: 
#>                          2.5%     25%     50%     75%   97.5%
#> (Intercept)            3.1631  3.7915  4.1292  4.5246  5.4862
#> scale(Elevation)      -1.0463 -0.6883 -0.5478 -0.4107 -0.1273
#> I(scale(Elevation)^2) -1.6625 -1.3297 -1.1768 -1.0550 -0.8096
#> 
#> Detection: 
#>                    2.5%     25%     50%     75%  97.5%
#> (Intercept)      0.4601  0.5906  0.6635  0.7402 0.8923
#> scale(day)       0.1574  0.2470  0.2900  0.3391 0.4312
#> scale(tod)      -0.1574 -0.0750 -0.0294  0.0203 0.1217
#> I(scale(day)^2) -0.2562 -0.1335 -0.0743 -0.0155 0.0963
#> 
#> Covariance: 
#>            2.5%    25%    50%    75%  97.5%
#> sigma.sq 0.4470 0.9317 1.3607 2.1120 4.1048
#> phi      0.0012 0.0038 0.0076 0.0148 0.0275
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
#> Total samples: 10000
#> Burn-in: 2000
#> Thin: 8
#> Total Posterior Samples: 1000
#> 
#> Bayesian p-value:  0.397 
#> Fit statistic:  freeman-tukey
```

### Model selection using WAIC and k-fold cross-validation

The `waicOcc` function computes the Widely Applicable InformatioN
Criterion (WAIC) for use in model selection and assessment.

``` r
waicOcc(out)
#>      elpd        pD      WAIC 
#> -675.0033   28.0920 1406.1906
```

Alternatively, we can perform k-fold cross-validation directly in our
call to `spPGOcc` using the `k.fold` argument and compare models using a
deviance scoring rule. We fit the model with `k.fold = 4` and so below
we access the devaince scoring rule from the 4-fold cross-validation.

``` r
out$k.fold.deviance
#> [1] 1509.138
```

### Prediction

Out-of-sample prediction is possible using the `predict` function, a set
of occurrence covariates at the new locations, and the spatial
coordinates of the new locations. The object `hbefElev` contains
elevation data across the entire Hubbard Brook Experimental Forest.

``` r
# First standardize elevation using mean and sd from fitted model
elev.pred <- (hbefElev$val - mean(btbwHBEF$occ.covs[, 1])) / sd(btbwHBEF$occ.covs[, 1])
coords.0 <- as.matrix(hbefElev[, c('Easting', 'Northing')])
X.0 <- cbind(1, elev.pred, elev.pred^2)
out.pred <- predict(out, X.0, coords.0, verbose = FALSE)
```

## Learn more

The package vignette “Fitting occupancy models with spOccupancy”
provides a more detailed description and tutorial of all functions in
`spOccupancy`. For full statistical details on the MCMC samplers used in
`spOccupancy`, see `MCMC samplers for models fit in spOccupancy`
vignette.
