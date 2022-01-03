
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spOccupancy <a href='https://www.jeffdoser.com/files/spoccupancy-web/'><img src="man/figures/logo.png" align="right" height="139" width="120"/></a>

[![](https://www.r-pkg.org/badges/version/spOccupancy?color=green)](https://cran.r-project.org/package=spOccupancy)
[![](http://cranlogs.r-pkg.org/badges/grand-total/spOccupancy?color=blue)](https://cran.r-project.org/package=spOccupancy)
[![](https://codecov.io/gh/doserjef/spOccupancy/branch/main/graph/badge.svg)](https://app.codecov.io/gh/doserjef/spOccupancy)

spOccupancy fits single-species, multi-species, and integrated spatial
occupancy models using Markov Chain Monte Carlo (MCMC). Models are fit
using Póly-Gamma data augmentation. Spatial models are fit using either
Gaussian processes or Nearest Neighbor Gaussian Processes (NNGP) for
large spatial datasets. The package provides functionality for data
integration of multiple single-species occupancy data sets using a joint
likelihood framework. Below we provide a very brief introduction to some
of the package’s functionality, and illustrate just one of the model
fitting funcitons. For more information, see the resources referenced at
the bottom of this page.

## Installation

You can install the released version of `spOccupancy` from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("spOccupancy")
```

## Functionality

| `spOccupancy` Function | Description                                                       |
| ---------------------- | ----------------------------------------------------------------- |
| `PGOcc`                | Single-species occupancy model                                    |
| `spPGOcc`              | Single-species spatial occupancy model                            |
| `intPGOcc`             | Single-species occupancy model with multiple data sources         |
| `spIntPGOcc`           | Single-species spatial occupancy model with multiple data sources |
| `msPGOcc`              | Multi-species occupancy model                                     |
| `spMsPGOcc`            | Multi-species spatial occupancy model                             |
| `ppcOcc`               | Posterior predictive check using Bayesian p-values                |
| `waicOcc`              | Compute Widely Applicable Information Criterion (WAIC)            |
| `simOcc`               | Simulate single-species occupancy data                            |
| `simMsOcc`             | Simulate multi-species occupancy data                             |
| `simIntOcc`            | Simulate single-species occupancy data from multiple data sources |

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

Below we fit a single-species spatial occupancy model to the BTBW data
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
rate of 0.43. We run 3 chains of the model each for 10,000 iterations
split into 400 batches each of length 25. For each chain, we discard the
first 6000 iterations as burn-in and use a thinning rate of 4 for a
resulting 3000 samples from the joint posterior. We fit the model using
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
#> Run Time (min): 1.9316
#> 
#> Occurrence (logit scale): 
#>                          Mean     SD    2.5%     50%   97.5%   Rhat      ESS
#> (Intercept)            4.1679 0.6201  3.0796  4.1129  5.5253 1.0281 247.0541
#> scale(Elevation)      -0.5374 0.2535 -1.0717 -0.5319 -0.0500 1.0330 558.6439
#> I(scale(Elevation)^2) -1.2247 0.2262 -1.7324 -1.2051 -0.8373 1.0363 251.1803
#> 
#> Detection (logit scale): 
#>                    Mean     SD    2.5%     50%  97.5%   Rhat      ESS
#> (Intercept)      0.6628 0.1127  0.4479  0.6615 0.8826 1.0003 5010.936
#> scale(day)       0.2900 0.0712  0.1511  0.2903 0.4301 1.0016 6000.000
#> scale(tod)      -0.0312 0.0702 -0.1680 -0.0310 0.1059 1.0011 6000.000
#> I(scale(day)^2) -0.0750 0.0858 -0.2424 -0.0762 0.0958 0.9999 6000.000
#> 
#> Spatial Covariance: 
#>            Mean     SD   2.5%    50%  97.5%   Rhat     ESS
#> sigma.sq 1.7061 1.2382 0.3851 1.3643 5.0950 1.0312 96.9992
#> phi      0.0064 0.0071 0.0006 0.0032 0.0268 1.1779 43.3338
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
#> Bayesian p-value:  0.4026667 
#> Fit statistic:  freeman-tukey
```

### Model selection using WAIC and k-fold cross-validation

The `waicOcc` function computes the Widely Applicable Information
Criterion (WAIC) for use in model selection and assessment (note that
due to Monte Carlo error your results will differ slightly).

``` r
waicOcc(out)
#>       elpd         pD       WAIC 
#> -676.57521   25.20569 1403.56181
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
#> [1] 1496.396
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
tutorial of all functions in `spOccupancy`. For full statistical details
on the MCMC samplers used in `spOccupancy`, see
`vignette("mcmcSamplers")`. In addition, see [our recent
paper](https://arxiv.org/abs/2111.12163) that describes the package in
more detail (Doser et al. 2021).

## References

Doser, J. W., Finley, A. O., Kéry, M., and Zipkin, E. F. (2021a).
spOccupancy: An R package for single-species, multi-species, and
integrated spatial occupancy models. arXiv preprint arxiv:2111.12163.
